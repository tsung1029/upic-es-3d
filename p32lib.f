c 3d parallel PIC library for MPI communications
c with 2D domain decomposition
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: june 16, 2008
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
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
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
c special case of only one grid per processor
      if (kyp.eq.1) then
         krr = kr + 1
         if (krr.ge.nvpy) krr = krr - nvpy
         krr = krr + kz
         kll = kl - 1
         if (kll.lt.0) kll = kll + nvpy
         kll = kll + kz
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
      call MPI_IRECV(scs(1,1,4,m),nxvz,mreal,kl-1,noff+3,lgrp,msid,ierr)
      call MPI_SEND(scs(1,1,1,m),nxvz,mreal,kr-1,noff+3,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scs(1,1,5,m),ngc*nxvz,mreal,kr-1,noff+4,lgrp,msid,i
     1err)
      call MPI_SEND(scs(1,1,2,m),ngc*nxvz,mreal,kl-1,noff+4,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      if (kyp.eq.1) then
         call MPI_IRECV(scs(1,1,6,m),ngc*nxvz,mreal,krr-1,noff+6,lgrp,ms
     1id,ierr)
         call MPI_SEND(scs(1,1,3,m),ngc*nxvz,mreal,kll-1,noff+6,lgrp,ier
     1r)
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
      call MPI_IRECV(f(1,1,1,m),nxvy,mreal,kl-1,noff+7,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,kzp+1,m),nxvy,mreal,kr-1,noff+7,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(f(1,1,kzp+2,m),ngc*nxvy,mreal,kr-1,noff+8,lgrp,msid
     1,ierr)
      call MPI_SEND(f(1,1,2,m),ngc*nxvy,mreal,kl-1,noff+8,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      if (kzp.eq.1) then
         call MPI_IRECV(f(1,1,kzp+3,m),ngc*nxvy,mreal,krr-1,noff+10,lgrp
     1,msid,ierr)
         call MPI_SEND(f(1,1,2,m),ngc*nxvy,mreal,kll-1,noff+10,lgrp,ierr
     1)
         call MPI_WAIT(msid,istatus,ierr)
      endif
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNCGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nzp
     1mx,mblok,nblok,ngds,idds,mter,nter)
c this subroutine copies data to guard cells in non-uniform partitions
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes three extra guard cells.
c scs/scr = scratch arrays for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c mter/nter = (0,1) = (no,yes) pass data to next processor only in y/z
c quadratic interpolation, for distributed data,
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, mblok, nblok, ngds
      integer idds, mter, nter
      integer nyzp
      real f, scs, scr
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx*ngds,2,mblok*nblok)
      dimension scr(nxv,nypmx*ngds,2,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, nps
      integer mnblok, nxvz, nxvzs, nyzp3, nxvy, nxvys, m, my, mz, j, k
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
      do 20 k = 1, nyzp(2,m)
      do 10 j = 1, nxv
      scs(j,k,1,m) = f(j,nyzp(1,m)+1,k+1,m)
      scs(j,k+nyzp(2,m),1,m) = f(j,2,k+1,m)
      scs(j,k+2*nyzp(2,m),1,m) = f(j,3,k+1,m)
   10 continue
   20 continue
   30 continue
c copy to guard cells in y
      do 120 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 110 my = 1, mblok
      m = my + moff
      nxvzs = nxv*nyzp(2,m)
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      kr = ky + 1
      if (kr.ge.nvpy) kr = kr - nvpy
      krr = kr + 1
      if (krr.ge.nvpy) krr = krr - nvpy
      krr = krr + kz
      kl = ky - 1
      if (kl.lt.0) kl = kl + nvpy
      kll = kl - 1
      if (kll.lt.0) kll = kll + nvpy
      kll = kll + kz
      kr = kr + kz
      kl = kl + kz
      ngc = 0
c special case of only one grid per processor
      if (nyzp(1,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (nyzp(1,kr).eq.1) then
c        do 50 k = 1, nyzp(2,m)
c        do 40 j = 1, nxv
c        scs(j,k,2,m) = scs(j,k,1,kl)
c        scs(j,k+nyzp(2,m),2,m) = scs(j,k+nyzp(2,m),1,kr)
c        scs(j,k+2*nyzp(2,m),2,m) = scs(j,k+nyzp(2,m),1,krr)
c  40    continue
c  50    continue
c     else
c        do 70 k = 1, nyzp(2,m)
c        do 60 j = 1, nxv
c        scs(j,k,2,m) = scs(j,k,1,kl)
c        scs(j,k+nyzp(2,m),2,m) = scs(j,k+nyzp(2,m),1,kr)
c        scs(j,k+2*nyzp(2,m),2,m) = scs(j,k+2*nyzp(2,m),1,kr)
c  60 continue
c  70 continue
c     endif
c this segment is used for mpi computers
      call MPI_IRECV(scs(1,1,2,m),nxvz,mreal,kl-1,noff+3,lgrp,msid,ierr)
      call MPI_SEND(scs(1,1,1,m),nxvzs,mreal,kr-1,noff+3,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scs(1,nyzp(2,m)+1,2,m),2*nxvz,mreal,kr-1,noff+4,lgr
     1p,msid,ierr)
      call MPI_SEND(scs(1,nyzp(2,m)+1,1,m),(2-ngc)*nxvzs,mreal,kl-1,noff
     1+4,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c special case of only one grid per processor in y
      if (mter.ge.1) go to 80
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      if (nps.eq.nxvzs) then
         call MPI_IRECV(scs(1,2*nyzp(2,m)+1,2,m),nxvz,mreal,krr-1,noff+6
     1,lgrp,msid,ierr)
      else
         call MPI_IRECV(scs(1,1,1,m),nxvz,mreal,krr-1,noff+6,lgrp,msid,i
     1err)
      endif
      call MPI_SEND(scs(1,nyzp(2,m)+1,1,m),nxvzs,mreal,kll-1,noff+6,lgrp
     1,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c copy guard cells
   80 do 100 k = 1, nyzp(2,m)
      do 90 j = 1, nxv
      f(j,1,k+1,m) = scs(j,k,2,m)
      f(j,nyzp(1,m)+2,k+1,m) = scs(j,k+nyzp(2,m),2,m)
      f(j,nyzp(1,m)+3,k+1,m) = scs(j,k+2*nyzp(2,m),2,m)
   90 continue
  100 continue
  110 continue
  120 continue
c buffer data in z
      do 150 m = 1, mnblok
      nyzp3 = nyzp(1,m) + 3
      do 140 k = 1, nyzp3
      do 130 j = 1, nxv
      scr(j,k,1,m) = f(j,k,nyzp(2,m)+1,m)
      scr(j,k+nyzp3,1,m) = f(j,k,2,m)
      scr(j,k+2*nyzp3,1,m) = f(j,k,3,m)
  130 continue
  140 continue
  150 continue
c copy to guard cells in z
      do 240 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 230 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(1,m) + 3
      nxvys = nxv*nyzp3
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      krr = kr + 1
      if (krr.ge.nvpz) krr = krr - nvpz
      krr = ky + nvpy*krr
      kr = ky + nvpy*kr
      kl = kz - 1
      if (kl.lt.0) kl = kl + nvpz
      kll = kl - 1
      if (kll.lt.0) kll = kll + nvpz
      kll = ky + nvpy*kll
      kl = ky + nvpy*kl
      ngc = 0
c special case of only one grid per processor
      if (nyzp(2,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (nyzp(2,kr).eq.1) then
c        do 170 k = 1, nyzp3
c        do 160 j = 1, nxv
c        scr(j,k,2,m) = scr(j,k,1,kl)
c        scr(j,k+nyzp3,2,m) = scr(j,k+nyzp3,1,kr)
c        scr(j,k+2*nyzp3,2,m) = scr(j,k+nyzp3,1,krr)
c 160    continue
c 170    continue
c     else
c        do 190 k = 1, nyzp3
c        do 180 j = 1, nxv
c        scr(j,k,2,m) = scr(j,k,1,kl)
c        scr(j,k+nyzp3,2,m) = scr(j,k+nyzp3,1,kr)
c        scr(j,k+2*nyzp3,2,m) = scr(j,k+2*nyzp3,1,kr)
c 180    continue
c 190    continue
c     endif
c this segment is used for mpi computers
      call MPI_IRECV(scr(1,1,2,m),nxvy,mreal,kl-1,noff+7,lgrp,msid,ierr)
      call MPI_SEND(scr(1,1,1,m),nxvys,mreal,kr-1,noff+7,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scr(1,nyzp3+1,2,m),2*nxvy,mreal,kr-1,noff+8,lgrp,ms
     1id,ierr)
      call MPI_SEND(scr(1,nyzp3+1,1,m),(2-ngc)*nxvys,mreal,kl-1,noff+8,l
     1grp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c special case of only one grid per processor in z
      if (nter.ge.1) go to 200
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      if (nps.eq.nxvys) then
         call MPI_IRECV(scr(1,2*nyzp3+1,2,m),nxvy,mreal,krr-1,noff+10,lg
     1rp,msid,ierr)
      else
         call MPI_IRECV(scr(1,1,1,m),nxvy,mreal,krr-1,noff+10,lgrp,msid,
     1ierr)
      endif
      call MPI_SEND(scr(1,nyzp3+1,1,m),nxvys,mreal,kll-1,noff+10,lgrp,ie
     1rr)
      call MPI_WAIT(msid,istatus,ierr)
c copy guard cells
  200 do 220 k = 1, nyzp3
      do 210 j = 1, nxv
      f(j,k,1,m) = scr(j,k,2,m)
      f(j,k,nyzp(2,m)+2,m) = scr(j,k+nyzp3,2,m)
      f(j,k,nyzp(2,m)+3,m) = scr(j,k+2*nyzp3,2,m)
  210 continue
  220 continue
  230 continue
  240 continue
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
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
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
      call MPI_IRECV(scs(1,1,2,m),nxvz,mreal,kr-1,noff+3,lgrp,msid,ierr)
      call MPI_SEND(scs(1,1,1,m),nxvz,mreal,kl-1,noff+3,lgrp,ierr)
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
      call MPI_IRECV(f(1,1,kzp+1,m),nxvy,mreal,kr-1,noff+4,lgrp,msid,ier
     1r)
      call MPI_SEND(f(1,1,1,m),nxvy,mreal,kl-1,noff+4,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,
     1mblok,nblok,ngds,idds)
c this subroutine copies data to guard cells in non-uniform partitions
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes one extra guard cell.
c scs(j,k,m) = scratch array for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c linear interpolation, for distributed data,
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, mblok, nblok, ngds
      integer idds
      integer nyzp
      real f, scs
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer ky, kz, js, ks, moff, noff, kr, kl, mnblok
      integer nxvz, nxvzs, nyzp1, nxvy, nxvys, m, my, mz, j, k
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
      do 20 k = 1, nyzp(2,m)
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
      nxvzs = nxv*nyzp(2,m)
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      kr = ky + 1
      if (kr.ge.nvpy) kr = kr - nvpy
      kl = ky - 1
      if (kl.lt.0) kl = kl + nvpy
      kr = kr + kz
      kl = kl + kz
c this segment is used for shared memory computers
c     do 50 k = 1, nyzp(2,m)
c     do 40 j = 1, nxv
c     scs(j,k,2,m) = scs(j,k,1,kr)
c  40 continue
c  50 continue
c this segment is used for mpi computers
      call MPI_IRECV(scs(1,1,2,m),nxvz,mreal,kr-1,noff+3,lgrp,msid,ierr)
      call MPI_SEND(scs(1,1,1,m),nxvzs,mreal,kl-1,noff+3,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c copy guard cells
      do 70 k = 1, nyzp(2,m)
      do 60 j = 1, nxv
      f(j,nyzp(1,m)+1,k,m) = scs(j,k,2,m)
   60 continue
   70 continue
   80 continue
   90 continue
c copy to guard cells in z
      do 130 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 120 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(1,m) + 1
      nxvys = nxv*nyzp1
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      kl = kz - 1
      if (kl.lt.0) kl = kl + nvpz
      kr = ky + nvpy*kr
      kl = ky + nvpy*kl
c this segment is used for shared memory computers
c     do 110 k = 1, nyzp1
c     do 100 j = 1, nxv
c     f(j,k,nyzp(2,m)+1,m) = f(j,k,1,kr)
c 100 continue
c 110 continue
c this segment is used for mpi computers
      call MPI_IRECV(f(1,1,nyzp(2,m)+1,m),nxvy,mreal,kr-1,noff+4,lgrp,ms
     1id,ierr)
      call MPI_SEND(f(1,1,1,m),nxvys,mreal,kl-1,noff+4,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PACGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx
     1,mblok,nblok,kyp,kzp,ngds)
c thus subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes three extra
c guard cells.
c scs/scr = scratch arrays for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, kyp, kzp
      real f, scs, scr
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(3,nxv,nypmx,ngds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, mnblok
      integer nx3, kyp3, kzp3, nxvz, nxvy, m, my, mz, j, k, n
      dimension istatus(lstat)
      nx3 = nx + 3
      kyp3 = kyp + 3
      kzp3 = kzp + 3
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c special case for one processor in y
      if (nvpy.eq.1) then
         do 40 m = 1, mnblok
         do 30 k = 1, kzp3
         do 20 j = 1, nx3
         do 10 n = 1, 3
         f(n,j,2,k,m) = f(n,j,2,k,m) + f(n,j,kyp+2,k,m)
         f(n,j,3,k,m) = f(n,j,3,k,m) + f(n,j,kyp+3,k,m)
         f(n,j,kyp+1,k,m) = f(n,j,kyp+1,k,m) + f(n,j,1,k,m)
   10    continue
   20    continue
   30    continue
   40    continue
         go to 170
      endif
c buffer data in y
      do 80 m = 1, mnblok
      do 70 k = 1, nzpmx
      do 60 j = 1, nxv
      do 50 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,kyp+2,k,m)
      scs(n,j,k,2,m) = f(n,j,kyp+3,k,m)
      scs(n,j,k,3,m) = f(n,j,1,k,m)
   50 continue
   60 continue
   70 continue
   80 continue
c add guard cells in y
      do 160 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 150 my = 1, mblok
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
c     do 110 k = 1, nzpmx
c     do 100 j = 1, nxv
c     do 90 n = 1, 3
c     scs(n,j,k,4,m) = scs(n,j,k,1,kl)
c     scs(n,j,k,5,m) = scs(n,j,k,2,kll)
c     scs(n,j,k,6,m) = scs(n,j,k,3,kr)
c  90 continue
c 100 continue
c 110 continue
c this segment is used for mpi computers
      call MPI_IRECV(scs(1,1,1,4,m),3*ngc*nxvz,mreal,kl-1,noff+1,lgrp,ms
     1id,ierr)
      call MPI_SEND(scs(1,1,1,1,m),3*ngc*nxvz,mreal,kr-1,noff+1,lgrp,ier
     1r)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scs(1,1,1,6,m),3*nxvz,mreal,kr-1,noff+2,lgrp,msid,i
     1err)
      call MPI_SEND(scs(1,1,1,3,m),3*nxvz,mreal,kl-1,noff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      if (kyp.eq.1) then
         call MPI_IRECV(scs(1,1,1,5,m),3*ngc*nxvz,mreal,kll-1,noff+5,lgr
     1p,msid,ierr)
         call MPI_SEND(scs(1,1,1,2,m),3*ngc*nxvz,mreal,krr-1,noff+5,lgrp
     1,ierr)
         call MPI_WAIT(msid,istatus,ierr)
      endif
      do 140 k = 1, kzp3
      do 130 j = 1, nx3
      do 120 n = 1, 3
      f(n,j,2,k,m) = f(n,j,2,k,m) + scs(n,j,k,4,m)
      f(n,j,ngc+1,k,m) = f(n,j,ngc+1,k,m) + scs(n,j,k,5,m)
      f(n,j,kyp+1,k,m) = f(n,j,kyp+1,k,m) + scs(n,j,k,6,m)
  120 continue
  130 continue
  140 continue
  150 continue
  160 continue
c special case for one processor in z
  170 if (nvpz.eq.1) then
         do 210 m = 1, mnblok
         do 200 k = 1, kyp3
         do 190 j = 1, nx3
         do 180 n = 1, 3
         f(n,j,k,2,m) = f(n,j,k,2,m) + f(n,j,k,kzp+2,m)
         f(n,j,k,3,m) = f(n,j,k,3,m) + f(n,j,k,kzp+3,m)
         f(n,j,k,kzp+1,m) = f(n,j,k,kzp+1,m) + f(n,j,k,1,m)
  180    continue
  190    continue
  200    continue
  210    continue
         return
      endif
c add guard cells in z
      do 290 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 280 my = 1, mblok
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
c     do 240 k = 1, nypmx
c     do 230 j = 1, nxv
c     do 220 n = 1, 3
c     scr(n,j,k,1,m) = f(n,j,k,kzp+2,kl)
c     scr(n,j,k,2,m) = f(n,j,k,kzp+3,kll)
c     scr(n,j,k,3,m) = f(n,j,k,1,kr)
c 220 continue
c 230 continue
c 240 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,3*ngc*nxvy,mreal,kl-1,noff+3,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,1,kzp+2,m),3*ngc*nxvy,mreal,kr-1,noff+3,lgrp,i
     1err)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scr(1,1,1,3,m),3*nxvy,mreal,kr-1,noff+4,lgrp,msid,i
     1err)
      call MPI_SEND(f(1,1,1,1,m),3*nxvy,mreal,kl-1,noff+4,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      if (kzp.eq.1) then
         call MPI_IRECV(scr(1,1,1,2,m),3*ngc*nxvy,mreal,kll-1,noff+6,lgr
     1p,msid,ierr)
         call MPI_SEND(f(1,1,1,kzp+3,m),3*ngc*nxvy,mreal,krr-1,noff+6,lg
     1rp,ierr)
         call MPI_WAIT(msid,istatus,ierr)
      endif
c add up the guard cells
      do 270 k = 1, kyp3
      do 260 j = 1, nx3
      do 250 n = 1, 3
      f(n,j,k,2,m) = f(n,j,k,2,m) + scr(n,j,k,1,m)
      f(n,j,k,ngc+1,m) = f(n,j,k,ngc+1,m) + scr(n,j,k,2,m)
      f(n,j,k,kzp+1,m) = f(n,j,k,kzp+1,m) + scr(n,j,k,3,m)
  250 continue
  260 continue
  270 continue
  280 continue
  290 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PAGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx,
     1mblok,nblok,kyp,kzp,ngds)
c thus subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes three extra
c guard cells.
c scs/scr = scratch arrays for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, kyp, kzp
      real f, scs, scr
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(nxv,nypmx,ngds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, mnblok
      integer nx3, kyp3, kzp3, nxvz, nxvy, m, my, mz, j, k
      dimension istatus(lstat)
      nx3 = nx + 3
      kyp3 = kyp + 3
      kzp3 = kzp + 3
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c special case for one processor in y
      if (nvpy.eq.1) then
         do 30 m = 1, mnblok
         do 20 k = 1, kzp3
         do 10 j = 1, nx3
         f(j,2,k,m) = f(j,2,k,m) + f(j,kyp+2,k,m)
         f(j,3,k,m) = f(j,3,k,m) + f(j,kyp+3,k,m)
         f(j,kyp+1,k,m) = f(j,kyp+1,k,m) + f(j,1,k,m)
   10    continue
   20    continue
   30    continue
         go to 130
      endif
c buffer data in y
      do 60 m = 1, mnblok
      do 50 k = 1, nzpmx
      do 40 j = 1, nxv
      scs(j,k,1,m) = f(j,kyp+2,k,m)
      scs(j,k,2,m) = f(j,kyp+3,k,m)
      scs(j,k,3,m) = f(j,1,k,m)
   40 continue
   50 continue
   60 continue
c add guard cells in y
      do 120 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 110 my = 1, mblok
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
c     do 80 k = 1, nzpmx
c     do 70 j = 1, nxv
c     scs(j,k,4,m) = scs(j,k,1,kl)
c     scs(j,k,5,m) = scs(j,k,2,kll)
c     scs(j,k,6,m) = scs(j,k,3,kr)
c  70 continue
c  80 continue
c this segment is used for mpi computers
      call MPI_IRECV(scs(1,1,4,m),ngc*nxvz,mreal,kl-1,noff+1,lgrp,msid,i
     1err)
      call MPI_SEND(scs(1,1,1,m),ngc*nxvz,mreal,kr-1,noff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scs(1,1,6,m),nxvz,mreal,kr-1,noff+2,lgrp,msid,ierr)
      call MPI_SEND(scs(1,1,3,m),nxvz,mreal,kl-1,noff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      if (kyp.eq.1) then
         call MPI_IRECV(scs(1,1,5,m),ngc*nxvz,mreal,kll-1,noff+5,lgrp,ms
     1id,ierr)
         call MPI_SEND(scs(1,1,2,m),ngc*nxvz,mreal,krr-1,noff+5,lgrp,ier
     1r)
         call MPI_WAIT(msid,istatus,ierr)
      endif
      do 100 k = 1, kzp3
      do 90 j = 1, nx3
      f(j,2,k,m) = f(j,2,k,m) + scs(j,k,4,m)
      f(j,ngc+1,k,m) = f(j,ngc+1,k,m) + scs(j,k,5,m)
      f(j,kyp+1,k,m) = f(j,kyp+1,k,m) + scs(j,k,6,m)
   90 continue
  100 continue
  110 continue
  120 continue
c special case for one processor in z
  130 if (nvpz.eq.1) then
         do 160 m = 1, mnblok
         do 150 k = 1, kyp3
         do 140 j = 1, nx3
         f(j,k,2,m) = f(j,k,2,m) + f(j,k,kzp+2,m)
         f(j,k,3,m) = f(j,k,3,m) + f(j,k,kzp+3,m)
         f(j,k,kzp+1,m) = f(j,k,kzp+1,m) + f(j,k,1,m)
  140    continue
  150    continue
  160    continue
         return
      endif
c add guard cells in z
      do 220 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 210 my = 1, mblok
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
c     do 180 k = 1, nypmx
c     do 170 j = 1, nxv
c     scr(j,k,1,m) = f(j,k,kzp+2,kl)
c     scr(j,k,2,m) = f(j,k,kzp+3,kll)
c     scr(j,k,3,m) = f(j,k,1,kr)
c 170 continue
c 180 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,ngc*nxvy,mreal,kl-1,noff+3,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,kzp+2,m),ngc*nxvy,mreal,kr-1,noff+3,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scr(1,1,3,m),nxvy,mreal,kr-1,noff+4,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,1,m),nxvy,mreal,kl-1,noff+4,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      if (kzp.eq.1) then
         call MPI_IRECV(scr(1,1,2,m),ngc*nxvy,mreal,kll-1,noff+6,lgrp,ms
     1id,ierr)
         call MPI_SEND(f(1,1,kzp+3,m),ngc*nxvy,mreal,krr-1,noff+6,lgrp,i
     1err)
         call MPI_WAIT(msid,istatus,ierr)
      endif
c add up the guard cells
      do 200 k = 1, kyp3
      do 190 j = 1, nx3
      f(j,k,2,m) = f(j,k,2,m) + scr(j,k,1,m)
      f(j,k,ngc+1,m) = f(j,k,ngc+1,m) + scr(j,k,2,m)
      f(j,k,kzp+1,m) = f(j,k,kzp+1,m) + scr(j,k,3,m)
  190 continue
  200 continue
  210 continue
  220 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNACGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypmx
     1,nzpmx,mblok,nblok,ngds,idds,mter,nter)
c this subroutine adds data from guard cells in non-uniform partitions
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes three extra guard cells.
c scs/scr = scratch arrays for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c mter/nter = (0,1) = (no,yes) pass data to next processor only in y/z
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds, mter, nter
      integer nyzp
      real f, scs, scr
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx*ngds,2,mblok*nblok)
      dimension scr(3,nxv,nypmx*ngds,2,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, mnblok
      integer nx3, nxvz, nxvzs, nyzp3, nxvy, nxvys, m, my, mz, j, k, n
      dimension istatus(lstat)
      nx3 = nx + 3
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c special case for one processor in y
      if (nvpy.eq.1) then
         do 40 m = 1, mnblok
         nyzp3 = nyzp(2,m) + 3
         do 30 k = 1, nyzp3
         do 20 j = 1, nx3
         do 10 n = 1, 3
         f(n,j,2,k,m) = f(n,j,2,k,m) + f(n,j,nyzp(1,m)+2,k,m)
         f(n,j,3,k,m) = f(n,j,3,k,m) + f(n,j,nyzp(1,m)+3,k,m)
         f(n,j,nyzp(1,m)+1,k,m) = f(n,j,nyzp(1,m)+1,k,m) + f(n,j,1,k,m)
         f(n,j,1,k,m) = 0.
         f(n,j,nyzp(1,m)+2,k,m) = 0.
         f(n,j,nyzp(1,m)+3,k,m) = 0.
   10    continue
   20    continue
   30    continue
   40    continue
         go to 240
      endif
c buffer data in y
      do 80 m = 1, mnblok
      nyzp3 = nyzp(2,m) + 3
      do 70 k = 1, nyzp3
      do 60 j = 1, nxv
      do 50 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,nyzp(1,m)+2,k,m)
      scs(n,j,k+nyzp3,1,m) = f(n,j,nyzp(1,m)+3,k,m)
      scs(n,j,k+2*nyzp3,1,m) = f(n,j,1,k,m)
   50 continue
   60 continue
   70 continue
   80 continue
c add guard cells in y
      do 230 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 220 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(2,m) + 3
      nxvzs = nxv*nyzp3
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      kr = ky + 1
      if (kr.ge.nvpy) kr = kr - nvpy
      krr = kr + 1
      if (krr.ge.nvpy) krr = krr - nvpy
      krr = krr + kz
      kl = ky - 1
      if (kl.lt.0) kl = kl + nvpy
      kll = kl - 1
      if (kll.lt.0) kll = kll + nvpy
      kll = kll + kz
      kr = kr + kz
      kl = kl + kz
      ngc = 0
c special case of only one grid per processor
      if (nyzp(1,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (nyzp(1,kl).eq.1) then
c        do 110 k = 1, nyzp3
c        do 100 j = 1, nxv
c        do 90 n = 1, 3
c        scs(n,j,k,2,m) = scs(n,j,k,1,kl) + scs(n,j,k+nyzp3,1,kll)
c        scs(n,j,k+nyzp3,2,m) = scs(n,j,k+nyzp3,1,kl)
c        scs(n,j,k+2*nyzp3,2,m) = scs(n,j,k+2*nyzp3,1,kr)
c  90    continue
c 100    continue
c 110    continue
c     else
c        do 140 k = 1, nyzp3
c        do 130 j = 1, nxv
c        do 120 n = 1, 3
c        scs(n,j,k,2,m) = scs(n,j,k,1,kl)
c        scs(n,j,k+nyzp3,2,m) = scs(n,j,k+nyzp3,1,kl)
c        scs(n,j,k+2*nyzp3,2,m) = scs(n,j,k+2*nyzp3,1,kr)
c 120    continue
c 130    continue
c 140    continue
c     endif
c this segment is used for mpi computers
      call MPI_IRECV(scs(1,1,1,2,m),6*nxvz,mreal,kl-1,noff+1,lgrp,msid,i
     1err)
      call MPI_SEND(scs(1,1,1,1,m),6*nxvzs,mreal,kr-1,noff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scs(1,1,2*nyzp3+1,2,m),3*nxvz,mreal,kr-1,noff+2,lgr
     1p,msid,ierr)
      call MPI_SEND(scs(1,1,2*nyzp3+1,1,m),3*nxvzs,mreal,kl-1,noff+2,lgr
     1p,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c special case of only one grid per processor in y
      if (mter.ge.1) go to 180
      call MPI_IRECV(ngc,1,mint,kl-1,noff+3,lgrp,msid,ierr)
      call MPI_IRECV(scs(1,1,1,1,m),3*nxvz,mreal,kll-1,noff+4,lgrp,nsid,
     1ierr)
      call MPI_SEND(nyzp(1,m),1,mint,kr-1,noff+3,lgrp,ierr)
      call MPI_SEND(scs(1,1,nyzp3+1,1,m),3*nxvzs,mreal,krr-1,noff+4,lgrp
     1,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_WAIT(nsid,istatus,ierr)
      if (ngc.eq.1) then
         do 170 k = 1, nyzp3
         do 160 j = 1, nx3
         do 150 n = 1, 3
         scs(n,j,k,2,m) = scs(n,j,k,2,m) + scs(n,j,k,1,m)
  150    continue
  160    continue
  170    continue
      endif
c add up the guard cells
  180 do 210 k = 1, nyzp3
      do 200 j = 1, nx3
      do 190 n = 1, 3
      f(n,j,2,k,m) = f(n,j,2,k,m) + scs(n,j,k,2,m)
      f(n,j,3,k,m) = f(n,j,3,k,m) + scs(n,j,k+nyzp3,2,m)
      f(n,j,nyzp(1,m)+1,k,m) = f(n,j,nyzp(1,m)+1,k,m) + scs(n,j,k+2*nyzp
     13,2,m)
      f(n,j,1,k,m) = 0.
      f(n,j,nyzp(1,m)+2,k,m) = 0.
      f(n,j,nyzp(1,m)+3,k,m) = 0.
  190 continue
  200 continue
  210 continue
  220 continue
  230 continue
c special case for one processor in z
  240 if (nvpz.eq.1) then
         do 280 m = 1, mnblok
         nyzp3 = nyzp(1,m) + 3
         do 270 k = 1, nyzp3
         do 260 j = 1, nx3
         do 250 n = 1, 3
         f(n,j,k,2,m) = f(n,j,k,2,m) + f(n,j,k,nyzp(2,m)+2,m)
         f(n,j,k,3,m) = f(n,j,k,3,m) + f(n,j,k,nyzp(2,m)+3,m)
         f(n,j,k,nyzp(2,m)+1,m) = f(n,j,k,nyzp(2,m)+1,m) + f(n,j,k,1,m)
         f(n,j,k,1,m) = 0.
         f(n,j,k,nyzp(2,m)+2,m) = 0.
         f(n,j,k,nyzp(2,m)+3,m) = 0.
  250    continue
  260    continue
  270    continue
  280    continue
         return
      endif
c buffer data in z
      do 320 m = 1, mnblok
      nyzp3 = nyzp(1,m) + 3
      do 310 k = 1, nyzp3
      do 300 j = 1, nxv
      do 290 n = 1, 3
      scr(n,j,k,1,m) = f(n,j,k,nyzp(2,m)+2,m)
      scr(n,j,k+nyzp3,1,m) = f(n,j,k,nyzp(2,m)+3,m)
      scr(n,j,k+2*nyzp3,1,m) = f(n,j,k,1,m)
  290 continue
  300 continue
  310 continue
  320 continue
c add guard cells in z
      do 470 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 460 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(1,m) + 3
      nxvys = nxv*nyzp3
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      krr = kr + 1
      if (krr.ge.nvpz) krr = krr - nvpz
      krr = ky + nvpy*krr
      kr = ky + nvpy*kr
      kl = kz - 1
      if (kl.lt.0) kl = kl + nvpz
      kll = kl - 1
      if (kll.lt.0) kll = kll + nvpz
      kll = ky + nvpy*kll
      kl = ky + nvpy*kl
      ngc = 0
c special case of only one grid per processor
      if (nyzp(2,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (nyzp(2,kl).eq.1) then
c        do 350 k = 1, nyzp3
c        do 340 j = 1, nxv
c        do 330 n = 1, 3
c        scr(n,j,k,2,m) = scr(n,j,k,1,kl) + scr(n,j,k+nyzp3,1,kll)
c        scr(n,j,k+nyzp3,2,m) = scr(n,j,k+nyzp3,1,kl)
c        scr(n,j,k+2*nyzp3,2,m) = scr(n,j,k+2*nyzp3,1,kr)
c 330    continue
c 340    continue
c 350    continue
c     else
c        do 380 k = 1, nyzp3
c        do 370 j = 1, nxv
c        do 360 n = 1, 3
c        scr(n,j,k,2,m) = scr(n,j,k,1,kl)
c        scr(n,j,k+nyzp3,2,m) = scr(n,j,k+nyzp3,1,kl)
c        scr(n,j,k+2*nyzp3,2,m) = scr(n,j,k+2*nyzp3,1,kr)
c 360    continue
c 370    continue
c 380    continue
c     endif
c this segment is used for mpi computers
      call MPI_IRECV(scr(1,1,1,2,m),6*nxvy,mreal,kl-1,noff+5,lgrp,msid,i
     1err)
      call MPI_SEND(scr(1,1,1,1,m),6*nxvys,mreal,kr-1,noff+5,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scr(1,1,2*nyzp3+1,2,m),3*nxvy,mreal,kr-1,noff+6,lgr
     1p,msid,i
     1err)
      call MPI_SEND(scr(1,1,2*nyzp3+1,1,m),3*nxvys,mreal,kl-1,noff+6,lgr
     1p,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c special case of only one grid per processor in z
      if (nter.ge.1) go to 420
      call MPI_IRECV(ngc,1,mint,kl-1,noff+7,lgrp,msid,ierr)
      call MPI_IRECV(scr(1,1,1,1,m),3*nxvy,mreal,kll-1,noff+8,lgrp,nsid,
     1ierr)
      call MPI_SEND(nyzp(2,m),1,mint,kr-1,noff+7,lgrp,ierr)
      call MPI_SEND(scr(1,1,nyzp3+1,1,m),3*nxvys,mreal,krr-1,noff+8,lgrp
     1,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_WAIT(nsid,istatus,ierr)
      if (ngc.eq.1) then
         do 410 k = 1, nyzp3
         do 400 j = 1, nx3
         do 390 n = 1, 3
         scr(n,j,k,2,m) = scr(n,j,k,2,m) + scr(n,j,k,1,m)
  390    continue
  400    continue
  410    continue
      endif
c add up the guard cells
  420 do 450 k = 1, nyzp3
      do 440 j = 1, nx3
      do 430 n = 1, 3
      f(n,j,k,2,m) = f(n,j,k,2,m) + scr(n,j,k,2,m)
      f(n,j,k,3,m) = f(n,j,k,3,m) + scr(n,j,k+nyzp3,2,m)
      f(n,j,k,nyzp(2,m)+1,m) = f(n,j,k,nyzp(2,m)+1,m) + scr(n,j,k+2*nyzp
     13,2,m)
      f(n,j,k,1,m) = 0.
      f(n,j,k,nyzp(2,m)+2,m) = 0.
      f(n,j,k,nyzp(2,m)+3,m) = 0.
  430 continue
  440 continue
  450 continue
  460 continue
  470 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNAGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypmx,
     1nzpmx,mblok,nblok,ngds,idds,mter,nter)
c this subroutine adds data from guard cells in non-uniform partitions
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes three extra guard cells.
c scs/scr = scratch arrays for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c mter/nter = (0,1) = (no,yes) pass data to next processor only in y/z
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds, mter, nter
      integer nyzp
      real f, scs, scr
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx*ngds,2,mblok*nblok)
      dimension scr(nxv,nypmx*ngds,2,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, mnblok
      integer nx3, nxvz, nxvzs, nyzp3, nxvy, nxvys, m, my, mz, j, k
      dimension istatus(lstat)
      nx3 = nx + 3
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c special case for one processor in y
      if (nvpy.eq.1) then
         do 30 m = 1, mnblok
         nyzp3 = nyzp(2,m) + 3
         do 20 k = 1, nyzp3
         do 10 j = 1, nx3
         f(j,2,k,m) = f(j,2,k,m) + f(j,nyzp(1,m)+2,k,m)
         f(j,3,k,m) = f(j,3,k,m) + f(j,nyzp(1,m)+3,k,m)
         f(j,nyzp(1,m)+1,k,m) = f(j,nyzp(1,m)+1,k,m) + f(j,1,k,m)
         f(j,1,k,m) = 0.
         f(j,nyzp(1,m)+2,k,m) = 0.
         f(j,nyzp(1,m)+3,k,m) = 0.
   10    continue
   20    continue
   30    continue
         go to 180
      endif
c buffer data in y
      do 60 m = 1, mnblok
      nyzp3 = nyzp(2,m) + 3
      do 50 k = 1, nyzp3
      do 40 j = 1, nxv
      scs(j,k,1,m) = f(j,nyzp(1,m)+2,k,m)
      scs(j,k+nyzp3,1,m) = f(j,nyzp(1,m)+3,k,m)
      scs(j,k+2*nyzp3,1,m) = f(j,1,k,m)
   40 continue
   50 continue
   60 continue
c add guard cells in y
      do 170 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 160 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(2,m) + 3
      nxvzs = nxv*nyzp3
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      kr = ky + 1
      if (kr.ge.nvpy) kr = kr - nvpy
      krr = kr + 1
      if (krr.ge.nvpy) krr = krr - nvpy
      krr = krr + kz
      kl = ky - 1
      if (kl.lt.0) kl = kl + nvpy
      kll = kl - 1
      if (kll.lt.0) kll = kll + nvpy
      kll = kll + kz
      kr = kr + kz
      kl = kl + kz
      ngc = 0
c special case of only one grid per processor
      if (nyzp(1,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (nyzp(1,kl).eq.1) then
c        do 80 k = 1, nyzp3
c        do 70 j = 1, nxv
c        scs(j,k,2,m) = scs(j,k,1,kl) + scs(j,k+nyzp3,1,kll)
c        scs(j,k+nyzp3,2,m) = scs(j,k+nyzp3,1,kl)
c        scs(j,k+2*nyzp3,2,m) = scs(j,k+2*nyzp3,1,kr)
c  70    continue
c  80    continue
c     else
c        do 100 k = 1, nyzp3
c        do 90 j = 1, nxv
c        scs(j,k,2,m) = scs(j,k,1,kl)
c        scs(j,k+nyzp3,2,m) = scs(j,k+nyzp3,1,kl)
c        scs(j,k+2*nyzp3,2,m) = scs(j,k+2*nyzp3,1,kr)
c  90    continue
c 100    continue
c     endif
c this segment is used for mpi computers
      call MPI_IRECV(scs(1,1,2,m),2*nxvz,mreal,kl-1,noff+1,lgrp,msid,ier
     1r)
      call MPI_SEND(scs(1,1,1,m),2*nxvzs,mreal,kr-1,noff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scs(1,2*nyzp3+1,2,m),nxvz,mreal,kr-1,noff+2,lgrp,ms
     1id,ierr)
      call MPI_SEND(scs(1,2*nyzp3+1,1,m),nxvzs,mreal,kl-1,noff+2,lgrp,ie
     1rr)
      call MPI_WAIT(msid,istatus,ierr)
c special case of only one grid per processor in y
      if (mter.ge.1) go to 130
      call MPI_IRECV(ngc,1,mint,kl-1,noff+3,lgrp,msid,ierr)
      call MPI_IRECV(scs(1,1,1,m),nxvz,mreal,kll-1,noff+4,lgrp,nsid,ierr
     1)
      call MPI_SEND(nyzp(1,m),1,mint,kr-1,noff+3,lgrp,ierr)
      call MPI_SEND(scs(1,nyzp3+1,1,m),nxvzs,mreal,krr-1,noff+4,lgrp,ier
     1r)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_WAIT(nsid,istatus,ierr)
      if (ngc.eq.1) then
         do 120 k = 1, nyzp3
         do 110 j = 1, nx3
         scs(j,k,2,m) = scs(j,k,2,m) + scs(j,k,1,m)
  110    continue
  120    continue
      endif
c add up the guard cells
  130 do 150 k = 1, nyzp3
      do 140 j = 1, nx3
      f(j,2,k,m) = f(j,2,k,m) + scs(j,k,2,m)
      f(j,3,k,m) = f(j,3,k,m) + scs(j,k+nyzp3,2,m)
      f(j,nyzp(1,m)+1,k,m) = f(j,nyzp(1,m)+1,k,m) + scs(j,k+2*nyzp3,2,m)
      f(j,1,k,m) = 0.
      f(j,nyzp(1,m)+2,k,m) = 0.
      f(j,nyzp(1,m)+3,k,m) = 0.
  140 continue
  150 continue
  160 continue
  170 continue
c special case for one processor in z
  180 if (nvpz.eq.1) then
         do 210 m = 1, mnblok
         nyzp3 = nyzp(1,m) + 3
         do 200 k = 1, nyzp3
         do 190 j = 1, nx3
         f(j,k,2,m) = f(j,k,2,m) + f(j,k,nyzp(2,m)+2,m)
         f(j,k,3,m) = f(j,k,3,m) + f(j,k,nyzp(2,m)+3,m)
         f(j,k,nyzp(2,m)+1,m) = f(j,k,nyzp(2,m)+1,m) + f(j,k,1,m)
         f(j,k,1,m) = 0.
         f(j,k,nyzp(2,m)+2,m) = 0.
         f(j,k,nyzp(2,m)+3,m) = 0.
  190    continue
  200    continue
  210    continue
         return
      endif
c buffer data in z
      do 240 m = 1, mnblok
      nyzp3 = nyzp(1,m) + 3
      do 230 k = 1, nyzp3
      do 220 j = 1, nxv
      scr(j,k,1,m) = f(j,k,nyzp(2,m)+2,m)
      scr(j,k+nyzp3,1,m) = f(j,k,nyzp(2,m)+3,m)
      scr(j,k+2*nyzp3,1,m) = f(j,k,1,m)
  220 continue
  230 continue
  240 continue
c add guard cells in z
      do 350 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 340 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(1,m) + 3
      nxvys = nxv*nyzp3
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      krr = kr + 1
      if (krr.ge.nvpz) krr = krr - nvpz
      krr = ky + nvpy*krr
      kr = ky + nvpy*kr
      kl = kz - 1
      if (kl.lt.0) kl = kl + nvpz
      kll = kl - 1
      if (kll.lt.0) kll = kll + nvpz
      kll = ky + nvpy*kll
      kl = ky + nvpy*kl
      ngc = 0
c special case of only one grid per processor
      if (nyzp(2,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (nyzp(2,kl).eq.1) then
c        do 260 k = 1, nyzp3
c        do 250 j = 1, nxv
c        scr(j,k,2,m) = scr(j,k,1,kl) + scr(j,k+nyzp3,1,kll)
c        scr(j,k+nyzp3,2,m) = scr(j,k+nyzp3,1,kl)
c        scr(j,k+2*nyzp3,2,m) = scr(j,k+2*nyzp3,1,kr)
c 250    continue
c 260    continue
c     else
c        do 280 k = 1, nyzp3
c        do 270 j = 1, nxv
c        scr(j,k,2,m) = scr(j,k,1,kl)
c        scr(j,k+nyzp3,2,m) = scr(j,k+nyzp3,1,kl)
c        scr(j,k+2*nyzp3,2,m) = scr(j,k+2*nyzp3,1,kr)
c 270    continue
c 280    continue
c     endif
c this segment is used for mpi computers
      call MPI_IRECV(scr(1,1,2,m),2*nxvy,mreal,kl-1,noff+5,lgrp,msid,ier
     1r)
      call MPI_SEND(scr(1,1,1,m),2*nxvys,mreal,kr-1,noff+5,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scr(1,2*nyzp3+1,2,m),nxvy,mreal,kr-1,noff+6,lgrp,ms
     1id,ierr)
      call MPI_SEND(scr(1,2*nyzp3+1,1,m),nxvys,mreal,kl-1,noff+6,lgrp,ie
     1rr)
      call MPI_WAIT(msid,istatus,ierr)
c special case of only one grid per processor in z
      if (nter.ge.1) go to 310
      call MPI_IRECV(ngc,1,mint,kl-1,noff+7,lgrp,msid,ierr)
      call MPI_IRECV(scr(1,1,1,m),nxvy,mreal,kll-1,noff+8,lgrp,nsid,ierr
     1)
      call MPI_SEND(nyzp(2,m),1,mint,kr-1,noff+7,lgrp,ierr)
      call MPI_SEND(scr(1,nyzp3+1,1,m),nxvys,mreal,krr-1,noff+8,lgrp,ier
     1r)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_WAIT(nsid,istatus,ierr)
      if (ngc.eq.1) then
         do 300 k = 1, nyzp3
         do 290 j = 1, nx3
         scr(j,k,2,m) = scr(j,k,2,m) + scr(j,k,1,m)
  290    continue
  300    continue
      endif
c add up the guard cells
  310 do 330 k = 1, nyzp3
      do 320 j = 1, nx3
      f(j,k,2,m) = f(j,k,2,m) + scr(j,k,2,m)
      f(j,k,3,m) = f(j,k,3,m) + scr(j,k+nyzp3,2,m)
      f(j,k,nyzp(2,m)+1,m) = f(j,k,nyzp(2,m)+1,m) + scr(j,k+2*nyzp3,2,m)
      f(j,k,1,m) = 0.
      f(j,k,nyzp(2,m)+2,m) = 0.
      f(j,k,nyzp(2,m)+3,m) = 0.
  320 continue
  330 continue
  340 continue
  350 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PACGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpm
     1x,mblok,nblok,kyp,kzp,ngds)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes one extra guard
c cell.
c scs/scr = scratch array for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c linear interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, kyp, kzp
      real f, scs, scr
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(3,nxv,nypmx,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer ky, kz, js, ks, moff, noff, kr, kl, mnblok
      integer nx1, kyp1, kzp1, nxvz, nxvy, m, my, mz, j, k, n
      dimension istatus(lstat)
      nx1 = nx + 1
      kyp1 = kyp + 1
      kzp1 = kzp + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c special case for one processor in y
      if (nvpy.eq.1) then
         do 40 m = 1, mnblok
         do 30 k = 1, kzp1
         do 20 j = 1, nx1
         do 10 n = 1, 3
         f(n,j,1,k,m) = f(n,j,1,k,m) + f(n,j,kyp+1,k,m)
   10    continue
   20    continue
   30    continue
   40    continue
         go to 170
      endif
c buffer data in y
      do 80 m = 1, mnblok
      do 70 k = 1, nzpmx
      do 60 j = 1, nxv
      do 50 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,kyp+1,k,m)
   50 continue
   60 continue
   70 continue
   80 continue
c add guard cells in y
      do 160 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 150 my = 1, mblok
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
c     do 110 k = 1, nzpmx
c     do 100 j = 1, nxv
c     do 90 n = 1, 3
c     scs(n,j,k,2,m) = scs(n,j,k,1,kl)
c  90 continue
c 100 continue
c 110 continue
c this segment is used for mpi computers
      call MPI_IRECV(scs(1,1,1,2,m),3*nxvz,mreal,kl-1,noff+1,lgrp,msid,i
     1err)
      call MPI_SEND(scs(1,1,1,1,m),3*nxvz,mreal,kr-1,noff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      do 140 k = 1, kzp1
      do 130 j = 1, nx1
      do 120 n = 1, 3
      f(n,j,1,k,m) = f(n,j,1,k,m) + scs(n,j,k,2,m)
  120 continue
  130 continue
  140 continue
  150 continue
  160 continue
c special case for one processor in z
  170 if (nvpz.eq.1) then
         do 210 m = 1, mnblok
         do 200 k = 1, kyp1
         do 190 j = 1, nx1
         do 180 n = 1, 3
         f(n,j,k,1,m) = f(n,j,k,1,m) + f(n,j,k,kzp+1,m)
  180    continue
  190    continue
  200    continue
  210    continue
         return
      endif
c add guard cells in z
      do 290 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 280 my = 1, mblok
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
c     do 240 k = 1, nypmx
c     do 230 j = 1, nxv
c     do 220 n = 1, 3
c     scr(n,j,k,m) = f(n,j,k,kzp+1,kl)
c 220 continue
c 230 continue
c 240 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,3*nxvy,mreal,kl-1,noff+2,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,1,kzp+1,m),3*nxvy,mreal,kr-1,noff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      do 270 k = 1, kyp1
      do 260 j = 1, nx1
      do 250 n = 1, 3
      f(n,j,k,1,m) = f(n,j,k,1,m) + scr(n,j,k,m)
  250 continue
  260 continue
  270 continue
  280 continue
  290 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PAGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx
     1,mblok,nblok,kyp,kzp,ngds)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes one extra guard
c cell.
c scs/scr = scratch array for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c linear interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, kyp, kzp
      real f, scs, scr
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(nxv,nypmx,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer ky, kz, js, ks, moff, noff, kr, kl, mnblok
      integer nx1, kyp1, kzp1, nxvz, nxvy, m, my, mz, j, k
      dimension istatus(lstat)
      nx1 = nx + 1
      kyp1 = kyp + 1
      kzp1 = kzp + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c special case for one processor in y
      if (nvpy.eq.1) then
         do 30 m = 1, mnblok
         do 20 k = 1, kzp1
         do 10 j = 1, nx1
         f(j,1,k,m) = f(j,1,k,m) + f(j,kyp+1,k,m)
   10    continue
   20    continue
   30    continue
         go to 130
      endif
c buffer data in y
      do 60 m = 1, mnblok
      do 50 k = 1, nzpmx
      do 40 j = 1, nxv
      scs(j,k,1,m) = f(j,kyp+1,k,m)
   40 continue
   50 continue
   60 continue
c add guard cells in y
      do 120 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 110 my = 1, mblok
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
c     do 80 k = 1, nzpmx
c     do 70 j = 1, nxv
c     scs(j,k,2,m) = scs(j,k,1,kl)
c  70 continue
c  80 continue
c this segment is used for mpi computers
      call MPI_IRECV(scs(1,1,2,m),nxvz,mreal,kl-1,noff+1,lgrp,msid,ierr)
      call MPI_SEND(scs(1,1,1,m),nxvz,mreal,kr-1,noff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      do 100 k = 1, kzp1
      do 90 j = 1, nx1
      f(j,1,k,m) = f(j,1,k,m) + scs(j,k,2,m)
   90 continue
  100 continue
  110 continue
  120 continue
c special case for one processor in z
  130 if (nvpz.eq.1) then
         do 160 m = 1, mnblok
         do 150 k = 1, kyp1
         do 140 j = 1, nx1
         f(j,k,1,m) = f(j,k,1,m) + f(j,k,kzp+1,m)
  140    continue
  150    continue
  160    continue
         return
      endif
c add guard cells in z
      do 220 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 210 my = 1, mblok
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
c     do 180 k = 1, nypmx
c     do 170 j = 1, nxv
c     scr(j,k,m) = f(j,k,kzp+1,kl)
c 170 continue
c 180 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,nxvy,mreal,kl-1,noff+2,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,kzp+1,m),nxvy,mreal,kr-1,noff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      do 200 k = 1, kyp1
      do 190 j = 1, nx1
      f(j,k,1,m) = f(j,k,1,m) + scr(j,k,m)
  190 continue
  200 continue
  210 continue
  220 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNACGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypm
     1x,nzpmx,mblok,nblok,ngds,idds)
c this subroutine adds data from guard cells in non-uniform partitions
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes one extra guard cell.
c scs/scr = scratch arrays for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c linear interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds
      integer nyzp
      real f, scs, scr
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(3,nxv,nypmx,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer ky, kz, js, ks, moff, noff, kr, kl, mnblok
      integer nx1, nxvz, nxvzs, nyzp1, nxvy, nxvys, m, my, mz, j, k, n
      dimension istatus(lstat)
      nx1 = nx + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c special case for one processor in y
      if (nvpy.eq.1) then
         do 40 m = 1, mnblok
         do 30 k = 1, nyzp(2,m)+1
         do 20 j = 1, nx1
         do 10 n = 1, 3
         f(n,j,1,k,m) = f(n,j,1,k,m) + f(n,j,nyzp(1,m)+1,k,m)
         f(n,j,nyzp(1,m)+1,k,m) = 0.
   10    continue
   20    continue
   30    continue
   40    continue
         go to 170
      endif
c buffer data in y
      do 80 m = 1, mnblok
      do 70 k = 1, nyzp(2,m)+1
      do 60 j = 1, nxv
      do 50 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,nyzp(1,m)+1,k,m)
   50 continue
   60 continue
   70 continue
   80 continue
c add guard cells in y
      do 160 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 150 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(2,m) + 1
      nxvzs = nxv*nyzp1
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      kr = ky + 1
      if (kr.ge.nvpy) kr = kr - nvpy
      kl = ky - 1
      if (kl.lt.0) kl = kl + nvpy
      kr = kr + kz
      kl = kl + kz
c this segment is used for shared memory computers
c     do 110 k = 1, nyzp1
c     do 100 j = 1, nxv
c     do 90 n = 1, 3
c     scs(n,j,k,2,m) = scs(n,j,k,1,kl)
c  90 continue
c 100 continue
c 110 continue
c this segment is used for mpi computers
      call MPI_IRECV(scs(1,1,1,2,m),3*nxvz,mreal,kl-1,noff+1,lgrp,msid,i
     1err)
      call MPI_SEND(scs(1,1,1,1,m),3*nxvzs,mreal,kr-1,noff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
      do 140 k = 1, nyzp1
      do 130 j = 1, nx1
      do 120 n = 1, 3
      f(n,j,1,k,m) = f(n,j,1,k,m) + scs(n,j,k,2,m)
      f(n,j,nyzp(1,m)+1,k,m) = 0.
  120 continue
  130 continue
  140 continue
  150 continue
  160 continue
c special case for one processor in z
  170 if (nvpz.eq.1) then
         do 210 m = 1, mnblok
         do 200 k = 1, nyzp(1,m)+1
         do 190 j = 1, nx1
         do 180 n = 1, 3
         f(n,j,k,1,m) = f(n,j,k,1,m) + f(n,j,k,nyzp(2,m)+1,m)
         f(n,j,k,nyzp(2,m)+1,m) = 0.
  180    continue
  190    continue
  200    continue
  210    continue
         return
      endif
c add guard cells in z
      do 290 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 280 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(1,m) + 1
      nxvys = nxv*nyzp1
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      kl = kz - 1
      if (kl.lt.0) kl = kl + nvpz
      kr = ky + nvpy*kr
      kl = ky + nvpy*kl
c this segment is used for shared memory computers
c     do 240 k = 1, nyzp1
c     do 230 j = 1, nxv
c     do 220 n = 1, 3
c     scr(n,j,k,m) = f(n,j,k,nyzp(2,m)+1,kl)
c 220 continue
c 230 continue
c 240 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,3*nxvy,mreal,kl-1,noff+2,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,1,nyzp(2,m)+1,m),3*nxvys,mreal,kr-1,noff+2,lgr
     1p,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
      do 270 k = 1, nyzp1
      do 260 j = 1, nx1
      do 250 n = 1, 3
      f(n,j,k,1,m) = f(n,j,k,1,m) + scr(n,j,k,m)
      f(n,j,k,nyzp(2,m)+1,m) = 0.
  250 continue
  260 continue
  270 continue
  280 continue
  290 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNAGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypmx
     1,nzpmx,mblok,nblok,ngds,idds)
c this subroutine adds data from guard cells in non-uniform partitions
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes one extra guard cell.
c scs/scr = scratch arrays for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c linear interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds
      integer nyzp
      real f, scs, scr
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(nxv,nypmx,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer ky, kz, js, ks, moff, noff, kr, kl, mnblok
      integer nx1, nxvz, nxvzs, nyzp1, nxvy, nxvys, m, my, mz, j, k
      dimension istatus(lstat)
      nx1 = nx + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c special case for one processor in y
      if (nvpy.eq.1) then
         do 30 m = 1, mnblok
         do 20 k = 1, nyzp(2,m)+1
         do 10 j = 1, nx1
         f(j,1,k,m) = f(j,1,k,m) + f(j,nyzp(1,m)+1,k,m)
         f(j,nyzp(1,m)+1,k,m) = 0.
   10    continue
   20    continue
   30    continue
         go to 130
      endif
c buffer data in y
      do 60 m = 1, mnblok
      do 50 k = 1, nyzp(2,m)+1
      do 40 j = 1, nxv
      scs(j,k,1,m) = f(j,nyzp(1,m)+1,k,m)
   40 continue
   50 continue
   60 continue
c add guard cells in y
      do 120 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 110 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(2,m) + 1
      nxvzs = nxv*nyzp1
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      kr = ky + 1
      if (kr.ge.nvpy) kr = kr - nvpy
      kl = ky - 1
      if (kl.lt.0) kl = kl + nvpy
      kr = kr + kz
      kl = kl + kz
c this segment is used for shared memory computers
c     do 80 k = 1, nyzp1
c     do 70 j = 1, nxv
c     scs(j,k,2,m) = scs(j,k,1,kl)
c  70 continue
c  80 continue
c this segment is used for mpi computers
      call MPI_IRECV(scs(1,1,2,m),nxvz,mreal,kl-1,noff+1,lgrp,msid,ierr)
      call MPI_SEND(scs(1,1,1,m),nxvzs,mreal,kr-1,noff+1,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
      do 100 k = 1, nyzp1
      do 90 j = 1, nx1
      f(j,1,k,m) = f(j,1,k,m) + scs(j,k,2,m)
      f(j,nyzp(1,m)+1,k,m) = 0.
   90 continue
  100 continue
  110 continue
  120 continue
c special case for one processor in z
  130 if (nvpz.eq.1) then
         do 160 m = 1, mnblok
         do 150 k = 1, nyzp(1,m)+1
         do 140 j = 1, nx1
         f(j,k,1,m) = f(j,k,1,m) + f(j,k,nyzp(2,m)+1,m)
         f(j,k,nyzp(2,m)+1,m) = 0.
  140    continue
  150    continue
  160    continue
         return
      endif
c add guard cells in z
      do 220 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 210 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(1,m) + 1
      nxvys = nxv*nyzp1
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      kl = kz - 1
      if (kl.lt.0) kl = kl + nvpz
      kr = ky + nvpy*kr
      kl = ky + nvpy*kl
c this segment is used for shared memory computers
c     do 180 k = 1, nyzp1
c     do 170 j = 1, nxv
c     scr(j,k,m) = f(j,k,nyzp(2,m)+1,kl)
c 170 continue
c 180 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,nxvy,mreal,kl-1,noff+2,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,nyzp(2,m)+1,m),nxvys,mreal,kr-1,noff+2,lgrp,ie
     1rr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
      do 200 k = 1, nyzp1
      do 190 j = 1, nx1
      f(j,k,1,m) = f(j,k,1,m) + scr(j,k,m)
      f(j,k,nyzp(2,m)+1,m) = 0.
  190 continue
  200 continue
  210 continue
  220 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PTRPSIN32C(cu,cu3,scb,scd,nx,ny,nz,kstrt,nvpy,nvpz,nxv,
     1kyp,kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
c this subroutine creates a tripled vector array cu3 from a vector array
c cu, so that various 3d sine/cosine transforms can be performed with a
c 3d real to complex fft.  the x component is an odd function in y,
c y component is an odd function in x, and the z component is an odd
c function in both x and y.  Asummes vector cu vanishes at end points
c linear interpolation for distributed data with 2D domain decomposition
c cu3 array may be modified
c scb/scd = scratch arrays
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = second dimension of input array cu, must be >= nx
c kyp/kzp = number of data values per block in cu in y/z
c kypd = third dimension of input array cu, must be >= kyp
c kzpd = fourth dimension of input array cu, must be >= kzp
c kyp2/kzp2 = number of data values per block in cu3 in y/z
c kblok/lblok = number of data blocks in y/z
c k2blok/l2blok = number of data blocks in y/z for tripled data
      implicit none
      real cu, cu3, scb, scd
      integer nx, ny, nz, kstrt, nvpy, nvpz, nxv, kyp, kzp, kypd, kzpd
      integer kyp2, kzp2, kblok, lblok, k2blok, l2blok
      dimension cu(3,nxv,kypd,kzpd,kblok*lblok)
      dimension cu3(3,2*nxv,2*kypd,kzp2,k2blok*l2blok)
      dimension scb(3,nxv,kypd,kzpd,kblok*lblok)
      dimension scd(3,nxv,kzpd,2,kblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      logical lt1
      integer istatus, lsid, msid, nsid, ierr
      integer i, j, k, l, m, my, mz, nxs, nys, nzs, ks, js, ny2, nz2, kr
      integer kyb, kyb2, kzb, kzb2, nxvy, noff, moff, loff, koff, joff
      integer jm, km, jl, kl, ky, kz, kk, ll, mm, lm, l1, l2, k0, k1, k2
      integer ii
      real at1
      dimension istatus(lstat)
      nxs = nx - 1
      nys = ny - 1
      nzs = nz - 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      kyb = ny/kyp
      ny2 = ny + ny
      kyb2 = ny2/kyp2
      kzb = nz/kzp
      nz2 = nz + nz
      kzb2 = nz2/kzp2
      nxvy = nxv*kypd
      noff = kypd*kzpd + kyb*kzb
c copy to triple array in y direction
      do 200 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 190 my = 1, k2blok
      m = my + moff
      loff = kyp2*(my + js)
      jm = loff/kyp + 1
      loff = kzp2*(mz + ks)
      km = loff/kzp + 1
      loff = kyp*(my + js)
      jl = loff/kyp2 + 1
      loff = kzp*(mz + ks)
      kl = loff - kzp2*(loff/kzp2)
      kz = nvpy*(mz + ks)
c special case for one processor
      if ((kyb2*kzb2).eq.1) then
         do 60 l = 1, nzs
         do 30 k = 1, nys
         do 10 j = 1, nxs
         cu3(1,j+1,k+1,l+1,m) = cu(1,j+1,k+1,l+1,m)
         cu3(2,j+1,k+1,l+1,m) = cu(2,j+1,k+1,l+1,m)
         cu3(3,j+1,k+1,l+1,m) = cu(3,j+1,k+1,l+1,m)
         cu3(1,nx+j+1,k+1,l+1,m) = cu(1,nx-j+1,k+1,l+1,m)
         cu3(2,nx+j+1,k+1,l+1,m) = -cu(2,nx-j+1,k+1,l+1,m)
         cu3(3,nx+j+1,k+1,l+1,m) = -cu(3,nx-j+1,k+1,l+1,m)
         cu3(1,j+1,ny+k+1,l+1,m) = -cu(1,j+1,ny-k+1,l+1,m)
         cu3(2,j+1,ny+k+1,l+1,m) = cu(2,j+1,ny-k+1,l+1,m)
         cu3(3,j+1,ny+k+1,l+1,m) = -cu(3,j+1,ny-k+1,l+1,m)
         cu3(1,nx+j+1,ny+k+1,l+1,m) = -cu(1,nx-j+1,ny-k+1,l+1,m)
         cu3(2,nx+j+1,ny+k+1,l+1,m) = -cu(2,nx-j+1,ny-k+1,l+1,m)
         cu3(3,nx+j+1,ny+k+1,l+1,m) = cu(3,nx-j+1,ny-k+1,l+1,m)
         cu3(1,j+1,k+1,nz+l+1,m) = -cu(1,j+1,k+1,nz-l+1,m)
         cu3(2,j+1,k+1,nz+l+1,m) = -cu(2,j+1,k+1,nz-l+1,m)
         cu3(3,j+1,k+1,nz+l+1,m) = cu(3,j+1,k+1,nz-l+1,m)
         cu3(1,nx+j+1,k+1,nz+l+1,m) = -cu(1,nx-j+1,k+1,nz-l+1,m)
         cu3(2,nx+j+1,k+1,nz+l+1,m) = cu(2,nx-j+1,k+1,nz-l+1,m)
         cu3(3,nx+j+1,k+1,nz+l+1,m) = -cu(3,nx-j+1,k+1,nz-l+1,m)
         cu3(1,j+1,ny+k+1,nz+l+1,m) = cu(1,j+1,ny-k+1,nz-l+1,m)
         cu3(2,j+1,ny+k+1,nz+l+1,m) = -cu(2,j+1,ny-k+1,nz-l+1,m)
         cu3(3,j+1,ny+k+1,nz+l+1,m) = -cu(3,j+1,ny-k+1,nz-l+1,m)
         cu3(1,nx+j+1,ny+k+1,nz+l+1,m) = cu(1,nx-j+1,ny-k+1,nz-l+1,m)
         cu3(2,nx+j+1,ny+k+1,nz+l+1,m) = cu(2,nx-j+1,ny-k+1,nz-l+1,m)
         cu3(3,nx+j+1,ny+k+1,nz+l+1,m) = cu(3,nx-j+1,ny-k+1,nz-l+1,m)
   10    continue
         do 20 i = 1, 3
         cu3(i,1,k+1,l+1,m) = 0.
         cu3(i,nx+1,k+1,l+1,m) = 0.
         cu3(i,1,k+ny+1,l+1,m) = 0.
         cu3(i,nx+1,k+ny+1,l+1,m) = 0.
         cu3(i,1,k+1,nz+l+1,m) = 0.
         cu3(i,nx+1,k+1,nz+l+1,m) = 0.
         cu3(i,1,k+ny+1,nz+l+1,m) = 0.
         cu3(i,nx+1,k+ny+1,nz+l+1,m) = 0.
   20    continue
   30    continue
         do 50 j = 1, nx
         do 40 i = 1, 3
         cu3(i,j,1,l+1,m) = 0.
         cu3(i,j+nx,1,l+1,m) = 0.
         cu3(i,j,ny+1,l+1,m) = 0.
         cu3(i,j+nx,ny+1,l+1,m) = 0.
         cu3(i,j,1,nz+l+1,m) = 0.
         cu3(i,j+nx,1,nz+l+1,m) = 0.
         cu3(i,j,ny+1,nz+l+1,m) = 0.
         cu3(i,j+nx,ny+1,nz+l+1,m) = 0.
   40    continue
   50    continue
   60    continue
         nxs = nx + nx
         nys = ny + ny
         do 90 k = 1, nys
         do 80 j = 1, nxs
         do 70 i = 1, 3
         cu3(i,j,k,1,m) = 0.
         cu3(i,j,k,nz+1,m) = 0.
   70    continue
   80    continue
   90    continue
         return
      endif
c copy main data in y direction
      do 180 ii = 1, 2
c for odd rows in z
      if (ii.eq.1) then
         mm = jm + 2*kz
         lm = jl + kz/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in z
      else if (ii.eq.2) then
         if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 180
         mm = mm + nvpy
         lm = lm - nvpy/2
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kyb).and.(km.le.kzb)) then
c        do 130 l = 1, kzp
c        do 120 k = 1, kyp
c        do 110 j = 1, nx
c        do 100 i = 1, 3
c        cu3(i,j,k,l+loff,m) = cu(i,j,k,l,mm)
c 100    continue
c 110    continue
c 120    continue
c 130    continue
c        if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
c           do 170 l = 1, kzp
c           do 160 k = 1, kyp
c           do 150 j = 1, nx
c           do 140 i = 1, 3
c           cu3(i,j,k+kyp,l+loff,m) = cu(i,j,k,l,mm+1)
c 140       continue
c 150       continue
c 160       continue
c 170       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kyb).and.(km.le.kzb)) then
         call MPI_IRECV(cu3(1,1,1,loff+1,m),3*kzp*nxvy,mreal,mm-1,noff+i
     1i,lgrp,msid,ierr)
         if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
            call MPI_IRECV(scb(1,1,1,1,m),3*kzp*nxvy,mreal,mm,noff+ii,lg
     1rp,nsid,ierr)
         endif
      endif
      if ((jl.le.((kyb2-1)/2+1)).and.lt1) then
         call MPI_SEND(cu(1,1,1,1,m),3*kzp*nxvy,mreal,lm-1,noff+ii,lgrp,
     1ierr)
      endif
c wait for data and unpack it
      if ((jm.le.kyb).and.(km.le.kzb)) then
         call MPI_WAIT(msid,istatus,ierr)
         do 130 l = 1, kzp
         l1 = kzp - l + 1
         l2 = (l1 - 1)/4 + 1
         koff = kypd*(l1 - 4*(l2 - 1) - 1)
         do 120 k = 1, kypd
         k1 = kypd - k + 1
         k0 = k1 - 1 + koff
         k2 = k0/2 + 1
         joff = nxv*(k0 - 2*(k2 - 1))
         do 110 j = 1, nxv
         do 100 i = 1, 3
         cu3(i,j,k1,l1+loff,m) = cu3(i,j+joff,k2,l2+loff,m)
  100    continue
  110    continue
  120    continue
  130    continue
         if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 170 l = 1, kzp
            do 160 k = 1, kypd
            do 150 j = 1, nxv
            do 140 i = 1, 3
            cu3(i,j,k+kyp,l+loff,m) = scb(i,j,k,l,m)
  140       continue
  150       continue
  160       continue
  170       continue
         endif
      endif
  180 continue
  190 continue
  200 continue
c copying in y direction not needed
      if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 380
c copy reflected data in y direction
      do 370 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 360 my = 1, k2blok
      m = my + moff
      loff = kyp2*(my + js)
      jm = (ny2 - loff - 1)/kyp + 1
      loff = kzp2*(mz + ks)
      km = loff/kzp + 1
      loff = kyp*(my + js)
      jl = (ny2 - loff - 1)/kyp2 + 1
      ky = loff + kyp2*jl - ny2
      loff = kzp*(mz + ks)
      kl = loff - kzp2*(loff/kzp2)
      kz = nvpy*(mz + ks)
      do 350 ii = 3, 4
c for odd rows in z
      if (ii.eq.3) then
         mm = jm + 2*kz
         lm = jl + kz/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in z
      else if (ii.eq.4) then
         if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 350
         mm = mm + nvpy
         lm = lm - nvpy/2
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kyb).and.(km.le.kzb)) then
c        if ((jm+1).le.kyb) then
c           do 230 l = 1, kzp
c           do 220 j = 1, nx
c           do 210 i = 1, 3
c           cu3(i,j,1,l+loff,m) = cu(i,j,1,l,mm+1)
c 210       continue
c 220       continue
c 230       continue
c        endif
c        if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
c           do 270 l = 1, kzp
c           do 260 k = 1, kyp
c           do 250 j = 1, nx
c           do 240 i = 1, 3
c           cu3(i,j,k+kyp,l+loff,m) = cu(i,j,k,l,mm)
c 240       continue
c 250       continue
c 260       continue
c 270       continue
c        endif
c        if (kyp.gt.1) then
c           do 310 l = 1, kzp
c           do 300 k = 2, kyp
c           do 290 j = 1, nx
c           do 280 i = 1, 3
c           cu3(i,j,k,l+loff,m) = cu(i,j,k,l,mm-1)
c 280       continue
c 290       continue
c 300       continue
c 310       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if ((jm+1).le.kyb) then
            call MPI_IRECV(scd(1,1,1,2,m),3*kzp*nxv,mreal,mm,noff+ii,lgr
     1p,lsid,ierr)
         endif
         if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
            call MPI_IRECV(scb(1,1,1,1,m),3*kzp*nxvy,mreal,mm-1,noff+ii,
     1lgrp,msid,ierr)
         endif
         if (kyp.gt.1) then
            call MPI_IRECV(cu3(1,1,1,loff+1,m),3*kzp*nxvy,mreal,mm-2,nof
     1f+ii,lgrp,nsid,ierr)
         endif
      endif
      if ((jl.gt.((kyb2-1)/2+1)).and.(jl.le.kyb2).and.lt1) then
         if (ky.eq.0) then
            if ((jl+1).le.kyb2) then
               do 230 l = 1, kzp
               do 220 j = 1, nxv
               do 210 i = 1, 3
               scd(i,j,l,1,m) = cu(i,j,1,l,m)
  210          continue
  220          continue
  230          continue
               call MPI_SEND(scd(1,1,1,1,m),3*kzp*nxv,mreal,lm,noff+ii,l
     1grp,ierr)
            endif
            if (kyp.gt.1) then
               call MPI_SEND(cu(1,1,1,1,m),3*kzp*nxvy,mreal,lm-1,noff+ii
     1,lgrp,ierr)
            endif
         else
            call MPI_SEND(cu(1,1,1,1,m),3*kzp*nxvy,mreal,lm-1,noff+ii,lg
     1rp,ierr)
         endif
      endif
c wait for data and unpack it
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if (kyp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 270 l = 1, kzp
            l1 = kzp - l + 1
            l2 = (l1 - 1)/4 + 1
            koff = kypd*(l1 - 4*(l2 - 1) - 1)
            do 260 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1 + koff
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 250 j = 1, nxv
            do 240 i = 1, 3
            cu3(i,j,k1,l1+loff,m) = cu3(i,j+joff,k2,l2+loff,m)
  240       continue
  250       continue
  260       continue
  270       continue
         endif
         if ((jm+1).le.kyb) then
            call MPI_WAIT(lsid,istatus,ierr)
            do 300 l = 1, kzp
            do 290 j = 1, nxv
            do 280 i = 1, 3
            cu3(i,j,1,l+loff,m) = scd(i,j,l,2,m)
  280       continue
  290       continue
  300       continue
         endif
         if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 340 l = 1, kzp
            do 330 k = 1, kypd
            do 320 j = 1, nxv
            do 310 i = 1, 3
            cu3(i,j,k+kyp,l+loff,m) = scb(i,j,k,l,m)
  310       continue
  320       continue
  330       continue
  340       continue
         endif
      endif
  350 continue
  360 continue
  370 continue
c copying in z direction not needed
  380 if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 860
c copy reflected data in z direction
      do 560 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 550 my = 1, k2blok
      m = my + moff
      loff = kzp2*(mz + ks)
      jm = (nz2 - loff - 1)/kzp + 1
      loff = kyp2*(my + js)
      km = loff/kyp + 1
      loff = kzp*(mz + ks)
      jl = (nz2 - loff - 1)/kzp2 + 1
      kz = loff + kzp2*jl - nz2
      loff = kyp*(my + js)
      kl = loff - kyp2*(loff/kyp2)
      kr = my + js
      do 500 ii = 5, 6
c for odd rows in y
      if (ii.eq.5) then
         mm = nvpy*jm + 2*kr
         lm = nvpy*jl + kr/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in y
      else if (ii.eq.6) then
         if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 500
         mm = mm + 1
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kzb).and.(km.le.kyb)) then
c        if (kzp.gt.1) then
c           do 420 l = 2, kzp
c           do 410 k = 1, kypd
c           do 400 j = 1, nx
c           do 390 i = 1, 3
c           cu3(i,j,k,l+loff,m) = cu(i,j,k,l,mm-2*nvpy+1)
c 390       continue
c 400       continue
c 410       continue
c 420       continue
c        endif
c        if ((jm+1).le.kzb) then
c           do 450 k = 1, kypd
c           do 440 j = 1, nx
c           do 430 i = 1, 3
c           cu3(i,j,k,loff+1,m) = cu(i,j,k,1,mm+1)
c 430       continue
c 440       continue
c 450       continue
c        endif
c        if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
c           do 490 l = 1, kzp
c           do 480 k = 1, kypd
c           do 470 j = 1, nx
c           do 460 i = 1, 3
c           cu3(i,j,k+kyp,l+loff,m) = cu(i,j,k,l,mm-nvpy+1)
c 460       continue
c 470       continue
c 480       continue
c 490       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kzb).and.(km.le.kyb)) then
         if ((jm+1).le.kzb) then
            call MPI_IRECV(cu3(1,1,1,loff+1,m),3*nxvy,mreal,mm,noff+ii,l
     1grp,lsid,ierr)
         endif
         if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
            call MPI_IRECV(scb(1,1,1,1,m),3*kzp*nxvy,mreal,mm-nvpy,noff+
     1ii,lgrp,msid,ierr)
         endif
         if (kzp.gt.1) then
            call MPI_IRECV(cu3(1,1,1,loff+2,m),3*(kzp-1)*nxvy,mreal,mm-2
     1*nvpy,noff+ii,lgrp,nsid,ierr)
         endif
      endif
      if ((jl.gt.((kzb2-1)/2+1)).and.(jl.le.kzb2).and.lt1) then
         if (kz.eq.0) then
            if ((jl+1).le.kzb2) then
               call MPI_SEND(cu(1,1,1,1,m),3*nxvy,mreal,lm,noff+ii,lgrp,
     1ierr)
            endif
            if (kzp.gt.1) then
               call MPI_SEND(cu(1,1,1,2,m),3*(kzp-1)*nxvy,mreal,lm-nvpy,
     1noff+ii,lgrp,ierr)
            endif
         else
            call MPI_SEND(cu(1,1,1,1,m),3*kzp*nxvy,mreal,lm-nvpy,noff+ii
     1,lgrp,ierr)
         endif
      endif
c wait for data and unpack it
      if ((jm.le.kzb).and.(km.le.kyb)) then
         if (kzp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 420 l = 2, kzp
            l1 = kzp - l + 2
            l2 = (l1 + 2)/4 + 1
            koff = kypd*(l1 - 4*(l2 - 1) + 2)
            do 410 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1 + koff
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 400 j = 1, nxv
            do 390 i = 1, 3
            cu3(i,j,k1,l1+loff,m) = cu3(i,j+joff,k2,l2+loff,m)
  390       continue
  400       continue
  410       continue
  420       continue
         endif
         if ((jm+1).le.kzb) then
            call MPI_WAIT(lsid,istatus,ierr)
            do 450 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 440 j = 1, nxv
            do 430 i = 1, 3
            cu3(i,j,k1,loff+1,m) = cu3(i,j+joff,k2,loff+1,m)
  430       continue
  440       continue
  450       continue
         endif
         if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 490 l = 1, kzp
            do 480 k = 1, kypd
            do 470 j = 1, nxv
            do 460 i = 1, 3
            cu3(i,j,k+kyp,l+loff,m) = scb(i,j,k,l,m)
  460       continue
  470       continue
  480       continue
  490       continue
         endif
      endif
  500 continue
c switch internal data
      if ((jm.le.kzb).and.(km.le.kyb)) then
         do 540 l = 1, kzp
         do 530 k = 1, kyp
         do 520 j = 1, nxv
         do 510 i = 1, 3
         at1 = cu3(i,j,k+kyp,l,m)
         cu3(i,j,k+kyp,l,m) = cu3(i,j,k,l+kzp,m)
         cu3(i,j,k,l+kzp,m) = at1
  510    continue
  520    continue
  530    continue
  540    continue
      endif
  550 continue
  560 continue
c copying in y direction not needed
      if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 860
c copy reflected data in y and z direction
      do 740 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 730 my = 1, k2blok
      m = my + moff
      loff = kzp2*(mz + ks)
      jm = (nz2 - loff - 1)/kzp + 1
      loff = kyp2*(my + js)
      km = (ny2 - loff - 1)/kyp + 1
      loff = kzp*(mz + ks)
      jl = (nz2 - loff - 1)/kzp2 + 1
      kz = loff + kzp2*jl - nz2
      loff = kyp*(my + js)
      kl = loff - kyp2*(loff/kyp2)
      kr = nvpy - (my + js) - 1
      do 680 ii = 7, 8
c for odd rows in y
      if (ii.eq.7) then
         mm = nvpy*jm + 2*kr
         lm = nvpy*jl + (nvpy + kr)/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in y
      else if (ii.eq.8) then
         if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 680
         mm = mm + 1
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kzb).and.(km.le.kyb)) then
c        if (kzp.gt.1) then
c           do 600 l = 2, kzp
c           do 590 k = 1, kypd
c           do 580 j = 1, nx
c           do 570 i = 1, 3
c           cu3(i,j,k,l+loff,m) = cu(i,j,k,l,mm-2*nvpy+1)
c 570       continue
c 580       continue
c 590       continue
c 600       continue
c        endif
c        if ((jm+1).le.kzb) then
c           do 630 k = 1, kypd
c           do 620 j = 1, nx
c           do 610 i = 1, 3
c           cu3(i,j,k,loff+1,m) = cu(i,j,k,1,mm+1)
c 610       continue
c 620       continue
c 630       continue
c        endif
c        if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
c           do 670 l = 1, kzp
c           do 660 k = 1, kypd
c           do 650 j = 1, nx
c           do 640 i = 1, 3
c           cu3(i,j,k+kyp,l+loff,m) = cu(i,j,k,l,mm-nvpy+1)
c 640       continue
c 650       continue
c 660       continue
c 670       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kzb).and.(km.le.kyb)) then
         if ((jm+1).le.kzb) then
            call MPI_IRECV(cu3(1,1,1,loff+1,m),3*nxvy,mreal,mm,noff+ii,l
     1grp,lsid,ierr)
         endif
         if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
            call MPI_IRECV(scb(1,1,1,1,m),3*kzp*nxvy,mreal,mm-nvpy,noff+
     1ii,lgrp,msid,ierr)
         endif
         if (kzp.gt.1) then
            call MPI_IRECV(cu3(1,1,1,loff+2,m),3*(kzp-1)*nxvy,mreal,mm-2
     1*nvpy,noff+ii,lgrp,nsid,ierr)
         endif
      endif
      if ((jl.gt.((kzb2-1)/2+1)).and.(jl.le.kzb2).and.lt1) then
         if (kz.eq.0) then
            if ((jl+1).le.kzb2) then
               call MPI_SEND(cu(1,1,1,1,m),3*nxvy,mreal,lm,noff+ii,lgrp,
     1ierr)
            endif
            if (kzp.gt.1) then
               call MPI_SEND(cu(1,1,1,2,m),3*(kzp-1)*nxvy,mreal,lm-nvpy,
     1noff+ii,lgrp,ierr)
            endif
         else
            call MPI_SEND(cu(1,1,1,1,m),3*kzp*nxvy,mreal,lm-nvpy,noff+ii
     1,lgrp,ierr)
         endif
      endif
c wait for data and unpack it
      if ((jm.le.kzb).and.(km.le.kyb)) then
         if (kzp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 600 l = 2, kzp
            l1 = kzp - l + 2
            l2 = (l1 + 2)/4 + 1
            koff = kypd*(l1 - 4*(l2 - 1) + 2)
            do 590 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1 + koff
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 580 j = 1, nxv
            do 570 i = 1, 3
            cu3(i,j,k1,l1+loff,m) = cu3(i,j+joff,k2,l2+loff,m)
  570       continue
  580       continue
  590       continue
  600       continue
         endif
         if ((jm+1).le.kzb) then
            call MPI_WAIT(lsid,istatus,ierr)
            do 630 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 620 j = 1, nxv
            do 610 i = 1, 3
            cu3(i,j,k1,loff+1,m) = cu3(i,j+joff,k2,loff+1,m)
  610       continue
  620       continue
  630       continue
         endif
         if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 670 l = 1, kzp
            do 660 k = 1, kypd
            do 650 j = 1, nxv
            do 640 i = 1, 3
            cu3(i,j,k+kyp,l+loff,m) = scb(i,j,k,l,m)
  640       continue
  650       continue
  660       continue
  670       continue
         endif
      endif
  680 continue
c switch internal data
      if ((jm.le.kzb).and.(km.le.kyb)) then
         do 720 l = 1, kzp
         do 710 k = 1, kyp
         do 700 j = 1, nxv
         do 690 i = 1, 3
         at1 = cu3(i,j,k+kyp,l,m)
         cu3(i,j,k+kyp,l,m) = cu3(i,j,k,l+kzp,m)
         cu3(i,j,k,l+kzp,m) = at1
  690    continue
  700    continue
  710    continue
  720    continue
      endif
  730 continue
  740 continue
c finish copy reflected data in y and z direction
      do 850 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 840 my = 1, k2blok
      m = my + moff
      loff = kyp2*(my + js)
      jm = (ny2 - loff - 1)/kyp + 1
      loff = kzp2*(mz + ks)
      km = (nz2 - loff - 1)/kzp + 1
      loff = kyp*(my + js)
      jl = (ny2 - loff - 1)/kyp2 + 1
      ky = loff + kyp2*jl - ny2
      loff = kzp*(mz + ks)
      kl = loff - kzp2*(loff/kzp2)
      kz = nvpy*(nvpz - (mz + ks) - 1)
      do 810 ii = 9, 10
c for odd rows in z
      if (ii.eq.9) then
         mm = jm + 2*kz
         lm = jl + (nvpy*(nvpz - 1) + kz)/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in z
      else if (ii.eq.10) then
         mm = mm + nvpy
         lm = lm + nvpy/2
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kyb).and.(km.le.kzb)) then
c        if ((jm+1).le.kyb) then
c           do 770 l = 1, kzp
c           do 760 j = 1, nx
c           do 750 i = 1, 3
c           cu3(i,j,1,l+loff,m) = cu(i,j,1,l,mm+1)
c 750       continue
c 760       continue
c 770       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if ((jm+1).le.kyb) then
            call MPI_IRECV(scd(1,1,1,2,m),3*kzp*nxv,mreal,mm,noff+ii,lgr
     1p,lsid,ierr)
         endif
      endif
      if ((jl.gt.((kyb2-1)/2+1)).and.(jl.le.kyb2).and.lt1) then
         if (ky.eq.0) then
            if ((jl+1).le.kyb2) then
               do 770 l = 1, kzp
               do 760 j = 1, nxv
               do 750 i = 1, 3
               scd(i,j,l,1,m) = cu(i,j,1,l,m)
  750          continue
  760          continue
  770          continue
               call MPI_SEND(scd(1,1,1,1,m),3*kzp*nxv,mreal,lm,noff+ii,l
     1grp,ierr)
            endif
         endif
      endif
c wait for data and unpack it
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if ((jm+1).le.kyb) then
            call MPI_WAIT(lsid,istatus,ierr)
            do 800 l = 1, kzp
            do 790 j = 1, nxv
            do 780 i = 1, 3
            cu3(i,j,1,l+loff,m) = scd(i,j,l,2,m)
  780       continue
  790       continue
  800       continue
         endif
      endif
  810 continue
c last point
      lm = jl + nvpy*(nvpz + mz + ks)/2
      mm = jm + nvpy*(2*(mz + ks) - nvpz)
      lt1 = (jl+1).le.kyb2
      loff = kzp*(mz + ks)
      jl = (nz2 - loff - 1)/kzp2 + 1
      kz = loff + kzp2*jl - nz2
      loff = kyp*(my + js)
      kl = loff - kyp2*(loff/kyp2)
c this segment is used for shared memory computers
c     if ((jm.le.kyb).and.(km.le.kzb)) then
c        if (((jm+1).le.kyb).and.((km+1).le.kzb)) then
c           do 830 j = 1, nx
c           do 820 i = 1, 3
c           cu3(i,j,1,1,m) = cu(i,j,1,1,mm+1)
c 820       continue
c 830       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if (((jm+1).le.kyb).and.((km+1).le.kzb)) then
            call MPI_IRECV(cu3(1,1,1,1,m),3*nxv,mreal,mm,noff+11,lgrp,ls
     1id,ierr)
         endif
      endif
      if ((jl.gt.((kzb2-1)/2+1)).and.(jl.le.kzb2).and.(kl.eq.0)) then
         if (kz.eq.0) then
            if (((jl+1).le.kzb2).and.lt1) then
               call MPI_SEND(cu(1,1,1,1,m),3*nxv,mreal,lm,noff+11,lgrp,i
     1err)
            endif
         endif
      endif
c wait for data
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if (((jm+1).le.kyb).and.((km+1).le.kzb)) then
            call MPI_WAIT(lsid,istatus,ierr)
         endif
      endif
  840 continue
  850 continue
c create odd array
  860 do 1130 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      loff = kzp2*(mz + ks)
      do 1120 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      do 1110 l = 1, kzp2
      ll = l + loff
      if ((ll.eq.1).or.(ll.eq.(nz+1))) then
         do 890 k = 1, kyp2
         do 880 j = 1, nx
         do 870 i = 1, 3
         cu3(i,j,k,l,m) = 0.
         cu3(i,j+nx,k,l,m) = 0.
  870    continue
  880    continue
  890    continue
      else if (ll.le.nz) then
         do 950 k = 1, kyp2
         kk = k + koff
         if ((kk.eq.1).or.(kk.eq.(ny+1))) then
            do 910 j = 1, nx
            do 900 i = 1, 3
            cu3(i,j,k,l,m) = 0.
            cu3(i,j+nx,k,l,m) = 0.
  900       continue
  910       continue
         else if (kk.le.ny) then
            do 920 j = 1, nxs
            cu3(1,nx+j+1,k,l,m) = cu3(1,nx-j+1,k,l,m)
            cu3(2,nx+j+1,k,l,m) = -cu3(2,nx-j+1,k,l,m)
            cu3(3,nx+j+1,k,l,m) = -cu3(3,nx-j+1,k,l,m)
  920       continue
            cu3(1,1,k,l,m) = 0.
            cu3(2,1,k,l,m) = 0.
            cu3(3,1,k,l,m) = 0.
            cu3(1,nx+1,k,l,m) = 0.
            cu3(2,nx+1,k,l,m) = 0.
            cu3(3,nx+1,k,l,m) = 0.
         else if (kk.gt.(ny+1)) then
            if (k.eq.1) then
               do 930 j = 1, nxs
               cu3(1,nx+j+1,k,l,m) = -cu3(1,nx-j+1,k,l,m)
               cu3(2,nx+j+1,k,l,m) = -cu3(2,nx-j+1,k,l,m)
               cu3(3,nx+j+1,k,l,m) = cu3(3,nx-j+1,k,l,m)
  930          continue
            else
               do 940 j = 1, nxs
               cu3(1,nx+j+1,k,l,m) = -cu3(1,nx-j+1,kyp2-k+2,l,m)
               cu3(2,nx+j+1,k,l,m) = -cu3(2,nx-j+1,kyp2-k+2,l,m)
               cu3(3,nx+j+1,k,l,m) = cu3(3,nx-j+1,kyp2-k+2,l,m)
  940          continue
            endif
            cu3(1,1,k,l,m) = 0.
            cu3(2,1,k,l,m) = 0.
            cu3(3,1,k,l,m) = 0.
            cu3(1,nx+1,k,l,m) = 0.
            cu3(2,nx+1,k,l,m) = 0.
            cu3(3,nx+1,k,l,m) = 0.
         endif
  950    continue
      else if (ll.gt.(nz+1)) then
         if (l.eq.1) then
            do 1010 k = 1, kyp2
            kk = k + koff
            if (kk.le.ny) then
               do 960 j = 1, nxs
               cu3(1,nx+j+1,k,l,m) = -cu3(1,nx-j+1,k,l,m)
               cu3(2,nx+j+1,k,l,m) = cu3(2,nx-j+1,k,l,m)
               cu3(3,nx+j+1,k,l,m) = -cu3(3,nx-j+1,k,l,m)
  960          continue
            else if (kk.gt.(ny+1)) then
               if (k.eq.1) then
                  do 980 j = 1, nxs
                  do 970 i = 1, 3
                  cu3(i,nx+j+1,k,l,m) = cu3(i,nx-j+1,k,l,m)
  970             continue
  980             continue
               else
                  do 1000 j = 1, nxs
                  do 990 i = 1, 3
                  cu3(i,nx+j+1,k,l,m) = cu3(i,nx-j+1,kyp2-k+2,l,m)
  990             continue
 1000             continue
               endif
            endif
 1010       continue
         else
            do 1070 k = 1, kyp2
            kk = k + koff
            if (kk.le.ny) then
               do 1020 j = 1, nxs
               cu3(1,nx+j+1,k,l,m) = -cu3(1,nx-j+1,k,kzp2-l+2,m)
               cu3(2,nx+j+1,k,l,m) = cu3(2,nx-j+1,k,kzp2-l+2,m)
               cu3(3,nx+j+1,k,l,m) = -cu3(3,nx-j+1,k,kzp2-l+2,m)
 1020          continue
            else if (kk.gt.(ny+1)) then
               if (k.eq.1) then
                  do 1040 j = 1, nxs
                  do 1030 i = 1, 3
                  cu3(i,nx+j+1,k,l,m) = cu3(i,nx-j+1,k,kzp2-l+2,m)
 1030             continue
 1040             continue
               else
                  do 1060 j = 1, nxs
                  do 1050 i = 1, 3
                  cu3(i,nx+j+1,k,l,m) = cu3(i,nx-j+1,kyp2-k+2,kzp2-l+2,m
     1)
 1050             continue
 1060             continue
               endif
            endif
 1070       continue
         endif
         do 1100 k = 1, kyp2
         kk = k + koff
         if ((kk.eq.1).or.(kk.eq.(ny+1))) then
            do 1090 j = 1, nx
            do 1080 i = 1, 3
            cu3(i,j,k,l,m) = 0.
            cu3(i,j+nx,k,l,m) = 0.
 1080       continue
 1090       continue
         else
            cu3(1,1,k,l,m) = 0.
            cu3(2,1,k,l,m) = 0.
            cu3(3,1,k,l,m) = 0.
            cu3(1,nx+1,k,l,m) = 0.
            cu3(2,nx+1,k,l,m) = 0.
            cu3(3,nx+1,k,l,m) = 0.
         endif
 1100    continue
      endif
 1110 continue
 1120 continue
 1130 continue
c finish odd array
      do 1210 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      loff = kzp2*(mz + ks)
      do 1200 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      do 1190 l = 1, kzp2
      ll = l + loff
      if (ll.le.nz) then
         do 1160 k = 1, kyp2
         kk = k + koff
         if (kk.gt.(ny+1)) then
            if (k.eq.1) then
               do 1140 j = 1, nxs
               cu3(1,nx-j+1,k,l,m) = cu3(1,nx+j+1,k,l,m)
               cu3(2,nx-j+1,k,l,m) = -cu3(2,nx+j+1,k,l,m)
               cu3(3,nx-j+1,k,l,m) = -cu3(3,nx+j+1,k,l,m)
 1140          continue
               cu3(1,nx+1,k,l,m) = cu3(1,nx+1,k,l,m)
               cu3(2,nx+1,k,l,m) = -cu3(2,nx+1,k,l,m)
               cu3(3,nx+1,k,l,m) = -cu3(3,nx+1,k,l,m)
            else
               do 1150 j = 1, nxs
               cu3(1,nx-j+1,k,l,m) = cu3(1,nx+j+1,k,l,m)
               cu3(2,nx-j+1,k,l,m) = -cu3(2,nx+j+1,k,l,m)
               cu3(3,nx-j+1,k,l,m) = -cu3(3,nx+j+1,k,l,m)
 1150          continue
               cu3(1,nx+1,k,l,m) = cu3(1,nx+1,k,l,m)
               cu3(2,nx+1,k,l,m) = -cu3(2,nx+1,k,l,m)
               cu3(3,nx+1,k,l,m) = -cu3(3,nx+1,k,l,m)
            endif
         endif
 1160    continue
      else if (ll.gt.(nz+1)) then
         do 1180 k = 1, kyp2
         kk = k + koff
         if ((kk.le.ny).or.(kk.gt.(ny+1))) then
            do 1170 j = 1, nxs
            cu3(1,nx-j+1,k,l,m) = cu3(1,nx+j+1,k,l,m)
            cu3(2,nx-j+1,k,l,m) = -cu3(2,nx+j+1,k,l,m)
            cu3(3,nx-j+1,k,l,m) = -cu3(3,nx+j+1,k,l,m)
 1170       continue
            cu3(1,nx+1,k,l,m) = cu3(1,nx+1,k,l,m)
            cu3(2,nx+1,k,l,m) = -cu3(2,nx+1,k,l,m)
            cu3(3,nx+1,k,l,m) = -cu3(3,nx+1,k,l,m)
         endif
 1180    continue
      endif
 1190 continue
 1200 continue
 1210 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PTRPSIN32D(q,q3,scb,scd,nx,ny,nz,kstrt,nvpy,nvpz,nxv,ky
     1p,kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
c this subroutine creates an odd array q3 from an array q, so that
c a 3d sine transform can be performed with a 3d real to complex fft.
c linear interpolation for distributed data with 2D domain decomposition
c q3 array may be modified
c scb/scd = scratch arrays
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = first dimension of input array q, must be >= nx
c kyp/kzp = number of data values per block in q in y/z
c kypd = second dimension of input array q, must be >= kyp
c kzpd = third dimension of input array q, must be >= kzp
c kyp2/kzp2 = number of data values per block in q3 in y/z
c kblok/lblok = number of data blocks in y/z
c k2blok/l2blok = number of data blocks in y/z for tripled data
      implicit none
      real q, q3, scb, scd
      integer nx, ny, nz, kstrt, nvpy, nvpz, nxv, kyp, kzp, kypd, kzpd
      integer kyp2, kzp2, kblok, lblok, k2blok, l2blok
      dimension q(nxv,kypd,kzpd,kblok*lblok)
      dimension q3(2*nxv,2*kypd,kzp2,k2blok*l2blok)
      dimension scb(nxv,kypd,kzpd,kblok*lblok)
      dimension scd(nxv,kzpd,2,kblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      logical lt1
      integer istatus, lsid, msid, nsid, ierr
      integer i, j, k, l, m, my, mz, nxs, nys, nzs, ks, js, ny2, nz2, kr
      integer kyb, kyb2, kzb, kzb2, nxvy, noff, moff, loff, koff, joff
      integer jm, km, jl, kl, ky, kz, kk, ll, mm, lm, l1, l2, k0, k1, k2
      real at1
      dimension istatus(lstat)
      nxs = nx - 1
      nys = ny - 1
      nzs = nz - 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      kyb = ny/kyp
      ny2 = ny + ny
      kyb2 = ny2/kyp2
      kzb = nz/kzp
      nz2 = nz + nz
      kzb2 = nz2/kzp2
      nxvy = nxv*kypd
      noff = kypd*kzpd + kyb*kzb
c copy to triple array in y direction
      do 150 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 140 my = 1, k2blok
      m = my + moff
      loff = kyp2*(my + js)
      jm = loff/kyp + 1
      loff = kzp2*(mz + ks)
      km = loff/kzp + 1
      loff = kyp*(my + js)
      jl = loff/kyp2 + 1
      loff = kzp*(mz + ks)
      kl = loff - kzp2*(loff/kzp2)
      kz = nvpy*(mz + ks)
c special case for one processor
      if ((kyb2*kzb2).eq.1) then
         do 40 l = 1, nzs
         do 20 k = 1, nys
         do 10 j = 1, nxs
         q3(j+1,k+1,l+1,m) = q(j+1,k+1,l+1,m)
         q3(nx+j+1,k+1,l+1,m) = -q(nx-j+1,k+1,l+1,m)
         q3(j+1,ny+k+1,l+1,m) = -q(j+1,ny-k+1,l+1,m)
         q3(nx+j+1,ny+k+1,l+1,m) = q(nx-j+1,ny-k+1,l+1,m)
         q3(j+1,k+1,nz+l+1,m) = -q(j+1,k+1,nz-l+1,m)
         q3(nx+j+1,k+1,nz+l+1,m) = q(nx-j+1,k+1,nz-l+1,m)
         q3(j+1,ny+k+1,nz+l+1,m) = q(j+1,ny-k+1,nz-l+1,m)
         q3(nx+j+1,ny+k+1,nz+l+1,m) = -q(nx-j+1,ny-k+1,nz-l+1,m)
   10    continue
         q3(1,k+1,l+1,m) = 0.
         q3(nx+1,k+1,l+1,m) = 0.
         q3(1,k+ny+1,l+1,m) = 0.
         q3(nx+1,k+ny+1,l+1,m) = 0.
         q3(1,k+1,nz+l+1,m) = 0.
         q3(nx+1,k+1,nz+l+1,m) = 0.
         q3(1,k+ny+1,nz+l+1,m) = 0.
         q3(nx+1,k+ny+1,nz+l+1,m) = 0.
   20    continue
         do 30 j = 1, nx
         q3(j,1,l+1,m) = 0.
         q3(j+nx,1,l+1,m) = 0.
         q3(j,ny+1,l+1,m) = 0.
         q3(j+nx,ny+1,l+1,m) = 0.
         q3(j,1,nz+l+1,m) = 0.
         q3(j+nx,1,nz+l+1,m) = 0.
         q3(j,ny+1,nz+l+1,m) = 0.
         q3(j+nx,ny+1,nz+l+1,m) = 0.
   30    continue
   40    continue
         nxs = nx + nx
         nys = ny + ny
         do 60 k = 1, nys
         do 50 j = 1, nxs
         q3(j,k,1,m) = 0.
         q3(j,k,nz+1,m) = 0.
   50    continue
   60    continue
         return
      endif
c copy main data in y direction
      do 130 i = 1, 2
c for odd rows in z
      if (i.eq.1) then
         mm = jm + 2*kz
         lm = jl + kz/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in z
      else if (i.eq.2) then
         if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 130
         mm = mm + nvpy
         lm = lm - nvpy/2
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kyb).and.(km.le.kzb)) then
c        do 90 l = 1, kzp
c        do 80 k = 1, kyp
c        do 70 j = 1, nx
c        q3(j,k,l+loff,m) = q(j,k,l,mm)
c  70    continue
c  80    continue
c  90    continue
c        if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
c           do 120 l = 1, kzp
c           do 110 k = 1, kyp
c           do 100 j = 1, nx
c           q3(j,k+kyp,l+loff,m) = q(j,k,l,mm+1)
c 100       continue
c 110       continue
c 120       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kyb).and.(km.le.kzb)) then
         call MPI_IRECV(q3(1,1,loff+1,m),kzp*nxvy,mreal,mm-1,noff+i,lgrp
     1,msid,ierr)
         if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
            call MPI_IRECV(scb(1,1,1,m),kzp*nxvy,mreal,mm,noff+i,lgrp,ns
     1id,ierr)
         endif
      endif
      if ((jl.le.((kyb2-1)/2+1)).and.lt1) then
         call MPI_SEND(q(1,1,1,m),kzp*nxvy,mreal,lm-1,noff+i,lgrp,ierr)
      endif
c wait for data and unpack it
      if ((jm.le.kyb).and.(km.le.kzb)) then
         call MPI_WAIT(msid,istatus,ierr)
         do 90 l = 1, kzp
         l1 = kzp - l + 1
         l2 = (l1 - 1)/4 + 1
         koff = kypd*(l1 - 4*(l2 - 1) - 1)
         do 80 k = 1, kypd
         k1 = kypd - k + 1
         k0 = k1 - 1 + koff
         k2 = k0/2 + 1
         joff = nxv*(k0 - 2*(k2 - 1))
         do 70 j = 1, nxv
         q3(j,k1,l1+loff,m) = q3(j+joff,k2,l2+loff,m)
   70    continue
   80    continue
   90    continue
         if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 120 l = 1, kzp
            do 110 k = 1, kypd
            do 100 j = 1, nxv
            q3(j,k+kyp,l+loff,m) = scb(j,k,l,m)
  100       continue
  110       continue
  120       continue
         endif
      endif
  130 continue
  140 continue
  150 continue
c copying in y direction not needed
      if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 290
c copy reflected data in y direction
      do 280 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 270 my = 1, k2blok
      m = my + moff
      loff = kyp2*(my + js)
      jm = (ny2 - loff - 1)/kyp + 1
      loff = kzp2*(mz + ks)
      km = loff/kzp + 1
      loff = kyp*(my + js)
      jl = (ny2 - loff - 1)/kyp2 + 1
      ky = loff + kyp2*jl - ny2
      loff = kzp*(mz + ks)
      kl = loff - kzp2*(loff/kzp2)
      kz = nvpy*(mz + ks)
      do 260 i = 3, 4
c for odd rows in z
      if (i.eq.3) then
         mm = jm + 2*kz
         lm = jl + kz/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in z
      else if (i.eq.4) then
         if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 260
         mm = mm + nvpy
         lm = lm - nvpy/2
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kyb).and.(km.le.kzb)) then
c        if ((jm+1).le.kyb) then
c           do 170 l = 1, kzp
c           do 160 j = 1, nx
c           q3(j,1,l+loff,m) = q(j,1,l,mm+1)
c 160       continue
c 170       continue
c        endif
c        if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
c           do 200 l = 1, kzp
c           do 190 k = 1, kyp
c           do 180 j = 1, nx
c           q3(j,k+kyp,l+loff,m) = q(j,k,l,mm)
c 180       continue
c 190       continue
c 200       continue
c        endif
c        if (kyp.gt.1) then
c           do 230 l = 1, kzp
c           do 220 k = 2, kyp
c           do 210 j = 1, nx
c           q3(j,k,l+loff,m) = q(j,k,l,mm-1)
c 210       continue
c 220       continue
c 230       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if ((jm+1).le.kyb) then
            call MPI_IRECV(scd(1,1,2,m),kzp*nxv,mreal,mm,noff+i,lgrp,lsi
     1d,ierr)
         endif
         if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
            call MPI_IRECV(scb(1,1,1,m),kzp*nxvy,mreal,mm-1,noff+i,lgrp,
     1msid,ierr)
         endif
         if (kyp.gt.1) then
            call MPI_IRECV(q3(1,1,loff+1,m),kzp*nxvy,mreal,mm-2,noff+i,l
     1grp,nsid,ierr)
         endif
      endif
      if ((jl.gt.((kyb2-1)/2+1)).and.(jl.le.kyb2).and.lt1) then
         if (ky.eq.0) then
            if ((jl+1).le.kyb2) then
               do 170 l = 1, kzp
               do 160 j = 1, nxv
               scd(j,l,1,m) = q(j,1,l,m)
  160          continue
  170          continue
               call MPI_SEND(scd(1,1,1,m),kzp*nxv,mreal,lm,noff+i,lgrp,i
     1err)
            endif
            if (kyp.gt.1) then
               call MPI_SEND(q(1,1,1,m),kzp*nxvy,mreal,lm-1,noff+i,lgrp,
     1ierr)
            endif
         else
            call MPI_SEND(q(1,1,1,m),kzp*nxvy,mreal,lm-1,noff+i,lgrp,ier
     1r)
         endif
      endif
c wait for data and unpack it
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if (kyp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 200 l = 1, kzp
            l1 = kzp - l + 1
            l2 = (l1 - 1)/4 + 1
            koff = kypd*(l1 - 4*(l2 - 1) - 1)
            do 190 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1 + koff
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 180 j = 1, nxv
            q3(j,k1,l1+loff,m) = q3(j+joff,k2,l2+loff,m)
  180       continue
  190       continue
  200       continue
         endif
         if ((jm+1).le.kyb) then
            call MPI_WAIT(lsid,istatus,ierr)
            do 220 l = 1, kzp
            do 210 j = 1, nxv
            q3(j,1,l+loff,m) = scd(j,l,2,m)
  210       continue
  220       continue
         endif
         if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 250 l = 1, kzp
            do 240 k = 1, kypd
            do 230 j = 1, nxv
            q3(j,k+kyp,l+loff,m) = scb(j,k,l,m)
  230       continue
  240       continue
  250       continue
         endif
      endif
  260 continue
  270 continue
  280 continue
c copying in z direction not needed
  290 if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 660
c copy reflected data in z direction
      do 430 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 420 my = 1, k2blok
      m = my + moff
      loff = kzp2*(mz + ks)
      jm = (nz2 - loff - 1)/kzp + 1
      loff = kyp2*(my + js)
      km = loff/kyp + 1
      loff = kzp*(mz + ks)
      jl = (nz2 - loff - 1)/kzp2 + 1
      kz = loff + kzp2*jl - nz2
      loff = kyp*(my + js)
      kl = loff - kyp2*(loff/kyp2)
      kr = my + js
      do 380 i = 5, 6
c for odd rows in y
      if (i.eq.5) then
         mm = nvpy*jm + 2*kr
         lm = nvpy*jl + kr/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in y
      else if (i.eq.6) then
         if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 380
         mm = mm + 1
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kzb).and.(km.le.kyb)) then
c        if (kzp.gt.1) then
c           do 320 l = 2, kzp
c           do 310 k = 1, kypd
c           do 300 j = 1, nx
c           q3(j,k,l+loff,m) = q(j,k,l,mm-2*nvpy+1)
c 300       continue
c 310       continue
c 320       continue
c        endif
c        if ((jm+1).le.kzb) then
c           do 340 k = 1, kypd
c           do 330 j = 1, nx
c           q3(j,k,loff+1,m) = q(j,k,1,mm+1)
c 330       continue
c 340       continue
c        endif
c        if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
c           do 370 l = 1, kzp
c           do 360 k = 1, kypd
c           do 350 j = 1, nx
c           q3(j,k+kyp,l+loff,m) = q(j,k,l,mm-nvpy+1)
c 350       continue
c 360       continue
c 370       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kzb).and.(km.le.kyb)) then
         if ((jm+1).le.kzb) then
            call MPI_IRECV(q3(1,1,loff+1,m),nxvy,mreal,mm,noff+i,lgrp,ls
     1id,ierr)
         endif
         if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
            call MPI_IRECV(scb(1,1,1,m),kzp*nxvy,mreal,mm-nvpy,noff+i,lg
     1rp,msid,ierr)
         endif
         if (kzp.gt.1) then
            call MPI_IRECV(q3(1,1,loff+2,m),(kzp-1)*nxvy,mreal,mm-2*nvpy
     1,noff+i,lgrp,nsid,ierr)
         endif
      endif
      if ((jl.gt.((kzb2-1)/2+1)).and.(jl.le.kzb2).and.lt1) then
         if (kz.eq.0) then
            if ((jl+1).le.kzb2) then
               call MPI_SEND(q(1,1,1,m),nxvy,mreal,lm,noff+i,lgrp,ierr)
            endif
            if (kzp.gt.1) then
               call MPI_SEND(q(1,1,2,m),(kzp-1)*nxvy,mreal,lm-nvpy,noff+
     1i,lgrp,ierr)
            endif
         else
            call MPI_SEND(q(1,1,1,m),kzp*nxvy,mreal,lm-nvpy,noff+i,lgrp,
     1ierr)
         endif
      endif
c wait for data and unpack it
      if ((jm.le.kzb).and.(km.le.kyb)) then
         if (kzp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 320 l = 2, kzp
            l1 = kzp - l + 2
            l2 = (l1 + 2)/4 + 1
            koff = kypd*(l1 - 4*(l2 - 1) + 2)
            do 310 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1 + koff
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 300 j = 1, nxv
            q3(j,k1,l1+loff,m) = q3(j+joff,k2,l2+loff,m)
  300       continue
  310       continue
  320       continue
         endif
         if ((jm+1).le.kzb) then
            call MPI_WAIT(lsid,istatus,ierr)
            do 340 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 330 j = 1, nxv
            q3(j,k1,loff+1,m) = q3(j+joff,k2,loff+1,m)
  330       continue
  340       continue
         endif
         if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 370 l = 1, kzp
            do 360 k = 1, kypd
            do 350 j = 1, nxv
            q3(j,k+kyp,l+loff,m) = scb(j,k,l,m)
  350       continue
  360       continue
  370       continue
         endif
      endif
  380 continue
c switch internal data
      if ((jm.le.kzb).and.(km.le.kyb)) then
         do 410 l = 1, kzp
         do 400 k = 1, kyp
         do 390 j = 1, nxv
         at1 = q3(j,k+kyp,l,m)
         q3(j,k+kyp,l,m) = q3(j,k,l+kzp,m)
         q3(j,k,l+kzp,m) = at1
  390    continue
  400    continue
  410    continue
      endif
  420 continue
  430 continue
c copying in y direction not needed
      if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 660
c copy reflected data in y and z direction
      do 570 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 560 my = 1, k2blok
      m = my + moff
      loff = kzp2*(mz + ks)
      jm = (nz2 - loff - 1)/kzp + 1
      loff = kyp2*(my + js)
      km = (ny2 - loff - 1)/kyp + 1
      loff = kzp*(mz + ks)
      jl = (nz2 - loff - 1)/kzp2 + 1
      kz = loff + kzp2*jl - nz2
      loff = kyp*(my + js)
      kl = loff - kyp2*(loff/kyp2)
      kr = nvpy - (my + js) - 1
      do 520 i = 7, 8
c for odd rows in y
      if (i.eq.7) then
         mm = nvpy*jm + 2*kr
         lm = nvpy*jl + (nvpy + kr)/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in y
      else if (i.eq.8) then
         if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 520
         mm = mm + 1
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kzb).and.(km.le.kyb)) then
c        if (kzp.gt.1) then
c           do 460 l = 2, kzp
c           do 450 k = 1, kypd
c           do 440 j = 1, nx
c           q3(j,k,l+loff,m) = q(j,k,l,mm-2*nvpy+1)
c 440       continue
c 450       continue
c 460       continue
c        endif
c        if ((jm+1).le.kzb) then
c           do 480 k = 1, kypd
c           do 470 j = 1, nx
c           q3(j,k,loff+1,m) = q(j,k,1,mm+1)
c 470       continue
c 480       continue
c        endif
c        if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
c           do 510 l = 1, kzp
c           do 500 k = 1, kypd
c           do 490 j = 1, nx
c           q3(j,k+kyp,l+loff,m) = q(j,k,l,mm-nvpy+1)
c 490       continue
c 500       continue
c 510       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kzb).and.(km.le.kyb)) then
         if ((jm+1).le.kzb) then
            call MPI_IRECV(q3(1,1,loff+1,m),nxvy,mreal,mm,noff+i,lgrp,ls
     1id,ierr)
         endif
         if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
            call MPI_IRECV(scb(1,1,1,m),kzp*nxvy,mreal,mm-nvpy,noff+i,lg
     1rp,msid,ierr)
         endif
         if (kzp.gt.1) then
            call MPI_IRECV(q3(1,1,loff+2,m),(kzp-1)*nxvy,mreal,mm-2*nvpy
     1,noff+i,lgrp,nsid,ierr)
         endif
      endif
      if ((jl.gt.((kzb2-1)/2+1)).and.(jl.le.kzb2).and.lt1) then
         if (kz.eq.0) then
            if ((jl+1).le.kzb2) then
               call MPI_SEND(q(1,1,1,m),nxvy,mreal,lm,noff+i,lgrp,ierr)
            endif
            if (kzp.gt.1) then
               call MPI_SEND(q(1,1,2,m),(kzp-1)*nxvy,mreal,lm-nvpy,noff+
     1i,lgrp,ierr)
            endif
         else
            call MPI_SEND(q(1,1,1,m),kzp*nxvy,mreal,lm-nvpy,noff+i,lgrp,
     1ierr)
         endif
      endif
c wait for data and unpack it
      if ((jm.le.kzb).and.(km.le.kyb)) then
         if (kzp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 460 l = 2, kzp
            l1 = kzp - l + 2
            l2 = (l1 + 2)/4 + 1
            koff = kypd*(l1 - 4*(l2 - 1) + 2)
            do 450 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1 + koff
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 440 j = 1, nxv
            q3(j,k1,l1+loff,m) = q3(j+joff,k2,l2+loff,m)
  440       continue
  450       continue
  460       continue
         endif
         if ((jm+1).le.kzb) then
            call MPI_WAIT(lsid,istatus,ierr)
            do 480 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 470 j = 1, nxv
            q3(j,k1,loff+1,m) = q3(j+joff,k2,loff+1,m)
  470       continue
  480       continue
         endif
         if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 510 l = 1, kzp
            do 500 k = 1, kypd
            do 490 j = 1, nxv
            q3(j,k+kyp,l+loff,m) = scb(j,k,l,m)
  490       continue
  500       continue
  510       continue
         endif
      endif
  520 continue
c switch internal data
      if ((jm.le.kzb).and.(km.le.kyb)) then
         do 550 l = 1, kzp
         do 540 k = 1, kyp
         do 530 j = 1, nxv
         at1 = q3(j,k+kyp,l,m)
         q3(j,k+kyp,l,m) = q3(j,k,l+kzp,m)
         q3(j,k,l+kzp,m) = at1
  530    continue
  540    continue
  550    continue
      endif
  560 continue
  570 continue
c finish copy reflected data in y and z direction
      do 650 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 640 my = 1, k2blok
      m = my + moff
      loff = kyp2*(my + js)
      jm = (ny2 - loff - 1)/kyp + 1
      loff = kzp2*(mz + ks)
      km = (nz2 - loff - 1)/kzp + 1
      loff = kyp*(my + js)
      jl = (ny2 - loff - 1)/kyp2 + 1
      ky = loff + kyp2*jl - ny2
      loff = kzp*(mz + ks)
      kl = loff - kzp2*(loff/kzp2)
      kz = nvpy*(nvpz - (mz + ks) - 1)
      do 620 i = 9, 10
c for odd rows in z
      if (i.eq.9) then
         mm = jm + 2*kz
         lm = jl + (nvpy*(nvpz - 1) + kz)/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in z
      else if (i.eq.10) then
         mm = mm + nvpy
         lm = lm + nvpy/2
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kyb).and.(km.le.kzb)) then
c        if ((jm+1).le.kyb) then
c           do 590 l = 1, kzp
c           do 580 j = 1, nx
c           q3(j,1,l+loff,m) = q(j,1,l,mm+1)
c 580       continue
c 590       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if ((jm+1).le.kyb) then
            call MPI_IRECV(scd(1,1,2,m),kzp*nxv,mreal,mm,noff+i,lgrp,lsi
     1d,ierr)
         endif
      endif
      if ((jl.gt.((kyb2-1)/2+1)).and.(jl.le.kyb2).and.lt1) then
         if (ky.eq.0) then
            if ((jl+1).le.kyb2) then
               do 590 l = 1, kzp
               do 580 j = 1, nxv
               scd(j,l,1,m) = q(j,1,l,m)
  580          continue
  590          continue
               call MPI_SEND(scd(1,1,1,m),kzp*nxv,mreal,lm,noff+i,lgrp,i
     1err)
            endif
         endif
      endif
c wait for data and unpack it
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if ((jm+1).le.kyb) then
            call MPI_WAIT(lsid,istatus,ierr)
            do 610 l = 1, kzp
            do 600 j = 1, nxv
            q3(j,1,l+loff,m) = scd(j,l,2,m)
  600       continue
  610       continue
         endif
      endif
  620 continue
c last point
      lm = jl + nvpy*(nvpz + mz + ks)/2
      mm = jm + nvpy*(2*(mz + ks) - nvpz)
      lt1 = (jl+1).le.kyb2
      loff = kzp*(mz + ks)
      jl = (nz2 - loff - 1)/kzp2 + 1
      kz = loff + kzp2*jl - nz2
      loff = kyp*(my + js)
      kl = loff - kyp2*(loff/kyp2)
c this segment is used for shared memory computers
c     if ((jm.le.kyb).and.(km.le.kzb)) then
c        if (((jm+1).le.kyb).and.((km+1).le.kzb)) then
c           do 630 j = 1, nx
c           q3(j,1,1,m) = q(j,1,1,mm+1)
c 630       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if (((jm+1).le.kyb).and.((km+1).le.kzb)) then
            call MPI_IRECV(q3(1,1,1,m),nxv,mreal,mm,noff+11,lgrp,lsid,ie
     1rr)
         endif
      endif
      if ((jl.gt.((kzb2-1)/2+1)).and.(jl.le.kzb2).and.(kl.eq.0)) then
         if (kz.eq.0) then
            if (((jl+1).le.kzb2).and.lt1) then
               call MPI_SEND(q(1,1,1,m),nxv,mreal,lm,noff+11,lgrp,ierr)
            endif
         endif
      endif
c wait for data
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if (((jm+1).le.kyb).and.((km+1).le.kzb)) then
            call MPI_WAIT(lsid,istatus,ierr)
         endif
      endif
  640 continue
  650 continue
c create odd array
  660 do 860 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      loff = kzp2*(mz + ks)
      do 850 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      do 840 l = 1, kzp2
      ll = l + loff
      if ((ll.eq.1).or.(ll.eq.(nz+1))) then
         do 680 k = 1, kyp2
         do 670 j = 1, nx
         q3(j,k,l,m) = 0.
         q3(j+nx,k,l,m) = 0.
  670    continue
  680    continue
      else if (ll.le.nz) then
         do 730 k = 1, kyp2
         kk = k + koff
         if ((kk.eq.1).or.(kk.eq.(ny+1))) then
            do 690 j = 1, nx
            q3(j,k,l,m) = 0.
            q3(j+nx,k,l,m) = 0.
  690       continue
         else if (kk.le.ny) then
            do 700 j = 1, nxs
            q3(nx+j+1,k,l,m) = -q3(nx-j+1,k,l,m)
  700       continue
            q3(1,k,l,m) = 0.
            q3(nx+1,k,l,m) = 0.
         else if (kk.gt.(ny+1)) then
            if (k.eq.1) then
               do 710 j = 1, nxs
               q3(nx+j+1,k,l,m) = q3(nx-j+1,k,l,m)
  710          continue
            else
               do 720 j = 1, nxs
               q3(nx+j+1,k,l,m) = q3(nx-j+1,kyp2-k+2,l,m)
  720          continue
            endif
            q3(1,k,l,m) = 0.
            q3(nx+1,k,l,m) = 0.
         endif
  730    continue
      else if (ll.gt.(nz+1)) then
         if (l.eq.1) then
            do 770 k = 1, kyp2
            kk = k + koff
            if (kk.le.ny) then
               do 740 j = 1, nxs
               q3(nx+j+1,k,l,m) = q3(nx-j+1,k,l,m)
  740          continue
            else if (kk.gt.(ny+1)) then
               if (k.eq.1) then
                  do 750 j = 1, nxs
                  q3(nx+j+1,k,l,m) = -q3(nx-j+1,k,l,m)
  750             continue
               else
                  do 760 j = 1, nxs
                  q3(nx+j+1,k,l,m) = -q3(nx-j+1,kyp2-k+2,l,m)
  760             continue
               endif
            endif
  770       continue
         else
            do 810 k = 1, kyp2
            kk = k + koff
            if (kk.le.ny) then
               do 780 j = 1, nxs
               q3(nx+j+1,k,l,m) = q3(nx-j+1,k,kzp2-l+2,m)
  780          continue
            else if (kk.gt.(ny+1)) then
               if (k.eq.1) then
                  do 790 j = 1, nxs
                  q3(nx+j+1,k,l,m) = -q3(nx-j+1,k,kzp2-l+2,m)
  790             continue
               else
                  do 800 j = 1, nxs
                  q3(nx+j+1,k,l,m) = -q3(nx-j+1,kyp2-k+2,kzp2-l+2,m)
  800             continue
               endif
            endif
  810       continue
         endif
         do 830 k = 1, kyp2
         kk = k + koff
         if ((kk.eq.1).or.(kk.eq.(ny+1))) then
            do 820 j = 1, nx
            q3(j,k,l,m) = 0.
            q3(j+nx,k,l,m) = 0.
  820       continue
         else
            q3(1,k,l,m) = 0.
            q3(nx+1,k,l,m) = 0.
         endif
  830    continue
      endif
  840 continue
  850 continue
  860 continue
c finish odd array
      do 940 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      loff = kzp2*(mz + ks)
      do 930 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      do 920 l = 1, kzp2
      ll = l + loff
      if (ll.le.nz) then
         do 890 k = 1, kyp2
         kk = k + koff
         if (kk.gt.(ny+1)) then
            if (k.eq.1) then
               do 870 j = 1, nxs
               q3(nx-j+1,k,l,m) = -q3(nx+j+1,k,l,m)
  870          continue
               q3(nx+1,k,l,m) = -q3(nx+1,k,l,m)
            else
               do 880 j = 1, nxs
               q3(nx-j+1,k,l,m) = -q3(nx+j+1,k,l,m)
  880          continue
               q3(nx+1,k,l,m) = -q3(nx+1,k,l,m)
            endif
         endif
  890    continue
      else if (ll.gt.(nz+1)) then
         do 910 k = 1, kyp2
         kk = k + koff
         if ((kk.le.ny).or.(kk.gt.(ny+1))) then
            do 900 j = 1, nxs
            q3(nx-j+1,k,l,m) = -q3(nx+j+1,k,l,m)
  900       continue
            q3(nx+1,k,l,m) = -q3(nx+1,k,l,m)
         endif
  910    continue
      endif
  920 continue
  930 continue
  940 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PHAFTRP32C(fxyz,fxyz3,scb,scd,nx,ny,nz,kstrt,nvpy,nxv,k
     1yp,kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
c this subroutine copies data from a triple array to regular array
c with guard cells for vector field and linear interpolation
c for distributed data with 2D domain decomposition
c fxyz array may be modified
c scb/scd = scratch arrays
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c nvpy = number of real or virtual processors in y
c nxv = second dimension of input array fxyz, must be >= nx
c kyp/kzp = number of data values per block in fxyz in y/z
c kypd = third dimension of input array fxyz, must be >= kyp
c kzpd = fourth dimension of input array fxyz, must be >= kzp
c kyp2/kzp2 = number of data values per block in fxyz3 in y/z
c kblok/lblok = number of data blocks in y/z
c k2blok/l2blok = number of data blocks in y/z for tripled data
      implicit none
      real fxyz, fxyz3, scb, scd
      integer nx, ny, nz, kstrt, nvpy, nxv, kyp, kzp, kypd, kzpd
      integer kyp2, kzp2, kblok, lblok, k2blok, l2blok
      dimension fxyz(3,nxv,kypd,kzpd,kblok*lblok)
      dimension fxyz3(3,2*nxv,2*kypd,kzp2,k2blok*l2blok)
      dimension scb(3,nxv,kypd,kzpd,kblok*lblok)
      dimension scd(3,nxv,kzpd,2,kblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer i, j, k, l, m, my, mz, nx1, ny1, nz1, kzp1, ks, js, kr, ls
      integer kyb, kyb2, kzb, kzb2, nxvy, noff, moff, loff, koff, jm, km
      integer jl, kl, kz, mm, lm, is, ii
      dimension istatus(lstat)
      nx1 = nx + 1
      ny1 = ny + 1
      nz1 = nz + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      kyb = ny/kyp
      kyb2 = (ny + ny)/kyp2
      kzb = nz/kzp
      kzb2 = (nz + nz)/kzp2
      kzp1 = kzp + 1
      nxvy = nxv*kypd
      noff = kypd*kzpd + kyb*kzb
c copy from triple array in y direction
      do 370 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 360 my = 1, k2blok
      m = my + moff
      loff = kyp2*(my + js)
      jl = loff/kyp + 1
      loff = kzp2*(mz + ks)
      kl = loff/kzp + 1
      loff = kyp*(my + js)
      jm = loff/kyp2 + 1
      koff = loff - kyp2*(jm - 1)
      loff = kzp*(mz + ks)
      kr = loff/kzp2
      km = loff - kzp2*kr
      mm = jm + nvpy*kr
      kz = nvpy*(mz + ks)
      lm = jl + 2*kz
c special case for one processor
      if ((kyb2*kzb2).eq.1) then
         do 40 l = 1, nz1
         do 30 k = 1, ny1
         do 20 j = 1, nx1
         do 10 i = 1, 3
         fxyz(i,j,k,l,m) = fxyz3(i,j,k,l,m)
   10    continue
   20    continue
   30    continue
   40    continue
         return
      endif
      if (km.eq.0) then
         ls = kzp1
         is = 1
      else
         ls = kzp
         is = 2
      endif
c this segment is used for shared memory computers
c     if (jm.le.kyb) then
c        if ((koff.eq.0).and.(kyp.lt.kyp2)) then
c           do 80 l = 1, kzp1
c           do 70 k = 1, kyp1
c           do 60 j = 1, nx1
c           do 50 i = 1, 3
c           fxyz(i,j,k,l,m) = fxyz3(i,j,k,l,mm)
c  50       continue
c  60       continue
c  70       continue
c  80       continue
c        else
c           do 120 l = 1, kzp1
c           do 110 k = 1, kyp
c           do 100 j = 1, nx1
c           do 90 i = 1, 3
c           fxyz(i,j,k,l,m) = fxyz3(i,j,k+koff,l,mm)
c  90       continue
c 100       continue
c 110       continue
c 120       continue
c           do 150 l = 1, kzp1
c           do 140 j = 1, nx1
c           do 130 i = 1, 3
c           fxyz(i,j,kyp+1,l,m) = fxyz3(i,j,1,l,mm+1)
c 130       continue
c 140       continue
c 150       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jm.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_IRECV(fxyz(1,1,1,1,m),3*ls*nxvy,mreal,mm-1,noff+is,
     1lgrp,msid,ierr)
         else
            call MPI_IRECV(fxyz(1,1,1,1,m),3*ls*nxvy,mreal,mm-1,noff+is,
     1lgrp,msid,ierr)
            call MPI_IRECV(scd(1,1,1,2,m),3*ls*nxv,mreal,mm,noff+is,lgrp
     1,nsid,ierr)
         endif
      endif
      do 230 ii = 1, 2
c for odd rows in z
      if (ii.eq.1) then
         loff = 0
         ls = kzp1
         lm = jl + 2*kz
c for even rows in z
      else if (ii.eq.2) then
         if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 230
         loff = kzp
         ls = kzp
         lm = lm + nvpy
      endif
c pack data and send it
      if ((jl.le.kyb).and.(kl.le.kzb)) then
         if (kyp.lt.kyp2) then
            do 80 l = 1, ls
            do 70 k = 1, kypd
            do 60 j = 1, nxv
            do 50 i = 1, 3
            scb(i,j,k,l,m) = fxyz3(i,j,k,l+loff,m)
   50       continue
   60       continue
   70       continue
   80       continue
            call MPI_SEND(scb(1,1,1,1,m),3*ls*nxvy,mreal,lm-1,noff+ii,lg
     1rp,ierr)
            if (kyb2.gt.1) then
               do 120 l = 1, ls
               do 110 k = 1, kypd
               do 100 j = 1, nxv
               do 90 i = 1, 3
               scb(i,j,k,l,m) = fxyz3(i,j,k+kyp,l+loff,m)
   90          continue
  100          continue
  110          continue
  120          continue
               call MPI_SEND(scb(1,1,1,1,m),3*ls*nxvy,mreal,lm,noff+ii,l
     1grp,ierr)
            endif
         else
            do 160 l = 1, ls
            do 150 k = 1, kypd
            do 140 j = 1, nxv
            do 130 i = 1, 3
            scb(i,j,k,l,m) = fxyz3(i,j,k,l+loff,m)
  130       continue
  140       continue
  150       continue
  160       continue
            call MPI_SEND(scb(1,1,1,1,m),3*ls*nxvy,mreal,lm-1,noff+ii,lg
     1rp,ierr)
         endif
         if (jl.gt.1) then
            do 190 l = 1, ls
            do 180 j = 1, nxv
            do 170 i = 1, 3
            scd(i,j,l,1,m) = fxyz3(i,j,1,l+loff,m)
  170       continue
  180       continue
  190       continue
            call MPI_SEND(scd(1,1,1,1,m),3*ls*nxv,mreal,lm-2,noff+ii,lgr
     1p,ierr)
         endif
      else if ((jl.eq.(kyb+1)).and.(kl.le.kzb)) then
         do 220 l = 1, ls
         do 210 j = 1, nxv
         do 200 i = 1, 3
         scd(i,j,l,1,m) = fxyz3(i,j,1,l+loff,m)
  200    continue
  210    continue
  220    continue
         call MPI_SEND(scd(1,1,1,1,m),3*ls*nxv,mreal,lm-2,noff+ii,lgrp,i
     1err)
      endif
  230 continue
c wait for data
      if (km.eq.0) then
         ls = kzp1
      else
         ls = kzp
      endif
      if (jm.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_WAIT(msid,istatus,ierr)
         else
            call MPI_WAIT(msid,istatus,ierr)
            call MPI_WAIT(nsid,istatus,ierr)
            do 260 l = 1, ls
            do 250 j = 1, nxv
            do 240 i = 1, 3
            fxyz(i,j,kyp+1,l,m) = scd(i,j,l,2,m)
  240       continue
  250       continue
  260       continue
         endif
      endif
c last point
      if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 360
      mm = mm + nvpy - 1
      if (koff.eq.0) then
         is = 3
      else
         is = 4
      endif
c this segment is used for shared memory computers
c     if ((kr.lt.kzb).and.(km.ne.0)) then
c        do 290 k = 1, kypd
c        do 280 j = 1, nx1
c        do 270 i = 1, 3
c        fxyz(i,j,k,kzp+1,m) = fxyz3(i,j,k,1,mm+1)
c 270    continue
c 280    continue
c 290    continue
c     endif
c this segment is used for mpi computers
      if ((kr.lt.kzb).and.(km.ne.0)) then
         call MPI_IRECV(fxyz(1,1,1,kzp+1,m),3*nxvy,mreal,mm,noff+is,lgrp
     1,nsid,ierr)
      endif
      do 330 ii = 3, 4
c for odd rows in z
      if (ii.eq.3) then
         loff = 0
         lm = jl + 2*kz - nvpy + 1
c for even rows in z
      else if (ii.eq.4) then
         if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 330
         loff = kyp
         lm = lm + 1
      endif
c pack data and send it
      if (jl.le.kyb) then
         if ((kl.le.kzb).and.(kl.gt.1)) then
            do 290 k = 1, kypd
            do 280 j = 1, nxv
            do 270 i = 1, 3
            scb(i,j,k,1,m) = fxyz3(i,j,k+loff,1,m)
  270       continue
  280       continue
  290       continue
            call MPI_SEND(scb(1,1,1,1,m),3*nxvy,mreal,lm-2,noff+ii,lgrp,
     1ierr)
         else if (kl.eq.(kzb+1)) then
            do 320 k = 1, kypd
            do 310 j = 1, nxv
            do 300 i = 1, 3
            scb(i,j,k,1,m) = fxyz3(i,j,k+loff,1,m)
  300       continue
  310       continue
  320       continue
            call MPI_SEND(scb(1,1,1,1,m),3*nxvy,mreal,lm-2,noff+ii,lgrp,
     1ierr)
         endif
      endif
  330 continue
c wait for data
      if ((kr.lt.kzb).and.(km.ne.0)) then
         call MPI_WAIT(nsid,istatus,ierr)
      endif
c last point
      if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 360
      mm = mm + 1
      lm = jl + 2*kz - nvpy
c this segment is used for shared memory computers
c     if ((kr.lt.kzb).and.(km.ne.0).and.(koff.ne.0)) then
c        do 350 j = 1, nx1
c        do 340 i = 1, 3
c        fxyz(i,j,kyp+1,kzp+1,m) = fxyz3(i,j,1,1,mm+1)
c 340    continue
c 350    continue
c     endif
c this segment is used for mpi computers
      if ((kr.lt.kzb).and.(km.ne.0).and.(koff.ne.0)) then
         call MPI_IRECV(fxyz(1,1,kyp+1,kzp+1,m),3*nxv,mreal,mm,noff+5,lg
     1rp,nsid,ierr)
      endif
c send data
      if ((jl.le.(kyb+1)).and.(jl.gt.1)) then
         if ((kl.le.(kzb+1)).and.(kl.gt.1)) then
            call MPI_SEND(fxyz3(1,1,1,1,m),3*nxv,mreal,lm-2,noff+5,lgrp,
     1ierr)
         endif
      endif
c wait for data
      if ((kr.lt.kzb).and.(km.ne.0).and.(koff.ne.0)) then
         call MPI_WAIT(nsid,istatus,ierr)
      endif
  360 continue
  370 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PHAFTRP32D(q,q3,scb,scd,nx,ny,nz,kstrt,nvpy,nxv,kyp,kzp
     1,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
c this subroutine copies data from a triple array to regular array
c with guard cells for scalar field and linear interpolation
c for distributed data with 2D domain decomposition
c q array may be modified
c scb/scd = scratch arrays
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c nvpy = number of real or virtual processors in y
c nxv = second dimension of input array q, must be >= nx
c kyp/kzp = number of data values per block in q in y/z
c kypd = third dimension of input array q, must be >= kyp
c kzpd = fourth dimension of input array q, must be >= kzp
c kyp2/kzp2 = number of data values per block in q3 in y/z
c kblok/lblok = number of data blocks in y/z
c k2blok/l2blok = number of data blocks in y/z for tripled data
      implicit none
      real q, q3, scb, scd
      integer nx, ny, nz, kstrt, nvpy, nxv, kyp, kzp, kypd, kzpd
      integer kyp2, kzp2, kblok, lblok, k2blok, l2blok
      dimension q(nxv,kypd,kzpd,kblok*lblok)
      dimension q3(2*nxv,2*kypd,kzp2,k2blok*l2blok)
      dimension scb(nxv,kypd,kzpd,kblok*lblok)
      dimension scd(nxv,kzpd,2,kblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer j, k, l, m, my, mz, nx1, ny1, nz1, kzp1, ks, js, kr, ls
      integer kyb, kyb2, kzb, kzb2, nxvy, noff, moff, loff, koff, jm, km
      integer jl, kl, kz, mm, lm, is, ii
      dimension istatus(lstat)
      nx1 = nx + 1
      ny1 = ny + 1
      nz1 = nz + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      kyb = ny/kyp
      kyb2 = (ny + ny)/kyp2
      kzb = nz/kzp
      kzb2 = (nz + nz)/kzp2
      kzp1 = kzp + 1
      nxvy = nxv*kypd
      noff = kypd*kzpd + kyb*kzb
c copy from triple array in y direction
      do 270 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 260 my = 1, k2blok
      m = my + moff
      loff = kyp2*(my + js)
      jl = loff/kyp + 1
      loff = kzp2*(mz + ks)
      kl = loff/kzp + 1
      loff = kyp*(my + js)
      jm = loff/kyp2 + 1
      koff = loff - kyp2*(jm - 1)
      loff = kzp*(mz + ks)
      kr = loff/kzp2
      km = loff - kzp2*kr
      mm = jm + nvpy*kr
      kz = nvpy*(mz + ks)
      lm = jl + 2*kz
c special case for one processor
      if ((kyb2*kzb2).eq.1) then
         do 30 l = 1, nz1
         do 20 k = 1, ny1
         do 10 j = 1, nx1
         q(j,k,l,m) = q3(j,k,l,m)
   10    continue
   20    continue
   30    continue
         return
      endif
      if (km.eq.0) then
         ls = kzp1
         is = 1
      else
         ls = kzp
         is = 2
      endif
c this segment is used for shared memory computers
c     if (jm.le.kyb) then
c        if ((koff.eq.0).and.(kyp.lt.kyp2)) then
c           do 60 l = 1, kzp1
c           do 50 k = 1, kyp1
c           do 40 j = 1, nx1
c           q(j,k,l,m) = q3(j,k,l,mm)
c  40       continue
c  50       continue
c  60       continue
c        else
c           do 90 l = 1, kzp1
c           do 80 k = 1, kyp
c           do 70 j = 1, nx1
c           q(j,k,l,m) = q3(j,k+koff,l,mm)
c  70       continue
c  80       continue
c  90       continue
c           do 110 l = 1, kzp1
c           do 100 j = 1, nx1
c           q(j,kyp+1,l,m) = q3(j,1,l,mm+1)
c 100       continue
c 110       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jm.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_IRECV(q(1,1,1,m),ls*nxvy,mreal,mm-1,noff+is,lgrp,ms
     1id,ierr)
         else
            call MPI_IRECV(q(1,1,1,m),ls*nxvy,mreal,mm-1,noff+is,lgrp,ms
     1id,ierr)
            call MPI_IRECV(scd(1,1,2,m),ls*nxv,mreal,mm,noff+is,lgrp,nsi
     1d,ierr)
         endif
      endif
      do 170 ii = 1, 2
c for odd rows in z
      if (ii.eq.1) then
         loff = 0
         ls = kzp1
         lm = jl + 2*kz
c for even rows in z
      else if (ii.eq.2) then
         if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 170
         loff = kzp
         ls = kzp
         lm = lm + nvpy
      endif
c pack data and send it
      if ((jl.le.kyb).and.(kl.le.kzb)) then
         if (kyp.lt.kyp2) then
            do 60 l = 1, ls
            do 50 k = 1, kypd
            do 40 j = 1, nxv
            scb(j,k,l,m) = q3(j,k,l+loff,m)
   40       continue
   50       continue
   60       continue
            call MPI_SEND(scb(1,1,1,m),ls*nxvy,mreal,lm-1,noff+ii,lgrp,i
     1err)
            if (kyb2.gt.1) then
               do 90 l = 1, ls
               do 80 k = 1, kypd
               do 70 j = 1, nxv
               scb(j,k,l,m) = q3(j,k+kyp,l+loff,m)
   70          continue
   80          continue
   90          continue
               call MPI_SEND(scb(1,1,1,m),ls*nxvy,mreal,lm,noff+ii,lgrp,
     1ierr)
            endif
         else
            do 120 l = 1, ls
            do 110 k = 1, kypd
            do 100 j = 1, nxv
            scb(j,k,l,m) = q3(j,k,l+loff,m)
  100       continue
  110       continue
  120       continue
            call MPI_SEND(scb(1,1,1,m),ls*nxvy,mreal,lm-1,noff+ii,lgrp,i
     1err)
         endif
         if (jl.gt.1) then
            do 140 l = 1, ls
            do 130 j = 1, nxv
            scd(j,l,1,m) = q3(j,1,l+loff,m)
  130       continue
  140       continue
            call MPI_SEND(scd(1,1,1,m),ls*nxv,mreal,lm-2,noff+ii,lgrp,ie
     1rr)
         endif
      else if ((jl.eq.(kyb+1)).and.(kl.le.kzb)) then
         do 160 l = 1, ls
         do 150 j = 1, nxv
         scd(j,l,1,m) = q3(j,1,l+loff,m)
  150    continue
  160    continue
         call MPI_SEND(scd(1,1,1,m),ls*nxv,mreal,lm-2,noff+ii,lgrp,ierr)
      endif
  170 continue
c wait for data
      if (km.eq.0) then
         ls = kzp1
      else
         ls = kzp
      endif
      if (jm.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_WAIT(msid,istatus,ierr)
         else
            call MPI_WAIT(msid,istatus,ierr)
            call MPI_WAIT(nsid,istatus,ierr)
            do 190 l = 1, ls
            do 180 j = 1, nxv
            q(j,kyp+1,l,m) = scd(j,l,2,m)
  180       continue
  190       continue
         endif
      endif
c last point
      if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 260
      mm = mm + nvpy - 1
      if (koff.eq.0) then
         is = 3
      else
         is = 4
      endif
c this segment is used for shared memory computers
c     if ((kr.lt.kzb).and.(km.ne.0)) then
c        do 210 k = 1, kypd
c        do 200 j = 1, nx1
c        q(j,k,kzp+1,m) = q3(j,k,1,mm+1)
c 200    continue
c 210    continue
c     endif
c this segment is used for mpi computers
      if ((kr.lt.kzb).and.(km.ne.0)) then
         call MPI_IRECV(q(1,1,kzp+1,m),nxvy,mreal,mm,noff+is,lgrp,nsid,i
     1err)
      endif
      do 240 ii = 3, 4
c for odd rows in z
      if (ii.eq.3) then
         loff = 0
         lm = jl + 2*kz - nvpy + 1
c for even rows in z
      else if (ii.eq.4) then
         if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 240
         loff = kyp
         lm = lm + 1
      endif
c pack data and send it
      if (jl.le.kyb) then
         if ((kl.le.kzb).and.(kl.gt.1)) then
            do 210 k = 1, kypd
            do 200 j = 1, nxv
            scb(j,k,1,m) = q3(j,k+loff,1,m)
  200       continue
  210       continue
            call MPI_SEND(scb(1,1,1,m),nxvy,mreal,lm-2,noff+ii,lgrp,ierr
     1)
         else if (kl.eq.(kzb+1)) then
            do 230 k = 1, kypd
            do 220 j = 1, nxv
            scb(j,k,1,m) = q3(j,k+loff,1,m)
  220       continue
  230       continue
            call MPI_SEND(scb(1,1,1,m),nxvy,mreal,lm-2,noff+ii,lgrp,ierr
     1)
         endif
      endif
  240 continue
c wait for data
      if ((kr.lt.kzb).and.(km.ne.0)) then
         call MPI_WAIT(nsid,istatus,ierr)
      endif
c last point
      if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 260
      mm = mm + 1
      lm = jl + 2*kz - nvpy
c this segment is used for shared memory computers
c     if ((kr.lt.kzb).and.(km.ne.0).and.(koff.ne.0)) then
c        do 250 j = 1, nx1
c        q(j,kyp+1,kzp+1,m) = q3(j,1,1,mm+1)
c 250    continue
c     endif
c this segment is used for mpi computers
      if ((kr.lt.kzb).and.(km.ne.0).and.(koff.ne.0)) then
         call MPI_IRECV(q(1,kyp+1,kzp+1,m),nxv,mreal,mm,noff+5,lgrp,nsid
     1,ierr)
      endif
c send data
      if ((jl.le.(kyb+1)).and.(jl.gt.1)) then
         if ((kl.le.(kzb+1)).and.(kl.gt.1)) then
            call MPI_SEND(q3(1,1,1,m),nxv,mreal,lm-2,noff+5,lgrp,ierr)
         endif
      endif
c wait for data
      if ((kr.lt.kzb).and.(km.ne.0).and.(koff.ne.0)) then
         call MPI_WAIT(nsid,istatus,ierr)
      endif
  260 continue
  270 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLCGUARD32(f,scs,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx,mbl
     1ok,nblok,kyp,kzp,ngds)
c this subroutine copies data from field to particle partitions, copying
c data to guard cells, where the field and particle partitions are 
c assumed to be the same.  for vector data
c the field is replicated so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y and z direction.
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes three extra
c guard cells.
c scs = scratch array for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer kyp, kzp, ngds
      real f, scs
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx,2*ngds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, mnblok
      integer nx2, nyp3, kzp1, nxvz, nxvy, m, my, mz, j, k, n
      integer jr, jrr, jl, jll
      dimension istatus(lstat)
      nx2 = nx + 2
      nyp3 = kyp + 3
      kzp1 = kzp + 1
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
      do 40 m = 1, mnblok
      do 30 k = 1, nzpmx
      do 20 j = 1, nxv
      do 10 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,kyp+1,k,m)
      scs(n,j,k,2,m) = f(n,j,2,k,m)
      scs(n,j,k,3,m) = f(n,j,ngc+1,k,m)
      scs(n,j,k,5,m) = f(n,j,kyp+2,k,m)
   10 continue
   20 continue
   30 continue
   40 continue
c copy to guard cells in y
      do 270 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 260 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr
      krr = jr + kz
      jl = ky - 1
      jll = jl
      kll = jl + kz
c special case of only one grid per processor
      if (kyp.eq.1) then
         jrr = jr + 1
         krr = jrr + kz
         jll = jl - 1
         kll = jll + kz
      endif
      kr = jr + kz
      kl = jl + kz
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 70 k = 1, nzpmx
c        do 60 j = 1, nxv
c        do 50 n = 1, 3
c        scs(n,j,k,4,m) = scs(n,j,k,1,kl)
c  50    continue
c  60    continue
c  70    continue
c     else
c        do 100 k = 1, nzpmx
c        do 90 j = 1, nxv
c        do 80 n = 1, 3
c        scs(n,j,k,4,m) = f(n,j,3,k,m)
c  80    continue
c  90    continue
c 100    continue
c     endif
c     if (jr.lt.nvpy) then
c        do 130 k = 1, nzpmx
c        do 120 j = 1, nxv
c        do 110 n = 1, 3
c        scs(n,j,k,5,m) = scs(n,j,k,2,kr)
c        scs(n,j,k,6,m) = scs(n,j,k,3,krr)
c 110    continue
c 120    continue
c 130    continue
c     else
c        do 160 k = 1, kzp1
c        do 150 j = 2, nx2
c        do 140 n = 1, 3
c        scs(n,j,k+1,6,m) = 2.*scs(n,j,k+1,5,m) - scs(n,j,k+1,1,m)
c 140    continue
c 150    continue
c 160    continue
c     endif
c     if (kyp.eq.1) then
c        if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
c           do 190 k = 1, nzpmx
c           do 180 j = 1, nxv
c           do 170 n = 1, 3
c           scs(n,j,k,4,m) = scs(n,j,k,2,kr)
c 170       continue
c 180       continue
c 190       continue
c        endif
c        if (jr.eq.(nvpy-1)) then
c           do 220 k = 1, nzpmx
c           do 210 j = 1, nxv
c           do 200 n = 1, 3
c           scs(n,j,k,6,m) = scs(n,j,k,5,kr)
c 200       continue
c 210       continue
c 220       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,1,4,m),3*nxvz,mreal,kl-1,noff+3,lgrp,msi
     1d,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,1,m),3*nxvz,mreal,kr-1,noff+3,lgrp,ierr
     1)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 70 k = 1, nzpmx
         do 60 j = 1, nxv
         do 50 n = 1, 3
         scs(n,j,k,4,m) = f(n,j,3,k,m)
   50    continue
   60    continue
   70    continue
      endif
      if (jr.lt.nvpy) then
         call MPI_IRECV(scs(1,1,1,5,m),3*ngc*nxvz,mreal,kr-1,noff+4,lgrp
     1,msid,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scs(1,1,1,2,m),3*ngc*nxvz,mreal,kl-1,noff+4,lgrp,
     1ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 100 k = 1, kzp1
         do 90 j = 2, nx2
         do 80 n = 1, 3
         scs(n,j,k+1,6,m) = 2.*scs(n,j,k+1,5,m) - scs(n,j,k+1,1,m)
   80    continue
   90    continue
  100    continue
      endif
c special case of only one grid per processor
      if (kyp.eq.1) then
         if (jrr.lt.nvpy) then
            call MPI_IRECV(scs(1,1,1,6,m),3*nxvz,mreal,krr-1,noff+5,lgrp
     1,msid,ierr)
         else if (jr.lt.nvpy) then
            call MPI_IRECV(scs(1,1,1,6,m),3*nxvz,mreal,kr-1,noff+5,lgrp,
     1msid,ierr)
         endif
         if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
            call MPI_IRECV(scs(1,1,1,4,m),3*nxvz,mreal,kr-1,noff+5,lgrp,
     1nsid,ierr)
         endif
         if (jll.ge.0) then
            call MPI_SEND(scs(1,1,1,3,m),3*nxvz,mreal,kll-1,noff+5,lgrp,
     1ierr)
         else if (jl.eq.0) then
            call MPI_SEND(scs(1,1,1,3,m),3*nxvz,mreal,kl-1,noff+5,lgrp,i
     1err)
         endif
         if ((jl.eq.(nvpy-2)).and.(jl.ge.0)) then
            call MPI_SEND(scs(1,1,1,5,m),3*nxvz,mreal,kl-1,noff+5,lgrp,i
     1err)
         endif
         if (jr.lt.nvpy) then
            call MPI_WAIT(msid,istatus,ierr)
         endif
         if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
            call MPI_WAIT(nsid,istatus,ierr)
         endif
      endif
      do 250 k = 1, nzpmx
      do 240 j = 1, nxv
      do 230 n = 1, 3
      f(n,j,1,k,m) = scs(n,j,k,4,m)
      f(n,j,kyp+2,k,m) = scs(n,j,k,5,m)
      f(n,j,kyp+3,k,m) = scs(n,j,k,6,m)
  230 continue
  240 continue
  250 continue
  260 continue
  270 continue
c fix left edge
      do 320 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 310 my = 1, mblok
      m = my + moff
      kl = my + js
      if (kl.eq.0) then
         do 300 k = 1, kzp1
         do 290 j = 2, nx2
         do 280 n = 1, 3
         f(n,j,1,k+1,m) = 2.*f(n,j,2,k+1,m) - f(n,j,1,k+1,m)
  280    continue
  290    continue
  300    continue
      endif
  310 continue
  320 continue
c copy to guard cells in z
      do 520 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 510 my = 1, mblok
      m = my + moff
      ky = my + js + 1
      kz = mz + ks
      jr = kz + 1
      jrr = jr
      krr = ky + nvpy*jr
      jl = kz - 1
      jll = jl
      kll = ky + nvpy*jl
      ngc = 2
c special case of only one grid per processor
      if (kzp.eq.1) then
         jrr = jr + 1
         krr = ky + nvpy*jrr
         jll = jl - 1
         kll = ky + nvpy*jll
         ngc = 1
      endif
      kr = ky + nvpy*jr
      kl = ky + nvpy*jl
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 350 k = 1, nypmx
c        do 340 j = 1, nxv
c        do 330 n = 1, 3
c        f(n,j,k,1,m) = f(n,j,k,kzp+1,kl)
c 330    continue
c 340    continue
c 350    continue
c     else
c        do 380 k = 1, nypmx
c        do 370 j = 1, nxv
c        do 360 n = 1, 3
c        f(n,j,k,1,m) = f(n,j,k,3,m)
c 360    continue
c 370    continue
c 380    continue
c     endif
c     if (jr.lt.nvpz) then
c        do 410 k = 1, nypmx
c        do 400 j = 1, nxv
c        do 390 n = 1, 3
c        f(n,j,k,kzp+2,m) = f(n,j,k,2,kr)
c        f(n,j,k,kzp+3,m) = f(n,j,k,ngc+1,krr)
c 390    continue
c 400    continue
c 410    continue
c     else
c        do 440 k = 1, nyp3
c        do 430 j = 2, nx2
c        do 420 n = 1, 3
c        f(n,j,k,kzp+3,m) = 2.*f(n,j,k,kzp+2,m) - f(n,j,k,kzp+1,m)
c 420    continue
c 430    continue
c 440    continue
c     endif
c     if (kzp.eq.1) then
c        if ((jl.eq.(-1)).and.(jr.lt.nvpz)) then
c           do 470 k = 1, nypmx
c           do 460 j = 1, nxv
c           do 450 n = 1, 3
c           f(n,j,k,1,m) = f(n,j,k,2,kr)
c 450       continue
c 460       continue
c 470       continue
c        endif
c        if (jr.eq.(nvpz-1)) then
c           do 500 k = 1, nypmx
c           do 490 j = 1, nxv
c           do 480 n = 1, 3
c           f(n,j,k,kzp+3,m) = f(n,j,k,kzp+2,kr)
c 480       continue
c 490       continue
c 500       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(f(1,1,1,1,m),3*nxvy,mreal,kl-1,noff+6,lgrp,msid,
     1ierr)
      endif
      if (jr.lt.nvpz) then
         call MPI_SEND(f(1,1,1,kzp+1,m),3*nxvy,mreal,kr-1,noff+6,lgrp,ie
     1rr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 350 k = 1, nypmx
         do 340 j = 1, nxv
         do 330 n = 1, 3
         f(n,j,k,1,m) = f(n,j,k,3,m)
  330    continue
  340    continue
  350    continue
      endif
      if (jr.lt.nvpz) then
         call MPI_IRECV(f(1,1,1,kzp+2,m),3*ngc*nxvy,mreal,kr-1,noff+8,lg
     1rp,msid,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(f(1,1,1,2,m),3*ngc*nxvy,mreal,kl-1,noff+8,lgrp,ie
     1rr)
      endif
      if (jr.lt.nvpz) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 380 k = 1, nyp3
         do 370 j = 2, nx2
         do 360 n = 1, 3
         f(n,j,k,kzp+3,m) = 2.*f(n,j,k,kzp+2,m) - f(n,j,k,kzp+1,m)
  360    continue
  370    continue
  380    continue
      endif
c special case of only one grid per processor
      if (kzp.eq.1) then
         if (jrr.lt.nvpz) then
            call MPI_IRECV(f(1,1,1,kzp+3,m),3*nxvy,mreal,krr-1,noff+9,lg
     1rp,msid,ierr)
         else if (jr.lt.nvpz) then
            call MPI_IRECV(f(1,1,1,kzp+3,m),3*nxvy,mreal,kr-1,noff+9,lgr
     1p,msid,ierr)
         endif
         if ((jl.eq.(-1)).and.(jr.lt.nvpz)) then
            call MPI_IRECV(f(1,1,1,1,m),3*nxvy,mreal,kr-1,noff+9,lgrp,ns
     1id,ierr)
         endif
         if (jll.ge.0) then
            call MPI_SEND(f(1,1,1,2,m),3*nxvy,mreal,kll-1,noff+9,lgrp,ie
     1rr)
         else if (jl.eq.0) then
            call MPI_SEND(f(1,1,1,2,m),3*nxvy,mreal,kl-1,noff+9,lgrp,ier
     1r)
         endif
         if ((jl.eq.(nvpz-2)).and.(jl.ge.0)) then
            call MPI_SEND(f(1,1,1,kzp+2,m),3*nxvy,mreal,kl-1,noff+9,lgrp
     1,ierr)
         endif
         if (jr.lt.nvpz) then
            call MPI_WAIT(msid,istatus,ierr)
         endif
         if ((jl.eq.(-1)).and.(jr.lt.nvpz)) then
            call MPI_WAIT(nsid,istatus,ierr)
         endif
      endif
  510 continue
  520 continue
c fix left edge
      do 570 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 560 my = 1, mblok
      m = my + moff
      kl = mz + ks
      if (kl.eq.0) then
         do 550 k = 1, nyp3
         do 540 j = 2, nx2
         do 530 n = 1, 3
         f(n,j,k,1,m) = 2.*f(n,j,k,2,m) - f(n,j,k,1,m)
  530    continue
  540    continue
  550    continue
      endif
  560 continue
  570 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLDGUARD32(f,scs,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx,mbl
     1ok,nblok,kyp,kzp,ngds)
c this subroutine copies data from field to particle partitions, copying
c data to guard cells, where the field and particle partitions are 
c assumed to be the same.
c the field is replicated so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y and z direction.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes three extra
c guard cells.
c scs = scratch array for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer kyp, kzp, ngds
      real f, scs
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, mnblok
      integer nx2, nyp3, kzp1, nxvz, nxvy, m, my, mz, j, k
      integer jr, jrr, jl, jll
      dimension istatus(lstat)
      nx2 = nx + 2
      nyp3 = kyp + 3
      kzp1 = kzp + 1
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
      scs(j,k,5,m) = f(j,kyp+2,k,m)
   10 continue
   20 continue
   30 continue
c copy to guard cells in y
      do 190 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 180 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr
      krr = jr + kz
      jl = ky - 1
      jll = jl
      kll = jl + kz
c special case of only one grid per processor
      if (kyp.eq.1) then
         jrr = jr + 1
         krr = jrr + kz
         jll = jl - 1
         kll = jll + kz
      endif
      kr = jr + kz
      kl = jl + kz
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 50 k = 1, nzpmx
c        do 40 j = 1, nxv
c        scs(j,k,4,m) = scs(j,k,1,kl)
c  40    continue
c  50    continue
c     else
c        do 70 k = 1, nzpmx
c        do 60 j = 1, nxv
c        scs(j,k,4,m) = f(j,3,k,m)
c  60    continue
c  70    continue
c     endif
c     if (jr.lt.nvpy) then
c        do 90 k = 1, nzpmx
c        do 80 j = 1, nxv
c        scs(j,k,5,m) = scs(j,k,2,kr)
c        scs(j,k,6,m) = scs(j,k,3,krr)
c  80    continue
c  90    continue
c     else
c        do 110 k = 1, kzp1
c        do 100 j = 2, nx2
c        scs(j,k+1,6,m) = 2.*scs(j,k+1,5,m) - scs(j,k+1,1,m)
c 100    continue
c 110    continue
c     endif
c     if (kyp.eq.1) then
c        if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
c           do 130 k = 1, nzpmx
c           do 120 j = 1, nxv
c           scs(j,k,4,m) = scs(j,k,2,kr)
c 120       continue
c 130       continue
c        endif
c        if (jr.eq.(nvpy-1)) then
c           do 150 k = 1, nzpmx
c           do 140 j = 1, nxv
c           scs(j,k,6,m) = scs(j,k,5,kr)
c 140       continue
c 150       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,4,m),nxvz,mreal,kl-1,noff+3,lgrp,msid,ie
     1rr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,m),nxvz,mreal,kr-1,noff+3,lgrp,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 50 k = 1, nzpmx
         do 40 j = 1, nxv
         scs(j,k,4,m) = f(j,3,k,m)
   40    continue
   50    continue
      endif
      if (jr.lt.nvpy) then
         call MPI_IRECV(scs(1,1,5,m),ngc*nxvz,mreal,kr-1,noff+4,lgrp,msi
     1d,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scs(1,1,2,m),ngc*nxvz,mreal,kl-1,noff+4,lgrp,ierr
     1)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 70 k = 1, kzp1
         do 60 j = 2, nx2
         scs(j,k+1,6,m) = 2.*scs(j,k+1,5,m) - scs(j,k+1,1,m)
   60    continue
   70    continue
      endif
c special case of only one grid per processor
      if (kyp.eq.1) then
         if (jrr.lt.nvpy) then
            call MPI_IRECV(scs(1,1,6,m),nxvz,mreal,krr-1,noff+5,lgrp,msi
     1d,ierr)
         else if (jr.lt.nvpy) then
            call MPI_IRECV(scs(1,1,6,m),nxvz,mreal,kr-1,noff+5,lgrp,msid
     1,ierr)
         endif
         if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
            call MPI_IRECV(scs(1,1,4,m),nxvz,mreal,kr-1,noff+5,lgrp,nsid
     1,ierr)
         endif
         if (jll.ge.0) then
            call MPI_SEND(scs(1,1,3,m),nxvz,mreal,kll-1,noff+5,lgrp,ierr
     1)
         else if (jl.eq.0) then
            call MPI_SEND(scs(1,1,3,m),nxvz,mreal,kl-1,noff+5,lgrp,ierr)
         endif
         if ((jl.eq.(nvpy-2)).and.(jl.ge.0)) then
            call MPI_SEND(scs(1,1,5,m),nxvz,mreal,kl-1,noff+5,lgrp,ierr)
         endif
         if (jr.lt.nvpy) then
            call MPI_WAIT(msid,istatus,ierr)
         endif
         if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
            call MPI_WAIT(nsid,istatus,ierr)
         endif
      endif
      do 170 k = 1, nzpmx
      do 160 j = 1, nxv
      f(j,1,k,m) = scs(j,k,4,m)
      f(j,kyp+2,k,m) = scs(j,k,5,m)
      f(j,kyp+3,k,m) = scs(j,k,6,m)
  160 continue
  170 continue
  180 continue
  190 continue
c fix left edge
      do 230 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 220 my = 1, mblok
      m = my + moff
      kl = my + js
      if (kl.eq.0) then
         do 210 k = 1, kzp1
         do 200 j = 2, nx2
         f(j,1,k+1,m) = 2.*f(j,2,k+1,m) - f(j,1,k+1,m)
  200    continue
  210    continue
      endif
  220 continue
  230 continue
c copy to guard cells in z
      do 370 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 360 my = 1, mblok
      m = my + moff
      ky = my + js + 1
      kz = mz + ks
      jr = kz + 1
      jrr = jr
      krr = ky + nvpy*jr
      jl = kz - 1
      jll = jl
      kll = ky + nvpy*jl
      ngc = 2
c special case of only one grid per processor
      if (kzp.eq.1) then
         jrr = jr + 1
         krr = ky + nvpy*jrr
         jll = jl - 1
         kll = ky + nvpy*jll
         ngc = 1
      endif
      kr = ky + nvpy*jr
      kl = ky + nvpy*jl
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 250 k = 1, nypmx
c        do 240 j = 1, nxv
c        f(j,k,1,m) = f(j,k,kzp+1,kl)
c 240    continue
c 250    continue
c     else
c        do 270 k = 1, nypmx
c        do 260 j = 1, nxv
c        f(j,k,1,m) = f(j,k,3,m)
c 260    continue
c 270    continue
c     endif
c     if (jr.lt.nvpz) then
c        do 290 k = 1, nypmx
c        do 280 j = 1, nxv
c        f(j,k,kzp+2,m) = f(j,k,2,kr)
c        f(j,k,kzp+3,m) = f(j,k,ngc+1,krr)
c 280    continue
c 290    continue
c     else
c        do 310 k = 1, nyp3
c        do 300 j = 2, nx2
c        f(j,k,kzp+3,m) = 2.*f(j,k,kzp+2,m) - f(j,k,kzp+1,m)
c 300    continue
c 310    continue
c     endif
c     if (kzp.eq.1) then
c        if ((jl.eq.(-1)).and.(jr.lt.nvpz)) then
c           do 330 k = 1, nypmx
c           do 320 j = 1, nxv
c           f(j,k,1,m) = f(j,k,2,kr)
c 320       continue
c 330       continue
c        endif
c        if (jr.eq.(nvpz-1)) then
c           do 350 k = 1, nypmx
c           do 340 j = 1, nxv
c           f(j,k,kzp+3,m) = f(j,k,kzp+2,kr)
c 340       continue
c 350       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(f(1,1,1,m),nxvy,mreal,kl-1,noff+6,lgrp,msid,ierr
     1)
      endif
      if (jr.lt.nvpz) then
         call MPI_SEND(f(1,1,kzp+1,m),nxvy,mreal,kr-1,noff+6,lgrp,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 250 k = 1, nypmx
         do 240 j = 1, nxv
         f(j,k,1,m) = f(j,k,3,m)
  240    continue
  250    continue
      endif
      if (jr.lt.nvpz) then
         call MPI_IRECV(f(1,1,kzp+2,m),ngc*nxvy,mreal,kr-1,noff+8,lgrp,m
     1sid,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(f(1,1,2,m),ngc*nxvy,mreal,kl-1,noff+8,lgrp,ierr)
      endif
      if (jr.lt.nvpz) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 270 k = 1, nyp3
         do 260 j = 2, nx2
         f(j,k,kzp+3,m) = 2.*f(j,k,kzp+2,m) - f(j,k,kzp+1,m)
  260    continue
  270    continue
      endif
c special case of only one grid per processor
      if (kzp.eq.1) then
         if (jrr.lt.nvpz) then
            call MPI_IRECV(f(1,1,kzp+3,m),nxvy,mreal,krr-1,noff+9,lgrp,m
     1sid,ierr)
         else if (jr.lt.nvpz) then
            call MPI_IRECV(f(1,1,kzp+3,m),nxvy,mreal,kr-1,noff+9,lgrp,ms
     1id,ierr)
         endif
         if ((jl.eq.(-1)).and.(jr.lt.nvpz)) then
            call MPI_IRECV(f(1,1,1,m),nxvy,mreal,kr-1,noff+9,lgrp,nsid,i
     1err)
         endif
         if (jll.ge.0) then
            call MPI_SEND(f(1,1,2,m),nxvy,mreal,kll-1,noff+9,lgrp,ierr)
         else if (jl.eq.0) then
            call MPI_SEND(f(1,1,2,m),nxvy,mreal,kl-1,noff+9,lgrp,ierr)
         endif
         if ((jl.eq.(nvpz-2)).and.(jl.ge.0)) then
            call MPI_SEND(f(1,1,kzp+2,m),nxvy,mreal,kl-1,noff+9,lgrp,ier
     1r)
         endif
         if (jr.lt.nvpz) then
            call MPI_WAIT(msid,istatus,ierr)
         endif
         if ((jl.eq.(-1)).and.(jr.lt.nvpz)) then
            call MPI_WAIT(nsid,istatus,ierr)
         endif
      endif
  360 continue
  370 continue
c fix left edge
      do 410 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 400 my = 1, mblok
      m = my + moff
      kl = mz + ks
      if (kl.eq.0) then
         do 390 k = 1, nyp3
         do 380 j = 2, nx2
         f(j,k,1,m) = 2.*f(j,k,2,m) - f(j,k,1,m)
  380    continue
  390    continue
      endif
  400 continue
  410 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLCGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypmx
     1,nzpmx,mblok,nblok,ngds,idds,mter,nter)
c this subroutine copies data to guard cells in non-uniform partitions
c for vector data the field is replicated so as to disable quadratic
c interpolation within half a cell of the edges, and reduce it to linear
c interpolation in the y and z direction.
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes three extra guard cells.
c scs/scr = scratch arrays for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c it is assumed that nyzp(n,m) > 0.
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c mter/nter = (0,1) = (no,yes) pass data to next processor only in y/z
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds, mter, nter
      integer nyzp
      real f, scs, scr
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx*ngds,2,mblok*nblok)
      dimension scr(3,nxv,nypmx*ngds,2,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, nps
      integer mnblok, nx2, nyzp1, nxvz, nxvzs, nyzp3, nxvy, nxvys
      integer m, my, mz, j, k, n, jr, jrr, jl, jll
      dimension istatus(lstat)
      nx2 = nx + 2
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c buffer data in y
      do 40 m = 1, mnblok
      ngc = 0
      if (nyzp(1,m).eq.1) ngc = 1
      nyzp1 = nyzp(2,m) + 1
      do 30 k = 1, nyzp1
      do 20 j = 1, nxv
      do 10 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,nyzp(1,m)+1,k+1,m)
      scs(n,j,k+nyzp1,1,m) = f(n,j,2,k+1,m)
      scs(n,j,k+2*nyzp1,1,m) = f(n,j,3-ngc,k+1,m)
      scs(n,j,k+nyzp1,2,m) = f(n,j,nyzp(1,m)+2,k+1,m)
   10 continue
   20 continue
   30 continue
   40 continue
c copy to guard cells in y
      do 310 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 300 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(2,m) + 1
      nxvzs = nxv*nyzp1
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr + 1
      krr = jrr + kz
      jl = ky - 1
      jll = jl - 1
      kll = jll + kz
      kr = jr + kz
      kl = jl + kz
      ngc = 0
c special case of only one grid per processor
      if (nyzp(1,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 70 k = 1, nyzp1
c        do 60 j = 1, nxv
c        do 50 n = 1, 3
c        scs(n,j,k,2,m) = scs(n,j,k,1,kl)
c  50    continue
c  60    continue
c  70    continue
c     else
c        do 100 k = 1, nyzp1
c        do 90 j = 1, nxv
c        do 80 n = 1, 3
c        scs(n,j,k,2,m) = f(n,j,3,k+1,m)
c  80    continue
c  90    continue
c 100    continue
c     endif
c     if (jr.lt.nvpy) then
c        if (nyzp(1,kr).eq.1) then
c           do 130 k = 1, nyzp1
c           do 120 j = 1, nxv
c           do 110 n = 1, 3
c           scs(n,j,k+nyzp1,2,m) = scs(n,j,k+nyzp1,1,kr)
c           scs(n,j,k+2*nyzp1,2,m) = scs(n,j,k+nyzp1,1,krr)
c 110       continue
c 120       continue
c 130       continue
c        else
c           do 160 k = 1, nyzp1
c           do 150 j = 1, nxv
c           do 140 n = 1, 3
c           scs(n,j,k+nyzp1,2,m) = scs(n,j,k+nyzp1,1,kr)
c           scs(n,j,k+2*nyzp1,2,m) = scs(n,j,k+2*nyzp1,1,kr)
c 140       continue
c 150       continue
c 160       continue
c        endif
c     else
c        do 190 k = 1, nyzp1
c        do 180 j = 2, nx2
c        do 170 n = 1, 3
c        scs(n,j,k+2*nyzp1,2,m) = 2.*scs(n,j,k+nyzp1,2,m) - scs(n,j,k,1,
c    1m)
c 170    continue
c 180    continue
c 190    continue
c     endif
c     if (nyzp(1,m).eq.1) then
c        if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
c           do 220 k = 1, nyzp1
c           do 210 j = 1, nxv
c           do 200 n = 1, 3
c           scs(n,j,k,2,m) = scs(n,j,k+nyzp1,1,kr)
c 200       continue
c 210       continue
c 220       continue
c        endif
c        if (jr.eq.(nvpy-1)) then
c           do 230 k = 1, nyzp1
c           do 240 j = 1, nxv
c           do 250 n = 1, 3
c           scs(n,j,k+2*nyzp1,2,m) = scs(n,j,k+nyzp1,2,kr)
c 230       continue
c 240       continue
c 250       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,1,2,m),3*nxvz,mreal,kl-1,noff+3,lgrp,msi
     1d,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,1,m),3*nxvzs,mreal,kr-1,noff+3,lgrp,ier
     1r)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 70 k = 1, nyzp1
         do 60 j = 1, nxv
         do 50 n = 1, 3
         scs(n,j,k,2,m) = f(n,j,3,k+1,m)
   50    continue
   60    continue
   70    continue
      endif
      if (jr.lt.nvpy) then
         call MPI_IRECV(scs(1,1,nyzp1+1,2,m),6*nxvz,mreal,kr-1,noff+4,lg
     1rp,msid,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scs(1,1,nyzp1+1,1,m),3*(2-ngc)*nxvzs,mreal,kl-1,n
     1off+4,lgrp,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 100 k = 1, nyzp1
         do 90 j = 2, nx2
         do 80 n = 1, 3
         scs(n,j,k+2*nyzp1,2,m) = 2.*scs(n,j,k+nyzp1,2,m) - scs(n,j,k,1,
     1m)
   80    continue
   90    continue
  100    continue
      endif
c special case of only one grid per processor in y
      if (mter.ge.1) go to 260
      if (jr.lt.nvpy) call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      if (jrr.lt.nvpy) then
         if (nps.eq.(3*nxvzs)) then
            call MPI_IRECV(scs(1,1,2*nyzp1+1,2,m),3*nxvz,mreal,krr-1,nof
     1f+5,lgrp,msid,ierr)
         else
            call MPI_IRECV(scs(1,1,1,1,m),3*nxvz,mreal,krr-1,noff+5,lgrp
     1,msid,ierr)
         endif
      else if (jr.lt.nvpy) then
         if (nps.eq.(3*nxvzs)) then
            call MPI_IRECV(scs(1,1,2*nyzp1+1,2,m),3*nxvz,mreal,kr-1,noff
     1+5,lgrp,msid,ierr)
         else
            call MPI_IRECV(scs(1,1,1,1,m),3*nxvz,mreal,kr-1,noff+5,lgrp,
     1msid,ierr)
         endif
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
         if (ngc.eq.1) then
            call MPI_IRECV(scs(1,1,1,2,m),3*nxvz,mreal,kr-1,noff+5,lgrp,
     1nsid,ierr)
         else
            call MPI_IRECV(scs(1,1,1,1,m),3*nxvz,mreal,kr-1,noff+5,lgrp,
     1nsid,ierr)
         endif
      endif
      if (jll.ge.0) then
         call MPI_SEND(scs(1,1,2*nyzp1+1,1,m),3*nxvzs,mreal,kll-1,noff+5
     1,lgrp,ierr)
      else if (jl.eq.0) then
         call MPI_SEND(scs(1,1,2*nyzp1+1,1,m),3*nxvzs,mreal,kl-1,noff+5,
     1lgrp,ierr)
      endif
      if ((jl.eq.(nvpy-2)).and.(jl.ge.0)) then
         call MPI_SEND(scs(1,1,nyzp1+1,2,m),3*nxvzs,mreal,kl-1,noff+5,lg
     1rp,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
         call MPI_WAIT(nsid,istatus,ierr)
      endif
c copy guard cells
  260 do 290 k = 1, nyzp1
      do 280 j = 1, nxv
      do 270 n = 1, 3
      f(n,j,1,k+1,m) = scs(n,j,k,2,m)
      f(n,j,nyzp(1,m)+2,k+1,m) = scs(n,j,k+nyzp1,2,m)
      f(n,j,nyzp(1,m)+3,k+1,m) = scs(n,j,k+2*nyzp1,2,m)
  270 continue
  280 continue
  290 continue
  300 continue
  310 continue
c fix left edge
      do 360 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 350 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(2,m) + 1
      kl = my + js
      if (kl.eq.0) then
         do 340 k = 1, nyzp1
         do 330 j = 2, nx2
         do 320 n = 1, 3
         f(n,j,1,k+1,m) = 2.*f(n,j,2,k+1,m) - f(n,j,1,k+1,m)
  320    continue
  330    continue
  340    continue
      endif
  350 continue
  360 continue
c buffer data in z
      do 400 m = 1, mnblok
      ngc = 0
      if (nyzp(2,m).eq.1) ngc = 1
      nyzp3 = nyzp(1,m) + 3
      do 390 k = 1, nyzp3
      do 380 j = 1, nxv
      do 370 n = 1, 3
      scr(n,j,k,1,m) = f(n,j,k,nyzp(2,m)+1,m)
      scr(n,j,k+nyzp3,1,m) = f(n,j,k,2,m)
      scr(n,j,k+2*nyzp3,1,m) = f(n,j,k,3-ngc,m)
      scr(n,j,k+nyzp3,2,m) = f(n,j,k,nyzp(2,m)+2,m)
  370 continue
  380 continue
  390 continue
  400 continue
c copy to guard cells in z
      do 670 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 660 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(1,m) + 3
      nxvys = nxv*nyzp3
      ky = my + js + 1
      kz = mz + ks
      jr = kz + 1
      jrr = jr + 1
      krr = ky + nvpy*jrr
      jl = kz - 1
      jll = jl - 1
      kll = ky + nvpy*jll
      kr = ky + nvpy*jr
      kl = ky + nvpy*jl
      ngc = 0
c special case of only one grid per processor
      if (nyzp(2,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 430 k = 1, nyzp3
c        do 420 j = 1, nxv
c        do 410 n = 1, 3
c        scr(n,j,k,2,m) = scr(n,j,k,1,kl)
c 410    continue
c 420    continue
c 430    continue
c     else
c        do 460 k = 1, nyzp3
c        do 450 j = 1, nxv
c        do 440 n = 1, 3
c        scr(n,j,k,2,m) = f(n,j,k,3,m)
c 440    continue
c 450    continue
c 460    continue
c     endif
c     if (jr.lt.nvpz) then
c        if (nyzp(2,kr).eq.1) then
c           do 490 k = 1, nyzp3
c           do 480 j = 1, nxv
c           do 470 n = 1, 3
c           scr(n,j,k+nyzp3,2,m) = scr(n,j,k+nyzp3,1,kr)
c           scr(n,j,k+2*nyzp3,2,m) = scr(n,j,k+nyzp3,1,krr)
c 470       continue
c 480       continue
c 490       continue
c        else
c           do 520 k = 1, nyzp3
c           do 510 j = 1, nxv
c           do 500 n = 1, 3
c           scr(n,j,k+nyzp3,2,m) = scr(n,j,k+nyzp3,1,kr)
c           scr(n,j,k+2*nyzp3,2,m) = scr(n,j,k+2*nyzp3,1,kr)
c 500       continue
c 510       continue
c 520       continue
c        endif
c     else
c        do 550 k = 1, nyzp3
c        do 540 j = 2, nx2
c        do 530 n = 1, 3
c        scr(n,j,k+2*nyzp3,2,m) = 2.*scr(n,j,k+nyzp3,2,m) - scr(n,j,k,1,
c    1m)
c 530    continue
c 540    continue
c 550    continue
c     endif
c     if (nyzp(2,m).eq.1) then
c        if ((jl.eq.(-1)).and.(jr.lt.nvpz)) then
c           do 580 k = 1, nyzp3
c           do 570 j = 1, nxv
c           do 560 n = 1, 3
c           scr(n,j,k,2,m) = scr(n,j,k+nyzp3,1,kr)
c 560       continue
c 570       continue
c 580       continue
c        endif
c        if (jr.eq.(nvpz-1)) then
c           do 610 k = 1, nyzp3
c           do 600 j = 1, nxv
c           do 590 n = 1, 3
c           scr(n,j,k+2*nyzp3,2,m) = scr(n,j,k+nyzp3,2,kr)
c 590       continue
c 600       continue
c 610       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scr(1,1,1,2,m),3*nxvy,mreal,kl-1,noff+6,lgrp,msi
     1d,ierr)
      endif
      if (jr.lt.nvpz) then
         call MPI_SEND(scr(1,1,1,1,m),3*nxvys,mreal,kr-1,noff+6,lgrp,ier
     1r)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 430 k = 1, nyzp3
         do 420 j = 1, nxv
         do 410 n = 1, 3
         scr(n,j,k,2,m) = f(n,j,k,3,m)
  410    continue
  420    continue
  430    continue
      endif
      if (jr.lt.nvpz) then
         call MPI_IRECV(scr(1,1,nyzp3+1,2,m),6*nxvy,mreal,kr-1,noff+8,lg
     1rp,msid,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scr(1,1,nyzp3+1,1,m),3*(2-ngc)*nxvys,mreal,kl-1,n
     1off+8,lgrp,ierr)
      endif
      if (jr.lt.nvpz) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 460 k = 1, nyzp3
         do 450 j = 2, nx2
         do 440 n = 1, 3
         scr(n,j,k+2*nyzp3,2,m) = 2.*scr(n,j,k+nyzp3,2,m) - scr(n,j,k,1,
     1m)
  440    continue
  450    continue
  460    continue
      endif
c special case of only one grid per processor in z
      if (nter.ge.1) go to 620
      if (jr.lt.nvpz) call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      if (jrr.lt.nvpz) then
         if (nps.eq.(3*nxvys)) then
            call MPI_IRECV(scr(1,1,2*nyzp3+1,2,m),3*nxvy,mreal,krr-1,nof
     1f+9,lgrp,msid,ierr)
         else
            call MPI_IRECV(scr(1,1,1,1,m),3*nxvy,mreal,krr-1,noff+9,lgrp
     1,msid,ierr)
         endif
      else if (jr.lt.nvpz) then
         if (nps.eq.(3*nxvys)) then
            call MPI_IRECV(scr(1,1,2*nyzp3+1,2,m),3*nxvy,mreal,kr-1,noff
     1+9,lgrp,msid,ierr)
         else
            call MPI_IRECV(scr(1,1,1,1,m),3*nxvy,mreal,kr-1,noff+9,lgrp,
     1msid,ierr)
         endif
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpz)) then
         if (ngc.eq.1) then
            call MPI_IRECV(scr(1,1,1,2,m),3*nxvy,mreal,kr-1,noff+9,lgrp,
     1nsid,ierr)
         else
            call MPI_IRECV(scr(1,1,1,1,m),3*nxvy,mreal,kr-1,noff+9,lgrp,
     1nsid,ierr)
         endif
      endif
      if (jll.ge.0) then
         call MPI_SEND(scr(1,1,2*nyzp3+1,1,m),3*nxvys,mreal,kll-1,noff+9
     1,lgrp,ierr)
      else if (jl.eq.0) then
         call MPI_SEND(scr(1,1,2*nyzp3+1,1,m),3*nxvys,mreal,kl-1,noff+9,
     1lgrp,ierr)
      endif
      if ((jl.eq.(nvpz-2)).and.(jl.ge.0)) then
         call MPI_SEND(scr(1,1,nyzp3+1,2,m),3*nxvys,mreal,kl-1,noff+9,lg
     1rp,ierr)
      endif
      if (jr.lt.nvpz) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpz)) then
         call MPI_WAIT(nsid,istatus,ierr)
      endif
c copy guard cells
  620 do 650 k = 1, nyzp3
      do 640 j = 1, nxv
      do 630 n = 1, 3
      f(n,j,k,1,m) = scr(n,j,k,2,m)
      f(n,j,k,nyzp(2,m)+2,m) = scr(n,j,k+nyzp3,2,m)
      f(n,j,k,nyzp(2,m)+3,m) = scr(n,j,k+2*nyzp3,2,m)
  630 continue
  640 continue
  650 continue
  660 continue
  670 continue
c fix left edge
      do 720 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 710 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(1,m) + 3
      kl = mz + ks
      if (kl.eq.0) then
         do 700 k = 1, nyzp3
         do 690 j = 2, nx2
         do 680 n = 1, 3
         f(n,j,k,1,m) = 2.*f(n,j,k,2,m) - f(n,j,k,1,m)
  680    continue
  690    continue
  700    continue
      endif
  710 continue
  720 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLDGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypmx
     1,nzpmx,mblok,nblok,ngds,idds,mter,nter)
c this subroutine copies data to guard cells in non-uniform partitions
c for scalar data the field is replicated so as to disable quadratic
c interpolation within half a cell of the edges, and reduce it to linear
c interpolation in the y and z direction.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes three extra guard cells.
c scs/scr = scratch arrays for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c it is assumed that nyzp(n,m) > 0.
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c mter/nter = (0,1) = (no,yes) pass data to next processor only in y/z
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds, mter, nter
      integer nyzp
      real f, scs, scr
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx*ngds,2,mblok*nblok)
      dimension scr(nxv,nypmx*ngds,2,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, nps
      integer mnblok, nx2, nyzp1, nxvz, nxvzs, nyzp3, nxvy, nxvys
      integer m, my, mz, j, k, jr, jrr, jl, jll
      dimension istatus(lstat)
      nx2 = nx + 2
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c buffer data in y
      do 30 m = 1, mnblok
      ngc = 0
      if (nyzp(1,m).eq.1) ngc = 1
      nyzp1 = nyzp(2,m) + 1
      do 20 k = 1, nyzp1
      do 10 j = 1, nxv
      scs(j,k,1,m) = f(j,nyzp(1,m)+1,k+1,m)
      scs(j,k+nyzp1,1,m) = f(j,2,k+1,m)
      scs(j,k+2*nyzp1,1,m) = f(j,3-ngc,k+1,m)
      scs(j,k+nyzp1,2,m) = f(j,nyzp(1,m)+2,k+1,m)
   10 continue
   20 continue
   30 continue
c copy to guard cells in y
      do 220 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 210 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(2,m) + 1
      nxvzs = nxv*nyzp1
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr + 1
      krr = jrr + kz
      jl = ky - 1
      jll = jl - 1
      kll = jll + kz
      kr = jr + kz
      kl = jl + kz
      ngc = 0
c special case of only one grid per processor
      if (nyzp(1,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 50 k = 1, nyzp1
c        do 40 j = 1, nxv
c        scs(j,k,2,m) = scs(j,k,1,kl)
c  40    continue
c  50    continue
c     else
c        do 70 k = 1, nyzp1
c        do 60 j = 1, nxv
c        scs(j,k,2,m) = f(j,3,k+1,m)
c  60    continue
c  70    continue
c     endif
c     if (jr.lt.nvpy) then
c        if (nyzp(1,kr).eq.1) then
c           do 90 k = 1, nyzp1
c           do 80 j = 1, nxv
c           scs(j,k+nyzp1,2,m) = scs(j,k+nyzp1,1,kr)
c           scs(j,k+2*nyzp1,2,m) = scs(j,k+nyzp1,1,krr)
c  80       continue
c  90       continue
c        else
c           do 110 k = 1, nyzp1
c           do 100 j = 1, nxv
c           scs(j,k+nyzp1,2,m) = scs(j,k+nyzp1,1,kr)
c           scs(j,k+2*nyzp1,2,m) = scs(j,k+2*nyzp1,1,kr)
c 100       continue
c 110       continue
c        endif
c     else
c        do 130 k = 1, nyzp1
c        do 120 j = 2, nx2
c        scs(j,k+2*nyzp1,2,m) = 2.*scs(j,k+nyzp1,2,m) - scs(j,k,1,m)
c 120    continue
c 130    continue
c     endif
c     if (nyzp(1,m).eq.1) then
c        if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
c           do 150 k = 1, nyzp1
c           do 140 j = 1, nxv
c           scs(j,k,2,m) = scs(j,k+nyzp1,1,kr)
c 140       continue
c 150       continue
c        endif
c        if (jr.eq.(nvpy-1)) then
c           do 170 k = 1, nyzp1
c           do 160 j = 1, nxv
c           scs(j,k+2*nyzp1,2,m) = scs(j,k+nyzp1,2,kr)
c 160       continue
c 170       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,2,m),nxvz,mreal,kl-1,noff+3,lgrp,msid,ie
     1rr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,m),nxvzs,mreal,kr-1,noff+3,lgrp,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 50 k = 1, nyzp1
         do 40 j = 1, nxv
         scs(j,k,2,m) = f(j,3,k+1,m)
   40    continue
   50    continue
      endif
      if (jr.lt.nvpy) then
         call MPI_IRECV(scs(1,nyzp1+1,2,m),2*nxvz,mreal,kr-1,noff+4,lgrp
     1,msid,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scs(1,nyzp1+1,1,m),(2-ngc)*nxvzs,mreal,kl-1,noff+
     14,lgrp,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 70 k = 1, nyzp1
         do 60 j = 2, nx2
         scs(j,k+2*nyzp1,2,m) = 2.*scs(j,k+nyzp1,2,m) - scs(j,k,1,m)
   60    continue
   70    continue
      endif
c special case of only one grid per processor in y
      if (mter.ge.1) go to 180
      if (jr.lt.nvpy) call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      if (jrr.lt.nvpy) then
         if (nps.eq.nxvzs) then
            call MPI_IRECV(scs(1,2*nyzp1+1,2,m),nxvz,mreal,krr-1,noff+5,
     1lgrp,msid,ierr)
         else
            call MPI_IRECV(scs(1,1,1,m),nxvz,mreal,krr-1,noff+5,lgrp,msi
     1d,ierr)
         endif
      else if (jr.lt.nvpy) then
         if (nps.eq.nxvzs) then
            call MPI_IRECV(scs(1,2*nyzp1+1,2,m),nxvz,mreal,kr-1,noff+5,l
     1grp,msid,ierr)
         else
            call MPI_IRECV(scs(1,1,1,m),nxvz,mreal,kr-1,noff+5,lgrp,msid
     1,ierr)
         endif
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
         if (ngc.eq.1) then
            call MPI_IRECV(scs(1,1,2,m),nxvz,mreal,kr-1,noff+5,lgrp,nsid
     1,ierr)
         else
            call MPI_IRECV(scs(1,1,1,m),nxvz,mreal,kr-1,noff+5,lgrp,nsid
     1,ierr)
         endif
      endif
      if (jll.ge.0) then
         call MPI_SEND(scs(1,2*nyzp1+1,1,m),nxvzs,mreal,kll-1,noff+5,lgr
     1p,ierr)
      else if (jl.eq.0) then
         call MPI_SEND(scs(1,2*nyzp1+1,1,m),nxvzs,mreal,kl-1,noff+5,lgrp
     1,ierr)
      endif
      if ((jl.eq.(nvpy-2)).and.(jl.ge.0)) then
         call MPI_SEND(scs(1,nyzp1+1,2,m),nxvzs,mreal,kl-1,noff+5,lgrp,i
     1err)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
         call MPI_WAIT(nsid,istatus,ierr)
      endif
c copy guard cells
  180 do 200 k = 1, nyzp1
      do 190 j = 1, nxv
      f(j,1,k+1,m) = scs(j,k,2,m)
      f(j,nyzp(1,m)+2,k+1,m) = scs(j,k+nyzp1,2,m)
      f(j,nyzp(1,m)+3,k+1,m) = scs(j,k+2*nyzp1,2,m)
  190 continue
  200 continue
  210 continue
  220 continue
c fix left edge
      do 260 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 250 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(2,m) + 1
      kl = my + js
      if (kl.eq.0) then
         do 240 k = 1, nyzp1
         do 230 j = 2, nx2
         f(j,1,k+1,m) = 2.*f(j,2,k+1,m) - f(j,1,k+1,m)
  230    continue
  240    continue
      endif
  250 continue
  260 continue
c buffer data in z
      do 290 m = 1, mnblok
      ngc = 0
      if (nyzp(2,m).eq.1) ngc = 1
      nyzp3 = nyzp(1,m) + 3
      do 280 k = 1, nyzp3
      do 270 j = 1, nxv
      scr(j,k,1,m) = f(j,k,nyzp(2,m)+1,m)
      scr(j,k+nyzp3,1,m) = f(j,k,2,m)
      scr(j,k+2*nyzp3,1,m) = f(j,k,3-ngc,m)
      scr(j,k+nyzp3,2,m) = f(j,k,nyzp(2,m)+2,m)
  270 continue
  280 continue
  290 continue
c copy to guard cells in z
      do 480 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 470 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(1,m) + 3
      nxvys = nxv*nyzp3
      ky = my + js + 1
      kz = mz + ks
      jr = kz + 1
      jrr = jr + 1
      krr = ky + nvpy*jrr
      jl = kz - 1
      jll = jl - 1
      kll = ky + nvpy*jll
      kr = ky + nvpy*jr
      kl = ky + nvpy*jl
      ngc = 0
c special case of only one grid per processor
      if (nyzp(2,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 310 k = 1, nyzp3
c        do 300 j = 1, nxv
c        scr(j,k,2,m) = scr(j,k,1,kl)
c 300    continue
c 310    continue
c     else
c        do 330 k = 1, nyzp3
c        do 320 j = 1, nxv
c        scr(j,k,2,m) = f(j,k,3,m)
c 320    continue
c 330    continue
c     endif
c     if (jr.lt.nvpz) then
c        if (nyzp(2,kr).eq.1) then
c           do 350 k = 1, nyzp3
c           do 340 j = 1, nxv
c           scr(j,k+nyzp3,2,m) = scr(j,k+nyzp3,1,kr)
c           scr(j,k+2*nyzp3,2,m) = scr(j,k+nyzp3,1,krr)
c 340       continue
c 350       continue
c        else
c           do 370 k = 1, nyzp3
c           do 360 j = 1, nxv
c           scr(j,k+nyzp3,2,m) = scr(j,k+nyzp3,1,kr)
c           scr(j,k+2*nyzp3,2,m) = scr(j,k+2*nyzp3,1,kr)
c 360       continue
c 370       continue
c        endif
c     else
c        do 390 k = 1, nyzp3
c        do 380 j = 2, nx2
c        scr(j,k+2*nyzp3,2,m) = 2.*scr(j,k+nyzp3,2,m) - scr(j,k,1,m)
c 380    continue
c 390    continue
c     endif
c     if (nyzp(2,m).eq.1) then
c        if ((jl.eq.(-1)).and.(jr.lt.nvpz)) then
c           do 410 k = 1, nyzp3
c           do 400 j = 1, nxv
c           scr(j,k,2,m) = scr(j,k+nyzp3,1,kr)
c 400       continue
c 410       continue
c        endif
c        if (jr.eq.(nvpz-1)) then
c           do 430 k = 1, nyzp3
c           do 420 j = 1, nxv
c           scr(j,k+2*nyzp3,2,m) = scr(j,k+nyzp3,2,kr)
c 420       continue
c 430       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scr(1,1,2,m),nxvy,mreal,kl-1,noff+6,lgrp,msid,ie
     1rr)
      endif
      if (jr.lt.nvpz) then
         call MPI_SEND(scr(1,1,1,m),nxvys,mreal,kr-1,noff+6,lgrp,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 310 k = 1, nyzp3
         do 300 j = 1, nxv
         scr(j,k,2,m) = f(j,k,3,m)
  300    continue
  310    continue
      endif
      if (jr.lt.nvpz) then
         call MPI_IRECV(scr(1,nyzp3+1,2,m),2*nxvy,mreal,kr-1,noff+8,lgrp
     1,msid,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scr(1,nyzp3+1,1,m),(2-ngc)*nxvys,mreal,kl-1,noff+
     18,lgrp,ierr)
      endif
      if (jr.lt.nvpz) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 330 k = 1, nyzp3
         do 320 j = 2, nx2
         scr(j,k+2*nyzp3,2,m) = 2.*scr(j,k+nyzp3,2,m) - scr(j,k,1,m)
  320    continue
  330    continue
      endif
c special case of only one grid per processor in z
      if (nter.ge.1) go to 440
      if (jr.lt.nvpz) call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      if (jrr.lt.nvpz) then
         if (nps.eq.nxvys) then
            call MPI_IRECV(scr(1,2*nyzp3+1,2,m),nxvy,mreal,krr-1,noff+9,
     1lgrp,msid,ierr)
         else
            call MPI_IRECV(scr(1,1,1,m),nxvy,mreal,krr-1,noff+9,lgrp,msi
     1d,ierr)
         endif
      else if (jr.lt.nvpz) then
         if (nps.eq.nxvys) then
            call MPI_IRECV(scr(1,2*nyzp3+1,2,m),nxvy,mreal,kr-1,noff+9,l
     1grp,msid,ierr)
         else
            call MPI_IRECV(scr(1,1,1,m),nxvy,mreal,kr-1,noff+9,lgrp,msid
     1,ierr)
         endif
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpz)) then
         if (ngc.eq.1) then
            call MPI_IRECV(scr(1,1,2,m),nxvy,mreal,kr-1,noff+9,lgrp,nsid
     1,ierr)
         else
            call MPI_IRECV(scr(1,1,1,m),nxvy,mreal,kr-1,noff+9,lgrp,nsid
     1,ierr)
         endif
      endif
      if (jll.ge.0) then
         call MPI_SEND(scr(1,2*nyzp3+1,1,m),nxvys,mreal,kll-1,noff+9,lgr
     1p,ierr)
      else if (jl.eq.0) then
         call MPI_SEND(scr(1,2*nyzp3+1,1,m),nxvys,mreal,kl-1,noff+9,lgrp
     1,ierr)
      endif
      if ((jl.eq.(nvpz-2)).and.(jl.ge.0)) then
         call MPI_SEND(scr(1,nyzp3+1,2,m),nxvys,mreal,kl-1,noff+9,lgrp,i
     1err)
      endif
      if (jr.lt.nvpz) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpz)) then
         call MPI_WAIT(nsid,istatus,ierr)
      endif
c copy guard cells
  440 do 460 k = 1, nyzp3
      do 450 j = 1, nxv
      f(j,k,1,m) = scr(j,k,2,m)
      f(j,k,nyzp(2,m)+2,m) = scr(j,k+nyzp3,2,m)
      f(j,k,nyzp(2,m)+3,m) = scr(j,k+2*nyzp3,2,m)
  450 continue
  460 continue
  470 continue
  480 continue
c fix left edge
      do 520 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 510 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(1,m) + 3
      kl = mz + ks
      if (kl.eq.0) then
         do 500 k = 1, nyzp3
         do 490 j = 2, nx2
         f(j,k,1,m) = 2.*f(j,k,2,m) - f(j,k,1,m)
  490    continue
  500    continue
      endif
  510 continue
  520 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx
     1,mblok,nblok,ngds,idds)
c this subroutine copies data to guard cells in non-uniform partitions
c guard cell on last processor is presumed already set.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes one extra guard cell.
c scs(j,k,m) = scratch array for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c it is assumed that nyzp(n,m) > 0.
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c linear interpolation, for distributed data,
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, mblok, nblok, ngds
      integer idds
      integer nyzp
      real f, scs
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer ky, kz, js, ks, moff, noff, jr, jl, kr, kl, mnblok
      integer nxvz, nxvzs, nyzp1, nxvy, nxvys, m, my, mz, j, k
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
      do 20 k = 1, nyzp(2,m)
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
      nxvzs = nxv*nyzp(2,m)
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jl = ky - 1
      kr = jr + kz
      kl = jl + kz
c this segment is used for shared memory computers
c     if (jr.lt.nvpy) then
c        do 50 k = 1, nyzp(2,m)
c        do 40 j = 1, nxv
c        scs(j,k,2,m) = scs(j,k,1,kr)
c  40    continue
c  50    continue
c     endif
c this segment is used for mpi computers
      if (jr.lt.nvpy) then
         call MPI_IRECV(scs(1,1,2,m),nxvz,mreal,kr-1,noff+3,lgrp,msid,ie
     1rr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scs(1,1,1,m),nxvzs,mreal,kl-1,noff+3,lgrp,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
c copy guard cells
         do 70 k = 1, nyzp(2,m)
         do 60 j = 1, nxv
         f(j,nyzp(1,m)+1,k,m) = scs(j,k,2,m)
   60    continue
   70    continue
      endif
   80 continue
   90 continue
c copy to guard cells in z
      do 130 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 120 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(1,m) + 1
      nxvys = nxv*nyzp1
      ky = my + js + 1
      kz = mz + ks
      jr = kz + 1
      jl = kz - 1
      kr = ky + nvpy*jr
      kl = ky + nvpy*jl
c this segment is used for shared memory computers
c     if (jr.lt.nvpz) then
c        do 110 k = 1, nyzp1
c        do 100 j = 1, nxv
c        f(j,k,nyzp(2,m)+1,m) = f(j,k,1,kr)
c 100    continue
c 110    continue
c     endif
c this segment is used for mpi computers
      if (jr.lt.nvpz) then
         call MPI_IRECV(f(1,1,nyzp(2,m)+1,m),nxvy,mreal,kr-1,noff+4,lgrp
     1,msid,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(f(1,1,1,m),nxvys,mreal,kl-1,noff+4,lgrp,ierr)
      endif
      if (jr.lt.nvpz) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLACGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpm
     1x,mblok,nblok,kyp,kzp,ngds)
c thus subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y and z direction.
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes three extra
c guard cells.
c scs/scr = scratch arrays for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer kyp, kzp, ngds
      real f, scs, scr
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(3,nxv,nypmx,ngds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, mnblok
      integer nx3, nyp3, nzp3, nxvz, nxvy, m, my, mz, j, k, n
      integer jr, jrr, jl, jll
      dimension istatus(lstat)
      nx3 = nx + 3
      nyp3 = kyp + 3
      nzp3 = kzp + 3
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c buffer data in y
      do 40 m = 1, mnblok
      do 30 k = 1, nzpmx
      do 20 j = 1, nxv
      do 10 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,kyp+2,k,m)
      scs(n,j,k,2,m) = f(n,j,kyp+3,k,m)
      scs(n,j,k,3,m) = f(n,j,1,k,m)
   10 continue
   20 continue
   30 continue
   40 continue
c add guard cells in y
      do 300 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 290 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr
      krr = jr + kz
      jl = ky - 1
      jll = jl
      kll = jl + kz
      ngc = 2
c special case of only one grid per processor
      if (kyp.eq.1) then
         jrr = jr + 1
         krr = jrr + kz
         jll = jl - 1
         kll = jll + kz
         ngc = 1
      endif
      kr = jr + kz
      kl = jl + kz
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 70 k = 1, nzpmx
c        do 60 j = 1, nxv
c        do 50 n = 1, 3
c        scs(n,j,k,4,m) = scs(n,j,k,1,kl)
c        scs(n,j,k,5,m) = scs(n,j,k,2,kll)
c  50    continue
c  60    continue
c  70    continue
c     else
c        do 100 k = 1, nzp3
c        do 90 j = 1, nx3
c        do 80 n = 1, 3
c        scs(n,j,k,4,m) = 2.*scs(n,j,k,3,m)
c        scs(n,j,k,5,m) = -scs(n,j,k,3,m)
c  80    continue
c  90    continue
c 100    continue
c     endif
c     if (jr.lt.nvpy) then
c        do 130 k = 1, nzpmx
c        do 120 j = 1, nxv
c        do 110 n = 1, 3
c        scs(n,j,k,6,m) = scs(n,j,k,3,kr)
c 110    continue
c 120    continue
c 130    continue
c     else
c        do 160 k = 1, nzp3
c        do 150 j = 1, nx3
c        do 140 n = 1, 3
c        scs(n,j,k,6,m) = -scs(n,j,k,2,m)
c        f(n,j,kyp+2,k,m) = f(n,j,kyp+2,k,m) + 2.*scs(n,j,k,2,m)
c        f(n,j,kyp+3,k,m) = 0.
c 140    continue
c 150    continue
c 160    continue
c     endif
c     if (kyp.eq.1) then
c        if (jl.eq.0) then
c           do 190 k = 1, nzp3
c           do 180 j = 1, nx3
c           do 170 n = 1, 3
c           scs(n,j,k,4,m) = scs(n,j,k,1,kl)
c           scs(n,j,k,5,m) = -scs(n,j,k,3,kl)
c 170       continue
c 180       continue
c 190       continue
c        else if (jl.lt.0) then
c           do 220 k = 1, nzp3
c           do 210 j = 1, nx3
c           do 200 n = 1, 3
c           scs(n,j,k,5,m) = 0.
c 200       continue
c 210       continue
c 220       continue
c        endif
c last point is special with only one grid
c        if ((jl.eq.(nvpy-2)).and.(jl.ge.0)) then
c           do 250 k = 1, nzp3
c           do 240 j = 1, nx3
c           do 230 n = 1, 3
c           f(n,j,kyp+2,k,m) = f(n,j,kyp+2,k,m) + scs(n,j,k,2,kl)
c 230       continue
c 240       continue
c 250       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,1,4,m),3*ngc*nxvz,mreal,kl-1,noff+1,lgrp
     1,msid,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,1,m),3*ngc*nxvz,mreal,kr-1,noff+1,lgrp,
     1ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 70 k = 1, nzp3
         do 60 j = 1, nx3
         do 50 n = 1, 3
         scs(n,j,k,4,m) = 2.*scs(n,j,k,3,m)
         scs(n,j,k,5,m) = -scs(n,j,k,3,m)
   50    continue
   60    continue
   70    continue
      endif
      if (jr.lt.nvpy) then
         call MPI_IRECV(scs(1,1,1,6,m),3*nxvz,mreal,kr-1,noff+2,lgrp,msi
     1d,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scs(1,1,1,3,m),3*nxvz,mreal,kl-1,noff+2,lgrp,ierr
     1)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 100 k = 1, nzp3
         do 90 j = 1, nx3
         do 80 n = 1, 3
         scs(n,j,k,6,m) = -scs(n,j,k,2,m)
         f(n,j,kyp+2,k,m) = f(n,j,kyp+2,k,m) + 2.*scs(n,j,k,2,m)
         f(n,j,kyp+3,k,m) = 0.
   80    continue
   90    continue
  100    continue
      endif
c special case of only one grid per processor
      if (kyp.eq.1) then
         if (jll.ge.0) then
            call MPI_IRECV(scs(1,1,1,5,m),3*nxvz,mreal,kll-1,noff+5,lgrp
     1,msid,ierr)
         else if (jl.eq.0) then
            call MPI_IRECV(scs(1,1,1,5,m),3*nxvz,mreal,kl-1,noff+5,lgrp,
     1msid,ierr)
         endif
         if (jrr.lt.nvpy) then
            call MPI_SEND(scs(1,1,1,2,m),3*nxvz,mreal,krr-1,noff+5,lgrp,
     1ierr)
         endif
         if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
            call MPI_SEND(scs(1,1,1,3,m),3*nxvz,mreal,kr-1,noff+5,lgrp,i
     1err)
         endif
         if (jl.ge.0) then
            call MPI_WAIT(msid,istatus,ierr)
            if (jl.eq.0) then
               do 130 k = 1, nzp3
               do 120 j = 1, nx3
               do 110 n = 1, 3
               scs(n,j,k,5,m) = -scs(n,j,k,5,m)
  110          continue
  120          continue
  130          continue
            endif
         else
            do 160 k = 1, nzp3
            do 150 j = 1, nx3
            do 140 n = 1, 3
            scs(n,j,k,5,m) = 0.
  140       continue
  150       continue
  160       continue
         endif
c last point is special with only one grid
         if ((jl.eq.(nvpy-2)).and.(jl.ge.0)) then
            call MPI_IRECV(scs(1,1,1,2,m),3*nxvz,mreal,kl-1,noff+6,lgrp,
     1msid,ierr)
         endif
         if (jr.eq.(nvpy-1)) then
            call MPI_SEND(scs(1,1,1,2,m),3*nxvz,mreal,kr-1,noff+6,lgrp,i
     1err)
         endif
         if ((jl.eq.(nvpy-2)).and.(jl.ge.0)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 190 k = 1, nzp3
            do 180 j = 1, nx3
            do 170 n = 1, 3
            f(n,j,kyp+2,k,m) = f(n,j,kyp+2,k,m) + scs(n,j,k,2,m)
            f(n,j,kzp+3,k,m) = 0.
  170       continue
  180       continue
  190       continue
         endif
      endif
c add up the guard cells
      do 280 k = 1, nzp3
      do 270 j = 1, nx3
      do 260 n = 1, 3
      f(n,j,2,k,m) = f(n,j,2,k,m) + scs(n,j,k,4,m)
      f(n,j,ngc+1,k,m) = f(n,j,ngc+1,k,m) + scs(n,j,k,5,m)
      f(n,j,kyp+1,k,m) = f(n,j,kyp+1,k,m) + scs(n,j,k,6,m)
  260 continue
  270 continue
  280 continue
  290 continue
  300 continue
c zero out the left edge
      do 350 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 340 my = 1, mblok
      m = my + moff
      kl = my + js
      if (kl.eq.0) then
         do 330 k = 1, nzp3
         do 320 j = 1, nx3
         do 310 n = 1, 3
         f(n,j,1,k,m) = 0.
  310    continue
  320    continue
  330    continue
      endif
  340 continue
  350 continue
c add guard cells in z
      do 610 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 600 my = 1, mblok
      m = my + moff
      ky = my + js + 1
      kz = mz + ks
      jr = kz + 1
      jrr = jr
      krr = ky + nvpy*jr
      jl = kz - 1
      jll = jl
      kll = ky + nvpy*jl
      ngc = 2
c special case of only one grid per processor
      if (kzp.eq.1) then
         jrr = jr + 1
         krr = ky + nvpy*jrr
         jll = jl - 1
         kll = ky + nvpy*jll
         ngc = 1
      endif
      kr = ky + nvpy*jr
      kl = ky + nvpy*jl
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 380 k = 1, nypmx
c        do 370 j = 1, nxv
c        do 360 n = 1, 3
c        scr(n,j,k,1,m) = f(n,j,k,kzp+2,kl)
c        scr(n,j,k,2,m) = f(n,j,k,kzp+3,kll)
c 360    continue
c 370    continue
c 380    continue
c     else
c        do 410 k = 1, nyp3
c        do 400 j = 1, nx3
c        do 390 n = 1, 3
c        scr(n,j,k,1,m) = 2.*f(n,j,k,1,m)
c        scr(n,j,k,2,m) = -f(n,j,k,1,m)
c 390    continue
c 400    continue
c 410    continue
c     endif
c     if (jr.lt.nvpz) then
c        do 440 k = 1, nypmx
c        do 430 j = 1, nxv
c        do 420 n = 1, 3
c        scr(n,j,k,3,m) = f(n,j,k,1,kr)
c 420    continue
c 430    continue
c 440    continue
c     else
c        do 470 k = 1, nyp3
c        do 460 j = 1, nx3
c        do 450 n = 1, 3
c        scr(n,j,k,3,m) = -f(n,j,k,kzp+3,m)
c        f(n,j,k,kzp+2,m) = f(n,j,k,kzp+2,m) + 2.*f(n,j,k,kzp+3,m)
c        f(n,j,k,kzp+3,m) = 0.
c 450    continue
c 460    continue
c 470    continue
c     endif
c     if (kzp.eq.1) then
c        if (jl.eq.0) then
c           do 500 k = 1, nyp3
c           do 490 j = 1, nx3
c           do 480 n = 1, 3
c           scr(n,j,k,1,m) = f(n,j,k,kzp+2,kl)
c           scr(n,j,k,2,m) = -f(n,j,k,1,kl)
c 480       continue
c 490       continue
c 500       continue
c        else if (jl.lt.0) then
c           do 530 k = 1, nyp3
c           do 520 j = 1, nx3
c           do 510 n = 1, 3
c           scr(n,j,k,2,m) = 0.
c 510       continue
c 520       continue
c 530       continue
c        endif
c last point is special with only one grid
c        if ((jl.eq.(nvpz-2)).and.(jl.ge.0)) then
c           do 560 k = 1, nyp3
c           do 550 j = 1, nx3
c           do 540 n = 1, 3
c           f(n,j,k,kzp+2,m) = f(n,j,k,kzp+2,m) + f(n,j,k,kzp+3,kl)
c 540       continue
c 550       continue
c 560       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scr,3*ngc*nxvy,mreal,kl-1,noff+3,lgrp,msid,ierr)
      endif
      if (jr.lt.nvpz) then
         call MPI_SEND(f(1,1,1,kzp+2,m),3*ngc*nxvy,mreal,kr-1,noff+3,lgr
     1p,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 380 k = 1, nyp3
         do 370 j = 1, nx3
         do 360 n = 1, 3
         scr(n,j,k,1,m) = 2.*f(n,j,k,1,m)
         scr(n,j,k,2,m) = -f(n,j,k,1,m)
  360    continue
  370    continue
  380    continue
      endif
      if (jr.lt.nvpz) then
         call MPI_IRECV(scr(1,1,1,3,m),3*nxvy,mreal,kr-1,noff+4,lgrp,msi
     1d,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(f(1,1,1,1,m),3*nxvy,mreal,kl-1,noff+4,lgrp,ierr)
      endif
      if (jr.lt.nvpz) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 410 k = 1, nyp3
         do 400 j = 1, nx3
         do 390 n = 1, 3
         scr(n,j,k,3,m) = -f(n,j,k,kzp+3,m)
         f(n,j,k,kzp+2,m) = f(n,j,k,kzp+2,m) + 2.*f(n,j,k,kzp+3,m)
         f(n,j,k,kzp+3,m) = 0.
  390    continue
  400    continue
  410    continue
      endif
c special case of only one grid per processor
      if (kzp.eq.1) then
         if (jll.ge.0) then
            call MPI_IRECV(scr(1,1,1,2,m),3*nxvy,mreal,kll-1,noff+7,lgrp
     1,msid,ierr)
         else if (jl.eq.0) then
            call MPI_IRECV(scr(1,1,1,2,m),3*nxvy,mreal,kl-1,noff+7,lgrp,
     1msid,ierr)
         endif
         if (jrr.lt.nvpz) then
            call MPI_SEND(f(1,1,1,kzp+3,m),3*nxvy,mreal,krr-1,noff+7,lgr
     1p,ierr)
         endif
         if ((jl.eq.(-1)).and.(jr.lt.nvpz)) then
            call MPI_SEND(f(1,1,1,1,m),3*nxvy,mreal,kr-1,noff+7,lgrp,ier
     1r)
         endif
         if (jl.ge.0) then
            call MPI_WAIT(msid,istatus,ierr)
            if (jl.eq.0) then
               do 440 k = 1, nyp3
               do 430 j = 1, nx3
               do 420 n = 1, 3
               scr(n,j,k,2,m) = -scr(n,j,k,2,m)
  420          continue
  430          continue
  440          continue
            endif
         else
            do 470 k = 1, nyp3
            do 460 j = 1, nx3
            do 450 n = 1, 3
            scr(n,j,k,2,m) = 0.
  450       continue
  460       continue
  470       continue
         endif
c last point is special with only one grid
         if ((jl.eq.(nvpz-2)).and.(jl.ge.0)) then
            call MPI_IRECV(f(1,1,1,kzp+3,m),3*nxvy,mreal,kl-1,noff+8,lgr
     1p,msid,ierr)
         endif
         if (jr.eq.(nvpz-1)) then
            call MPI_SEND(f(1,1,1,kzp+3,m),3*nxvy,mreal,kr-1,noff+8,lgrp
     1,ierr)
         endif
         if ((jl.eq.(nvpz-2)).and.(jl.ge.0)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 500 k = 1, nyp3
            do 490 j = 1, nx3
            do 480 n = 1, 3
            f(n,j,k,kzp+2,m) = f(n,j,k,kzp+2,m) + f(n,j,k,kzp+3,m)
            f(n,j,k,kzp+3,m) = 0.
  480       continue
  490       continue
  500       continue
         endif
      endif
c add up the guard cells
      do 590 k = 1, nyp3
      do 580 j = 1, nx3
      do 570 n = 1, 3
      f(n,j,k,2,m) = f(n,j,k,2,m) + scr(n,j,k,1,m)
      f(n,j,k,ngc+1,m) = f(n,j,k,ngc+1,m) + scr(n,j,k,2,m)
      f(n,j,k,kzp+1,m) = f(n,j,k,kzp+1,m) + scr(n,j,k,3,m)
  570 continue
  580 continue
  590 continue
  600 continue
  610 continue
c zero out the front edge
      do 660 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 650 my = 1, mblok
      m = my + moff
      kl = mz + ks
      if (kl.eq.0) then
         do 640 k = 1, nyp3
         do 630 j = 1, nx3
         do 620 n = 1, 3
         f(n,j,k,1,m) = 0.
  620    continue
  630    continue
  640    continue
      endif
  650 continue
  660 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLAGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx
     1,mblok,nblok,kyp,kzp,ngds)
c thus subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y and z direction.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes three extra
c guard cells.
c scs/scr = scratch arrays for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer kyp, kzp, ngds
      real f, scs, scr
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(nxv,nypmx,ngds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, mnblok
      integer nx3, nyp3, nzp3, nxvz, nxvy, m, my, mz, j, k
      integer jr, jrr, jl, jll
      dimension istatus(lstat)
      nx3 = nx + 3
      nyp3 = kyp + 3
      nzp3 = kzp + 3
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
      scs(j,k,1,m) = f(j,kyp+2,k,m)
      scs(j,k,2,m) = f(j,kyp+3,k,m)
      scs(j,k,3,m) = f(j,1,k,m)
   10 continue
   20 continue
   30 continue
c add guard cells in y
      do 210 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 200 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr
      krr = jr + kz
      jl = ky - 1
      jll = jl
      kll = jl + kz
      ngc = 2
c special case of only one grid per processor
      if (kyp.eq.1) then
         jrr = jr + 1
         krr = jrr + kz
         jll = jl - 1
         kll = jll + kz
         ngc = 1
      endif
      kr = jr + kz
      kl = jl + kz
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 50 k = 1, nzpmx
c        do 40 j = 1, nxv
c        scs(j,k,4,m) = scs(j,k,1,kl)
c        scs(j,k,5,m) = scs(j,k,2,kll)
c  40    continue
c  50    continue
c     else
c        do 70 k = 1, nzp3
c        do 60 j = 1, nx3
c        scs(j,k,4,m) = 2.*scs(j,k,3,m)
c        scs(j,k,5,m) = -scs(j,k,3,m)
c  60    continue
c  70    continue
c     endif
c     if (jr.lt.nvpy) then
c        do 90 k = 1, nzpmx
c        do 80 j = 1, nxv
c        scs(j,k,6,m) = scs(j,k,3,kr)
c  80    continue
c  90    continue
c     else
c        do 110 k = 1, nzp3
c        do 100 j = 1, nx3
c        scs(j,k,6,m) = -scs(j,k,2,m)
c        f(j,kyp+2,k,m) = f(j,kyp+2,k,m) + 2.*scs(j,k,2,m)
c        f(j,kyp+3,k,m) = 0.
c 100    continue
c 110    continue
c     endif
c     if (kyp.eq.1) then
c        if (jl.eq.0) then
c           do 130 k = 1, nzp3
c           do 120 j = 1, nx3
c           scs(j,k,4,m) = scs(j,k,1,kl)
c           scs(j,k,5,m) = -scs(j,k,3,kl)
c 120       continue
c 130       continue
c        else if (jl.lt.0) then
c           do 150 k = 1, nzp3
c           do 140 j = 1, nx3
c           scs(j,k,5,m) = 0.
c 140       continue
c 150       continue
c        endif
c last point is special with only one grid
c        if ((jl.eq.(nvpy-2)).and.(jl.ge.0)) then
c           do 170 k = 1, nzp3
c           do 160 j = 1, nx3
c           f(j,kyp+2,k,m) = f(j,kyp+2,k,m) + scs(j,k,2,kl)
c 160       continue
c 170       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,4,m),ngc*nxvz,mreal,kl-1,noff+1,lgrp,msi
     1d,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,m),ngc*nxvz,mreal,kr-1,noff+1,lgrp,ierr
     1)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 50 k = 1, nzp3
         do 40 j = 1, nx3
         scs(j,k,4,m) = 2.*scs(j,k,3,m)
         scs(j,k,5,m) = -scs(j,k,3,m)
   40    continue
   50    continue
      endif
      if (jr.lt.nvpy) then
         call MPI_IRECV(scs(1,1,6,m),nxvz,mreal,kr-1,noff+2,lgrp,msid,ie
     1rr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scs(1,1,3,m),nxvz,mreal,kl-1,noff+2,lgrp,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 70 k = 1, nzp3
         do 60 j = 1, nx3
         scs(j,k,6,m) = -scs(j,k,2,m)
         f(j,kyp+2,k,m) = f(j,kyp+2,k,m) + 2.*scs(j,k,2,m)
         f(j,kyp+3,k,m) = 0.
   60    continue
   70    continue
      endif
c special case of only one grid per processor
      if (kyp.eq.1) then
         if (jll.ge.0) then
            call MPI_IRECV(scs(1,1,5,m),nxvz,mreal,kll-1,noff+5,lgrp,msi
     1d,ierr)
         else if (jl.eq.0) then
            call MPI_IRECV(scs(1,1,5,m),nxvz,mreal,kl-1,noff+5,lgrp,msid
     1,ierr)
         endif
         if (jrr.lt.nvpy) then
            call MPI_SEND(scs(1,1,2,m),nxvz,mreal,krr-1,noff+5,lgrp,ierr
     1)
         endif
         if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
            call MPI_SEND(scs(1,1,3,m),nxvz,mreal,kr-1,noff+5,lgrp,ierr)
         endif
         if (jl.ge.0) then
            call MPI_WAIT(msid,istatus,ierr)
            if (jl.eq.0) then
               do 90 k = 1, nzp3
               do 80 j = 1, nx3
               scs(j,k,5,m) = -scs(j,k,5,m)
   80          continue
   90          continue
            endif
         else
            do 110 k = 1, nzp3
            do 100 j = 1, nx3
            scs(j,k,5,m) = 0.
  100       continue
  110       continue
         endif
c last point is special with only one grid
         if ((jl.eq.(nvpy-2)).and.(jl.ge.0)) then
            call MPI_IRECV(scs(1,1,2,m),nxvz,mreal,kl-1,noff+6,lgrp,msid
     1,ierr)
         endif
         if (jr.eq.(nvpy-1)) then
            call MPI_SEND(scs(1,1,2,m),nxvz,mreal,kr-1,noff+6,lgrp,ierr)
         endif
         if ((jl.eq.(nvpy-2)).and.(jl.ge.0)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 130 k = 1, nzp3
            do 120 j = 1, nx3
            f(j,kyp+2,k,m) = f(j,kyp+2,k,m) + scs(j,k,2,m)
            f(j,kzp+3,k,m) = 0.
  120       continue
  130       continue
         endif
      endif
c add up the guard cells
      do 190 k = 1, nzp3
      do 180 j = 1, nx3
      f(j,2,k,m) = f(j,2,k,m) + scs(j,k,4,m)
      f(j,ngc+1,k,m) = f(j,ngc+1,k,m) + scs(j,k,5,m)
      f(j,kyp+1,k,m) = f(j,kyp+1,k,m) + scs(j,k,6,m)
  180 continue
  190 continue
  200 continue
  210 continue
c zero out the left edge
      do 250 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 240 my = 1, mblok
      m = my + moff
      kl = my + js
      if (kl.eq.0) then
         do 230 k = 1, nzp3
         do 220 j = 1, nx3
         f(j,1,k,m) = 0.
  220    continue
  230    continue
      endif
  240 continue
  250 continue
c add guard cells in z
      do 430 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 420 my = 1, mblok
      m = my + moff
      ky = my + js + 1
      kz = mz + ks
      jr = kz + 1
      jrr = jr
      krr = ky + nvpy*jr
      jl = kz - 1
      jll = jl
      kll = ky + nvpy*jl
      ngc = 2
c special case of only one grid per processor
      if (kzp.eq.1) then
         jrr = jr + 1
         krr = ky + nvpy*jrr
         jll = jl - 1
         kll = ky + nvpy*jll
         ngc = 1
      endif
      kr = ky + nvpy*jr
      kl = ky + nvpy*jl
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 270 k = 1, nypmx
c        do 260 j = 1, nxv
c        scr(j,k,1,m) = f(j,k,kzp+2,kl)
c        scr(j,k,2,m) = f(j,k,kzp+3,kll)
c 260    continue
c 270    continue
c     else
c        do 290 k = 1, nyp3
c        do 280 j = 1, nx3
c        scr(j,k,1,m) = 2.*f(j,k,1,m)
c        scr(j,k,2,m) = -f(j,k,1,m)
c 280    continue
c 290    continue
c     endif
c     if (jr.lt.nvpz) then
c        do 310 k = 1, nypmx
c        do 300 j = 1, nxv
c        scr(j,k,3,m) = f(j,k,1,kr)
c 300    continue
c 310    continue
c     else
c        do 330 k = 1, nyp3
c        do 320 j = 1, nx3
c        scr(j,k,3,m) = -f(j,k,kzp+3,m)
c        f(j,k,kzp+2,m) = f(j,k,kzp+2,m) + 2.*f(j,k,kzp+3,m)
c        f(j,k,kzp+3,m) = 0.
c 320    continue
c 330    continue
c     endif
c     if (kzp.eq.1) then
c        if (jl.eq.0) then
c           do 350 k = 1, nyp3
c           do 340 j = 1, nx3
c           scr(j,k,1,m) = f(j,k,kzp+2,kl)
c           scr(j,k,2,m) = -f(j,k,1,kl)
c 340       continue
c 350       continue
c        else if (jl.lt.0) then
c           do 370 k = 1, nyp3
c           do 360 j = 1, nx3
c           scr(j,k,2,m) = 0.
c 360       continue
c 370       continue
c        endif
c last point is special with only one grid
c        if ((jl.eq.(nvpz-2)).and.(jl.ge.0)) then
c           do 390 k = 1, nyp3
c           do 380 j = 1, nx3
c           f(j,k,kzp+2,m) = f(j,k,kzp+2,m) + f(j,k,kzp+3,kl)
c 380       continue
c 390       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scr,ngc*nxvy,mreal,kl-1,noff+3,lgrp,msid,ierr)
      endif
      if (jr.lt.nvpz) then
         call MPI_SEND(f(1,1,kzp+2,m),ngc*nxvy,mreal,kr-1,noff+3,lgrp,ie
     1rr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 270 k = 1, nyp3
         do 260 j = 1, nx3
         scr(j,k,1,m) = 2.*f(j,k,1,m)
         scr(j,k,2,m) = -f(j,k,1,m)
  260    continue
  270    continue
      endif
      if (jr.lt.nvpz) then
         call MPI_IRECV(scr(1,1,3,m),nxvy,mreal,kr-1,noff+4,lgrp,msid,ie
     1rr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(f(1,1,1,m),nxvy,mreal,kl-1,noff+4,lgrp,ierr)
      endif
      if (jr.lt.nvpz) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 290 k = 1, nyp3
         do 280 j = 1, nx3
         scr(j,k,3,m) = -f(j,k,kzp+3,m)
         f(j,k,kzp+2,m) = f(j,k,kzp+2,m) + 2.*f(j,k,kzp+3,m)
         f(j,k,kzp+3,m) = 0.
  280    continue
  290    continue
      endif
c special case of only one grid per processor
      if (kzp.eq.1) then
         if (jll.ge.0) then
            call MPI_IRECV(scr(1,1,2,m),nxvy,mreal,kll-1,noff+7,lgrp,msi
     1d,ierr)
         else if (jl.eq.0) then
            call MPI_IRECV(scr(1,1,2,m),nxvy,mreal,kl-1,noff+7,lgrp,msid
     1,ierr)
         endif
         if (jrr.lt.nvpz) then
            call MPI_SEND(f(1,1,kzp+3,m),nxvy,mreal,krr-1,noff+7,lgrp,ie
     1rr)
         endif
         if ((jl.eq.(-1)).and.(jr.lt.nvpz)) then
            call MPI_SEND(f(1,1,1,m),nxvy,mreal,kr-1,noff+7,lgrp,ierr)
         endif
         if (jl.ge.0) then
            call MPI_WAIT(msid,istatus,ierr)
            if (jl.eq.0) then
               do 310 k = 1, nyp3
               do 300 j = 1, nx3
               scr(j,k,2,m) = -scr(j,k,2,m)
  300          continue
  310          continue
            endif
         else
            do 330 k = 1, nyp3
            do 320 j = 1, nx3
            scr(j,k,2,m) = 0.
  320       continue
  330       continue
         endif
c last point is special with only one grid
         if ((jl.eq.(nvpz-2)).and.(jl.ge.0)) then
            call MPI_IRECV(f(1,1,kzp+3,m),nxvy,mreal,kl-1,noff+8,lgrp,ms
     1id,ierr)
         endif
         if (jr.eq.(nvpz-1)) then
            call MPI_SEND(f(1,1,kzp+3,m),nxvy,mreal,kr-1,noff+8,lgrp,ier
     1r)
         endif
         if ((jl.eq.(nvpz-2)).and.(jl.ge.0)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 350 k = 1, nyp3
            do 340 j = 1, nx3
            f(j,k,kzp+2,m) = f(j,k,kzp+2,m) + f(j,k,kzp+3,m)
            f(j,k,kzp+3,m) = 0.
  340       continue
  350       continue
         endif
      endif
c add up the guard cells
      do 410 k = 1, nyp3
      do 400 j = 1, nx3
      f(j,k,2,m) = f(j,k,2,m) + scr(j,k,1,m)
      f(j,k,ngc+1,m) = f(j,k,ngc+1,m) + scr(j,k,2,m)
      f(j,k,kzp+1,m) = f(j,k,kzp+1,m) + scr(j,k,3,m)
  400 continue
  410 continue
  420 continue
  430 continue
c zero out the front edge
      do 470 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 460 my = 1, mblok
      m = my + moff
      kl = mz + ks
      if (kl.eq.0) then
         do 450 k = 1, nyp3
         do 440 j = 1, nx3
         f(j,k,1,m) = 0.
  440    continue
  450    continue
      endif
  460 continue
  470 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLACGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypm
     1x,nzpmx,mblok,nblok,ngds,idds,mter,nter)
c this subroutine adds data from guard cells in non-uniform partitions
c for vector data, the field is added up so as to disable quadratic
c interpolation within half a cell of the edges, and reduce it to linear
c interpolation in the y and z direction.
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes three extra guard cells.
c scs/scr = scratch arrays for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c it is assumed the nyzp(n,m) > 0.
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c mter/nter = (0,1) = (no,yes) pass data to next processor only in y/z
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds, mter, nter
      integer nyzp
      real f, scs, scr
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx*ngds,2,mblok*nblok)
      dimension scr(3,nxv,nypmx*ngds,2,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, mnblok
      integer nx3, nxvz, nxvzs, nyzp3, nxvy, nxvys, m, my, mz, j, k, n
      integer jr, jrr, jl, jll
      dimension istatus(lstat)
      nx3 = nx + 3
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c buffer data in y
      do 40 m = 1, mnblok
      nyzp3 = nyzp(2,m) + 3
      do 30 k = 1, nyzp3
      do 20 j = 1, nxv
      do 10 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,nyzp(1,m)+2,k,m)
      scs(n,j,k+nyzp3,1,m) = f(n,j,nyzp(1,m)+3,k,m)
      scs(n,j,k+2*nyzp3,1,m) = f(n,j,1,k,m)
   10 continue
   20 continue
   30 continue
   40 continue
c add guard cells in y
      do 340 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 330 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(2,m) + 3
      nxvzs = nxv*nyzp3
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr + 1
      krr = jrr + kz
      jl = ky - 1
      jll = jl - 1
      kll = jll + kz
      kr = jr + kz
      kl = jl + kz
      ngc = 0
c special case of only one grid per processor
      if (nyzp(1,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        if (nyzp(1,kl).eq.1) then
c           do 70 k = 1, nyzp3
c           do 60 j = 1, nxv
c           do 50 n = 1, 3
c           scs(n,j,k,2,m) = scs(n,j,k,1,kl) + scs(n,j,k,1,kll)
c           scs(n,j,k+nyzp3,2,m) = scs(n,j,k+nyzp3,1,kl)
c  50       continue
c  60       continue
c  70       continue
c        else
c           do 100 k = 1, nyzp3
c           do 90 j = 1, nxv
c           do 80 n = 1, 3
c           scs(n,j,k,2,m) = scs(n,j,k,1,kl)
c           scs(n,j,k+nyzp3,2,m) = scs(n,j,k+nyzp3,1,kl)
c  80       continue
c  90       continue
c 100       continue
c        endif
c     else
c        do 130 k = 1, nyzp3
c        do 120 j = 1, nx3
c        do 110 n = 1, 3
c        scs(n,j,k,2,m) = 2.*scs(n,j,k+2*nyzp3,1,m)
c        scs(n,j,k+nyzp3,2,m) = -scs(n,j,k+2*nyzp3,1,m)
c 110    continue
c 120    continue
c 130    continue
c     endif
c     if (jr.lt.nvpy) then
c        do 160 k = 1, nyzp3
c        do 150 j = 1, nxv
c        do 140 n = 1, 3
c        scs(n,j,k+2*nyzp3,2,m) = scs(n,j,k+2*nyzp3,1,kr)
c 140    continue
c 150    continue
c 160    continue
c     else
c        do 190 k = 1, nyzp3
c        do 180 j = 1, nx3
c        do 170 n = 1, 3
c        scs(n,j,k+2*nyzp3,2,m) = -scs(n,j,k+nyzp3,1,m)
c        f(n,j,nyzp(1,m)+2,k,m) = f(n,j,nyzp(1,m)+2,k,m) + 2.*scs(n,j,k+
c    1nyzp3,1,m)
c        f(n,j,nyzp(1,m)+3,k,m) = 0.
c 170    continue
c 180    continue
c 190    continue
c     endif
c     if (jl.ge.0) then
c        if ((nyzp(1,kl).eq.1).and.(jl.eq.0)) then
c           do 220 k = 1, nyzp3
c           do 210 j = 1, nx3
c           do 200 n = 1, 3
c           scs(n,j,k,2,m) = scs(n,j,k,2,m) - scs(n,j,k+2*nyzp3,1,kl)
c 200       continue
c 210       continue
c 220       continue
c        endif
c     else if (nyzp(1,m).eq.1) then
c        do 250 k = 1, nyzp3
c        do 240 j = 1, nx3
c        do 230 n = 1, 3
c        scs(n,j,k+nyzp3,2,m) = 0.
c 230    continue
c 240    continue
c 250    continue
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,1,2,m),6*nxvz,mreal,kl-1,noff+1,lgrp,msi
     1d,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,1,m),6*nxvzs,mreal,kr-1,noff+1,lgrp,ier
     1r)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 70 k = 1, nyzp3
         do 60 j = 1, nx3
         do 50 n = 1, 3
         scs(n,j,k,2,m) = 2.*scs(n,j,k+2*nyzp3,1,m)
         scs(n,j,k+nyzp3,2,m) = -scs(n,j,k+2*nyzp3,1,m)
   50    continue
   60    continue
   70    continue
      endif
      if (jr.lt.nvpy) then
         call MPI_IRECV(scs(1,1,2*nyzp3+1,2,m),3*nxvz,mreal,kr-1,noff+2,
     1lgrp,msid,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scs(1,1,2*nyzp3+1,1,m),3*nxvzs,mreal,kl-1,noff+2,
     1lgrp,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 100 k = 1, nyzp3
         do 90 j = 1, nx3
         do 80 n = 1, 3
         scs(n,j,k+2*nyzp3,2,m) = -scs(n,j,k+nyzp3,1,m)
         f(n,j,nyzp(1,m)+2,k,m) = f(n,j,nyzp(1,m)+2,k,m) + 2.*scs(n,j,k+
     1nyzp3,1,m)
         f(n,j,nyzp(1,m)+3,k,m) = 0.
   80    continue
   90    continue
  100    continue
      endif
c special case of only one grid per processor in y
      if (mter.ge.1) go to 260
      if (jl.ge.0) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+3,lgrp,msid,ierr)
      endif
      if (jll.ge.0) then
         call MPI_IRECV(scs(1,1,1,1,m),3*nxvz,mreal,kll-1,noff+5,lgrp,ns
     1id,ierr)
      else if (jl.eq.0) then
         call MPI_IRECV(scs(1,1,1,1,m),3*nxvz,mreal,kl-1,noff+5,lgrp,nsi
     1d,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(nyzp(1,m),1,mint,kr-1,moff+3,lgrp,ierr)
      endif
      if (jrr.lt.nvpy) then
         call MPI_SEND(scs(1,1,nyzp3+1,1,m),3*nxvzs,mreal,krr-1,noff+5,l
     1grp,ierr)
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
         call MPI_SEND(scs(1,1,2*nyzp3+1,1,m),3*nxvzs,mreal,kr-1,noff+5,
     1lgrp,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         if (ngc.eq.1) then
            if (jl.eq.0) then
               do 130 k = 1, nyzp3
               do 120 j = 1, nx3
               do 110 n = 1, 3
               scs(n,j,k,2,m) = scs(n,j,k,2,m) - scs(n,j,k,1,m)
  110          continue
  120          continue
  130          continue
            else
               do 160 k = 1, nyzp3
               do 150 j = 1, nx3
               do 140 n = 1, 3
               scs(n,j,k,2,m) = scs(n,j,k,2,m) + scs(n,j,k,1,m)
  140          continue
  150          continue
  160          continue
            endif
         else if (nyzp(1,m).eq.1) then
            do 190 k = 1, nyzp3
            do 180 j = 1, nx3
            do 170 n = 1, 3
            scs(n,j,k+nyzp3,2,m) = 0.
  170       continue
  180       continue
  190       continue
         endif
      endif
c add up the guard cells
  260 do 290 k = 1, nyzp3
      do 280 j = 1, nx3
      do 270 n = 1, 3
      f(n,j,2,k,m) = f(n,j,2,k,m) + scs(n,j,k,2,m)
      f(n,j,3,k,m) = f(n,j,3,k,m) + scs(n,j,k+nyzp3,2,m)
      f(n,j,nyzp(1,m)+1,k,m) = f(n,j,nyzp(1,m)+1,k,m) + scs(n,j,k+2*nyzp
     13,2,m)
      f(n,j,1,k,m) = 0.
      f(n,j,nyzp(1,m)+3,k,m) = 0.
  270 continue
  280 continue
  290 continue
      if (jr.lt.nvpy) then
         do 320 k = 1, nyzp3
         do 310 j = 1, nx3
         do 300 n = 1, 3
         f(n,j,nyzp(1,m)+2,k,m) = 0.
  300    continue
  310    continue
  320    continue
      endif
  330 continue
  340 continue
c zero out the left edge
      do 390 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 380 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(2,m) + 3
      kl = my + js
      if (kl.eq.0) then
         do 370 k = 1, nyzp3
         do 360 j = 1, nx3
         do 350 n = 1, 3
         f(n,j,1,k,m) = 0.
  350    continue
  360    continue
  370    continue
      endif
  380 continue
  390 continue
c buffer data in z
      do 430 m = 1, mnblok
      nyzp3 = nyzp(1,m) + 3
      do 420 k = 1, nyzp3
      do 410 j = 1, nxv
      do 400 n = 1, 3
      scr(n,j,k,1,m) = f(n,j,k,nyzp(2,m)+2,m)
      scr(n,j,k+nyzp3,1,m) = f(n,j,k,nyzp(2,m)+3,m)
      scr(n,j,k+2*nyzp3,1,m) = f(n,j,k,1,m)
  400 continue
  410 continue
  420 continue
  430 continue
c add guard cells in z
      do 730 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 720 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(1,m) + 3
      nxvys = nxv*nyzp3
      ky = my + js + 1
      kz = mz + ks
      jr = kz + 1
      jrr = jr + 1
      krr = ky + nvpy*jrr
      jl = kz - 1
      jll = jl - 1
      kll = ky + nvpy*jll
      kr = ky + nvpy*jr
      kl = ky + nvpy*jl
      ngc = 0
c special case of only one grid per processor
      if (nyzp(2,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        if (nyzp(2,kl).eq.1) then
c           do 460 k = 1, nyzp3
c           do 450 j = 1, nxv
c           do 440 n = 1, 3
c           scr(n,j,k,2,m) = scr(n,j,k,1,kl) + scr(n,j,k,1,kll)
c           scr(n,j,k+nyzp3,2,m) = scr(n,j,k+nyzp3,1,kl)
c 440       continue
c 450       continue
c 460       continue
c       else
c           do 490 k = 1, nyzp3
c           do 480 j = 1, nxv
c           do 470 n = 1, 3
c           scr(n,j,k,2,m) = scr(n,j,k,1,kl)
c           scr(n,j,k+nyzp3,2,m) = scr(n,j,k+nyzp3,1,kl)
c 470       continue
c 480       continue
c 490       continue
c        endif
c     else
c        do 520 k = 1, nyzp3
c        do 510 j = 1, nx3
c        do 500 n = 1, 3
c        scr(n,j,k,2,m) = 2.*scr(n,j,k+2*nyzp3,1,m)
c        scr(n,j,k+nyzp3,2,m) = -scr(n,j,k+2*nyzp3,1,m)
c 500    continue
c 510    continue
c 520    continue
c     endif
c     if (jr.lt.nvpz) then
c        do 550 k = 1, nyzp3
c        do 540 j = 1, nxv
c        do 530 n = 1, 3
c        scr(n,j,k+2*nyzp3,2,m) = scr(n,j,k+2*nyzp3,1,kr)
c 530    continue
c 540    continue
c 550    continue
c     else
c        do 580 k = 1, nyzp3
c        do 570 j = 1, nx3
c        do 560 n = 1, 3
c        scr(n,j,k+2*nyzp3,2,m) = -scr(n,j,k+nyzp3,1,m)
c        f(n,j,k,nyzp(2,m)+2,m) = f(n,j,k,nyzp(2,m)+2,m) + 2.*scr(n,j,k+
c    1nyzp3,1,m)
c        f(n,j,k,nyzp(2,m)+3,m) = 0.
c 560    continue
c 570    continue
c 580    continue
c     endif
c     if (jl.ge.0) then
c        if ((nyzp(2,kl).eq.1).and.(jl.eq.0)) then
c           do 610 k = 1, nyzp3
c           do 600 j = 1, nx3
c           do 590 n = 1, 3
c           scr(n,j,k,2,m) = scr(n,j,k,2,m) - scr(n,j,k+2*nyzp3,1,kl)
c 590       continue
c 600       continue
c 610       continue
c        endif
c     else if (nyzp(2,m).eq.1) then
c        do 640 k = 1, nyzp3
c        do 630 j = 1, nx3
c        do 620 n = 1, 3
c        scr(n,j,k+nyzp3,2,m) = 0.
c 620    continue
c 630    continue
c 640    continue
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scr(1,1,1,2,m),6*nxvy,mreal,kl-1,noff+3,lgrp,msi
     1d,ierr)
      endif
      if (jr.lt.nvpz) then
         call MPI_SEND(scr(1,1,1,1,m),6*nxvys,mreal,kr-1,noff+3,lgrp,ier
     1r)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 460 k = 1, nyzp3
         do 450 j = 1, nx3
         do 440 n = 1, 3
         scr(n,j,k,2,m) = 2.*scr(n,j,k+2*nyzp3,1,m)
         scr(n,j,k+nyzp3,2,m) = -scr(n,j,k+2*nyzp3,1,m)
  440    continue
  450    continue
  460    continue
      endif
      if (jr.lt.nvpz) then
         call MPI_IRECV(scr(1,1,2*nyzp3+1,2,m),3*nxvy,mreal,kr-1,noff+4,
     1lgrp,msid,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scr(1,1,2*nyzp3+1,1,m),3*nxvys,mreal,kl-1,noff+4,
     1lgrp,ierr)
      endif
      if (jr.lt.nvpz) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 490 k = 1, nyzp3
         do 480 j = 1, nx3
         do 470 n = 1, 3
         scr(n,j,k+2*nyzp3,2,m) = -scr(n,j,k+nyzp3,1,m)
         f(n,j,k,nyzp(2,m)+2,m) = f(n,j,k,nyzp(2,m)+2,m) + 2.*scr(n,j,k+
     1nyzp3,1,m)
         f(n,j,k,nyzp(2,m)+3,m) = 0.
  470    continue
  480    continue
  490    continue
      endif
c special case of only one grid per processor in z
      if (nter.ge.1) go to 650
      if (jl.ge.0) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+3,lgrp,msid,ierr)
      endif
      if (jll.ge.0) then
         call MPI_IRECV(scr(1,1,1,1,m),3*nxvy,mreal,kll-1,noff+7,lgrp,ns
     1id,ierr)
      else if (jl.eq.0) then
         call MPI_IRECV(scr(1,1,1,1,m),3*nxvy,mreal,kl-1,noff+7,lgrp,nsi
     1d,ierr)
      endif
      if (jr.lt.nvpz) then
         call MPI_SEND(nyzp(2,m),1,mint,kr-1,moff+3,lgrp,ierr)
      endif
      if (jrr.lt.nvpz) then
         call MPI_SEND(scr(1,1,nyzp3+1,1,m),3*nxvys,mreal,krr-1,noff+7,l
     1grp,ierr)
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpz)) then
         call MPI_SEND(scr(1,1,2*nyzp3+1,1,m),3*nxvys,mreal,kr-1,noff+7,
     1lgrp,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         if (ngc.eq.1) then
            if (jl.eq.0) then
               do 520 k = 1, nyzp3
               do 510 j = 1, nx3
               do 500 n = 1, 3
               scr(n,j,k,2,m) = scr(n,j,k,2,m) - scr(n,j,k,1,m)
  500          continue
  510          continue
  520          continue
            else
               do 550 k = 1, nyzp3
               do 540 j = 1, nx3
               do 530 n = 1, 3
               scr(n,j,k,2,m) = scr(n,j,k,2,m) + scr(n,j,k,1,m)
  530          continue
  540          continue
  550          continue
            endif
         else if (nyzp(2,m).eq.1) then
            do 580 k = 1, nyzp3
            do 570 j = 1, nx3
            do 560 n = 1, 3
            scr(n,j,k+nyzp3,2,m) = 0.
  560       continue
  570       continue
  580       continue
         endif
      endif
c add up the guard cells
  650 do 680 k = 1, nyzp3
      do 670 j = 1, nx3
      do 660 n = 1, 3
      f(n,j,k,2,m) = f(n,j,k,2,m) + scr(n,j,k,2,m)
      f(n,j,k,3,m) = f(n,j,k,3,m) + scr(n,j,k+nyzp3,2,m)
      f(n,j,k,nyzp(2,m)+1,m) = f(n,j,k,nyzp(2,m)+1,m) + scr(n,j,k+2*nyzp
     13,2,m)
      f(n,j,k,1,m) = 0.
      f(n,j,k,nyzp(2,m)+3,m) = 0.
  660 continue
  670 continue
  680 continue
      if (jr.lt.nvpz) then
         do 710 k = 1, nyzp3
         do 700 j = 1, nx3
         do 690 n = 1, 3
         f(n,j,k,nyzp(2,m)+2,m) = 0.
  690    continue
  700    continue
  710    continue
      endif
  720 continue
  730 continue
c zero out the front edge
      do 780 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 770 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(1,m) + 3
      kl = mz + ks
      if (kl.eq.0) then
         do 760 k = 1, nyzp3
         do 750 j = 1, nx3
         do 740 n = 1, 3
         f(n,j,k,1,m) = 0.
  740    continue
  750    continue
  760    continue
      endif
  770 continue
  780 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLAGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypmx
     1,nzpmx,mblok,nblok,ngds,idds,mter,nter)
c this subroutine adds data from guard cells in non-uniform partitions
c for scalar data, the field is added up so as to disable quadratic
c interpolation within half a cell of the edges, and reduce it to linear
c interpolation in the y and z direction.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes three extra guard cells.
c scs/scr = scratch arrays for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c it is assumed the nyzp(n,m) > 0.
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c mter/nter = (0,1) = (no,yes) pass data to next processor only in y/z
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds, mter, nter
      integer nyzp
      real f, scs, scr
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx*ngds,2,mblok*nblok)
      dimension scr(nxv,nypmx*ngds,2,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, mnblok
      integer nx3, nxvz, nxvzs, nyzp3, nxvy, nxvys, m, my, mz, j, k
      integer jr, jrr, jl, jll
      dimension istatus(lstat)
      nx3 = nx + 3
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c buffer data in y
      do 30 m = 1, mnblok
      nyzp3 = nyzp(2,m) + 3
      do 20 k = 1, nyzp3
      do 10 j = 1, nxv
      scs(j,k,1,m) = f(j,nyzp(1,m)+2,k,m)
      scs(j,k+nyzp3,1,m) = f(j,nyzp(1,m)+3,k,m)
      scs(j,k+2*nyzp3,1,m) = f(j,1,k,m)
   10 continue
   20 continue
   30 continue
c add guard cells in y
      do 240 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 230 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(2,m) + 3
      nxvzs = nxv*nyzp3
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr + 1
      krr = jrr + kz
      jl = ky - 1
      jll = jl - 1
      kll = jll + kz
      kr = jr + kz
      kl = jl + kz
      ngc = 0
c special case of only one grid per processor
      if (nyzp(1,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        if (nyzp(1,kl).eq.1) then
c           do 50 k = 1, nyzp3
c           do 40 j = 1, nxv
c           scs(j,k,2,m) = scs(j,k,1,kl) + scs(j,k,1,kll)
c           scs(j,k+nyzp3,2,m) = scs(j,k+nyzp3,1,kl)
c  40       continue
c  50       continue
c        else
c           do 70 k = 1, nyzp3
c           do 60 j = 1, nxv
c           scs(j,k,2,m) = scs(j,k,1,kl)
c           scs(j,k+nyzp3,2,m) = scs(j,k+nyzp3,1,kl)=
c  60       continue
c  70       continue
c        endif
c     else
c        do 90 k = 1, nyzp3
c        do 80 j = 1, nx3
c        scs(j,k,2,m) = 2.*scs(j,k+2*nyzp3,1,m)
c        scs(j,k+nyzp3,2,m) = -scs(j,k+2*nyzp3,1,m)
c  80    continue
c  90    continue
c     endif
c     if (jr.lt.nvpy) then
c        do 110 k = 1, nyzp3
c        do 100 j = 1, nxv
c        scs(j,k+2*nyzp3,2,m) = scs(j,k+2*nyzp3,1,kr)
c 100    continue
c 110    continue
c     else
c        do 130 k = 1, nyzp3
c        do 120 j = 1, nx3
c        scs(j,k+2*nyzp3,2,m) = -scs(j,k+nyzp3,1,m)
c        f(j,nyzp(1,m)+2,k,m) = f(j,nyzp(1,m)+2,k,m) + 2.*scs(j,k+nyzp3,
c    11,m)
c        f(j,nyzp(1,m)+3,k,m) = 0.
c 120    continue
c 130    continue
c     endif
c     if (jl.ge.0) then
c        if ((nyzp(1,kl).eq.1).and.(jl.eq.0)) then
c           do 150 k = 1, nyzp3
c           do 140 j = 1, nx3
c           scs(j,k,2,m) = scs(j,k,2,m) - scs(j,k+2*nyzp3,1,kl)
c 140       continue
c 150       continue
c        endif
c     else if (nyzp(1,m).eq.1) then
c        do 170 k = 1, nyzp3
c        do 160 j = 1, nx3
c        scs(j,k+nyzp3,2,m) = 0.
c 160    continue
c 170    continue
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,2,m),2*nxvz,mreal,kl-1,noff+1,lgrp,msid,
     1ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,m),2*nxvzs,mreal,kr-1,noff+1,lgrp,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 50 k = 1, nyzp3
         do 40 j = 1, nx3
         scs(j,k,2,m) = 2.*scs(j,k+2*nyzp3,1,m)
         scs(j,k+nyzp3,2,m) = -scs(j,k+2*nyzp3,1,m)
   40    continue
   50    continue
      endif
      if (jr.lt.nvpy) then
         call MPI_IRECV(scs(1,2*nyzp3+1,2,m),nxvz,mreal,kr-1,noff+2,lgrp
     1,msid,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scs(1,2*nyzp3+1,1,m),nxvzs,mreal,kl-1,noff+2,lgrp
     1,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 70 k = 1, nyzp3
         do 60 j = 1, nx3
         scs(j,k+2*nyzp3,2,m) = -scs(j,k+nyzp3,1,m)
         f(j,nyzp(1,m)+2,k,m) = f(j,nyzp(1,m)+2,k,m) + 2.*scs(j,k+nyzp3,
     11,m)
         f(j,nyzp(1,m)+3,k,m) = 0.
   60    continue
   70    continue
      endif
c special case of only one grid per processor in y
      if (mter.ge.1) go to 180
      if (jl.ge.0) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+3,lgrp,msid,ierr)
      endif
      if (jll.ge.0) then
         call MPI_IRECV(scs(1,1,1,m),nxvz,mreal,kll-1,noff+5,lgrp,nsid,i
     1err)
      else if (jl.eq.0) then
         call MPI_IRECV(scs(1,1,1,m),nxvz,mreal,kl-1,noff+5,lgrp,nsid,ie
     1rr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(nyzp(1,m),1,mint,kr-1,moff+3,lgrp,ierr)
      endif
      if (jrr.lt.nvpy) then
         call MPI_SEND(scs(1,nyzp3+1,1,m),nxvzs,mreal,krr-1,noff+5,lgrp,
     1ierr)
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
         call MPI_SEND(scs(1,2*nyzp3+1,1,m),nxvzs,mreal,kr-1,noff+5,lgrp
     1,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         if (ngc.eq.1) then
            if (jl.eq.0) then
               do 90 k = 1, nyzp3
               do 80 j = 1, nx3
               scs(j,k,2,m) = scs(j,k,2,m) - scs(j,k,1,m)
   80          continue
   90          continue
            else
               do 110 k = 1, nyzp3
               do 100 j = 1, nx3
               scs(j,k,2,m) = scs(j,k,2,m) + scs(j,k,1,m)
  100          continue
  110          continue
            endif
         else if (nyzp(1,m).eq.1) then
            do 130 k = 1, nyzp3
            do 120 j = 1, nx3
            scs(j,k+nyzp3,2,m) = 0.
  120       continue
  130       continue
         endif
      endif
c add up the guard cells
  180 do 200 k = 1, nyzp3
      do 190 j = 1, nx3
      f(j,2,k,m) = f(j,2,k,m) + scs(j,k,2,m)
      f(j,3,k,m) = f(j,3,k,m) + scs(j,k+nyzp3,2,m)
      f(j,nyzp(1,m)+1,k,m) = f(j,nyzp(1,m)+1,k,m) + scs(j,k+2*nyzp3,2,m)
      f(j,1,k,m) = 0.
      f(j,nyzp(1,m)+3,k,m) = 0.
  190 continue
  200 continue
      if (jr.lt.nvpy) then
         do 220 k = 1, nyzp3
         do 210 j = 1, nx3
         f(j,nyzp(1,m)+2,k,m) = 0.
  210    continue
  220    continue
      endif
  230 continue
  240 continue
c zero out the left edge
      do 280 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 270 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(2,m) + 3
      kl = my + js
      if (kl.eq.0) then
         do 260 k = 1, nyzp3
         do 250 j = 1, nx3
         f(j,1,k,m) = 0.
  250    continue
  260    continue
      endif
  270 continue
  280 continue
c buffer data in z
      do 310 m = 1, mnblok
      nyzp3 = nyzp(1,m) + 3
      do 300 k = 1, nyzp3
      do 290 j = 1, nxv
      scr(j,k,1,m) = f(j,k,nyzp(2,m)+2,m)
      scr(j,k+nyzp3,1,m) = f(j,k,nyzp(2,m)+3,m)
      scr(j,k+2*nyzp3,1,m) = f(j,k,1,m)
  290 continue
  300 continue
  310 continue
c add guard cells in z
      do 520 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 510 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(1,m) + 3
      nxvys = nxv*nyzp3
      ky = my + js + 1
      kz = mz + ks
      jr = kz + 1
      jrr = jr + 1
      krr = ky + nvpy*jrr
      jl = kz - 1
      jll = jl - 1
      kll = ky + nvpy*jll
      kr = ky + nvpy*jr
      kl = ky + nvpy*jl
      ngc = 0
c special case of only one grid per processor
      if (nyzp(2,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        if (nyzp(2,kl).eq.1) then
c           do 330 k = 1, nyzp3
c           do 320 j = 1, nxv
c           scr(j,k,2,m) = scr(j,k,1,kl) + scr(j,k,1,kll)
c           scr(j,k+nyzp3,2,m) = scr(j,k+nyzp3,1,kl)
c 320       continue
c 330       continue
c       else
c           do 350 k = 1, nyzp3
c           do 340 j = 1, nxv
c           scr(j,k,2,m) = scr(j,k,1,kl)
c           scr(j,k+nyzp3,2,m) = scr(j,k+nyzp3,1,kl)
c 340       continue
c 350       continue
c        endif
c     else
c        do 370 k = 1, nyzp3
c        do 360 j = 1, nx3
c        scr(j,k,2,m) = 2.*scr(j,k+2*nyzp3,1,m)
c        scr(j,k+nyzp3,2,m) = -scr(j,k+2*nyzp3,1,m)
c 360    continue
c 370    continue
c     endif
c     if (jr.lt.nvpz) then
c        do 390 k = 1, nyzp3
c        do 380 j = 1, nxv
c        scr(j,k+2*nyzp3,2,m) = scr(j,k+2*nyzp3,1,kr)
c 380    continue
c 390    continue
c     else
c        do 410 k = 1, nyzp3
c        do 400 j = 1, nx3
c        scr(j,k+2*nyzp3,2,m) = -scr(j,k+nyzp3,1,m)
c        f(j,k,nyzp(2,m)+2,m) = f(j,k,nyzp(2,m)+2,m) + 2.*scr(j,k+nyzp3,
c    11,m)
c        f(j,k,nyzp(2,m)+3,m) = 0.
c 400    continue
c 410    continue
c     endif
c     if (jl.ge.0) then
c        if ((nyzp(2,kl).eq.1).and.(jl.eq.0)) then
c           do 430 k = 1, nyzp3
c           do 420 j = 1, nx3
c           scr(j,k,2,m) = scr(j,k,2,m) - scr(j,k+2*nyzp3,1,kl)
c 420       continue
c 430       continue
c        endif
c     else if (nyzp(2,m).eq.1) then
c        do 450 k = 1, nyzp3
c        do 440 j = 1, nx3
c        scr(j,k+nyzp3,2,m) = 0.
c 440    continue
c 450    continue
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scr(1,1,2,m),2*nxvy,mreal,kl-1,noff+3,lgrp,msid,
     1ierr)
      endif
      if (jr.lt.nvpz) then
         call MPI_SEND(scr(1,1,1,m),2*nxvys,mreal,kr-1,noff+3,lgrp,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 330 k = 1, nyzp3
         do 320 j = 1, nx3
         scr(j,k,2,m) = 2.*scr(j,k+2*nyzp3,1,m)
         scr(j,k+nyzp3,2,m) = -scr(j,k+2*nyzp3,1,m)
  320    continue
  330    continue
      endif
      if (jr.lt.nvpz) then
         call MPI_IRECV(scr(1,2*nyzp3+1,2,m),nxvy,mreal,kr-1,noff+4,lgrp
     1,msid,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scr(1,2*nyzp3+1,1,m),nxvys,mreal,kl-1,noff+4,lgrp
     1,ierr)
      endif
      if (jr.lt.nvpz) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 350 k = 1, nyzp3
         do 340 j = 1, nx3
         scr(j,k+2*nyzp3,2,m) = -scr(j,k+nyzp3,1,m)
         f(j,k,nyzp(2,m)+2,m) = f(j,k,nyzp(2,m)+2,m) + 2.*scr(j,k+nyzp3,
     11,m)
         f(j,k,nyzp(2,m)+3,m) = 0.
  340    continue
  350    continue
      endif
c special case of only one grid per processor in z
      if (nter.ge.1) go to 460
      if (jl.ge.0) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+3,lgrp,msid,ierr)
      endif
      if (jll.ge.0) then
         call MPI_IRECV(scr(1,1,1,m),nxvy,mreal,kll-1,noff+7,lgrp,nsid,i
     1err)
      else if (jl.eq.0) then
         call MPI_IRECV(scr(1,1,1,m),nxvy,mreal,kl-1,noff+7,lgrp,nsid,ie
     1rr)
      endif
      if (jr.lt.nvpz) then
         call MPI_SEND(nyzp(2,m),1,mint,kr-1,moff+3,lgrp,ierr)
      endif
      if (jrr.lt.nvpz) then
         call MPI_SEND(scr(1,nyzp3+1,1,m),nxvys,mreal,krr-1,noff+7,lgrp,
     1ierr)
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpz)) then
         call MPI_SEND(scr(1,2*nyzp3+1,1,m),nxvys,mreal,kr-1,noff+7,lgrp
     1,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         if (ngc.eq.1) then
            if (jl.eq.0) then
               do 370 k = 1, nyzp3
               do 360 j = 1, nx3
               scr(j,k,2,m) = scr(j,k,2,m) - scr(j,k,1,m)
  360          continue
  370          continue
            else
               do 390 k = 1, nyzp3
               do 380 j = 1, nx3
               scr(j,k,2,m) = scr(j,k,2,m) + scr(j,k,1,m)
  380          continue
  390          continue
            endif
         else if (nyzp(2,m).eq.1) then
            do 410 k = 1, nyzp3
            do 400 j = 1, nx3
            scr(j,k+nyzp3,2,m) = 0.
  400       continue
  410       continue
         endif
      endif
c add up the guard cells
  460 do 480 k = 1, nyzp3
      do 470 j = 1, nx3
      f(j,k,2,m) = f(j,k,2,m) + scr(j,k,2,m)
      f(j,k,3,m) = f(j,k,3,m) + scr(j,k+nyzp3,2,m)
      f(j,k,nyzp(2,m)+1,m) = f(j,k,nyzp(2,m)+1,m) + scr(j,k+2*nyzp3,2,m)
      f(j,k,1,m) = 0.
      f(j,k,nyzp(2,m)+3,m) = 0.
  470 continue
  480 continue
      if (jr.lt.nvpz) then
         do 500 k = 1, nyzp3
         do 490 j = 1, nx3
         f(j,k,nyzp(2,m)+2,m) = 0.
  490    continue
  500    continue
      endif
  510 continue
  520 continue
c zero out the front edge
      do 560 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 550 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(1,m) + 3
      kl = mz + ks
      if (kl.eq.0) then
         do 540 k = 1, nyzp3
         do 530 j = 1, nx3
         f(j,k,1,m) = 0.
  530    continue
  540    continue
      endif
  550 continue
  560 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLACGUARDS32(f,scs,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx,m
     1blok,nblok,kyp,kzp,ngds)
c this subroutine corrects the charge density data for particle boundary
c conditions which keep particles one grid away from the edges
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y and z direction.
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes three extra
c guard cells.
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      real f, scs
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer kyp, kzp, ngds
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx,ngds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, mnblok
      integer kzp1, nxvz, nxvy, m, my, mz, j, k, n
      integer jr, jrr, jl, jll
      dimension istatus(lstat)
      kzp1 = kzp + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c buffer data in y
      do 40 m = 1, mnblok
      do 30 k = 1, nzpmx
      do 20 j = 1, nxv
      do 10 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,kyp+2,k,m)
      scs(n,j,k,2,m) = f(n,j,2,k,m)
   10 continue
   20 continue
   30 continue
   40 continue
c add guard cells in y
      do 360 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 350 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr
      krr = jr + kz
      jl = ky - 1
      jll = jl
      kll = jl + kz
c special case of only one grid per processor
      if (kyp.eq.1) then
         jrr = jr + 1
         krr = jrr + kz
         jll = jl - 1
         kll = jll + kz
      endif
      kr = jr + kz
      kl = jl + kz
c fix edges if all points are on the same processor
      if (jl.eq.(-1)) then
         if (kyp.gt.2) then
            do 70 k = 1, kzp1
            do 60 j = 2, nx
            do 50 n = 1, 3
            f(n,j+1,3,k+1,m) = f(n,j+1,3,k+1,m) + 2.*f(n,j+1,2,k+1,m)
            f(n,j+1,4,k+1,m) = f(n,j+1,4,k+1,m) - f(n,j+1,2,k+1,m)
            f(n,j+1,2,k+1,m) = 0.
   50       continue
   60       continue
   70       continue
         else if (kyp.eq.2) then
            do 100 k = 1, kzp1
            do 90 j = 2, nx
            do 80 n = 1, 3
            f(n,j+1,3,k+1,m) = f(n,j+1,3,k+1,m) + 2.*f(n,j+1,2,k+1,m)
   80       continue
   90       continue
  100       continue
         endif
      endif
      if (jr.eq.nvpy) then
         if (kyp.gt.1) then
            do 130 k = 1, kzp1
            do 120 j = 2, nx
            do 110 n = 1, 3
            f(n,j+1,kyp,k+1,m) = f(n,j+1,kyp,k+1,m) - f(n,j+1,kyp+2,k+1,
     1m)
            f(n,j+1,kyp+1,k+1,m) = f(n,j+1,kyp+1,k+1,m) + 2.*f(n,j+1,kyp
     1+2,k+1,m)
            f(n,j+1,kyp+2,k+1,m) = 0.
  110       continue
  120       continue
  130       continue
         else if (kyp.eq.1) then
            do 160 k = 1, kzp1
            do 150 j = 2, nx
            do 140 n = 1, 3
            f(n,j+1,kyp+1,k+1,m) = f(n,j+1,kyp+1,k+1,m) + 2.*f(n,j+1,kyp
     1+2,k+1,m)
            f(n,j+1,kyp+2,k+1,m) = 0.
  140       continue
  150       continue
  160       continue
         endif
      endif
c this segment is used for shared memory computers
c     if (kyp.eq.2) then
c        if (jl.eq.0) then
c           do 190 k = 1, kzp1
c           do 180 j = 2, nx
c           do 170 n = 1, 3
c           f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) - scs(n,j+1,k+1,2,kl)
c 170       continue
c 180       continue
c 190       continue
c        endif
c     else if (kyp.eq.1) then
c        if (jl.eq.0) then
c           do 220 k = 1, kzp1
c           do 210 j = 2, nx
c           do 200 n = 1, 3
c           f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) + 2.*scs(n,j+1,k+1,2,kl)
c 200       continue
c 210       continue
c 220       continue
c        endif
c        if (jll.eq.0) then
c           do 250 k = 1, kzp1
c           do 240 j = 2, nx
c           do 230 n = 1, 3
c           f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) - scs(n,j+1,k+1,2,kll)
c 230       continue
c 240       continue
c 250       continue
c        endif
c        if (jr.eq.(nvpy-1)) then
c           do 280 k = 1, kzp1
c           do 270 j = 2, nx
c           do 260 n = 1, 3
c           f(n,j+1,kyp+1,k+1,m) = f(n,j+1,kyp+1,k+1,m) - scs(n,j+1,k+1,
c    11,kr)
c 260       continue
c 270       continue
c 280       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (kyp.eq.2) then
         if (jl.eq.0) then
            call MPI_IRECV(scs(1,1,1,3,m),3*nxvz,mreal,kl-1,noff+1,lgrp,
     1msid,ierr)
         endif
         if (jl.eq.(-1)) then
            call MPI_SEND(scs(1,1,1,2,m),3*nxvz,mreal,kr-1,noff+1,lgrp,i
     1err)
            do 190 k = 1, kzp1
            do 180 j = 2, nx
            do 170 n = 1, 3
            f(n,j+1,2,k+1,m) = 0.
  170       continue
  180       continue
  190       continue
         endif
         if (jl.eq.0) then
            call MPI_WAIT(msid,istatus,ierr)
            do 220 k = 1, kzp1
            do 210 j = 2, nx
            do 200 n = 1, 3
            f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) - scs(n,j+1,k+1,3,m)
            f(n,j+1,1,k+1,m) = 0.
  200       continue
  210       continue
  220       continue
         endif
      else if (kyp.eq.1) then
         if (jl.eq.0) then
            call MPI_IRECV(scs(1,1,1,3,m),3*nxvz,mreal,kl-1,noff+1,lgrp,
     1msid,ierr)
         endif
         if (jll.eq.0) then
            call MPI_IRECV(scs(1,1,1,3,m),3*nxvz,mreal,kll-1,noff+1,lgrp
     1,nsid,ierr)
         endif
         if (jl.eq.(-1)) then
            call MPI_SEND(scs(1,1,1,2,m),3*nxvz,mreal,kr-1,noff+1,lgrp,i
     1err)
            call MPI_SEND(scs(1,1,1,2,m),3*nxvz,mreal,krr-1,noff+1,lgrp,
     1ierr)
            do 250 k = 1, kzp1
            do 240 j = 2, nx
            do 230 n = 1, 3
            f(n,j+1,2,k+1,m) = 0.
  230       continue
  240       continue
  250       continue
         endif
         if (jl.eq.0) then
            call MPI_WAIT(msid,istatus,ierr)
            do 280 k = 1, kzp1
            do 270 j = 2, nx
            do 260 n = 1, 3
            f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) + 2.*scs(n,j+1,k+1,3,m)
            f(n,j+1,1,k+1,m) = 0.
  260       continue
  270       continue
  280       continue
         endif
         if (jll.eq.0) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 310 k = 1, kzp1
            do 300 j = 2, nx
            do 290 n = 1, 3
            f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) - scs(n,j+1,k+1,3,m)
            f(n,j+1,1,k+1,m) = 0.
  290       continue
  300       continue
  310       continue
         endif
         if (jr.eq.(nvpy-1)) then
            call MPI_IRECV(scs(1,1,1,3,m),3*nxvz,mreal,kr-1,noff+2,lgrp,
     1msid,ierr)
         endif
         if (jr.eq.nvpy) then
            call MPI_SEND(scs(1,1,1,1,m),3*nxvz,mreal,kl-1,noff+2,lgrp,i
     1err)
         endif
         if (jr.eq.(nvpy-1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 340 k = 1, kzp1
            do 330 j = 2, nx
            do 320 n = 1, 3
            f(n,j+1,kyp+1,k+1,m) = f(n,j+1,kyp+1,k+1,m) - scs(n,j+1,k+1,
     13,m)
            f(n,j+1,1,k+1,m) = 0.
  320       continue
  330       continue
  340       continue
         endif
      endif
  350 continue
  360 continue
c add guard cells in z
      do 710 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 700 my = 1, mblok
      m = my + moff
      ky = my + js + 1
      kz = mz + ks
      jr = kz + 1
      jrr = jr
      krr = ky + nvpy*jr
      jl = kz - 1
      jll = jl
      kll = ky + nvpy*jl
c special case of only one grid per processor
      if (kzp.eq.1) then
         jrr = jr + 1
         krr = ky + nvpy*jrr
         jll = jl - 1
         kll = ky + nvpy*jll
      endif
      kr = ky + nvpy*jr
      kl = ky + nvpy*jl
c fix edges if all points are on the same processor
      if (jl.eq.(-1)) then
         if (kzp.gt.2) then
            do 390 k = 1, kyp
            do 380 j = 2, nx
            do 370 n = 1, 3
            f(n,j+1,k+1,3,m) = f(n,j+1,k+1,3,m) + 2.*f(n,j+1,k+1,2,m)
            f(n,j+1,k+1,4,m) = f(n,j+1,k+1,4,m) - f(n,j+1,k+1,2,m)
            f(n,j+1,k+1,2,m) = 0.
  370       continue
  380       continue
  390       continue
         else if (kzp.eq.2) then
            do 420 k = 1, kyp
            do 410 j = 2, nx
            do 400 n = 1, 3
            f(n,j+1,k+1,3,m) = f(n,j+1,k+1,3,m) + 2.*f(n,j+1,k+1,2,m)
  400       continue
  410       continue
  420       continue
         endif
      endif
      if (jr.eq.nvpz) then
         if (kzp.gt.1) then
            do 450 k = 1, kyp
            do 440 j = 2, nx
            do 430 n = 1, 3
            f(n,j+1,k+1,kzp,m) = f(n,j+1,k+1,kzp,m) - f(n,j+1,k+1,kzp+2,
     1m)
            f(n,j+1,k+1,kzp+1,m) = f(n,j+1,k+1,kzp+1,m) + 2.*f(n,j+1,k+1
     1,kzp+2,m)
            f(n,j+1,k+1,kzp+2,m) = 0.
  430       continue
  440       continue
  450       continue
         else if (kzp.eq.1) then
            do 480 k = 1, kyp
            do 470 j = 2, nx
            do 460 n = 1, 3
            f(n,j+1,k+1,kzp+1,m) = f(n,j+1,k+1,kzp+1,m) + 2.*f(n,j+1,k+1
     1,kzp+2,m)
  460       continue
  470       continue
  480       continue
         endif
      endif
c this segment is used for shared memory computers
c     if (kzp.eq.2) then
c        if (jl.eq.0) then
c           do 510 k = 1, kyp
c           do 500 j = 2, nx
c           do 490 n = 1, 3
c           f(n,j+1,k+1,2,m) = f(n,j+1,k+1,2,m) - f(n,j+1,k+1,2,kl)
c 490       continue
c 500       continue
c 510       continue
c        endif
c     else if (kzp.eq.1) then
c        if (jl.eq.0) then
c           do 540 k = 1, kyp
c           do 530 j = 2, nx
c           do 520 n = 1, 3
c           f(n,j+1,k+1,2,m) = f(n,j+1,k+1,2,m) + 2.*f(n,j+1,k+1,2,kl)
c 520       continue
c 530       continue
c 540       continue
c        endif
c        if (jll.eq.0) then
c           do 570 k = 1, kyp
c           do 560 j = 2, nx
c           do 550 n = 1, 3
c           f(n,j+1,k+1,2,m) = f(n,j+1,k+1,2,m) - f(n,j+1,k+1,2,kll)
c 550       continue
c 560       continue
c 570       continue
c        endif
c        if (jr.eq.(nvpz-1)) then
c           do 600 k = 1, kyp
c           do 590 j = 2, nx
c           do 580 n = 1, 3
c           f((n,j+1,k+1,kzp+1,m) = f((n,j+1,k+1,kzp+1,m) - f((n,j+1,k+1
c    1,kzp+2,kr)
c 580       continue
c 590       continue
c 600       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (kzp.eq.2) then
         if (jl.eq.0) then
            call MPI_IRECV(f(1,1,1,1,m),3*nxvy,mreal,kl-1,noff+3,lgrp,ms
     1id,ierr)
         endif
         if (jl.eq.(-1)) then
            call MPI_SEND(f(1,1,1,2,m),3*nxvy,mreal,kr-1,noff+3,lgrp,ier
     1r)
            do 510 k = 1, kyp
            do 500 j = 2, nx
            do 490 n = 1, 3
            f(n,j+1,k+1,2,m) = 0.
  490       continue
  500       continue
  510       continue
         endif
         if (jl.eq.0) then
            call MPI_WAIT(msid,istatus,ierr)
            do 540 k = 1, kyp
            do 530 j = 2, nx
            do 520 n = 1, 3
            f(n,j+1,k+1,2,m) = f(n,j+1,k+1,2,m) - f(n,j+1,k+1,1,m)
            f(n,j+1,k+1,1,m) = 0.
  520       continue
  530       continue
  540       continue
         endif
      else if (kzp.eq.1) then
         if (jl.eq.0) then
            call MPI_IRECV(f(1,1,1,1,m),3*nxvy,mreal,kl-1,noff+3,lgrp,ms
     1id,ierr)
         endif
         if (jll.eq.0) then
            call MPI_IRECV(f(1,1,1,1,m),3*nxvy,mreal,kll-1,noff+3,lgrp,n
     1sid,ierr)
         endif
         if (jl.eq.(-1)) then
            call MPI_SEND(f(1,1,1,2,m),3*nxvy,mreal,kr-1,noff+3,lgrp,ier
     1r)
            call MPI_SEND(f(1,1,1,2,m),3*nxvy,mreal,krr-1,noff+3,lgrp,ie
     1rr)
            do 570 k = 1, kyp
            do 560 j = 2, nx
            do 550 n = 1, 3
            f(n,j+1,k+1,2,m) = 0.
  550       continue
  560       continue
  570       continue
         endif
         if (jl.eq.0) then
            call MPI_WAIT(msid,istatus,ierr)
            do 600 k = 1, kyp
            do 590 j = 2, nx
            do 580 n = 1, 3
            f(n,j+1,k+1,2,m) = f(n,j+1,k+1,2,m) + 2.*f(n,j+1,k+1,1,m)
            f(n,j+1,k+1,1,m) = 0.
  580       continue
  590       continue
  600       continue
         endif
         if (jll.eq.0) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 630 k = 1, kyp
            do 620 j = 2, nx
            do 610 n = 1, 3
            f(n,j+1,k+1,2,m) = f(n,j+1,k+1,2,m) - f(n,j+1,k+1,1,m)
            f(n,j+1,k+1,1,m) = 0.
  610       continue
  620       continue
  630       continue
         endif
         if (jr.eq.(nvpz-1)) then
            call MPI_IRECV(f(1,1,1,1,m),3*nxvy,mreal,kr-1,noff+4,lgrp,ms
     1id,ierr)
         endif
         if (jr.eq.nvpz) then
            call MPI_SEND(f(1,1,1,kzp+2,m),3*nxvy,mreal,kl-1,noff+4,lgrp
     1,ierr)
            do 660 k = 1, kyp
            do 650 j = 2, nx
            do 640 n = 1, 3
            f(n,j+1,k+1,kzp+2,m) = 0.
  640       continue
  650       continue
  660       continue
         endif
         if (jr.eq.(nvpz-1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 690 k = 1, kyp
            do 680 j = 2, nx
            do 670 n = 1, 3
            f(n,j+1,k+1,kzp+1,m) = f(n,j+1,k+1,kzp+1,m) - f(n,j+1,k+1,1,
     1m)
            f(n,j+1,k+1,1,m) = 0.
  670       continue
  680       continue
  690       continue
         endif
      endif
  700 continue
  710 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLAGUARDS32(f,scs,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx,mb
     1lok,nblok,kyp,kzp,ngds)
c this subroutine corrects the charge density data for particle boundary
c conditions which keep particles one grid away from the edges
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y and z direction.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes three extra
c guard cells.
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      real f, scs
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer kyp, kzp, ngds
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,ngds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, mnblok
      integer kzp1, nxvz, nxvy, m, my, mz, j, k
      integer jr, jrr, jl, jll
      dimension istatus(lstat)
      kzp1 = kzp + 1
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
      scs(j,k,1,m) = f(j,kyp+2,k,m)
      scs(j,k,2,m) = f(j,2,k,m)
   10 continue
   20 continue
   30 continue
c add guard cells in y
      do 250 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 240 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr
      krr = jr + kz
      jl = ky - 1
      jll = jl
      kll = jl + kz
c special case of only one grid per processor
      if (kyp.eq.1) then
         jrr = jr + 1
         krr = jrr + kz
         jll = jl - 1
         kll = jll + kz
      endif
      kr = jr + kz
      kl = jl + kz
c fix edges if all points are on the same processor
      if (jl.eq.(-1)) then
         if (kyp.gt.2) then
            do 50 k = 1, kzp1
            do 40 j = 2, nx
            f(j+1,3,k+1,m) = f(j+1,3,k+1,m) + 2.*f(j+1,2,k+1,m)
            f(j+1,4,k+1,m) = f(j+1,4,k+1,m) - f(j+1,2,k+1,m)
            f(j+1,2,k+1,m) = 0.
   40       continue
   50       continue
         else if (kyp.eq.2) then
            do 70 k = 1, kzp1
            do 60 j = 2, nx
            f(j+1,3,k+1,m) = f(j+1,3,k+1,m) + 2.*f(j+1,2,k+1,m)
   60       continue
   70       continue
         endif
      endif
      if (jr.eq.nvpy) then
         if (kyp.gt.1) then
            do 90 k = 1, kzp1
            do 80 j = 2, nx
            f(j+1,kyp,k+1,m) = f(j+1,kyp,k+1,m) - f(j+1,kyp+2,k+1,m)
            f(j+1,kyp+1,k+1,m) = f(j+1,kyp+1,k+1,m) + 2.*f(j+1,kyp+2,k+1
     1,m)
            f(j+1,kyp+2,k+1,m) = 0.
   80       continue
   90       continue
         else if (kyp.eq.1) then
            do 110 k = 1, kzp1
            do 100 j = 2, nx
            f(j+1,kyp+1,k+1,m) = f(j+1,kyp+1,k+1,m) + 2.*f(j+1,kyp+2,k+1
     1,m)
            f(j+1,kyp+2,k+1,m) = 0.
  100       continue
  110       continue
         endif
      endif
c this segment is used for shared memory computers
c     if (kyp.eq.2) then
c        if (jl.eq.0) then
c           do 130 k = 1, kzp1
c           do 120 j = 2, nx
c           f(j+1,2,k+1,m) = f(j+1,2,k+1,m) - scs(j+1,k+1,2,kl)
c 120       continue
c 130       continue
c        endif
c     else if (kyp.eq.1) then
c        if (jl.eq.0) then
c           do 150 k = 1, kzp1
c           do 140 j = 2, nx
c           f(j+1,2,k+1,m) = f(j+1,2,k+1,m) + 2.*scs(j+1,k+1,2,kl)
c 140       continue
c 150       continue
c        endif
c        if (jll.eq.0) then
c           do 170 k = 1, kzp1
c           do 160 j = 2, nx
c           f(j+1,2,k+1,m) = f(j+1,2,k+1,m) - scs(j+1,k+1,2,kll)
c 160       continue
c 170       continue
c        endif
c        if (jr.eq.(nvpy-1)) then
c           do 190 k = 1, kzp1
c           do 180 j = 2, nx
c           f(j+1,kyp+1,k+1,m) = f(j+1,kyp+1,k+1,m) - scs(j+1,k+1,1,kr)
c 180       continue
c 190       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (kyp.eq.2) then
         if (jl.eq.0) then
            call MPI_IRECV(scs(1,1,3,m),nxvz,mreal,kl-1,noff+1,lgrp,msid
     1,ierr)
         endif
         if (jl.eq.(-1)) then
            call MPI_SEND(scs(1,1,2,m),nxvz,mreal,kr-1,noff+1,lgrp,ierr)
            do 130 k = 1, kzp1
            do 120 j = 2, nx
            f(j+1,2,k+1,m) = 0.
  120       continue
  130       continue
         endif
         if (jl.eq.0) then
            call MPI_WAIT(msid,istatus,ierr)
            do 150 k = 1, kzp1
            do 140 j = 2, nx
            f(j+1,2,k+1,m) = f(j+1,2,k+1,m) - scs(j+1,k+1,3,m)
            f(j+1,1,k+1,m) = 0.
  140       continue
  150       continue
         endif
      else if (kyp.eq.1) then
         if (jl.eq.0) then
            call MPI_IRECV(scs(1,1,3,m),nxvz,mreal,kl-1,noff+1,lgrp,msid
     1,ierr)
         endif
         if (jll.eq.0) then
            call MPI_IRECV(scs(1,1,3,m),nxvz,mreal,kll-1,noff+1,lgrp,nsi
     1d,ierr)
         endif
         if (jl.eq.(-1)) then
            call MPI_SEND(scs(1,1,2,m),nxvz,mreal,kr-1,noff+1,lgrp,ierr)
            call MPI_SEND(scs(1,1,2,m),nxvz,mreal,krr-1,noff+1,lgrp,ierr
     1)
            do 170 k = 1, kzp1
            do 160 j = 2, nx
            f(j+1,2,k+1,m) = 0.
  160       continue
  170       continue
         endif
         if (jl.eq.0) then
            call MPI_WAIT(msid,istatus,ierr)
            do 190 k = 1, kzp1
            do 180 j = 2, nx
            f(j+1,2,k+1,m) = f(j+1,2,k+1,m) + 2.*scs(j+1,k+1,3,m)
            f(j+1,1,k+1,m) = 0.
  180       continue
  190       continue
         endif
         if (jll.eq.0) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 210 k = 1, kzp1
            do 200 j = 2, nx
            f(j+1,2,k+1,m) = f(j+1,2,k+1,m) - scs(j+1,k+1,3,m)
            f(j+1,1,k+1,m) = 0.
  200       continue
  210       continue
         endif
         if (jr.eq.(nvpy-1)) then
            call MPI_IRECV(scs(1,1,3,m),nxvz,mreal,kr-1,noff+2,lgrp,msid
     1,ierr)
         endif
         if (jr.eq.nvpy) then
            call MPI_SEND(scs(1,1,1,m),nxvz,mreal,kl-1,noff+2,lgrp,ierr)
         endif
         if (jr.eq.(nvpy-1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 230 k = 1, kzp1
            do 220 j = 2, nx
            f(j+1,kyp+1,k+1,m) = f(j+1,kyp+1,k+1,m) - scs(j+1,k+1,3,m)
            f(j+1,1,k+1,m) = 0.
  220       continue
  230       continue
         endif
      endif
  240 continue
  250 continue
c add guard cells in z
      do 490 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 480 my = 1, mblok
      m = my + moff
      ky = my + js + 1
      kz = mz + ks
      jr = kz + 1
      jrr = jr
      krr = ky + nvpy*jr
      jl = kz - 1
      jll = jl
      kll = ky + nvpy*jl
c special case of only one grid per processor
      if (kzp.eq.1) then
         jrr = jr + 1
         krr = ky + nvpy*jrr
         jll = jl - 1
         kll = ky + nvpy*jll
      endif
      kr = ky + nvpy*jr
      kl = ky + nvpy*jl
c fix edges if all points are on the same processor
      if (jl.eq.(-1)) then
         if (kzp.gt.2) then
            do 270 k = 1, kyp
            do 260 j = 2, nx
            f(j+1,k+1,3,m) = f(j+1,k+1,3,m) + 2.*f(j+1,k+1,2,m)
            f(j+1,k+1,4,m) = f(j+1,k+1,4,m) - f(j+1,k+1,2,m)
            f(j+1,k+1,2,m) = 0.
  260       continue
  270       continue
         else if (kzp.eq.2) then
            do 290 k = 1, kyp
            do 280 j = 2, nx
            f(j+1,k+1,3,m) = f(j+1,k+1,3,m) + 2.*f(j+1,k+1,2,m)
  280       continue
  290       continue
         endif
      endif
      if (jr.eq.nvpz) then
         if (kzp.gt.1) then
            do 310 k = 1, kyp
            do 300 j = 2, nx
            f(j+1,k+1,kzp,m) = f(j+1,k+1,kzp,m) - f(j+1,k+1,kzp+2,m)
            f(j+1,k+1,kzp+1,m) = f(j+1,k+1,kzp+1,m) + 2.*f(j+1,k+1,kzp+2
     1,m)
            f(j+1,k+1,kzp+2,m) = 0.
  300       continue
  310       continue
         else if (kzp.eq.1) then
            do 330 k = 1, kyp
            do 320 j = 2, nx
            f(j+1,k+1,kzp+1,m) = f(j+1,k+1,kzp+1,m) + 2.*f(j+1,k+1,kzp+2
     1,m)
  320       continue
  330       continue
         endif
      endif
c this segment is used for shared memory computers
c     if (kzp.eq.2) then
c        if (jl.eq.0) then
c           do 350 k = 1, kyp
c           do 340 j = 2, nx
c           f(j+1,k+1,2,m) = f(j+1,k+1,2,m) - f(j+1,k+1,2,kl)
c 340       continue
c 350       continue
c        endif
c     else if (kzp.eq.1) then
c        if (jl.eq.0) then
c           do 370 k = 1, kyp
c           do 360 j = 2, nx
c           f(j+1,k+1,2,m) = f(j+1,k+1,2,m) + 2.*f(j+1,k+1,2,kl)
c 360       continue
c 370       continue
c        endif
c        if (jll.eq.0) then
c           do 390 k = 1, kyp
c           do 380 j = 2, nx
c           f(j+1,k+1,2,m) = f(j+1,k+1,2,m) - f(j+1,k+1,2,kll)
c 380       continue
c 390       continue
c        endif
c        if (jr.eq.(nvpz-1)) then
c           do 410 k = 1, kyp
c           do 400 j = 2, nx
c           f(j+1,k+1,kzp+1,m) = f(j+1,k+1,kzp+1,m) - f(j+1,k+1,kzp+2,kr
c    1)
c 400       continue
c 410       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (kzp.eq.2) then
         if (jl.eq.0) then
            call MPI_IRECV(f(1,1,1,m),nxvy,mreal,kl-1,noff+3,lgrp,msid,i
     1err)
         endif
         if (jl.eq.(-1)) then
            call MPI_SEND(f(1,1,2,m),nxvy,mreal,kr-1,noff+3,lgrp,ierr)
            do 350 k = 1, kyp
            do 340 j = 2, nx
            f(j+1,k+1,2,m) = 0.
  340       continue
  350       continue
         endif
         if (jl.eq.0) then
            call MPI_WAIT(msid,istatus,ierr)
            do 370 k = 1, kyp
            do 360 j = 2, nx
            f(j+1,k+1,2,m) = f(j+1,k+1,2,m) - f(j+1,k+1,1,m)
            f(j+1,k+1,1,m) = 0.
  360       continue
  370       continue
         endif
      else if (kzp.eq.1) then
         if (jl.eq.0) then
            call MPI_IRECV(f(1,1,1,m),nxvy,mreal,kl-1,noff+3,lgrp,msid,i
     1err)
         endif
         if (jll.eq.0) then
            call MPI_IRECV(f(1,1,1,m),nxvy,mreal,kll-1,noff+3,lgrp,nsid,
     1ierr)
         endif
         if (jl.eq.(-1)) then
            call MPI_SEND(f(1,1,2,m),nxvy,mreal,kr-1,noff+3,lgrp,ierr)
            call MPI_SEND(f(1,1,2,m),nxvy,mreal,krr-1,noff+3,lgrp,ierr)
            do 390 k = 1, kyp
            do 380 j = 2, nx
            f(j+1,k+1,2,m) = 0.
  380       continue
  390       continue
         endif
         if (jl.eq.0) then
            call MPI_WAIT(msid,istatus,ierr)
            do 410 k = 1, kyp
            do 400 j = 2, nx
            f(j+1,k+1,2,m) = f(j+1,k+1,2,m) + 2.*f(j+1,k+1,1,m)
            f(j+1,k+1,1,m) = 0.
  400       continue
  410       continue
         endif
         if (jll.eq.0) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 430 k = 1, kyp
            do 420 j = 2, nx
            f(j+1,k+1,2,m) = f(j+1,k+1,2,m) - f(j+1,k+1,1,m)
            f(j+1,k+1,1,m) = 0.
  420       continue
  430       continue
         endif
         if (jr.eq.(nvpz-1)) then
            call MPI_IRECV(f(1,1,1,m),nxvy,mreal,kr-1,noff+4,lgrp,msid,i
     1err)
         endif
         if (jr.eq.nvpz) then
            call MPI_SEND(f(1,1,kzp+2,m),nxvy,mreal,kl-1,noff+4,lgrp,ier
     1r)
            do 450 k = 1, kyp
            do 440 j = 2, nx
            f(j+1,k+1,kzp+2,m) = 0.
  440       continue
  450       continue
         endif
         if (jr.eq.(nvpz-1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 470 k = 1, kyp
            do 460 j = 2, nx
            f(j+1,k+1,kzp+1,m) = f(j+1,k+1,kzp+1,m) - f(j+1,k+1,1,m)
            f(j+1,k+1,1,m) = 0.
  460       continue
  470       continue
         endif
      endif
  480 continue
  490 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLACGUARDS32(f,scs,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypmx,n
     1zpmx,mblok,nblok,ngds,idds,mter,nter)
c this subroutine corrects current density data for particle boundary
c conditions which keep particles one grid away from the edges
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y and z direction.
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes three extra guard cells.
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c it is assumed the nyzp(n,m) > 0.
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c mter/nter = (0,1) = (no,yes) pass data to next processor only in y/z
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      real f, scs
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds, mter, nter
      integer nyzp
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx,ngds,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, nps
      integer mnblok, kyp1, kzp1, nxvz, nxvzs, nxvy, nxvys
      integer m, my, mz, j, k, n, jr, jrr, jl, jll
      dimension istatus(lstat)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c buffer data in y
      do 40 m = 1, mnblok
      kzp1 = nyzp(2,m) + 1
      do 30 k = 1, kzp1 + 1
      do 20 j = 1, nxv
      do 10 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,nyzp(1,m)+2,k,m)
      scs(n,j,k,2,m) = f(n,j,2,k,m)
   10 continue
   20 continue
   30 continue
   40 continue
c add guard cells in y
      do 400 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 390 my = 1, mblok
      m = my + moff
      kzp1 = nyzp(2,m) + 1
      nxvzs = nxv*(kzp1 + 1)
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr + 1
      krr = jrr + kz
      jl = ky - 1
      jll = jl - 1
      kll = jll + kz
      kr = jr + kz
      kl = jl + kz
c fix edges if all points are on the same processor
      if (jl.eq.(-1)) then
         if (nyzp(1,m).gt.2) then
            do 70 k = 1, kzp1
            do 60 j = 2, nx
            do 50 n = 1, 3
            f(n,j+1,3,k+1,m) = f(n,j+1,3,k+1,m) + 2.*f(n,j+1,2,k+1,m)
            f(n,j+1,4,k+1,m) = f(n,j+1,4,k+1,m) - f(n,j+1,2,k+1,m)
            f(n,j+1,2,k+1,m) = 0.
   50       continue
   60       continue
   70       continue
         else if (nyzp(1,m).eq.2) then
            do 100 k = 1, kzp1
            do 90 j = 2, nx
            do 80 n = 1, 3
            f(n,j+1,3,k+1,m) = f(n,j+1,3,k+1,m) + 2.*f(n,j+1,2,k+1,m)
   80       continue
   90       continue
  100       continue
         endif
      endif
      if (jr.eq.nvpy) then
         if (nyzp(1,m).gt.1) then
            do 130 k = 1, kzp1
            do 120 j = 2, nx
            do 110 n = 1, 3
            f(n,j+1,nyzp(1,m),k+1,m) = f(n,j+1,nyzp(1,m),k+1,m) - f(n,j+
     11,nyzp(1,m)+2,k+1,m)
            f(n,j+1,nyzp(1,m)+1,k+1,m) = f(n,j+1,nyzp(1,m)+1,k+1,m) + 2.
     1*f(n,j+1,nyzp(1,m)+2,k+1,m)
            f(n,j+1,nyzp(1,m)+2,k+1,m) = 0.
  110       continue
  120       continue
  130       continue
         else if (nyzp(1,m).eq.1) then
            do 160 k = 1, kzp1
            do 150 j = 2, nx
            do 140 n = 1, 3
            f(n,j+1,nyzp(1,m)+1,k+1,m) = f(n,j+1,nyzp(1,m)+1,k+1,m) + 2.
     1*f(n,j+1,nyzp(1,m)+2,k+1,m)
            f(n,j+1,nyzp(1,m)+2,k+1,m) = 0.
  140       continue
  150       continue
  160       continue
         endif
      endif
c this segment is used for shared memory computers
c     if (mter.ge.2) go to 390
c     if (mter.eq.1) go to 260
c     if (jll.eq.0) then
c        if (nyzp(1,kll).eq.1).and.(nyzp(1,kl).eq.1)) then
c           do 190 k = 1, kzp1
c           do 180 j = 2, nx
c           do 170 n = 1, 3
c           f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) - scs(n,j+1,k+1,2,kll)
c 170       continue
c 180       continue
c 190       continue
c        endif
c     endif
c     if (jr.eq.(nvpy-1)) then
c        if (nyzp(1,kr).eq.1) then
c           do 250 k = 1, kzp1
c           do 240 j = 2, nx
c           do 230 n = 1, 3
c           f(n,j+1,nyzp(1,m)+1,k+1,m) = f(n,j+1,nyzp(1,m)+1,k+1,m) - sc
c    1s(n,j+1,k+1,1,kr)
c 230       continue
c 240       continue
c 250       continue
c        endif
c     endif
c 260 if (jl.eq.0) then
c        if (nyzp(1,kl).le.2) then
c           if (nyzp(1,kl).eq.1) then
c              if (nyzp(1,m).gt.1) then
c                 do 320 k = 1, kzp1
c                 do 310 j = 2, nx
c                 do 300 n = 1, 3
c                 f(n,j+1,3,k+1,m) = f(n,j+1,3,k+1,m) - scs(n,j+1,k+1,2,
c    1kl)
c 300             continue
c 310             continue
c 320             continue
c              endif
c              do 350 k = 1, kzp1
c              do 340 j = 2, nx
c              do 330 n = 1, 3
c              f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) + 2.*scs(n,j+1,k+1,2,
c    1kl)
c 330          continue
c 340          continue
c 350          continue
c           else
c              do 380 k = 1, kzp1
c              do 370 j = 2, nx
c              do 360 n = 1, 3
c              f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) - scs(n,j+1,k+1,2,kl)
c 360          continue
c 370          continue
c 380          continue
c           endif
c        endif
c     endif
c this segment is used for mpi computers
      if (mter.ge.2) go to 390
      if (mter.eq.1) go to 260
      if (jll.eq.0) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+4,lgrp,msid,ierr)
         call MPI_IRECV(scs(1,1,1,3,m),3*nxvz,mreal,kll-1,moff+1,lgrp,ns
     1id,ierr)
      endif
      if (jl.eq.0) then
         call MPI_SEND(nyzp(1,m),1,mint,kr-1,moff+4,lgrp,ierr)
      endif
      if (jl.eq.(-1)) then
         if (nyzp(1,m).eq.1) then
            call MPI_SEND(scs(1,1,1,2,m),3*nxvzs,mreal,krr-1,moff+1,lgrp
     1,ierr)
         else
            nps = 0
            call MPI_SEND(scs(1,1,1,2,m),nps,mreal,krr-1,moff+1,lgrp,ier
     1r)
         endif
      endif
      if (jll.eq.0) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            if (ngc.eq.1) then
               do 190 k = 1, kzp1
               do 180 j = 2, nx
               do 170 n = 1, 3
               f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) - scs(n,j+1,k+1,3,m)
  170          continue
  180          continue
  190          continue
            endif
            do 220 k = 1, kzp1
            do 210 j = 2, nx
            do 200 n = 1, 3
            f(n,j+1,1,k+1,m) = 0.
  200       continue
  210       continue
  220       continue
         endif
      endif
      if (jr.eq.(nvpy-1)) then
         call MPI_IRECV(scs(1,1,1,3,m),3*nxvz,mreal,kr-1,moff+2,lgrp,msi
     1d,ierr)
      endif
      if (jr.eq.nvpy) then
         if (nyzp(1,m).eq.1) then
            call MPI_SEND(scs(1,1,1,1,m),3*nxvzs,mreal,kl-1,moff+2,lgrp,
     1ierr)
         else
            nps = 0
            call MPI_SEND(scs(1,1,1,1,m),nps,mreal,kl-1,moff+2,lgrp,ierr
     1)
         endif
      endif
      if (jr.eq.(nvpy-1)) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            do 250 k = 1, kzp1
            do 240 j = 2, nx
            do 230 n = 1, 3
            f(n,j+1,nyzp(1,m)+1,k+1,m) = f(n,j+1,nyzp(1,m)+1,k+1,m) - sc
     1s(n,j+1,k+1,3,m)
            f(n,j+1,1,k+1,m) = 0.
  230       continue
  240       continue
  250       continue
         endif
      endif
  260 if (jl.eq.0) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+5,lgrp,msid,ierr)
         call MPI_IRECV(scs(1,1,1,3,m),3*nxvz,mreal,kl-1,moff+1,lgrp,nsi
     1d,ierr)
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
         call MPI_SEND(nyzp(1,m),1,mint,kr-1,moff+5,lgrp,ierr)
         if (nyzp(1,m).le.2) then
            call MPI_SEND(scs(1,1,1,2,m),3*nxvzs,mreal,kr-1,moff+1,lgrp,
     1ierr)
            do 290 k = 1, kzp1
            do 280 j = 2, nx
            do 270 n = 1, 3
            f(n,j+1,2,k+1,m) = 0.
  270       continue
  280       continue
  290       continue
         else
            nps = 0
            call MPI_SEND(scs(1,1,1,2,m),nps,mreal,kr-1,moff+1,lgrp,ierr
     1)
         endif
      endif
      if (jl.eq.0) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            if (ngc.eq.1) then
               if (nyzp(1,m).gt.1) then
                  do 320 k = 1, kzp1
                  do 310 j = 2, nx
                  do 300 n = 1, 3
                  f(n,j+1,3,k+1,m) = f(n,j+1,3,k+1,m) - scs(n,j+1,k+1,3,
     1m)
  300             continue
  310             continue
  320             continue
               endif
               do 350 k = 1, kzp1
               do 340 j = 2, nx
               do 330 n = 1, 3
               f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) + 2.*scs(n,j+1,k+1,3,
     1m)
               f(n,j+1,1,k+1,m) = 0.
  330          continue
  340          continue
  350          continue
            else
               do 380 k = 1, kzp1
               do 370 j = 2, nx
               do 360 n = 1, 3
               f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) - scs(n,j+1,k+1,3,m)
               f(n,j+1,1,k+1,m) = 0.
  360          continue
  370          continue
  380          continue
            endif
         endif
      endif
  390 continue
  400 continue
c add guard cells in z
      do 760 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 750 my = 1, mblok
      m = my + moff
      kyp1 = nyzp(1,m)
      nxvys = nxv*(kyp1 + 1)
      ky = my + js + 1
      kz = mz + ks
      jr = kz + 1
      jrr = jr + 1
      krr = ky + nvpy*jrr
      jl = kz - 1
      jll = jl - 1
      kll = ky + nvpy*jll
      kr = ky + nvpy*jr
      kl = ky + nvpy*jl
c fix edges if all points are on the same processor
      if (jl.eq.(-1)) then
         if (nyzp(2,m).gt.2) then
            do 430 k = 1, kyp1
            do 420 j = 2, nx
            do 410 n = 1, 3
            f(n,j+1,k+1,3,m) = f(n,j+1,k+1,3,m) + 2.*f(n,j+1,k+1,2,m)
            f(n,j+1,k+1,4,m) = f(n,j+1,k+1,4,m) - f(n,j+1,k+1,2,m)
            f(n,j+1,k+1,2,m) = 0.
  410       continue
  420       continue
  430       continue
         else if (nyzp(2,m).eq.2) then
            do 460 k = 1, kyp1
            do 450 j = 2, nx
            do 440 n = 1, 3
            f(n,j+1,k+1,3,m) = f(n,j+1,k+1,3,m) + 2.*f(n,j+1,k+1,2,m)
  440       continue
  450       continue
  460       continue
         endif
      endif
      if (jr.eq.nvpz) then
         if (nyzp(2,m).gt.1) then
            do 490 k = 1, kyp1
            do 480 j = 2, nx
            do 470 n = 1, 3
            f(n,j+1,k+1,nyzp(2,m),m) = f(n,j+1,k+1,nyzp(2,m),m) - f(n,j+
     11,k+1,nyzp(2,m)+2,m)
            f(n,j+1,k+1,nyzp(2,m)+1,m) = f(n,j+1,k+1,nyzp(2,m)+1,m) + 2.
     1*f(n,j+1,k+1,nyzp(2,m)+2,m)
            f(n,j+1,k+1,nyzp(2,m)+2,m) = 0.
  470       continue
  480       continue
  490       continue
         else if (nyzp(2,m).eq.1) then
            do 520 k = 1, kyp1
            do 510 j = 2, nx
            do 500 n = 1, 3
            f(n,j+1,k+1,nyzp(2,m)+1,m) = f(n,j+1,k+1,nyzp(2,m)+1,m) + 2.
     1*f(n,j+1,k+1,nyzp(2,m)+2,m)
  500       continue
  510       continue
  520       continue
         endif
      endif
c this segment is used for shared memory computers
c     if (mter.ge.2) go to 750
c     if (mter.eq.1) go to 620
c     if (jll.eq.0) then
c        if (nyzp(2,kll).eq.1).and.(nyzp(2,kl).eq.1)) then
c           do 550 k = 1, kyp1
c           do 540 j = 2, nx
c           do 530 n = 1, 3
c           f(n,j+1,k+1,2,m) = f(n,j+1,k+1,2,m) - f(n,j+1,k+1,2,kll)
c 530       continue
c 540       continue
c 550       continue
c        endif
c     endif
c     if (jr.eq.(nvpz-1)) then
c        if (nyzp(2,kr).eq.1) then
c           do 610 k = 1, kyp1
c           do 600 j = 2, nx
c           do 590 n = 1, 3
c           f(n,j+1,k+1,nyzp(2,m)+1,m) = f(n,j+1,k+1,nyzp(2,m)+1,m) - f(
c    1n,j+1,k+1,nyzp(2,kr)+2,kr)
c 590       continue
c 600       continue
c 610       continue
c        endif
c     endif
c 620 if (jl.eq.0) then
c        if (nyzp(kl).le.2) then
c           if (nzp(2,kl).eq.1) then
c              if (nyzp(2,m).gt.1) then
c                 do 680 k = 1, kyp1
c                 do 670 j = 2, nx
c                 do 660 n = 1, 3
c                 f(n,j+1,k+1,3,m) = f(n,j+1,k+1,3,m) - f(n,j+1,k+1,2,kl
c    1)
c 660             continue
c 670             continue
c 680             continue
c              endif
c              do 710 k = 1, kyp1
c              do 700 j = 2, nx
c              do 690 n = 1, 3
c              f(n,j+1,k+1,2,m) = f(n,j+1,k+1,2,m) + 2.*f(n,j+1,k+1,2,kl)
c 690          continue
c 700          continue
c 710          continue
c           else
c              do 740 k = 1, kyp1
c              do 730 j = 2, nx
c              do 720 n = 1, 3
c              f(n,j+1,k+1,2,m) = f(n,j+1,k+1,2,m) - f(n,j+1,k+1,2,kl)
c 720          continue
c 730          continue
c 740          continue
c           endif
c        endif
c     endif
c this segment is used for mpi computers
      if (nter.ge.2) go to 750
      if (nter.eq.1) go to 620
      if (jll.eq.0) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+4,lgrp,msid,ierr)
         call MPI_IRECV(f(1,1,1,1,m),3*nxvy,mreal,kll-1,moff+1,lgrp,nsid
     1,ierr)
      endif
      if (jl.eq.0) then
         call MPI_SEND(nyzp(2,m),1,mint,kr-1,moff+4,lgrp,ierr)
      endif
      if (jl.eq.(-1)) then
         if (nyzp(2,m).eq.1) then
            call MPI_SEND(f(1,1,1,2,m),3*nxvys,mreal,krr-1,moff+1,lgrp,i
     1err)
         else
            nps = 0
            call MPI_SEND(f(1,1,1,2,m),nps,mreal,krr-1,moff+1,lgrp,ierr)
         endif
      endif
      if (jll.eq.0) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            if (ngc.eq.1) then
               do 550 k = 1, kyp1
               do 540 j = 2, nx
               do 530 n = 1, 3
               f(n,j+1,k+1,2,m) = f(n,j+1,k+1,2,m) - f(n,j+1,k+1,1,m)
  530          continue
  540          continue
  550          continue
            endif
            do 580 k = 1, kyp1
            do 570 j = 2, nx
            do 560 n = 1, 3
            f(n,j+1,k+1,1,m) = 0.
  560       continue
  570       continue
  580       continue
         endif
      endif
      if (jr.eq.(nvpz-1)) then
         call MPI_IRECV(f(1,1,1,1,m),3*nxvy,mreal,kr-1,moff+2,lgrp,msid,
     1ierr)
      endif
      if (jr.eq.nvpz) then
         if (nyzp(2,m).eq.1) then
            call MPI_SEND(f(1,1,1,nyzp(2,m)+2,m),3*nxvys,mreal,kl-1,moff
     1+2,lgrp,ierr)
         else
            nps = 0
            call MPI_SEND(f(1,1,1,nyzp(2,m)+2,m),nps,mreal,kl-1,moff+2,l
     1grp,ierr)
         endif
      endif
      if (jr.eq.(nvpz-1)) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            do 610 k = 1, kyp1
            do 600 j = 2, nx
            do 590 n = 1, 3
            f(n,j+1,k+1,nyzp(2,m)+1,m) = f(n,j+1,k+1,nyzp(2,m)+1,m) - f(
     1n,j+1,k+1,1,m)
            f(n,j+1,k+1,1,m) = 0.
  590       continue
  600       continue
  610       continue
         endif
      endif
  620 if (jl.eq.0) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+5,lgrp,msid,ierr)
         call MPI_IRECV(f(1,1,1,1,m),3*nxvy,mreal,kl-1,moff+1,lgrp,nsid,
     1ierr)
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpz)) then
         call MPI_SEND(nyzp(2,m),1,mint,kr-1,moff+5,lgrp,ierr)
         if (nyzp(2,m).le.2) then
            call MPI_SEND(f(1,1,1,2,m),3*nxvys,mreal,kr-1,moff+1,lgrp,ie
     1rr)
            do 650 k = 1, kyp1
            do 640 j = 2, nx
            do 630 n = 1, 3
            f(n,j+1,k+1,2,m) = 0.
  630       continue
  640       continue
  650       continue
         else
            nps = 0
            call MPI_SEND(f(1,1,1,2,m),nps,mreal,kr-1,moff+1,lgrp,ierr)
         endif
      endif
      if (jl.eq.0) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            if (ngc.eq.1) then
               if (nyzp(2,m).gt.1) then
                  do 680 k = 1, kyp1
                  do 670 j = 2, nx
                  do 660 n = 1, 3
                  f(n,j+1,k+1,3,m) = f(n,j+1,k+1,3,m) - f(n,j+1,k+1,1,m)
  660             continue
  670             continue
  680             continue
               endif
               do 710 k = 1, kyp1
               do 700 j = 2, nx
               do 690 n = 1, 3
               f(n,j+1,k+1,2,m) = f(n,j+1,k+1,2,m) + 2.*f(n,j+1,k+1,1,m)
               f(n,j+1,k+1,1,m) = 0.
  690          continue
  700          continue
  710          continue
            else
               do 740 k = 1, kyp1
               do 730 j = 2, nx
               do 720 n = 1, 3
               f(n,j+1,k+1,2,m) = f(n,j+1,k+1,2,m) - f(n,j+1,k+1,1,m)
               f(n,j+1,k+1,1,m) = 0.
  720          continue
  730          continue
  740          continue
            endif
         endif
      endif
  750 continue
  760 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLAGUARDS32(f,scs,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypmx,nz
     1pmx,mblok,nblok,ngds,idds,mter,nter)
c this subroutine corrects current density data for particle boundary
c conditions which keep particles one grid away from the edges
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y and z direction.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes three extra guard cells.
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c it is assumed the nyzp(n,m) > 0.
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c mter/nter = (0,1) = (no,yes) pass data to next processor only in y/z
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      real f, scs
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds, mter, nter
      integer nyzp
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,ngds,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, nps
      integer mnblok, kyp1, kzp1, nxvz, nxvzs, nxvy, nxvys, m, my, mz
      integer j, k, jr, jrr, jl, jll
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
      kzp1 = nyzp(2,m) + 1
      do 20 k = 1, kzp1 + 1
      do 10 j = 1, nxv
      scs(j,k,1,m) = f(j,nyzp(1,m)+2,k,m)
      scs(j,k,2,m) = f(j,2,k,m)
   10 continue
   20 continue
   30 continue
c add guard cells in y
      do 280 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 270 my = 1, mblok
      m = my + moff
      kzp1 = nyzp(2,m) + 1
      nxvzs = nxv*(kzp1 + 1)
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr + 1
      krr = jrr + kz
      jl = ky - 1
      jll = jl - 1
      kll = jll + kz
      kr = jr + kz
      kl = jl + kz
c fix edges if all points are on the same processor
      if (jl.eq.(-1)) then
         if (nyzp(1,m).gt.2) then
            do 50 k = 1, kzp1
            do 40 j = 2, nx
            f(j+1,3,k+1,m) = f(j+1,3,k+1,m) + 2.*f(j+1,2,k+1,m)
            f(j+1,4,k+1,m) = f(j+1,4,k+1,m) - f(j+1,2,k+1,m)
            f(j+1,2,k+1,m) = 0.
   40       continue
   50       continue
         else if (nyzp(1,m).eq.2) then
            do 70 k = 1, kzp1
            do 60 j = 2, nx
            f(j+1,3,k+1,m) = f(j+1,3,k+1,m) + 2.*f(j+1,2,k+1,m)
   60       continue
   70       continue
         endif
      endif
      if (jr.eq.nvpy) then
         if (nyzp(1,m).gt.1) then
            do 90 k = 1, kzp1
            do 80 j = 2, nx
            f(j+1,nyzp(1,m),k+1,m) = f(j+1,nyzp(1,m),k+1,m) - f(j+1,nyzp
     1(1,m)+2,k+1,m)
            f(j+1,nyzp(1,m)+1,k+1,m) = f(j+1,nyzp(1,m)+1,k+1,m) + 2.*f(j
     1+1,nyzp(1,m)+2,k+1,m)
            f(j+1,nyzp(1,m)+2,k+1,m) = 0.
   80       continue
   90       continue
         else if (nyzp(1,m).eq.1) then
            do 110 k = 1, kzp1
            do 100 j = 2, nx
            f(j+1,nyzp(1,m)+1,k+1,m) = f(j+1,nyzp(1,m)+1,k+1,m) + 2.*f(j
     1+1,nyzp(1,m)+2,k+1,m)
            f(j+1,nyzp(1,m)+2,k+1,m) = 0.
  100       continue
  110       continue
         endif
      endif
c this segment is used for shared memory computers
c     if (mter.ge.2) go to 270
c     if (mter.eq.1) go to 180
c     if (jll.eq.0) then
c        if (nyzp(1,kll).eq.1).and.(nyzp(1,kl).eq.1)) then
c           do 130 k = 1, kzp1
c           do 120 j = 2, nx
c           f(j+1,2,k+1,m) = f(j+1,2,k+1,m) - scs(j+1,k+1,2,kll)
c 120       continue
c 130       continue
c        endif
c     endif
c     if (jr.eq.(nvpy-1)) then
c        if (nyzp(1,kr).eq.1) then
c           do 170 k = 1, kzp1
c           do 160 j = 2, nx
c           f(j+1,nyzp(1,m)+1,k+1,m) = f(j+1,nyzp(1,m)+1,k+1,m) - scs(j+
c    11,k+1,1,kr)
c 160       continue
c 170       continue
c        endif
c     endif
c 180 if (jl.eq.0) then
c        if (nyzp(1,kl).le.2) then
c           if (nyzp(1,kl).eq.1) then
c              if (nyzp(1,m).gt.1) then
c                 do 220 k = 1, kzp1
c                 do 210 j = 2, nx
c                 f(j+1,3,k+1,m) = f(j+1,3,k+1,m) - scs(j+1,k+1,2,kl)
c 210             continue
c 220             continue
c              endif
c              do 240 k = 1, kzp1
c              do 230 j = 2, nx
c              f(j+1,2,k+1,m) = f(j+1,2,k+1,m) + 2.*scs(j+1,k+1,2,kl)
c 230          continue
c 240          continue
c           else
c              do 260 k = 1, kzp1
c              do 250 j = 2, nx
c              f(j+1,2,k+1,m) = f(j+1,2,k+1,m) - scs(j+1,k+1,2,kl)
c 250          continue
c 260          continue
c           endif
c        endif
c     endif
c this segment is used for mpi computers
      if (mter.ge.2) go to 270
      if (mter.eq.1) go to 180
      if (jll.eq.0) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+4,lgrp,msid,ierr)
         call MPI_IRECV(scs(1,1,3,m),nxvz,mreal,kll-1,moff+1,lgrp,nsid,i
     1err)
      endif
      if (jl.eq.0) then
         call MPI_SEND(nyzp(1,m),1,mint,kr-1,moff+4,lgrp,ierr)
      endif
      if (jl.eq.(-1)) then
         if (nyzp(1,m).eq.1) then
            call MPI_SEND(scs(1,1,2,m),nxvzs,mreal,krr-1,moff+1,lgrp,ier
     1r)
         else
            nps = 0
            call MPI_SEND(scs(1,1,2,m),nps,mreal,krr-1,moff+1,lgrp,ierr)
         endif
      endif
      if (jll.eq.0) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            if (ngc.eq.1) then
               do 130 k = 1, kzp1
               do 120 j = 2, nx
               f(j+1,2,k+1,m) = f(j+1,2,k+1,m) - scs(j+1,k+1,3,m)
  120          continue
  130          continue
            endif
            do 150 k = 1, kzp1
            do 140 j = 2, nx
            f(j+1,1,k+1,m) = 0.
  140       continue
  150       continue
         endif
      endif
      if (jr.eq.(nvpy-1)) then
         call MPI_IRECV(scs(1,1,3,m),nxvz,mreal,kr-1,moff+2,lgrp,msid,ie
     1rr)
      endif
      if (jr.eq.nvpy) then
         if (nyzp(1,m).eq.1) then
            call MPI_SEND(scs(1,1,1,m),nxvzs,mreal,kl-1,moff+2,lgrp,ierr
     1)
         else
            nps = 0
            call MPI_SEND(scs(1,1,1,m),nps,mreal,kl-1,moff+2,lgrp,ierr)
         endif
      endif
      if (jr.eq.(nvpy-1)) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            do 170 k = 1, kzp1
            do 160 j = 2, nx
            f(j+1,nyzp(1,m)+1,k+1,m) = f(j+1,nyzp(1,m)+1,k+1,m) - scs(j+
     11,k+1,3,m)
            f(j+1,1,k+1,m) = 0.
  160       continue
  170       continue
         endif
      endif
  180 if (jl.eq.0) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+5,lgrp,msid,ierr)
         call MPI_IRECV(scs(1,1,3,m),nxvz,mreal,kl-1,moff+1,lgrp,nsid,ie
     1rr)
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
         call MPI_SEND(nyzp(1,m),1,mint,kr-1,moff+5,lgrp,ierr)
         if (nyzp(1,m).le.2) then
            call MPI_SEND(scs(1,1,2,m),nxvzs,mreal,kr-1,moff+1,lgrp,ierr
     1)
            do 200 k = 1, kzp1
            do 190 j = 2, nx
            f(j+1,2,k+1,m) = 0.
  190       continue
  200       continue
         else
            nps = 0
            call MPI_SEND(scs(1,1,2,m),nps,mreal,kr-1,moff+1,lgrp,ierr)
         endif
      endif
      if (jl.eq.0) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            if (ngc.eq.1) then
               if (nyzp(1,m).gt.1) then
                  do 220 k = 1, kzp1
                  do 210 j = 2, nx
                  f(j+1,3,k+1,m) = f(j+1,3,k+1,m) - scs(j+1,k+1,3,m)
  210             continue
  220             continue
               endif
               do 240 k = 1, kzp1
               do 230 j = 2, nx
               f(j+1,2,k+1,m) = f(j+1,2,k+1,m) + 2.*scs(j+1,k+1,3,m)
               f(j+1,1,k+1,m) = 0.
  230          continue
  240          continue
            else
               do 260 k = 1, kzp1
               do 250 j = 2, nx
               f(j+1,2,k+1,m) = f(j+1,2,k+1,m) - scs(j+1,k+1,3,m)
               f(j+1,1,k+1,m) = 0.
  250          continue
  260          continue
            endif
         endif
      endif
  270 continue
  280 continue
c add guard cells in z
      do 530 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 520 my = 1, mblok
      m = my + moff
      kyp1 = nyzp(1,m)
      nxvys = nxv*(kyp1 + 1)
      ky = my + js + 1
      kz = mz + ks
      jr = kz + 1
      jrr = jr + 1
      krr = ky + nvpy*jrr
      jl = kz - 1
      jll = jl - 1
      kll = ky + nvpy*jll
      kr = ky + nvpy*jr
      kl = ky + nvpy*jl
c fix edges if all points are on the same processor
      if (jl.eq.(-1)) then
         if (nyzp(2,m).gt.2) then
            do 300 k = 1, kyp1
            do 290 j = 2, nx
            f(j+1,k+1,3,m) = f(j+1,k+1,3,m) + 2.*f(j+1,k+1,2,m)
            f(j+1,k+1,4,m) = f(j+1,k+1,4,m) - f(j+1,k+1,2,m)
            f(j+1,k+1,2,m) = 0.
  290       continue
  300       continue
         else if (nyzp(2,m).eq.2) then
            do 320 k = 1, kyp1
            do 310 j = 2, nx
            f(j+1,k+1,3,m) = f(j+1,k+1,3,m) + 2.*f(j+1,k+1,2,m)
  310       continue
  320       continue
         endif
      endif
      if (jr.eq.nvpz) then
         if (nyzp(2,m).gt.1) then
            do 340 k = 1, kyp1
            do 330 j = 2, nx
            f(j+1,k+1,nyzp(2,m),m) = f(j+1,k+1,nyzp(2,m),m) - f(j+1,k+1,
     1nyzp(2,m)+2,m)
            f(j+1,k+1,nyzp(2,m)+1,m) = f(j+1,k+1,nyzp(2,m)+1,m) + 2.*f(j
     1+1,k+1,nyzp(2,m)+2,m)
            f(j+1,k+1,nyzp(2,m)+2,m) = 0.
  330       continue
  340       continue
         else if (nyzp(2,m).eq.1) then
            do 360 k = 1, kyp1
            do 350 j = 2, nx
            f(j+1,k+1,nyzp(2,m)+1,m) = f(j+1,k+1,nyzp(2,m)+1,m) + 2.*f(j
     1+1,k+1,nyzp(2,m)+2,m)
  350       continue
  360       continue
         endif
      endif
c this segment is used for shared memory computers
c     if (mter.ge.2) go to 520
c     if (mter.eq.1) go to 430
c     if (jll.eq.0) then
c        if (nyzp(2,kll).eq.1).and.(nyzp(2,kl).eq.1)) then
c           do 380 k = 1, kyp1
c           do 370 j = 2, nx
c           f(j+1,k+1,2,m) = f(j+1,k+1,2,m) - f(j+1,k+1,2,kll)
c 370       continue
c 380       continue
c        endif
c     endif
c     if (jr.eq.(nvpz-1)) then
c        if (nyzp(2,kr).eq.1) then
c           do 420 k = 1, kyp1
c           do 410 j = 2, nx
c           f(j+1,k+1,nyzp(2,m)+1,m) = f(j+1,k+1,nyzp(2,m)+1,m) - f(j+1,
c    1k+1,nyzp(2,kr)+2,kr)
c 410       continue
c 420       continue
c        endif
c     endif
c 430 if (jl.eq.0) then
c        if (nyzp(kl).le.2) then
c           if (nzp(2,kl).eq.1) then
c              if (nyzp(2,m).gt.1) then
c                 do 470 k = 1, kyp1
c                 do 460 j = 2, nx
c                 f(j+1,k+1,3,m) = f(j+1,k+1,3,m) - f(j+1,k+1,2,kl)
c 460             continue
c 470             continue
c              endif
c              do 490 k = 1, kyp1
c              do 480 j = 2, nx
c              f(j+1,k+1,2,m) = f(j+1,k+1,2,m) + 2.*f(j+1,k+1,2,kl)
c 480          continue
c 490          continue
c           else
c              do 510 k = 1, kyp1
c              do 500 j = 2, nx
c              f(j+1,k+1,2,m) = f(j+1,k+1,2,m) - f(j+1,k+1,2,kl)
c 500          continue
c 510          continue
c           endif
c        endif
c     endif
c this segment is used for mpi computers
      if (nter.ge.2) go to 520
      if (nter.eq.1) go to 430
      if (jll.eq.0) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+4,lgrp,msid,ierr)
         call MPI_IRECV(f(1,1,1,m),nxvy,mreal,kll-1,moff+1,lgrp,nsid,ier
     1r)
      endif
      if (jl.eq.0) then
         call MPI_SEND(nyzp(2,m),1,mint,kr-1,moff+4,lgrp,ierr)
      endif
      if (jl.eq.(-1)) then
         if (nyzp(2,m).eq.1) then
            call MPI_SEND(f(1,1,2,m),nxvys,mreal,krr-1,moff+1,lgrp,ierr)
         else
            nps = 0
            call MPI_SEND(f(1,1,2,m),nps,mreal,krr-1,moff+1,lgrp,ierr)
         endif
      endif
      if (jll.eq.0) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            if (ngc.eq.1) then
               do 380 k = 1, kyp1
               do 370 j = 2, nx
               f(j+1,k+1,2,m) = f(j+1,k+1,2,m) - f(j+1,k+1,1,m)
  370          continue
  380          continue
            endif
            do 400 k = 1, kyp1
            do 390 j = 2, nx
            f(j+1,k+1,1,m) = 0.
  390       continue
  400       continue
         endif
      endif
      if (jr.eq.(nvpz-1)) then
         call MPI_IRECV(f(1,1,1,m),nxvy,mreal,kr-1,moff+2,lgrp,msid,ierr
     1)
      endif
      if (jr.eq.nvpz) then
         if (nyzp(2,m).eq.1) then
            call MPI_SEND(f(1,1,nyzp(2,m)+2,m),nxvys,mreal,kl-1,moff+2,l
     1grp,ierr)
         else
            nps = 0
            call MPI_SEND(f(1,1,nyzp(2,m)+2,m),nps,mreal,kl-1,moff+2,lgr
     1p,ierr)
         endif
      endif
      if (jr.eq.(nvpz-1)) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            do 420 k = 1, kyp1
            do 410 j = 2, nx
            f(j+1,k+1,nyzp(2,m)+1,m) = f(j+1,k+1,nyzp(2,m)+1,m) - f(j+1,
     1k+1,1,m)
            f(j+1,k+1,1,m) = 0.
  410       continue
  420       continue
         endif
      endif
  430 if (jl.eq.0) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+5,lgrp,msid,ierr)
         call MPI_IRECV(f(1,1,1,m),nxvy,mreal,kl-1,moff+1,lgrp,nsid,ierr
     1)
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpz)) then
         call MPI_SEND(nyzp(2,m),1,mint,kr-1,moff+5,lgrp,ierr)
         if (nyzp(2,m).le.2) then
            call MPI_SEND(f(1,1,2,m),nxvys,mreal,kr-1,moff+1,lgrp,ierr)
            do 450 k = 1, kyp1
            do 440 j = 2, nx
            f(j+1,k+1,2,m) = 0.
  440       continue
  450       continue
         else
            nps = 0
            call MPI_SEND(f(1,1,2,m),nps,mreal,kr-1,moff+1,lgrp,ierr)
         endif
      endif
      if (jl.eq.0) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            if (ngc.eq.1) then
               if (nyzp(2,m).gt.1) then
                  do 470 k = 1, kyp1
                  do 460 j = 2, nx
                  f(j+1,k+1,3,m) = f(j+1,k+1,3,m) - f(j+1,k+1,1,m)
  460             continue
  470             continue
               endif
               do 490 k = 1, kyp1
               do 480 j = 2, nx
               f(j+1,k+1,2,m) = f(j+1,k+1,2,m) + 2.*f(j+1,k+1,1,m)
               f(j+1,k+1,1,m) = 0.
  480          continue
  490          continue
            else
               do 510 k = 1, kyp1
               do 500 j = 2, nx
               f(j+1,k+1,2,m) = f(j+1,k+1,2,m) - f(j+1,k+1,1,m)
               f(j+1,k+1,1,m) = 0.
  500          continue
  510          continue
            endif
         endif
      endif
  520 continue
  530 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLACGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzp
     1mx,mblok,nblok,kyp,kzp,ngds)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c no copying is done at the boundary edges.
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes one extra guard
c cell.
c scs/scr = scratch array for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c linear interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, kyp, kzp
      real f, scs, scr
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(3,nxv,nypmx,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer ky, kz, js, ks, moff, noff, jr, jl, kr, kl, mnblok
      integer nx1, nyp1, nzp1, nxvz, nxvy, m, my, mz, j, k, n
      dimension istatus(lstat)
      nx1 = nx + 1
      nyp1 = kyp + 1
      nzp1 = kzp + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c buffer data in y
      do 40 m = 1, mnblok
      do 30 k = 1, nzpmx
      do 20 j = 1, nxv
      do 10 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,kyp+1,k,m)
   10 continue
   20 continue
   30 continue
   40 continue
c add guard cells in y
      do 150 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 140 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jl = ky - 1
      kr = jr + kz
      kl = jl + kz
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 70 k = 1, nzpmx
c        do 60 j = 1, nxv
c        do 50 n = 1, 3
c        scs(n,j,k,2,m) = scs(n,j,k,1,kl)
c  50    continue
c  60    continue
c  70    continue
c     else
c        do 100 k = 1, nzp1
c        do 90 j = 1, nx1
c        do 80 n = 1, 3
c        scs(n,j,k,2,m) = 0.
c  80    continue
c  90    continue
c 100    continue
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,1,2,m),3*nxvz,mreal,kl-1,noff+1,lgrp,msi
     1d,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,1,m),3*nxvz,mreal,kr-1,noff+1,lgrp,ierr
     1)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 70 k = 1, nzp1
         do 60 j = 1, nx1
         do 50 n = 1, 3
         scs(n,j,k,2,m) = 0.
   50    continue
   60    continue
   70    continue
      endif
c add up the guard cells
      do 130 k = 1, nzp1
      do 120 j = 1, nx1
      do 110 n = 1, 3
      f(n,j,1,k,m) = f(n,j,1,k,m) + scs(n,j,k,2,m)
  110 continue
  120 continue
  130 continue
  140 continue
  150 continue
c add guard cells in z
      do 260 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 250 my = 1, mblok
      m = my + moff
      ky = my + js + 1
      kz = mz + ks
      jr = kz + 1
      jl = kz - 1
      kr = ky + nvpy*jr
      kl = ky + nvpy*jl
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 180 k = 1, nypmx
c        do 170 j = 1, nxv
c        do 160 n = 1, 3
c        scr(n,j,k,m) = f(n,j,k,kzp+1,kl)
c 160    continue
c 170    continue
c 180    continue
c     else
c        do 210 k = 1, nyp1
c        do 200 j = 1, nx1
c        do 190 n = 1, 3
c        scr(n,j,k,m) = 0.
c 190    continue
c 200    continue
c 210    continue
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scr,3*nxvy,mreal,kl-1,noff+2,lgrp,msid,ierr)
      endif
      if (jr.lt.nvpz) then
         call MPI_SEND(f(1,1,1,kzp+1,m),3*nxvy,mreal,kr-1,noff+2,lgrp,ie
     1rr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 180 k = 1, nyp1
         do 170 j = 1, nx1
         do 160 n = 1, 3
         scr(n,j,k,m) = 0.
  160    continue
  170    continue
  180    continue
      endif
c add up the guard cells
      do 240 k = 1, nyp1
      do 230 j = 1, nx1
      do 220 n = 1, 3
      f(n,j,k,1,m) = f(n,j,k,1,m) + scr(n,j,k,m)
  220 continue
  230 continue
  240 continue
  250 continue
  260 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLAGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpm
     1x,mblok,nblok,kyp,kzp,ngds)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c no copying is done at the boundary edges.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes one extra guard
c cell.
c scs/scr = scratch array for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c linear interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, kyp, kzp
      real f, scs, scr
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(nxv,nypmx,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer ky, kz, js, ks, moff, noff, jr, jl, kr, kl, mnblok
      integer nx1, nyp1, nzp1, nxvz, nxvy, m, my, mz, j, k
      dimension istatus(lstat)
      nx1 = nx + 1
      nyp1 = kyp + 1
      nzp1 = kzp + 1
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
      do 110 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 100 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jl = ky - 1
      kr = jr + kz
      kl = jl + kz
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 50 k = 1, nzpmx
c        do 40 j = 1, nxv
c        scs(j,k,2,m) = scs(j,k,1,kl)
c  40    continue
c  50    continue
c     else
c        do 70 k = 1, nzp1
c        do 60 j = 1, nx1
c        scs(j,k,2,m) = 0.
c  60    continue
c  70    continue
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,2,m),nxvz,mreal,kl-1,noff+1,lgrp,msid,ie
     1rr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,m),nxvz,mreal,kr-1,noff+1,lgrp,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 50 k = 1, nzp1
         do 40 j = 1, nx1
         scs(j,k,2,m) = 0.
   40    continue
   50    continue
      endif
c add up the guard cells
      do 90 k = 1, nzp1
      do 80 j = 1, nx1
      f(j,1,k,m) = f(j,1,k,m) + scs(j,k,2,m)
   80 continue
   90 continue
  100 continue
  110 continue
c add guard cells in z
      do 190 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 180 my = 1, mblok
      m = my + moff
      ky = my + js + 1
      kz = mz + ks
      jr = kz + 1
      jl = kz - 1
      kr = ky + nvpy*jr
      kl = ky + nvpy*jl
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 130 k = 1, nypmx
c        do 120 j = 1, nxv
c        scr(j,k,m) = f(j,k,kzp+1,kl)
c 120    continue
c 130    continue
c     else
c        do 150 k = 1, nyp1
c        do 140 j = 1, nx1
c        scr(j,k,m) = 0.
c 140    continue
c 150    continue
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scr,nxvy,mreal,kl-1,noff+2,lgrp,msid,ierr)
      endif
      if (jr.lt.nvpz) then
         call MPI_SEND(f(1,1,kzp+1,m),nxvy,mreal,kr-1,noff+2,lgrp,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 130 k = 1, nyp1
         do 120 j = 1, nx1
         scr(j,k,m) = 0.
  120    continue
  130    continue
      endif
c add up the guard cells
      do 170 k = 1, nyp1
      do 160 j = 1, nx1
      f(j,k,1,m) = f(j,k,1,m) + scr(j,k,m)
  160 continue
  170 continue
  180 continue
  190 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLACGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nyp
     1mx,nzpmx,mblok,nblok,ngds,idds)
c this subroutine adds data from guard cells in non-uniform partitions
c for vector data.  no copying is done at the boundary edges.
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes one extra guard cell.
c scs/scr = scratch arrays for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c linear interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds
      integer nyzp
      real f, scs, scr
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(3,nxv,nypmx,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer ky, kz, js, ks, moff, noff, jr, jl, kr, kl, mnblok
      integer nx1, nxvz, nxvzs, nyzp1, nxvy, nxvys, m, my, mz, j, k, n
      dimension istatus(lstat)
      nx1 = nx + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c special case for one processor in y
      if (nvpy.eq.1) go to 190
c buffer data in y
      do 40 m = 1, mnblok
      do 30 k = 1, nyzp(2,m)+1
      do 20 j = 1, nxv
      do 10 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,nyzp(1,m)+1,k,m)
   10 continue
   20 continue
   30 continue
   40 continue
c add guard cells in y
      do 180 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 170 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(2,m) + 1
      nxvzs = nxv*nyzp1
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jl = ky - 1
      kr = jr + kz
      kl = jl + kz
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 70 k = 1, nyzp1
c        do 60 j = 1, nxv
c        do 50 n = 1, 3
c        scs(n,j,k,2,m) = scs(n,j,k,1,kl)
c  50    continue
c  60    continue
c  70    continue
c     else
c        do 100 k = 1, nyzp1
c        do 90 j = 1, nxv
c        do 80 n = 1, 3
c        scs(n,j,k,2,m) = 0.
c  80    continue
c  90    continue
c 100    continue
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,1,2,m),3*nxvz,mreal,kl-1,noff+1,lgrp,msi
     1d,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,1,m),3*nxvzs,mreal,kr-1,noff+1,lgrp,ier
     1r)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 100 k = 1, nyzp1
         do 90 j = 1, nxv
         do 80 n = 1, 3
         scs(n,j,k,2,m) = 0.
   80    continue
   90    continue
  100    continue
      endif
c add up the guard cells
      do 130 k = 1, nyzp1
      do 120 j = 1, nx1
      do 110 n = 1, 3
      f(n,j,1,k,m) = f(n,j,1,k,m) + scs(n,j,k,2,m)
  110 continue
  120 continue
  130 continue
      if (jr.lt.nvpy) then
         do 160 k = 1, nyzp1
         do 150 j = 1, nx1
         do 140 n = 1, 3
         f(n,j,nyzp(1,m)+1,k,m) = 0.
  140    continue
  150    continue
  160    continue
      endif
  170 continue
  180 continue
c special case for one processor in z
  190 if (nvpz.eq.1) return
c add guard cells in z
      do 330 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 320 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(1,m) + 1
      nxvys = nxv*nyzp1
      ky = my + js + 1
      kz = mz + ks
      jr = kz + 1
      jl = kz - 1
      kr = ky + nvpy*jr
      kl = ky + nvpy*jl
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 220 k = 1, nyzp1
c        do 210 j = 1, nxv
c        do 200 n = 1, 3
c        scr(n,j,k,m) = f(n,j,k,nyzp(2,m)+1,kl)
c 200    continue
c 210    continue
c 220    continue
c     else
c        do 250 k = 1, nyzp1
c        do 240 j = 1, nxv
c        do 230 n = 1, 3
c        scr(n,j,k,m) = 0.
c 230    continue
c 240    continue
c 250    continue
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scr,3*nxvy,mreal,kl-1,noff+2,lgrp,msid,ierr)
      endif
      if (jr.lt.nvpz) then
         call MPI_SEND(f(1,1,1,nyzp(2,m)+1,m),3*nxvys,mreal,kr-1,noff+2,
     1lgrp,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 250 k = 1, nyzp1
         do 240 j = 1, nxv
         do 230 n = 1, 3
         scr(n,j,k,m) = 0.
  230    continue
  240    continue
  250    continue
      endif
c add up the guard cells
      do 280 k = 1, nyzp1
      do 270 j = 1, nx1
      do 260 n = 1, 3
      f(n,j,k,1,m) = f(n,j,k,1,m) + scr(n,j,k,m)
  260 continue
  270 continue
  280 continue
      if (jr.lt.nvpz) then
         do 310 k = 1, nyzp1
         do 300 j = 1, nx1
         do 290 n = 1, 3
         f(n,j,k,nyzp(2,m)+1,m) = 0.
  290    continue
  300    continue
  310    continue
      endif
  320 continue
  330 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNLAGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypm
     1x,nzpmx,mblok,nblok,ngds,idds)
c this subroutine adds data from guard cells in non-uniform partitions
c for scalar data.  no copying is done at the boundary edges.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes one extra guard cell.
c scs/scr = scratch arrays for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c linear interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds
      integer nyzp
      real f, scs, scr
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(nxv,nypmx,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer ky, kz, js, ks, moff, noff, jr, jl, kr, kl, mnblok
      integer nx1, nxvz, nxvzs, nyzp1, nxvy, nxvys, m, my, mz, j, k
      dimension istatus(lstat)
      nx1 = nx + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c special case for one processor in y
      if (nvpy.eq.1) go to 140
c buffer data in y
      do 30 m = 1, mnblok
      do 20 k = 1, nyzp(2,m)+1
      do 10 j = 1, nxv
      scs(j,k,1,m) = f(j,nyzp(1,m)+1,k,m)
   10 continue
   20 continue
   30 continue
c add guard cells in y
      do 130 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 120 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(2,m) + 1
      nxvzs = nxv*nyzp1
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jl = ky - 1
      kr = jr + kz
      kl = jl + kz
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 50 k = 1, nyzp1
c        do 40 j = 1, nxv
c        scs(j,k,2,m) = scs(j,k,1,kl)
c  40    continue
c  50    continue
c     else
c        do 70 k = 1, nyzp1
c        do 60 j = 1, nxv
c        scs(j,k,2,m) = 0.
c  60    continue
c  70    continue
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,2,m),nxvz,mreal,kl-1,noff+1,lgrp,msid,ie
     1rr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,m),nxvzs,mreal,kr-1,noff+1,lgrp,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 70 k = 1, nyzp1
         do 60 j = 1, nxv
         scs(j,k,2,m) = 0.
   60    continue
   70    continue
      endif
c add up the guard cells
      do 90 k = 1, nyzp1
      do 80 j = 1, nx1
      f(j,1,k,m) = f(j,1,k,m) + scs(j,k,2,m)
   80 continue
   90 continue
      if (jr.lt.nvpy) then
         do 110 k = 1, nyzp1
         do 100 j = 1, nx1
         f(j,nyzp(1,m)+1,k,m) = 0.
  100    continue
  110    continue
      endif
  120 continue
  130 continue
c special case for one processor in z
  140 if (nvpz.eq.1) return
c add guard cells in z
      do 240 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 230 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(1,m) + 1
      nxvys = nxv*nyzp1
      ky = my + js + 1
      kz = mz + ks
      jr = kz + 1
      jl = kz - 1
      kr = ky + nvpy*jr
      kl = ky + nvpy*jl
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 160 k = 1, nyzp1
c        do 150 j = 1, nxv
c        scr(j,k,m) = f(j,k,nyzp(2,m)+1,kl)
c 150    continue
c 160    continue
c     else
c        do 180 k = 1, nyzp1
c        do 170 j = 1, nxv
c        scr(j,k,m) = 0.
c 170    continue
c 180    continue
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scr,nxvy,mreal,kl-1,noff+2,lgrp,msid,ierr)
      endif
      if (jr.lt.nvpz) then
         call MPI_SEND(f(1,1,nyzp(2,m)+1,m),nxvys,mreal,kr-1,noff+2,lgrp
     1,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 180 k = 1, nyzp1
         do 170 j = 1, nxv
         scr(j,k,m) = 0.
  170    continue
  180    continue
      endif
c add up the guard cells
      do 200 k = 1, nyzp1
      do 190 j = 1, nx1
      f(j,k,1,m) = f(j,k,1,m) + scr(j,k,m)
  190 continue
  200 continue
      if (jr.lt.nvpz) then
         do 220 k = 1, nyzp1
         do 210 j = 1, nx1
         f(j,k,nyzp(2,m)+1,m) = 0.
  210    continue
  220    continue
      endif
  230 continue
  240 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMCGUARD32(f,scs,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx,mbl
     1ok,nblok,kyp,kzp,ngds)
c this subroutine copies data from field to particle partitions, copying
c data to guard cells, where the field and particle partitions are 
c assumed to be the same.  for vector data
c the field is replicated so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y direction, and periodic in the z direction.
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes three extra
c guard cells.
c scs = scratch array for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer kyp, kzp, ngds
      real f, scs
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx,2*ngds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, mnblok
      integer nx2, nyp3, kzp1, nxvz, nxvy, m, my, mz, j, k, n
      integer jr, jrr, jl, jll
      dimension istatus(lstat)
      nx2 = nx + 2
      nyp3 = kyp + 3
      kzp1 = kzp + 1
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
      do 40 m = 1, mnblok
      do 30 k = 1, nzpmx
      do 20 j = 1, nxv
      do 10 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,kyp+1,k,m)
      scs(n,j,k,2,m) = f(n,j,2,k,m)
      scs(n,j,k,3,m) = f(n,j,ngc+1,k,m)
      scs(n,j,k,5,m) = f(n,j,kyp+2,k,m)
   10 continue
   20 continue
   30 continue
   40 continue
c copy to guard cells in y
      do 270 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 260 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr
      krr = jr + kz
      jl = ky - 1
      jll = jl
      kll = jl + kz
c special case of only one grid per processor
      if (kyp.eq.1) then
         jrr = jr + 1
         krr = jrr + kz
         jll = jl - 1
         kll = jll + kz
      endif
      kr = jr + kz
      kl = jl + kz
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 70 k = 1, nzpmx
c        do 60 j = 1, nxv
c        do 50 n = 1, 3
c        scs(n,j,k,4,m) = scs(n,j,k,1,kl)
c  50    continue
c  60    continue
c  70    continue
c     else
c        do 100 k = 1, nzpmx
c        do 90 j = 1, nxv
c        do 80 n = 1, 3
c        scs(n,j,k,4,m) = f(n,j,3,k,m)
c  80    continue
c  90    continue
c 100    continue
c     endif
c     if (jr.lt.nvpy) then
c        do 130 k = 1, nzpmx
c        do 120 j = 1, nxv
c        do 110 n = 1, 3
c        scs(n,j,k,5,m) = scs(n,j,k,2,kr)
c        scs(n,j,k,6,m) = scs(n,j,k,3,krr)
c 110    continue
c 120    continue
c 130    continue
c     else
c        do 160 k = 1, kzp1
c        do 150 j = 2, nx2
c        do 140 n = 1, 3
c        scs(n,j,k+1,6,m) = 2.*scs(n,j,k+1,5,m) - scs(n,j,k+1,1,m)
c 140    continue
c 150    continue
c 160    continue
c     endif
c     if (kyp.eq.1) then
c        if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
c           do 190 k = 1, nzpmx
c           do 180 j = 1, nxv
c           do 170 n = 1, 3
c           scs(n,j,k,4,m) = scs(n,j,k,2,kr)
c 170       continue
c 180       continue
c 190       continue
c        endif
c        if (jr.eq.(nvpy-1)) then
c           do 220 k = 1, nzpmx
c           do 210 j = 1, nxv
c           do 200 n = 1, 3
c           scs(n,j,k,6,m) = scs(n,j,k,5,kr)
c 200       continue
c 210       continue
c 220       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,1,4,m),3*nxvz,mreal,kl-1,noff+3,lgrp,msi
     1d,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,1,m),3*nxvz,mreal,kr-1,noff+3,lgrp,ierr
     1)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 70 k = 1, nzpmx
         do 60 j = 1, nxv
         do 50 n = 1, 3
         scs(n,j,k,4,m) = f(n,j,3,k,m)
   50    continue
   60    continue
   70    continue
      endif
      if (jr.lt.nvpy) then
         call MPI_IRECV(scs(1,1,1,5,m),3*ngc*nxvz,mreal,kr-1,noff+4,lgrp
     1,msid,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scs(1,1,1,2,m),3*ngc*nxvz,mreal,kl-1,noff+4,lgrp,
     1ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 100 k = 1, kzp1
         do 90 j = 2, nx2
         do 80 n = 1, 3
         scs(n,j,k+1,6,m) = 2.*scs(n,j,k+1,5,m) - scs(n,j,k+1,1,m)
   80    continue
   90    continue
  100    continue
      endif
c special case of only one grid per processor
      if (kyp.eq.1) then
         if (jrr.lt.nvpy) then
            call MPI_IRECV(scs(1,1,1,6,m),3*nxvz,mreal,krr-1,noff+5,lgrp
     1,msid,ierr)
         else if (jr.lt.nvpy) then
            call MPI_IRECV(scs(1,1,1,6,m),3*nxvz,mreal,kr-1,noff+5,lgrp,
     1msid,ierr)
         endif
         if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
            call MPI_IRECV(scs(1,1,1,4,m),3*nxvz,mreal,kr-1,noff+5,lgrp,
     1nsid,ierr)
         endif
         if (jll.ge.0) then
            call MPI_SEND(scs(1,1,1,3,m),3*nxvz,mreal,kll-1,noff+5,lgrp,
     1ierr)
         else if (jl.eq.0) then
            call MPI_SEND(scs(1,1,1,3,m),3*nxvz,mreal,kl-1,noff+5,lgrp,i
     1err)
         endif
         if ((jl.eq.(nvpy-2)).and.(jl.ge.0)) then
            call MPI_SEND(scs(1,1,1,5,m),3*nxvz,mreal,kl-1,noff+5,lgrp,i
     1err)
         endif
         if (jr.lt.nvpy) then
            call MPI_WAIT(msid,istatus,ierr)
         endif
         if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
            call MPI_WAIT(nsid,istatus,ierr)
         endif
      endif
      do 250 k = 1, nzpmx
      do 240 j = 1, nxv
      do 230 n = 1, 3
      f(n,j,1,k,m) = scs(n,j,k,4,m)
      f(n,j,kyp+2,k,m) = scs(n,j,k,5,m)
      f(n,j,kyp+3,k,m) = scs(n,j,k,6,m)
  230 continue
  240 continue
  250 continue
  260 continue
  270 continue
c fix left edge
      do 320 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 310 my = 1, mblok
      m = my + moff
      kl = my + js
      if (kl.eq.0) then
         do 300 k = 1, kzp1
         do 290 j = 2, nx2
         do 280 n = 1, 3
         f(n,j,1,k+1,m) = 2.*f(n,j,2,k+1,m) - f(n,j,1,k+1,m)
  280    continue
  290    continue
  300    continue
      endif
  310 continue
  320 continue
c copy to guard cells in z
      do 370 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 360 my = 1, mblok
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
c     do 350 k = 1, nypmx
c     do 340 j = 1, nxv
c     do 330 n = 1, 3
c     f(n,j,k,1,m) = f(n,j,k,kzp+1,kl)
c     f(n,j,k,kzp+2,m) = f(n,j,k,2,kr)
c     f(n,j,k,kzp+3,m) = f(n,j,k,ngc+1,krr)
c 330 continue
c 340 continue
c 350 continue
c this segment is used for mpi computers
      call MPI_IRECV(f(1,1,1,1,m),3*nxvy,mreal,kl-1,noff+6,lgrp,msid,ier
     1r)
      call MPI_SEND(f(1,1,1,kzp+1,m),3*nxvy,mreal,kr-1,noff+6,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(f(1,1,1,kzp+2,m),3*ngc*nxvy,mreal,kr-1,noff+8,lgrp,
     1msid,ierr)
      call MPI_SEND(f(1,1,1,2,m),3*ngc*nxvy,mreal,kl-1,noff+8,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      if (kzp.eq.1) then
         call MPI_IRECV(f(1,1,1,kzp+3,m),3*ngc*nxvy,mreal,krr-1,noff+9,l
     1grp,msid,ierr)
         call MPI_SEND(f(1,1,1,2,m),3*ngc*nxvy,mreal,kll-1,noff+9,lgrp,i
     1err)
         call MPI_WAIT(msid,istatus,ierr)
      endif
  360 continue
  370 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMDGUARD32(f,scs,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx,mbl
     1ok,nblok,kyp,kzp,ngds)
c this subroutine copies data from field to particle partitions, copying
c data to guard cells, where the field and particle partitions are 
c assumed to be the same.
c the field is replicated so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y direction, and periodic in the z direction.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes three extra
c guard cells.
c scs = scratch array for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer kyp, kzp, ngds
      real f, scs
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, mnblok
      integer nx2, nyp3, kzp1, nxvz, nxvy, m, my, mz, j, k
      integer jr, jrr, jl, jll
      dimension istatus(lstat)
      nx2 = nx + 2
      nyp3 = kyp + 3
      kzp1 = kzp + 1
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
      scs(j,k,5,m) = f(j,kyp+2,k,m)
   10 continue
   20 continue
   30 continue
c copy to guard cells in y
      do 190 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 180 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr
      krr = jr + kz
      jl = ky - 1
      jll = jl
      kll = jl + kz
c special case of only one grid per processor
      if (kyp.eq.1) then
         jrr = jr + 1
         krr = jrr + kz
         jll = jl - 1
         kll = jll + kz
      endif
      kr = jr + kz
      kl = jl + kz
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 50 k = 1, nzpmx
c        do 40 j = 1, nxv
c        scs(j,k,4,m) = scs(j,k,1,kl)
c  40    continue
c  50    continue
c     else
c        do 70 k = 1, nzpmx
c        do 60 j = 1, nxv
c        scs(j,k,4,m) = f(j,3,k,m)
c  60    continue
c  70    continue
c     endif
c     if (jr.lt.nvpy) then
c        do 90 k = 1, nzpmx
c        do 80 j = 1, nxv
c        scs(j,k,5,m) = scs(j,k,2,kr)
c        scs(j,k,6,m) = scs(j,k,3,krr)
c  80    continue
c  90    continue
c     else
c        do 110 k = 1, kzp1
c        do 100 j = 2, nx2
c        scs(j,k+1,6,m) = 2.*scs(j,k+1,5,m) - scs(j,k+1,1,m)
c 100    continue
c 110    continue
c     endif
c     if (kyp.eq.1) then
c        if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
c           do 130 k = 1, nzpmx
c           do 120 j = 1, nxv
c           scs(j,k,4,m) = scs(j,k,2,kr)
c 120       continue
c 130       continue
c        endif
c        if (jr.eq.(nvpy-1)) then
c           do 150 k = 1, nzpmx
c           do 140 j = 1, nxv
c           scs(j,k,6,m) = scs(j,k,5,kr)
c 140       continue
c 150       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,4,m),nxvz,mreal,kl-1,noff+3,lgrp,msid,ie
     1rr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,m),nxvz,mreal,kr-1,noff+3,lgrp,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 50 k = 1, nzpmx
         do 40 j = 1, nxv
         scs(j,k,4,m) = f(j,3,k,m)
   40    continue
   50    continue
      endif
      if (jr.lt.nvpy) then
         call MPI_IRECV(scs(1,1,5,m),ngc*nxvz,mreal,kr-1,noff+4,lgrp,msi
     1d,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scs(1,1,2,m),ngc*nxvz,mreal,kl-1,noff+4,lgrp,ierr
     1)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 70 k = 1, kzp1
         do 60 j = 2, nx2
         scs(j,k+1,6,m) = 2.*scs(j,k+1,5,m) - scs(j,k+1,1,m)
   60    continue
   70    continue
      endif
c special case of only one grid per processor
      if (kyp.eq.1) then
         if (jrr.lt.nvpy) then
            call MPI_IRECV(scs(1,1,6,m),nxvz,mreal,krr-1,noff+5,lgrp,msi
     1d,ierr)
         else if (jr.lt.nvpy) then
            call MPI_IRECV(scs(1,1,6,m),nxvz,mreal,kr-1,noff+5,lgrp,msid
     1,ierr)
         endif
         if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
            call MPI_IRECV(scs(1,1,4,m),nxvz,mreal,kr-1,noff+5,lgrp,nsid
     1,ierr)
         endif
         if (jll.ge.0) then
            call MPI_SEND(scs(1,1,3,m),nxvz,mreal,kll-1,noff+5,lgrp,ierr
     1)
         else if (jl.eq.0) then
            call MPI_SEND(scs(1,1,3,m),nxvz,mreal,kl-1,noff+5,lgrp,ierr)
         endif
         if ((jl.eq.(nvpy-2)).and.(jl.ge.0)) then
            call MPI_SEND(scs(1,1,5,m),nxvz,mreal,kl-1,noff+5,lgrp,ierr)
         endif
         if (jr.lt.nvpy) then
            call MPI_WAIT(msid,istatus,ierr)
         endif
         if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
            call MPI_WAIT(nsid,istatus,ierr)
         endif
      endif
      do 170 k = 1, nzpmx
      do 160 j = 1, nxv
      f(j,1,k,m) = scs(j,k,4,m)
      f(j,kyp+2,k,m) = scs(j,k,5,m)
      f(j,kyp+3,k,m) = scs(j,k,6,m)
  160 continue
  170 continue
  180 continue
  190 continue
c fix left edge
      do 230 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 220 my = 1, mblok
      m = my + moff
      kl = my + js
      if (kl.eq.0) then
         do 210 k = 1, kzp1
         do 200 j = 2, nx2
         f(j,1,k+1,m) = 2.*f(j,2,k+1,m) - f(j,1,k+1,m)
  200    continue
  210    continue
      endif
  220 continue
  230 continue
c copy to guard cells in z
      do 270 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 260 my = 1, mblok
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
c     do 250 k = 1, nypmx
c     do 240 j = 1, nxv
c     f(j,k,1,m) = f(j,k,kzp+1,kl)
c     f(j,k,kzp+2,m) = f(j,k,2,kr)
c     f(j,k,kzp+3,m) = f(j,k,ngc+1,krr)
c 240 continue
c 250 continue
c this segment is used for mpi computers
      call MPI_IRECV(f(1,1,1,m),nxvy,mreal,kl-1,noff+6,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,kzp+1,m),nxvy,mreal,kr-1,noff+6,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(f(1,1,kzp+2,m),ngc*nxvy,mreal,kr-1,noff+8,lgrp,msid
     1,ierr)
      call MPI_SEND(f(1,1,2,m),ngc*nxvy,mreal,kl-1,noff+8,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      if (kzp.eq.1) then
         call MPI_IRECV(f(1,1,kzp+3,m),ngc*nxvy,mreal,krr-1,noff+9,lgrp,
     1msid,ierr)
         call MPI_SEND(f(1,1,2,m),ngc*nxvy,mreal,kll-1,noff+9,lgrp,ierr)
         call MPI_WAIT(msid,istatus,ierr)
      endif
  260 continue
  270 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNMCGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypmx
     1,nzpmx,mblok,nblok,ngds,idds,mter,nter)
c this subroutine copies data to guard cells in non-uniform partitions
c for vector data the field is replicated so as to disable quadratic
c interpolation within half a cell of the edges, and reduce it to linear
c interpolation in the y direction, and periodic in the z direction.
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes three extra guard cells.
c scs/scr = scratch arrays for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c it is assumed that nyzp(n,m) > 0.
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c mter/nter = (0,1) = (no,yes) pass data to next processor only in y/z
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds, mter, nter
      integer nyzp
      real f, scs, scr
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx*ngds,2,mblok*nblok)
      dimension scr(3,nxv,nypmx*ngds,2,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, nps
      integer mnblok, nx2, nyzp1, nxvz, nxvzs, nyzp3, nxvy, nxvys
      integer m, my, mz, j, k, n, jr, jrr, jl, jll
      dimension istatus(lstat)
      nx2 = nx + 2
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c buffer data in y
      do 40 m = 1, mnblok
      ngc = 0
      if (nyzp(1,m).eq.1) ngc = 1
      nyzp1 = nyzp(2,m) + 1
      do 30 k = 1, nyzp1
      do 20 j = 1, nxv
      do 10 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,nyzp(1,m)+1,k+1,m)
      scs(n,j,k+nyzp1,1,m) = f(n,j,2,k+1,m)
      scs(n,j,k+2*nyzp1,1,m) = f(n,j,3-ngc,k+1,m)
      scs(n,j,k+nyzp1,2,m) = f(n,j,nyzp(1,m)+2,k+1,m)
   10 continue
   20 continue
   30 continue
   40 continue
c copy to guard cells in y
      do 310 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 300 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(2,m) + 1
      nxvzs = nxv*nyzp1
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr + 1
      krr = jrr + kz
      jl = ky - 1
      jll = jl - 1
      kll = jll + kz
      kr = jr + kz
      kl = jl + kz
      ngc = 0
c special case of only one grid per processor
      if (nyzp(1,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 70 k = 1, nyzp1
c        do 60 j = 1, nxv
c        do 50 n = 1, 3
c        scs(n,j,k,2,m) = scs(n,j,k,1,kl)
c  50    continue
c  60    continue
c  70    continue
c     else
c        do 100 k = 1, nyzp1
c        do 90 j = 1, nxv
c        do 80 n = 1, 3
c        scs(n,j,k,2,m) = f(n,j,3,k+1,m)
c  80    continue
c  90    continue
c 100    continue
c     endif
c     if (jr.lt.nvpy) then
c        if (nyzp(1,kr).eq.1) then
c           do 130 k = 1, nyzp1
c           do 120 j = 1, nxv
c           do 110 n = 1, 3
c           scs(n,j,k+nyzp1,2,m) = scs(n,j,k+nyzp1,1,kr)
c           scs(n,j,k+2*nyzp1,2,m) = scs(n,j,k+nyzp1,1,krr)
c 110       continue
c 120       continue
c 130       continue
c        else
c           do 160 k = 1, nyzp1
c           do 150 j = 1, nxv
c           do 140 n = 1, 3
c           scs(n,j,k+nyzp1,2,m) = scs(n,j,k+nyzp1,1,kr)
c           scs(n,j,k+2*nyzp1,2,m) = scs(n,j,k+2*nyzp1,1,kr)
c 140       continue
c 150       continue
c 160       continue
c        endif
c     else
c        do 190 k = 1, nyzp1
c        do 180 j = 2, nx2
c        do 170 n = 1, 3
c        scs(n,j,k+2*nyzp1,2,m) = 2.*scs(n,j,k+nyzp1,2,m) - scs(n,j,k,1,
c    1m)
c 170    continue
c 180    continue
c 190    continue
c     endif
c     if (nyzp(1,m).eq.1) then
c        if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
c           do 220 k = 1, nyzp1
c           do 210 j = 1, nxv
c           do 200 n = 1, 3
c           scs(n,j,k,2,m) = scs(n,j,k+nyzp1,1,kr)
c 200       continue
c 210       continue
c 220       continue
c        endif
c        if (jr.eq.(nvpy-1)) then
c           do 230 k = 1, nyzp1
c           do 240 j = 1, nxv
c           do 250 n = 1, 3
c           scs(n,j,k+2*nyzp1,2,m) = scs(n,j,k+nyzp1,2,kr)
c 230       continue
c 240       continue
c 250       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,1,2,m),3*nxvz,mreal,kl-1,noff+3,lgrp,msi
     1d,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,1,m),3*nxvzs,mreal,kr-1,noff+3,lgrp,ier
     1r)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 70 k = 1, nyzp1
         do 60 j = 1, nxv
         do 50 n = 1, 3
         scs(n,j,k,2,m) = f(n,j,3,k+1,m)
   50    continue
   60    continue
   70    continue
      endif
      if (jr.lt.nvpy) then
         call MPI_IRECV(scs(1,1,nyzp1+1,2,m),6*nxvz,mreal,kr-1,noff+4,lg
     1rp,msid,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scs(1,1,nyzp1+1,1,m),3*(2-ngc)*nxvzs,mreal,kl-1,n
     1off+4,lgrp,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 100 k = 1, nyzp1
         do 90 j = 2, nx2
         do 80 n = 1, 3
         scs(n,j,k+2*nyzp1,2,m) = 2.*scs(n,j,k+nyzp1,2,m) - scs(n,j,k,1,
     1m)
   80    continue
   90    continue
  100    continue
      endif
c special case of only one grid per processor in y
      if (mter.ge.1) go to 260
      if (jr.lt.nvpy) call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      if (jrr.lt.nvpy) then
         if (nps.eq.(3*nxvzs)) then
            call MPI_IRECV(scs(1,1,2*nyzp1+1,2,m),3*nxvz,mreal,krr-1,nof
     1f+5,lgrp,msid,ierr)
         else
            call MPI_IRECV(scs(1,1,1,1,m),3*nxvz,mreal,krr-1,noff+5,lgrp
     1,msid,ierr)
         endif
      else if (jr.lt.nvpy) then
         if (nps.eq.(3*nxvzs)) then
            call MPI_IRECV(scs(1,1,2*nyzp1+1,2,m),3*nxvz,mreal,kr-1,noff
     1+5,lgrp,msid,ierr)
         else
            call MPI_IRECV(scs(1,1,1,1,m),3*nxvz,mreal,kr-1,noff+5,lgrp,
     1msid,ierr)
         endif
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
         if (ngc.eq.1) then
            call MPI_IRECV(scs(1,1,1,2,m),3*nxvz,mreal,kr-1,noff+5,lgrp,
     1nsid,ierr)
         else
            call MPI_IRECV(scs(1,1,1,1,m),3*nxvz,mreal,kr-1,noff+5,lgrp,
     1nsid,ierr)
         endif
      endif
      if (jll.ge.0) then
         call MPI_SEND(scs(1,1,2*nyzp1+1,1,m),3*nxvzs,mreal,kll-1,noff+5
     1,lgrp,ierr)
      else if (jl.eq.0) then
         call MPI_SEND(scs(1,1,2*nyzp1+1,1,m),3*nxvzs,mreal,kl-1,noff+5,
     1lgrp,ierr)
      endif
      if ((jl.eq.(nvpy-2)).and.(jl.ge.0)) then
         call MPI_SEND(scs(1,1,nyzp1+1,2,m),3*nxvzs,mreal,kl-1,noff+5,lg
     1rp,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
         call MPI_WAIT(nsid,istatus,ierr)
      endif
c copy guard cells
  260 do 290 k = 1, nyzp1
      do 280 j = 1, nxv
      do 270 n = 1, 3
      f(n,j,1,k+1,m) = scs(n,j,k,2,m)
      f(n,j,nyzp(1,m)+2,k+1,m) = scs(n,j,k+nyzp1,2,m)
      f(n,j,nyzp(1,m)+3,k+1,m) = scs(n,j,k+2*nyzp1,2,m)
  270 continue
  280 continue
  290 continue
  300 continue
  310 continue
c fix left edge
      do 360 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 350 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(2,m) + 1
      kl = my + js
      if (kl.eq.0) then
         do 340 k = 1, nyzp1
         do 330 j = 2, nx2
         do 320 n = 1, 3
         f(n,j,1,k+1,m) = 2.*f(n,j,2,k+1,m) - f(n,j,1,k+1,m)
  320    continue
  330    continue
  340    continue
      endif
  350 continue
  360 continue
c buffer data in z
      do 400 m = 1, mnblok
      nyzp3 = nyzp(1,m) + 3
      do 390 k = 1, nyzp3
      do 380 j = 1, nxv
      do 370 n = 1, 3
      scr(n,j,k,1,m) = f(n,j,k,nyzp(2,m)+1,m)
      scr(n,j,k+nyzp3,1,m) = f(n,j,k,2,m)
      scr(n,j,k+2*nyzp3,1,m) = f(n,j,k,3,m)
  370 continue
  380 continue
  390 continue
  400 continue
c copy to guard cells in z
      do 520 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 510 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(1,m) + 3
      nxvys = nxv*nyzp3
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      krr = kr + 1
      if (krr.ge.nvpz) krr = krr - nvpz
      krr = ky + nvpy*krr
      kr = ky + nvpy*kr
      kl = kz - 1
      if (kl.lt.0) kl = kl + nvpz
      kll = kl - 1
      if (kll.lt.0) kll = kll + nvpz
      kll = ky + nvpy*kll
      kl = ky + nvpy*kl
      ngc = 0
c special case of only one grid per processor
      if (nyzp(2,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (nyzp(2,kr).eq.1) then
c        do 430 k = 1, nyzp3
c        do 420 j = 1, nxv
c        do 410 n = 1, 3
c        scr(n,j,k,2,m) = scr(n,j,k,1,kl)
c        scr(n,j,k+nyzp3,2,m) = scr(n,j,k+nyzp3,1,kr)
c        scr(n,j,k+2*nyzp3,2,m) = scr(n,j,k+nyzp3,1,krr)
c 410    continue
c 420    continue
c 430    continue
c     else
c        do 460 k = 1, nyzp3
c        do 450 j = 1, nxv
c        do 440 n = 1, 3
c        scr(n,j,k,2,m) = scr(n,j,k,1,kl)
c        scr(n,j,k+nyzp3,2,m) = scr(n,j,k+nyzp3,1,kr)
c        scr(n,j,k+2*nyzp3,2,m) = scr(n,j,k+2*nyzp3,1,kr)
c 440    continue
c 450    continue
c 460    continue
c     endif
c this segment is used for mpi computers
      call MPI_IRECV(scr(1,1,1,2,m),3*nxvy,mreal,kl-1,noff+7,lgrp,msid,i
     1err)
      call MPI_SEND(scr(1,1,1,1,m),3*nxvys,mreal,kr-1,noff+7,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scr(1,1,nyzp3+1,2,m),6*nxvy,mreal,kr-1,noff+8,lgrp,
     1msid,ierr)
      call MPI_SEND(scr(1,1,nyzp3+1,1,m),3*(2-ngc)*nxvys,mreal,kl-1,noff
     1+8,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c special case of only one grid per processor in z
      if (nter.ge.1) go to 470
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      if (nps.eq.(3*nxvys)) then
         call MPI_IRECV(scr(1,1,2*nyzp3+1,2,m),3*nxvy,mreal,krr-1,noff+1
     10,lgrp,msid,ierr)
      else
         call MPI_IRECV(scr(1,1,1,1,m),3*nxvy,mreal,krr-1,noff+10,lgrp,m
     1sid,ierr)
      endif
      call MPI_SEND(scr(1,1,nyzp3+1,1,m),3*nxvys,mreal,kll-1,noff+10,lgr
     1p,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c copy guard cells
  470 do 500 k = 1, nyzp3
      do 490 j = 1, nxv
      do 480 n = 1, 3
      f(n,j,k,1,m) = scr(n,j,k,2,m)
      f(n,j,k,nyzp(2,m)+2,m) = scr(n,j,k+nyzp3,2,m)
      f(n,j,k,nyzp(2,m)+3,m) = scr(n,j,k+2*nyzp3,2,m)
  480 continue
  490 continue
  500 continue
  510 continue
  520 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNMDGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypmx
     1,nzpmx,mblok,nblok,ngds,idds,mter,nter)
c this subroutine copies data to guard cells in non-uniform partitions
c for scalar data the field is replicated so as to disable quadratic
c interpolation within half a cell of the edges, and reduce it to linear
c interpolation in the y direction, and periodic in the z direction.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes three extra guard cells.
c scs/scr = scratch arrays for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c it is assumed that nyzp(n,m) > 0.
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c mter/nter = (0,1) = (no,yes) pass data to next processor only in y/z
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds, mter, nter
      integer nyzp
      real f, scs, scr
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx*ngds,2,mblok*nblok)
      dimension scr(nxv,nypmx*ngds,2,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, nps
      integer mnblok, nx2, nyzp1, nxvz, nxvzs, nyzp3, nxvy, nxvys
      integer m, my, mz, j, k, jr, jrr, jl, jll
      dimension istatus(lstat)
      nx2 = nx + 2
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c buffer data in y
      do 30 m = 1, mnblok
      ngc = 0
      if (nyzp(1,m).eq.1) ngc = 1
      nyzp1 = nyzp(2,m) + 1
      do 20 k = 1, nyzp1
      do 10 j = 1, nxv
      scs(j,k,1,m) = f(j,nyzp(1,m)+1,k+1,m)
      scs(j,k+nyzp1,1,m) = f(j,2,k+1,m)
      scs(j,k+2*nyzp1,1,m) = f(j,3-ngc,k+1,m)
      scs(j,k+nyzp1,2,m) = f(j,nyzp(1,m)+2,k+1,m)
   10 continue
   20 continue
   30 continue
c copy to guard cells in y
      do 220 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 210 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(2,m) + 1
      nxvzs = nxv*nyzp1
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr + 1
      krr = jrr + kz
      jl = ky - 1
      jll = jl - 1
      kll = jll + kz
      kr = jr + kz
      kl = jl + kz
      ngc = 0
c special case of only one grid per processor
      if (nyzp(1,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 50 k = 1, nyzp1
c        do 40 j = 1, nxv
c        scs(j,k,2,m) = scs(j,k,1,kl)
c  40    continue
c  50    continue
c     else
c        do 70 k = 1, nyzp1
c        do 60 j = 1, nxv
c        scs(j,k,2,m) = f(j,3,k+1,m)
c  60    continue
c  70    continue
c     endif
c     if (jr.lt.nvpy) then
c        if (nyzp(1,kr).eq.1) then
c           do 90 k = 1, nyzp1
c           do 80 j = 1, nxv
c           scs(j,k+nyzp1,2,m) = scs(j,k+nyzp1,1,kr)
c           scs(j,k+2*nyzp1,2,m) = scs(j,k+nyzp1,1,krr)
c  80       continue
c  90       continue
c        else
c           do 110 k = 1, nyzp1
c           do 100 j = 1, nxv
c           scs(j,k+nyzp1,2,m) = scs(j,k+nyzp1,1,kr)
c           scs(j,k+2*nyzp1,2,m) = scs(j,k+2*nyzp1,1,kr)
c 100       continue
c 110       continue
c        endif
c     else
c        do 130 k = 1, nyzp1
c        do 120 j = 2, nx2
c        scs(j,k+2*nyzp1,2,m) = 2.*scs(j,k+nyzp1,2,m) - scs(j,k,1,m)
c 120    continue
c 130    continue
c     endif
c     if (nyzp(1,m).eq.1) then
c        if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
c           do 150 k = 1, nyzp1
c           do 140 j = 1, nxv
c           scs(j,k,2,m) = scs(j,k+nyzp1,1,kr)
c 140       continue
c 150       continue
c        endif
c        if (jr.eq.(nvpy-1)) then
c           do 170 k = 1, nyzp1
c           do 160 j = 1, nxv
c           scs(j,k+2*nyzp1,2,m) = scs(j,k+nyzp1,2,kr)
c 160       continue
c 170       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,2,m),nxvz,mreal,kl-1,noff+3,lgrp,msid,ie
     1rr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,m),nxvzs,mreal,kr-1,noff+3,lgrp,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 50 k = 1, nyzp1
         do 40 j = 1, nxv
         scs(j,k,2,m) = f(j,3,k+1,m)
   40    continue
   50    continue
      endif
      if (jr.lt.nvpy) then
         call MPI_IRECV(scs(1,nyzp1+1,2,m),2*nxvz,mreal,kr-1,noff+4,lgrp
     1,msid,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scs(1,nyzp1+1,1,m),(2-ngc)*nxvzs,mreal,kl-1,noff+
     14,lgrp,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 70 k = 1, nyzp1
         do 60 j = 2, nx2
         scs(j,k+2*nyzp1,2,m) = 2.*scs(j,k+nyzp1,2,m) - scs(j,k,1,m)
   60    continue
   70    continue
      endif
c special case of only one grid per processor in y
      if (mter.ge.1) go to 180
      if (jr.lt.nvpy) call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      if (jrr.lt.nvpy) then
         if (nps.eq.nxvzs) then
            call MPI_IRECV(scs(1,2*nyzp1+1,2,m),nxvz,mreal,krr-1,noff+5,
     1lgrp,msid,ierr)
         else
            call MPI_IRECV(scs(1,1,1,m),nxvz,mreal,krr-1,noff+5,lgrp,msi
     1d,ierr)
         endif
      else if (jr.lt.nvpy) then
         if (nps.eq.nxvzs) then
            call MPI_IRECV(scs(1,2*nyzp1+1,2,m),nxvz,mreal,kr-1,noff+5,l
     1grp,msid,ierr)
         else
            call MPI_IRECV(scs(1,1,1,m),nxvz,mreal,kr-1,noff+5,lgrp,msid
     1,ierr)
         endif
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
         if (ngc.eq.1) then
            call MPI_IRECV(scs(1,1,2,m),nxvz,mreal,kr-1,noff+5,lgrp,nsid
     1,ierr)
         else
            call MPI_IRECV(scs(1,1,1,m),nxvz,mreal,kr-1,noff+5,lgrp,nsid
     1,ierr)
         endif
      endif
      if (jll.ge.0) then
         call MPI_SEND(scs(1,2*nyzp1+1,1,m),nxvzs,mreal,kll-1,noff+5,lgr
     1p,ierr)
      else if (jl.eq.0) then
         call MPI_SEND(scs(1,2*nyzp1+1,1,m),nxvzs,mreal,kl-1,noff+5,lgrp
     1,ierr)
      endif
      if ((jl.eq.(nvpy-2)).and.(jl.ge.0)) then
         call MPI_SEND(scs(1,nyzp1+1,2,m),nxvzs,mreal,kl-1,noff+5,lgrp,i
     1err)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
         call MPI_WAIT(nsid,istatus,ierr)
      endif
c copy guard cells
  180 do 200 k = 1, nyzp1
      do 190 j = 1, nxv
      f(j,1,k+1,m) = scs(j,k,2,m)
      f(j,nyzp(1,m)+2,k+1,m) = scs(j,k+nyzp1,2,m)
      f(j,nyzp(1,m)+3,k+1,m) = scs(j,k+2*nyzp1,2,m)
  190 continue
  200 continue
  210 continue
  220 continue
c fix left edge
      do 260 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 250 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(2,m) + 1
      kl = my + js
      if (kl.eq.0) then
         do 240 k = 1, nyzp1
         do 230 j = 2, nx2
         f(j,1,k+1,m) = 2.*f(j,2,k+1,m) - f(j,1,k+1,m)
  230    continue
  240    continue
      endif
  250 continue
  260 continue
c buffer data in z
      do 290 m = 1, mnblok
      nyzp3 = nyzp(1,m) + 3
      do 280 k = 1, nyzp3
      do 270 j = 1, nxv
      scr(j,k,1,m) = f(j,k,nyzp(2,m)+1,m)
      scr(j,k+nyzp3,1,m) = f(j,k,2,m)
      scr(j,k+2*nyzp3,1,m) = f(j,k,3,m)
  270 continue
  280 continue
  290 continue
c copy to guard cells in z
      do 380 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 370 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(1,m) + 3
      nxvys = nxv*nyzp3
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      krr = kr + 1
      if (krr.ge.nvpz) krr = krr - nvpz
      krr = ky + nvpy*krr
      kr = ky + nvpy*kr
      kl = kz - 1
      if (kl.lt.0) kl = kl + nvpz
      kll = kl - 1
      if (kll.lt.0) kll = kll + nvpz
      kll = ky + nvpy*kll
      kl = ky + nvpy*kl
      ngc = 0
c special case of only one grid per processor
      if (nyzp(2,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (nyzp(2,kr).eq.1) then
c        do 310 k = 1, nyzp3
c        do 300 j = 1, nxv
c        scr(j,k,2,m) = scr(j,k,1,kl)
c        scr(j,k+nyzp3,2,m) = scr(j,k+nyzp3,1,kr)
c        scr(j,k+2*nyzp3,2,m) = scr(j,k+nyzp3,1,krr)
c 300    continue
c 310    continue
c     else
c        do 330 k = 1, nyzp3
c        do 320 j = 1, nxv
c        scr(j,k,2,m) = scr(j,k,1,kl)
c        scr(j,k+nyzp3,2,m) = scr(j,k+nyzp3,1,kr)
c        scr(j,k+2*nyzp3,2,m) = scr(j,k+2*nyzp3,1,kr)
c 320    continue
c 330    continue
c     endif
c this segment is used for mpi computers
      call MPI_IRECV(scr(1,1,2,m),nxvy,mreal,kl-1,noff+7,lgrp,msid,ierr)
      call MPI_SEND(scr(1,1,1,m),nxvys,mreal,kr-1,noff+7,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scr(1,nyzp3+1,2,m),2*nxvy,mreal,kr-1,noff+8,lgrp,ms
     1id,ierr)
      call MPI_SEND(scr(1,nyzp3+1,1,m),(2-ngc)*nxvys,mreal,kl-1,noff+8,l
     1grp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c special case of only one grid per processor in z
      if (nter.ge.1) go to 340
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      if (nps.eq.nxvys) then
         call MPI_IRECV(scr(1,2*nyzp3+1,2,m),nxvy,mreal,krr-1,noff+10,lg
     1rp,msid,ierr)
      else
         call MPI_IRECV(scr(1,1,1,m),nxvy,mreal,krr-1,noff+10,lgrp,msid,
     1ierr)
      endif
      call MPI_SEND(scr(1,nyzp3+1,1,m),nxvys,mreal,kll-1,noff+10,lgrp,ie
     1rr)
      call MPI_WAIT(msid,istatus,ierr)
c copy guard cells
  340 do 360 k = 1, nyzp3
      do 350 j = 1, nxv
      f(j,k,1,m) = scr(j,k,2,m)
      f(j,k,nyzp(2,m)+2,m) = scr(j,k+nyzp3,2,m)
      f(j,k,nyzp(2,m)+3,m) = scr(j,k+2*nyzp3,2,m)
  350 continue
  360 continue
  370 continue
  380 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMCGUARD32L(f,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mblok,nbl
     1ok,kzp)
c this subroutine copies data from field to particle partitions, copying
c data to guard cells, where the field and particle partitions are 
c assumed to be the same.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes one extra guard
c cell.
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = first dimension of f, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kzp = number of complex grids in z for each field partition.
c linear interpolation, for distributed data,
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, mblok, nblok, kzp
      real f
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer ky, kz, js, ks, moff, noff, kr, kl, nxvy, m, my, mz
      integer j, k
      dimension istatus(lstat)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      nxvy = nxv*nypmx
c copy to guard cells in z
      do 40 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 30 my = 1, mblok
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
c     do 20 k = 1, nypmx
c     do 10 j = 1, nxv
c     f(j,k,kzp+1,m) = f(j,k,1,kr)
c  10 continue
c  20 continue
c this segment is used for mpi computers
      call MPI_IRECV(f(1,1,kzp+1,m),nxvy,mreal,kr-1,noff+4,lgrp,msid,ier
     1r)
      call MPI_SEND(f(1,1,1,m),nxvy,mreal,kl-1,noff+4,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNMCGUARD32L(f,scs,nyzp,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx
     1,mblok,nblok,ngds,idds)
c this subroutine copies data to guard cells in non-uniform partitions
c guard cell on last processor is presumed already set in y.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes one extra guard cell.
c scs(j,k,m) = scratch array for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c it is assumed that nyzp(n,m) > 0.
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c linear interpolation, for distributed data,
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nxv, nypmx, nzpmx, mblok, nblok, ngds
      integer idds
      integer nyzp
      real f, scs
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer ky, kz, js, ks, moff, noff, jr, jl, kr, kl, mnblok
      integer nxvz, nxvzs, nyzp1, nxvy, nxvys, m, my, mz, j, k
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
      do 20 k = 1, nyzp(2,m)
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
      nxvzs = nxv*nyzp(2,m)
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jl = ky - 1
      kr = jr + kz
      kl = jl + kz
c this segment is used for shared memory computers
c     if (jr.lt.nvpy) then
c        do 50 k = 1, nyzp(2,m)
c        do 40 j = 1, nxv
c        scs(j,k,2,m) = scs(j,k,1,kr)
c  40    continue
c  50    continue
c     endif
c this segment is used for mpi computers
      if (jr.lt.nvpy) then
         call MPI_IRECV(scs(1,1,2,m),nxvz,mreal,kr-1,noff+3,lgrp,msid,ie
     1rr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scs(1,1,1,m),nxvzs,mreal,kl-1,noff+3,lgrp,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
c copy guard cells
         do 70 k = 1, nyzp(2,m)
         do 60 j = 1, nxv
         f(j,nyzp(1,m)+1,k,m) = scs(j,k,2,m)
   60    continue
   70    continue
      endif
   80 continue
   90 continue
c copy to guard cells in z
      do 130 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 120 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(1,m) + 1
      nxvys = nxv*nyzp1
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      kl = kz - 1
      if (kl.lt.0) kl = kl + nvpz
      kr = ky + nvpy*kr
      kl = ky + nvpy*kl
c this segment is used for shared memory computers
c     do 110 k = 1, nyzp1
c     do 100 j = 1, nxv
c     f(j,k,nyzp(2,m)+1,m) = f(j,k,1,kr)
c 100 continue
c 110 continue
c this segment is used for mpi computers
      call MPI_IRECV(f(1,1,nyzp(2,m)+1,m),nxvy,mreal,kr-1,noff+4,lgrp,ms
     1id,ierr)
      call MPI_SEND(f(1,1,1,m),nxvys,mreal,kl-1,noff+4,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMACGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpm
     1x,mblok,nblok,kyp,kzp,ngds)
c thus subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y direction, and periodic in the z direction.
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes three extra
c guard cells.
c scs/scr = scratch arrays for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer kyp, kzp, ngds
      real f, scs, scr
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(3,nxv,nypmx,ngds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, mnblok
      integer nx3, nyp3, nzp3, nxvz, nxvy, m, my, mz, j, k, n
      integer jr, jrr, jl, jll
      dimension istatus(lstat)
      nx3 = nx + 3
      nyp3 = kyp + 3
      nzp3 = kzp + 3
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c buffer data in y
      do 40 m = 1, mnblok
      do 30 k = 1, nzpmx
      do 20 j = 1, nxv
      do 10 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,kyp+2,k,m)
      scs(n,j,k,2,m) = f(n,j,kyp+3,k,m)
      scs(n,j,k,3,m) = f(n,j,1,k,m)
   10 continue
   20 continue
   30 continue
   40 continue
c add guard cells in y
      do 300 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 290 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr
      krr = jr + kz
      jl = ky - 1
      jll = jl
      kll = jl + kz
      ngc = 2
c special case of only one grid per processor
      if (kyp.eq.1) then
         jrr = jr + 1
         krr = jrr + kz
         jll = jl - 1
         kll = jll + kz
         ngc = 1
      endif
      kr = jr + kz
      kl = jl + kz
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 70 k = 1, nzpmx
c        do 60 j = 1, nxv
c        do 50 n = 1, 3
c        scs(n,j,k,4,m) = scs(n,j,k,1,kl)
c        scs(n,j,k,5,m) = scs(n,j,k,2,kll)
c  50    continue
c  60    continue
c  70    continue
c     else
c        do 100 k = 1, nzp3
c        do 90 j = 1, nx3
c        do 80 n = 1, 3
c        scs(n,j,k,4,m) = 2.*scs(n,j,k,3,m)
c        scs(n,j,k,5,m) = -scs(n,j,k,3,m)
c  80    continue
c  90    continue
c 100    continue
c     endif
c     if (jr.lt.nvpy) then
c        do 130 k = 1, nzpmx
c        do 120 j = 1, nxv
c        do 110 n = 1, 3
c        scs(n,j,k,6,m) = scs(n,j,k,3,kr)
c 110    continue
c 120    continue
c 130    continue
c     else
c        do 160 k = 1, nzp3
c        do 150 j = 1, nx3
c        do 140 n = 1, 3
c        scs(n,j,k,6,m) = -scs(n,j,k,2,m)
c        f(n,j,kyp+2,k,m) = f(n,j,kyp+2,k,m) + 2.*scs(n,j,k,2,m)
c        f(n,j,kyp+3,k,m) = 0.
c 140    continue
c 150    continue
c 160    continue
c     endif
c     if (kyp.eq.1) then
c        if (jl.eq.0) then
c           do 190 k = 1, nzp3
c           do 180 j = 1, nx3
c           do 170 n = 1, 3
c           scs(n,j,k,4,m) = scs(n,j,k,1,kl)
c           scs(n,j,k,5,m) = -scs(n,j,k,3,kl)
c 170       continue
c 180       continue
c 190       continue
c        else if (jl.lt.0) then
c           do 220 k = 1, nzp3
c           do 210 j = 1, nx3
c           do 200 n = 1, 3
c           scs(n,j,k,5,m) = 0.
c 200       continue
c 210       continue
c 220       continue
c        endif
c last point is special with only one grid
c        if ((jl.eq.(nvpy-2)).and.(jl.ge.0)) then
c           do 250 k = 1, nzp3
c           do 240 j = 1, nx3
c           do 230 n = 1, 3
c           f(n,j,kyp+2,k,m) = f(n,j,kyp+2,k,m) + scs(n,j,k,2,kl)
c 230       continue
c 240       continue
c 250       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,1,4,m),3*ngc*nxvz,mreal,kl-1,noff+1,lgrp
     1,msid,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,1,m),3*ngc*nxvz,mreal,kr-1,noff+1,lgrp,
     1ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 70 k = 1, nzp3
         do 60 j = 1, nx3
         do 50 n = 1, 3
         scs(n,j,k,4,m) = 2.*scs(n,j,k,3,m)
         scs(n,j,k,5,m) = -scs(n,j,k,3,m)
   50    continue
   60    continue
   70    continue
      endif
      if (jr.lt.nvpy) then
         call MPI_IRECV(scs(1,1,1,6,m),3*nxvz,mreal,kr-1,noff+2,lgrp,msi
     1d,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scs(1,1,1,3,m),3*nxvz,mreal,kl-1,noff+2,lgrp,ierr
     1)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 100 k = 1, nzp3
         do 90 j = 1, nx3
         do 80 n = 1, 3
         scs(n,j,k,6,m) = -scs(n,j,k,2,m)
         f(n,j,kyp+2,k,m) = f(n,j,kyp+2,k,m) + 2.*scs(n,j,k,2,m)
         f(n,j,kyp+3,k,m) = 0.
   80    continue
   90    continue
  100    continue
      endif
c special case of only one grid per processor
      if (kyp.eq.1) then
         if (jll.ge.0) then
            call MPI_IRECV(scs(1,1,1,5,m),3*nxvz,mreal,kll-1,noff+5,lgrp
     1,msid,ierr)
         else if (jl.eq.0) then
            call MPI_IRECV(scs(1,1,1,5,m),3*nxvz,mreal,kl-1,noff+5,lgrp,
     1msid,ierr)
         endif
         if (jrr.lt.nvpy) then
            call MPI_SEND(scs(1,1,1,2,m),3*nxvz,mreal,krr-1,noff+5,lgrp,
     1ierr)
         endif
         if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
            call MPI_SEND(scs(1,1,1,3,m),3*nxvz,mreal,kr-1,noff+5,lgrp,i
     1err)
         endif
         if (jl.ge.0) then
            call MPI_WAIT(msid,istatus,ierr)
            if (jl.eq.0) then
               do 130 k = 1, nzp3
               do 120 j = 1, nx3
               do 110 n = 1, 3
               scs(n,j,k,5,m) = -scs(n,j,k,5,m)
  110          continue
  120          continue
  130          continue
            endif
         else
            do 160 k = 1, nzp3
            do 150 j = 1, nx3
            do 140 n = 1, 3
            scs(n,j,k,5,m) = 0.
  140       continue
  150       continue
  160       continue
         endif
c last point is special with only one grid
         if ((jl.eq.(nvpy-2)).and.(jl.ge.0)) then
            call MPI_IRECV(scs(1,1,1,2,m),3*nxvz,mreal,kl-1,noff+6,lgrp,
     1msid,ierr)
         endif
         if (jr.eq.(nvpy-1)) then
            call MPI_SEND(scs(1,1,1,2,m),3*nxvz,mreal,kr-1,noff+6,lgrp,i
     1err)
         endif
         if ((jl.eq.(nvpy-2)).and.(jl.ge.0)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 190 k = 1, nzp3
            do 180 j = 1, nx3
            do 170 n = 1, 3
            f(n,j,kyp+2,k,m) = f(n,j,kyp+2,k,m) + scs(n,j,k,2,m)
            f(n,j,kzp+3,k,m) = 0.
  170       continue
  180       continue
  190       continue
         endif
      endif
c add up the guard cells
      do 280 k = 1, nzp3
      do 270 j = 1, nx3
      do 260 n = 1, 3
      f(n,j,2,k,m) = f(n,j,2,k,m) + scs(n,j,k,4,m)
      f(n,j,ngc+1,k,m) = f(n,j,ngc+1,k,m) + scs(n,j,k,5,m)
      f(n,j,kyp+1,k,m) = f(n,j,kyp+1,k,m) + scs(n,j,k,6,m)
  260 continue
  270 continue
  280 continue
  290 continue
  300 continue
c zero out the left edge
      do 350 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 340 my = 1, mblok
      m = my + moff
      kl = my + js
      if (kl.eq.0) then
         do 330 k = 1, nzp3
         do 320 j = 1, nx3
         do 310 n = 1, 3
         f(n,j,1,k,m) = 0.
  310    continue
  320    continue
  330    continue
      endif
  340 continue
  350 continue
c add guard cells in z
      do 400 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 390 my = 1, mblok
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
c     do 280 k = 1, nypmx
c     do 270 j = 1, nxv
c     do 360 n = 1, 3
c     scr(n,j,k,1,m) = f(n,j,k,kzp+2,kl)
c     scr(n,j,k,2,m) = f(n,j,k,kzp+3,kll)
c     scr(n,j,k,3,m) = f(n,j,k,1,kr)
c 360 continue
c 370 continue
c 380 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,3*ngc*nxvy,mreal,kl-1,noff+7,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,1,kzp+2,m),3*ngc*nxvy,mreal,kr-1,noff+7,lgrp,i
     1err)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scr(1,1,1,3,m),3*nxvy,mreal,kr-1,noff+8,lgrp,msid,i
     1err)
      call MPI_SEND(f(1,1,1,1,m),3*nxvy,mreal,kl-1,noff+8,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      if (kzp.eq.1) then
         call MPI_IRECV(scr(1,1,1,2,m),3*ngc*nxvy,mreal,kll-1,noff+9,lgr
     1p,msid,ierr)
         call MPI_SEND(f(1,1,1,kzp+3,m),3*ngc*nxvy,mreal,krr-1,noff+9,lg
     1rp,ierr)
         call MPI_WAIT(msid,istatus,ierr)
      endif
c add up the guard cells
      do 380 k = 1, nypmx
      do 370 j = 1, nxv
      do 360 n = 1, 3
      f(n,j,k,2,m) = f(n,j,k,2,m) + scr(n,j,k,1,m)
      f(n,j,k,ngc+1,m) = f(n,j,k,ngc+1,m) + scr(n,j,k,2,m)
      f(n,j,k,kzp+1,m) = f(n,j,k,kzp+1,m) + scr(n,j,k,3,m)
  360 continue
  370 continue
  380 continue
  390 continue
  400 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMAGUARD32(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpmx
     1,mblok,nblok,kyp,kzp,ngds)
c thus subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y direction, and periodic in the z direction.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes three extra
c guard cells.
c scs/scr = scratch arrays for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer kyp, kzp, ngds
      real f, scs, scr
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(nxv,nypmx,ngds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, mnblok
      integer nx3, nyp3, nzp3, nxvz, nxvy, m, my, mz, j, k
      integer jr, jrr, jl, jll
      dimension istatus(lstat)
      nx3 = nx + 3
      nyp3 = kyp + 3
      nzp3 = kzp + 3
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
      scs(j,k,1,m) = f(j,kyp+2,k,m)
      scs(j,k,2,m) = f(j,kyp+3,k,m)
      scs(j,k,3,m) = f(j,1,k,m)
   10 continue
   20 continue
   30 continue
c add guard cells in y
      do 210 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 200 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr
      krr = jr + kz
      jl = ky - 1
      jll = jl
      kll = jl + kz
      ngc = 2
c special case of only one grid per processor
      if (kyp.eq.1) then
         jrr = jr + 1
         krr = jrr + kz
         jll = jl - 1
         kll = jll + kz
         ngc = 1
      endif
      kr = jr + kz
      kl = jl + kz
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 50 k = 1, nzpmx
c        do 40 j = 1, nxv
c        scs(j,k,4,m) = scs(j,k,1,kl)
c        scs(j,k,5,m) = scs(j,k,2,kll)
c  40    continue
c  50    continue
c     else
c        do 70 k = 1, nzp3
c        do 60 j = 1, nx3
c        scs(j,k,4,m) = 2.*scs(j,k,3,m)
c        scs(j,k,5,m) = -scs(j,k,3,m)
c  60    continue
c  70    continue
c     endif
c     if (jr.lt.nvpy) then
c        do 90 k = 1, nzpmx
c        do 80 j = 1, nxv
c        scs(j,k,6,m) = scs(j,k,3,kr)
c  80    continue
c  90    continue
c     else
c        do 110 k = 1, nzp3
c        do 100 j = 1, nx3
c        scs(j,k,6,m) = -scs(j,k,2,m)
c        f(j,kyp+2,k,m) = f(j,kyp+2,k,m) + 2.*scs(j,k,2,m)
c        f(j,kyp+3,k,m) = 0.
c 100    continue
c 110    continue
c     endif
c     if (kyp.eq.1) then
c        if (jl.eq.0) then
c           do 130 k = 1, nzp3
c           do 120 j = 1, nx3
c           scs(j,k,4,m) = scs(j,k,1,kl)
c           scs(j,k,5,m) = -scs(j,k,3,kl)
c 120       continue
c 130       continue
c        else if (jl.lt.0) then
c           do 150 k = 1, nzp3
c           do 140 j = 1, nx3
c           scs(j,k,5,m) = 0.
c 140       continue
c 150       continue
c        endif
c last point is special with only one grid
c        if ((jl.eq.(nvpy-2)).and.(jl.ge.0)) then
c           do 170 k = 1, nzp3
c           do 160 j = 1, nx3
c           f(j,kyp+2,k,m) = f(j,kyp+2,k,m) + scs(j,k,2,kl)
c 160       continue
c 170       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,4,m),ngc*nxvz,mreal,kl-1,noff+1,lgrp,msi
     1d,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,m),ngc*nxvz,mreal,kr-1,noff+1,lgrp,ierr
     1)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 50 k = 1, nzp3
         do 40 j = 1, nx3
         scs(j,k,4,m) = 2.*scs(j,k,3,m)
         scs(j,k,5,m) = -scs(j,k,3,m)
   40    continue
   50    continue
      endif
      if (jr.lt.nvpy) then
         call MPI_IRECV(scs(1,1,6,m),nxvz,mreal,kr-1,noff+2,lgrp,msid,ie
     1rr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scs(1,1,3,m),nxvz,mreal,kl-1,noff+2,lgrp,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 70 k = 1, nzp3
         do 60 j = 1, nx3
         scs(j,k,6,m) = -scs(j,k,2,m)
         f(j,kyp+2,k,m) = f(j,kyp+2,k,m) + 2.*scs(j,k,2,m)
         f(j,kyp+3,k,m) = 0.
   60    continue
   70    continue
      endif
c special case of only one grid per processor
      if (kyp.eq.1) then
         if (jll.ge.0) then
            call MPI_IRECV(scs(1,1,5,m),nxvz,mreal,kll-1,noff+5,lgrp,msi
     1d,ierr)
         else if (jl.eq.0) then
            call MPI_IRECV(scs(1,1,5,m),nxvz,mreal,kl-1,noff+5,lgrp,msid
     1,ierr)
         endif
         if (jrr.lt.nvpy) then
            call MPI_SEND(scs(1,1,2,m),nxvz,mreal,krr-1,noff+5,lgrp,ierr
     1)
         endif
         if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
            call MPI_SEND(scs(1,1,3,m),nxvz,mreal,kr-1,noff+5,lgrp,ierr)
         endif
         if (jl.ge.0) then
            call MPI_WAIT(msid,istatus,ierr)
            if (jl.eq.0) then
               do 90 k = 1, nzp3
               do 80 j = 1, nx3
               scs(j,k,5,m) = -scs(j,k,5,m)
   80          continue
   90          continue
            endif
         else
            do 110 k = 1, nzp3
            do 100 j = 1, nx3
            scs(j,k,5,m) = 0.
  100       continue
  110       continue
         endif
c last point is special with only one grid
         if ((jl.eq.(nvpy-2)).and.(jl.ge.0)) then
            call MPI_IRECV(scs(1,1,2,m),nxvz,mreal,kl-1,noff+6,lgrp,msid
     1,ierr)
         endif
         if (jr.eq.(nvpy-1)) then
            call MPI_SEND(scs(1,1,2,m),nxvz,mreal,kr-1,noff+6,lgrp,ierr)
         endif
         if ((jl.eq.(nvpy-2)).and.(jl.ge.0)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 130 k = 1, nzp3
            do 120 j = 1, nx3
            f(j,kyp+2,k,m) = f(j,kyp+2,k,m) + scs(j,k,2,m)
            f(j,kzp+3,k,m) = 0.
  120       continue
  130       continue
         endif
      endif
c add up the guard cells
      do 190 k = 1, nzp3
      do 180 j = 1, nx3
      f(j,2,k,m) = f(j,2,k,m) + scs(j,k,4,m)
      f(j,ngc+1,k,m) = f(j,ngc+1,k,m) + scs(j,k,5,m)
      f(j,kyp+1,k,m) = f(j,kyp+1,k,m) + scs(j,k,6,m)
  180 continue
  190 continue
  200 continue
  210 continue
c zero out the left edge
      do 250 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 240 my = 1, mblok
      m = my + moff
      kl = my + js
      if (kl.eq.0) then
         do 230 k = 1, nzp3
         do 220 j = 1, nx3
         f(j,1,k,m) = 0.
  220    continue
  230    continue
      endif
  240 continue
  250 continue
c add guard cells in z
      do 310 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 300 my = 1, mblok
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
c     do 270 k = 1, nypmx
c     do 260 j = 1, nxv
c     scr(j,k,1,m) = f(j,k,kzp+2,kl)
c     scr(j,k,2,m) = f(j,k,kzp+3,kll)
c     scr(j,k,3,m) = f(j,k,1,kr)
c 260 continue
c 270 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,ngc*nxvy,mreal,kl-1,noff+7,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,kzp+2,m),ngc*nxvy,mreal,kr-1,noff+7,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scr(1,1,3,m),nxvy,mreal,kr-1,noff+8,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,1,m),nxvy,mreal,kl-1,noff+8,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      if (kzp.eq.1) then
         call MPI_IRECV(scr(1,1,2,m),ngc*nxvy,mreal,kll-1,noff+9,lgrp,ms
     1id,ierr)
         call MPI_SEND(f(1,1,kzp+3,m),ngc*nxvy,mreal,krr-1,noff+9,lgrp,i
     1err)
         call MPI_WAIT(msid,istatus,ierr)
      endif
c add up the guard cells
      do 290 k = 1, nypmx
      do 280 j = 1, nxv
      f(j,k,2,m) = f(j,k,2,m) + scr(j,k,1,m)
      f(j,k,ngc+1,m) = f(j,k,ngc+1,m) + scr(j,k,2,m)
      f(j,k,kzp+1,m) = f(j,k,kzp+1,m) + scr(j,k,3,m)
  280 continue
  290 continue
  300 continue
  310 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNMACGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypm
     1x,nzpmx,mblok,nblok,ngds,idds,mter,nter)
c this subroutine adds data from guard cells in non-uniform partitions
c for vector data, the field is added up so as to disable quadratic
c interpolation within half a cell of the edges, and reduce it to linear
c interpolation in the y direction, and periodic in z direction.
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes three extra guard cells.
c scs/scr = scratch arrays for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c it is assumed the nyzp(n,m) > 0.
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c mter/nter = (0,1) = (no,yes) pass data to next processor only in y/z
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds, mter, nter
      integer nyzp
      real f, scs, scr
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx*ngds,2,mblok*nblok)
      dimension scr(3,nxv,nypmx*ngds,2,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, mnblok
      integer nx3, nxvz, nxvzs, nyzp3, nxvy, nxvys, m, my, mz, j, k, n
      integer jr, jrr, jl, jll
      dimension istatus(lstat)
      nx3 = nx + 3
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c buffer data in y
      do 40 m = 1, mnblok
      nyzp3 = nyzp(2,m) + 3
      do 30 k = 1, nyzp3
      do 20 j = 1, nxv
      do 10 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,nyzp(1,m)+2,k,m)
      scs(n,j,k+nyzp3,1,m) = f(n,j,nyzp(1,m)+3,k,m)
      scs(n,j,k+2*nyzp3,1,m) = f(n,j,1,k,m)
   10 continue
   20 continue
   30 continue
   40 continue
c add guard cells in y
      do 340 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 330 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(2,m) + 3
      nxvzs = nxv*nyzp3
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr + 1
      krr = jrr + kz
      jl = ky - 1
      jll = jl - 1
      kll = jll + kz
      kr = jr + kz
      kl = jl + kz
      ngc = 0
c special case of only one grid per processor
      if (nyzp(1,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        if (nyzp(1,kl).eq.1) then
c           do 70 k = 1, nyzp3
c           do 60 j = 1, nxv
c           do 50 n = 1, 3
c           scs(n,j,k,2,m) = scs(n,j,k,1,kl) + scs(n,j,k,1,kll)
c           scs(n,j,k+nyzp3,2,m) = scs(n,j,k+nyzp3,1,kl)
c  50       continue
c  60       continue
c  70       continue
c        else
c           do 100 k = 1, nyzp3
c           do 90 j = 1, nxv
c           do 80 n = 1, 3
c           scs(n,j,k,2,m) = scs(n,j,k,1,kl)
c           scs(n,j,k+nyzp3,2,m) = scs(n,j,k+nyzp3,1,kl)
c  80       continue
c  90       continue
c 100       continue
c        endif
c     else
c        do 130 k = 1, nyzp3
c        do 120 j = 1, nx3
c        do 110 n = 1, 3
c        scs(n,j,k,2,m) = 2.*scs(n,j,k+2*nyzp3,1,m)
c        scs(n,j,k+nyzp3,2,m) = -scs(n,j,k+2*nyzp3,1,m)
c 110    continue
c 120    continue
c 130    continue
c     endif
c     if (jr.lt.nvpy) then
c        do 160 k = 1, nyzp3
c        do 150 j = 1, nxv
c        do 140 n = 1, 3
c        scs(n,j,k+2*nyzp3,2,m) = scs(n,j,k+2*nyzp3,1,kr)
c 140    continue
c 150    continue
c 160    continue
c     else
c        do 190 k = 1, nyzp3
c        do 180 j = 1, nx3
c        do 170 n = 1, 3
c        scs(n,j,k+2*nyzp3,2,m) = -scs(n,j,k+nyzp3,1,m)
c        f(n,j,nyzp(1,m)+2,k,m) = f(n,j,nyzp(1,m)+2,k,m) + 2.*scs(n,j,k+
c    1nyzp3,1,m)
c        f(n,j,nyzp(1,m)+3,k,m) = 0.
c 170    continue
c 180    continue
c 190    continue
c     endif
c     if (jl.ge.0) then
c        if ((nyzp(1,kl).eq.1).and.(jl.eq.0)) then
c           do 220 k = 1, nyzp3
c           do 210 j = 1, nx3
c           do 200 n = 1, 3
c           scs(n,j,k,2,m) = scs(n,j,k,2,m) - scs(n,j,k+2*nyzp3,1,kl)
c 200       continue
c 210       continue
c 220       continue
c        endif
c     else if (nyzp(1,m).eq.1) then
c        do 250 k = 1, nyzp3
c        do 240 j = 1, nx3
c        do 230 n = 1, 3
c        scs(n,j,k+nyzp3,2,m) = 0.
c 230    continue
c 240    continue
c 250    continue
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,1,2,m),6*nxvz,mreal,kl-1,noff+1,lgrp,msi
     1d,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,1,m),6*nxvzs,mreal,kr-1,noff+1,lgrp,ier
     1r)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 70 k = 1, nyzp3
         do 60 j = 1, nx3
         do 50 n = 1, 3
         scs(n,j,k,2,m) = 2.*scs(n,j,k+2*nyzp3,1,m)
         scs(n,j,k+nyzp3,2,m) = -scs(n,j,k+2*nyzp3,1,m)
   50    continue
   60    continue
   70    continue
      endif
      if (jr.lt.nvpy) then
         call MPI_IRECV(scs(1,1,2*nyzp3+1,2,m),3*nxvz,mreal,kr-1,noff+2,
     1lgrp,msid,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scs(1,1,2*nyzp3+1,1,m),3*nxvzs,mreal,kl-1,noff+2,
     1lgrp,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 100 k = 1, nyzp3
         do 90 j = 1, nx3
         do 80 n = 1, 3
         scs(n,j,k+2*nyzp3,2,m) = -scs(n,j,k+nyzp3,1,m)
         f(n,j,nyzp(1,m)+2,k,m) = f(n,j,nyzp(1,m)+2,k,m) + 2.*scs(n,j,k+
     1nyzp3,1,m)
         f(n,j,nyzp(1,m)+3,k,m) = 0.
   80    continue
   90    continue
  100    continue
      endif
c special case of only one grid per processor in y
      if (mter.ge.1) go to 260
      if (jl.ge.0) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+3,lgrp,msid,ierr)
      endif
      if (jll.ge.0) then
         call MPI_IRECV(scs(1,1,1,1,m),3*nxvz,mreal,kll-1,noff+5,lgrp,ns
     1id,ierr)
      else if (jl.eq.0) then
         call MPI_IRECV(scs(1,1,1,1,m),3*nxvz,mreal,kl-1,noff+5,lgrp,nsi
     1d,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(nyzp(1,m),1,mint,kr-1,moff+3,lgrp,ierr)
      endif
      if (jrr.lt.nvpy) then
         call MPI_SEND(scs(1,1,nyzp3+1,1,m),3*nxvzs,mreal,krr-1,noff+5,l
     1grp,ierr)
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
         call MPI_SEND(scs(1,1,2*nyzp3+1,1,m),3*nxvzs,mreal,kr-1,noff+5,
     1lgrp,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         if (ngc.eq.1) then
            if (jl.eq.0) then
               do 130 k = 1, nyzp3
               do 120 j = 1, nx3
               do 110 n = 1, 3
               scs(n,j,k,2,m) = scs(n,j,k,2,m) - scs(n,j,k,1,m)
  110          continue
  120          continue
  130          continue
            else
               do 160 k = 1, nyzp3
               do 150 j = 1, nx3
               do 140 n = 1, 3
               scs(n,j,k,2,m) = scs(n,j,k,2,m) + scs(n,j,k,1,m)
  140          continue
  150          continue
  160          continue
            endif
         else if (nyzp(1,m).eq.1) then
            do 190 k = 1, nyzp3
            do 180 j = 1, nx3
            do 170 n = 1, 3
            scs(n,j,k+nyzp3,2,m) = 0.
  170       continue
  180       continue
  190       continue
         endif
      endif
c add up the guard cells
  260 do 290 k = 1, nyzp3
      do 280 j = 1, nx3
      do 270 n = 1, 3
      f(n,j,2,k,m) = f(n,j,2,k,m) + scs(n,j,k,2,m)
      f(n,j,3,k,m) = f(n,j,3,k,m) + scs(n,j,k+nyzp3,2,m)
      f(n,j,nyzp(1,m)+1,k,m) = f(n,j,nyzp(1,m)+1,k,m) + scs(n,j,k+2*nyzp
     13,2,m)
      f(n,j,1,k,m) = 0.
      f(n,j,nyzp(1,m)+3,k,m) = 0.
  270 continue
  280 continue
  290 continue
      if (jr.lt.nvpy) then
         do 320 k = 1, nyzp3
         do 310 j = 1, nx3
         do 300 n = 1, 3
         f(n,j,nyzp(1,m)+2,k,m) = 0.
  300    continue
  310    continue
  320    continue
      endif
  330 continue
  340 continue
c zero out the left edge
      do 390 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 380 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(2,m) + 3
      kl = my + js
      if (kl.eq.0) then
         do 370 k = 1, nyzp3
         do 360 j = 1, nx3
         do 350 n = 1, 3
         f(n,j,1,k,m) = 0.
  350    continue
  360    continue
  370    continue
      endif
  380 continue
  390 continue
c special case for one processor in z
      if (nvpz.eq.1) then
         do 430 m = 1, mnblok
         nyzp3 = nyzp(1,m) + 3
         do 420 k = 1, nyzp3
         do 410 j = 1, nx3
         do 400 n = 1, 3
         f(n,j,k,2,m) = f(n,j,k,2,m) + f(n,j,k,nyzp(2,m)+2,m)
         f(n,j,k,3,m) = f(n,j,k,3,m) + f(n,j,k,nyzp(2,m)+3,m)
         f(n,j,k,nyzp(2,m)+1,m) = f(n,j,k,nyzp(2,m)+1,m) + f(n,j,k,1,m)
         f(n,j,k,1,m) = 0.
         f(n,j,k,nyzp(2,m)+2,m) = 0.
         f(n,j,k,nyzp(2,m)+3,m) = 0.
  400    continue
  410    continue
  420    continue
  430    continue
         return
      endif
c buffer data in z
      do 470 m = 1, mnblok
      nyzp3 = nyzp(1,m) + 3
      do 460 k = 1, nyzp3
      do 450 j = 1, nxv
      do 440 n = 1, 3
      scr(n,j,k,1,m) = f(n,j,k,nyzp(2,m)+2,m)
      scr(n,j,k+nyzp3,1,m) = f(n,j,k,nyzp(2,m)+3,m)
      scr(n,j,k+2*nyzp3,1,m) = f(n,j,k,1,m)
  440 continue
  450 continue
  460 continue
  470 continue
c add guard cells in z
      do 620 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 610 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(1,m) + 3
      nxvys = nxv*nyzp3
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      krr = kr + 1
      if (krr.ge.nvpz) krr = krr - nvpz
      krr = ky + nvpy*krr
      kr = ky + nvpy*kr
      kl = kz - 1
      if (kl.lt.0) kl = kl + nvpz
      kll = kl - 1
      if (kll.lt.0) kll = kll + nvpz
      kll = ky + nvpy*kll
      kl = ky + nvpy*kl
      ngc = 0
c special case of only one grid per processor
      if (nyzp(2,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (nyzp(2,kl).eq.1) then
c        do 400 k = 1, nyzp3
c        do 490 j = 1, nxv
c        do 480 n = 1, 3
c        scr(n,j,k,2,m) = scr(n,j,k,1,kl) + scr(n,j,k+nyzp3,1,kll)
c        scr(n,j,k+nyzp3,2,m) = scr(n,j,k+nyzp3,1,kl)
c        scr(n,j,k+2*nyzp3,2,m) = scr(n,j,k+2*nyzp3,1,kr)
c 480    continue
c 490    continue
c 500    continue
c     else
c        do 530 k = 1, nyzp3
c        do 520 j = 1, nxv
c        do 510 n = 1, 3
c        scr(n,j,k,2,m) = scr(n,j,k,1,kl)
c        scr(n,j,k+nyzp3,2,m) = scr(n,j,k+nyzp3,1,kl)
c        scr(n,j,k+2*nyzp3,2,m) = scr(n,j,k+2*nyzp3,1,kr)
c 510    continue
c 520    continue
c 530    continue
c     endif
c this segment is used for mpi computers
      call MPI_IRECV(scr(1,1,1,2,m),6*nxvy,mreal,kl-1,noff+5,lgrp,msid,i
     1err)
      call MPI_SEND(scr(1,1,1,1,m),6*nxvys,mreal,kr-1,noff+5,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scr(1,1,2*nyzp3+1,2,m),3*nxvy,mreal,kr-1,noff+6,lgr
     1p,msid,i
     1err)
      call MPI_SEND(scr(1,1,2*nyzp3+1,1,m),3*nxvys,mreal,kl-1,noff+6,lgr
     1p,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c special case of only one grid per processor in z
      if (nter.ge.1) go to 570
      call MPI_IRECV(ngc,1,mint,kl-1,noff+7,lgrp,msid,ierr)
      call MPI_IRECV(scr(1,1,1,1,m),3*nxvy,mreal,kll-1,noff+8,lgrp,nsid,
     1ierr)
      call MPI_SEND(nyzp(2,m),1,mint,kr-1,noff+7,lgrp,ierr)
      call MPI_SEND(scr(1,1,nyzp3+1,1,m),3*nxvys,mreal,krr-1,noff+8,lgrp
     1,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_WAIT(nsid,istatus,ierr)
      if (ngc.eq.1) then
         do 560 k = 1, nyzp3
         do 550 j = 1, nx3
         do 540 n = 1, 3
         scr(n,j,k,2,m) = scr(n,j,k,2,m) + scr(n,j,k,1,m)
  540    continue
  550    continue
  560    continue
      endif
c add up the guard cells
  570 do 600 k = 1, nyzp3
      do 590 j = 1, nx3
      do 580 n = 1, 3
      f(n,j,k,2,m) = f(n,j,k,2,m) + scr(n,j,k,2,m)
      f(n,j,k,3,m) = f(n,j,k,3,m) + scr(n,j,k+nyzp3,2,m)
      f(n,j,k,nyzp(2,m)+1,m) = f(n,j,k,nyzp(2,m)+1,m) + scr(n,j,k+2*nyzp
     13,2,m)
      f(n,j,k,1,m) = 0.
      f(n,j,k,nyzp(2,m)+2,m) = 0.
      f(n,j,k,nyzp(2,m)+3,m) = 0.
  580 continue
  590 continue
  600 continue
  610 continue
  620 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNMAGUARD32(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypmx
     1,nzpmx,mblok,nblok,ngds,idds,mter,nter)
c this subroutine adds data from guard cells in non-uniform partitions
c for scalar data, the field is added up so as to disable quadratic
c interpolation within half a cell of the edges, and reduce it to linear
c interpolation in the y direction, and periodic in z direction.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes three extra guard cells.
c scs/scr = scratch arrays for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c it is assumed the nyzp(n,m) > 0.
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c mter/nter = (0,1) = (no,yes) pass data to next processor only in y/z
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds, mter, nter
      integer nyzp
      real f, scs, scr
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx*ngds,2,mblok*nblok)
      dimension scr(nxv,nypmx*ngds,2,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, mnblok
      integer nx3, nxvz, nxvzs, nyzp3, nxvy, nxvys, m, my, mz, j, k
      integer jr, jrr, jl, jll
      dimension istatus(lstat)
      nx3 = nx + 3
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c buffer data in y
      do 30 m = 1, mnblok
      nyzp3 = nyzp(2,m) + 3
      do 20 k = 1, nyzp3
      do 10 j = 1, nxv
      scs(j,k,1,m) = f(j,nyzp(1,m)+2,k,m)
      scs(j,k+nyzp3,1,m) = f(j,nyzp(1,m)+3,k,m)
      scs(j,k+2*nyzp3,1,m) = f(j,1,k,m)
   10 continue
   20 continue
   30 continue
c add guard cells in y
      do 240 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 230 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(2,m) + 3
      nxvzs = nxv*nyzp3
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr + 1
      krr = jrr + kz
      jl = ky - 1
      jll = jl - 1
      kll = jll + kz
      kr = jr + kz
      kl = jl + kz
      ngc = 0
c special case of only one grid per processor
      if (nyzp(1,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        if (nyzp(1,kl).eq.1) then
c           do 50 k = 1, nyzp3
c           do 40 j = 1, nxv
c           scs(j,k,2,m) = scs(j,k,1,kl) + scs(j,k,1,kll)
c           scs(j,k+nyzp3,2,m) = scs(j,k+nyzp3,1,kl)
c  40       continue
c  50       continue
c        else
c           do 70 k = 1, nyzp3
c           do 60 j = 1, nxv
c           scs(j,k,2,m) = scs(j,k,1,kl)
c           scs(j,k+nyzp3,2,m) = scs(j,k+nyzp3,1,kl)=
c  60       continue
c  70       continue
c        endif
c     else
c        do 90 k = 1, nyzp3
c        do 80 j = 1, nx3
c        scs(j,k,2,m) = 2.*scs(j,k+2*nyzp3,1,m)
c        scs(j,k+nyzp3,2,m) = -scs(j,k+2*nyzp3,1,m)
c  80    continue
c  90    continue
c     endif
c     if (jr.lt.nvpy) then
c        do 110 k = 1, nyzp3
c        do 100 j = 1, nxv
c        scs(j,k+2*nyzp3,2,m) = scs(j,k+2*nyzp3,1,kr)
c 100    continue
c 110    continue
c     else
c        do 130 k = 1, nyzp3
c        do 120 j = 1, nx3
c        scs(j,k+2*nyzp3,2,m) = -scs(j,k+nyzp3,1,m)
c        f(j,nyzp(1,m)+2,k,m) = f(j,nyzp(1,m)+2,k,m) + 2.*scs(j,k+nyzp3,
c    11,m)
c        f(j,nyzp(1,m)+3,k,m) = 0.
c 120    continue
c 130    continue
c     endif
c     if (jl.ge.0) then
c        if ((nyzp(1,kl).eq.1).and.(jl.eq.0)) then
c           do 150 k = 1, nyzp3
c           do 140 j = 1, nx3
c           scs(j,k,2,m) = scs(j,k,2,m) - scs(j,k+2*nyzp3,1,kl)
c 140       continue
c 150       continue
c        endif
c     else if (nyzp(1,m).eq.1) then
c        do 170 k = 1, nyzp3
c        do 160 j = 1, nx3
c        scs(j,k+nyzp3,2,m) = 0.
c 160    continue
c 170    continue
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,2,m),2*nxvz,mreal,kl-1,noff+1,lgrp,msid,
     1ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,m),2*nxvzs,mreal,kr-1,noff+1,lgrp,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 50 k = 1, nyzp3
         do 40 j = 1, nx3
         scs(j,k,2,m) = 2.*scs(j,k+2*nyzp3,1,m)
         scs(j,k+nyzp3,2,m) = -scs(j,k+2*nyzp3,1,m)
   40    continue
   50    continue
      endif
      if (jr.lt.nvpy) then
         call MPI_IRECV(scs(1,2*nyzp3+1,2,m),nxvz,mreal,kr-1,noff+2,lgrp
     1,msid,ierr)
      endif
      if (jl.ge.0) then
         call MPI_SEND(scs(1,2*nyzp3+1,1,m),nxvzs,mreal,kl-1,noff+2,lgrp
     1,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 70 k = 1, nyzp3
         do 60 j = 1, nx3
         scs(j,k+2*nyzp3,2,m) = -scs(j,k+nyzp3,1,m)
         f(j,nyzp(1,m)+2,k,m) = f(j,nyzp(1,m)+2,k,m) + 2.*scs(j,k+nyzp3,
     11,m)
         f(j,nyzp(1,m)+3,k,m) = 0.
   60    continue
   70    continue
      endif
c special case of only one grid per processor in y
      if (mter.ge.1) go to 180
      if (jl.ge.0) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+3,lgrp,msid,ierr)
      endif
      if (jll.ge.0) then
         call MPI_IRECV(scs(1,1,1,m),nxvz,mreal,kll-1,noff+5,lgrp,nsid,i
     1err)
      else if (jl.eq.0) then
         call MPI_IRECV(scs(1,1,1,m),nxvz,mreal,kl-1,noff+5,lgrp,nsid,ie
     1rr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(nyzp(1,m),1,mint,kr-1,moff+3,lgrp,ierr)
      endif
      if (jrr.lt.nvpy) then
         call MPI_SEND(scs(1,nyzp3+1,1,m),nxvzs,mreal,krr-1,noff+5,lgrp,
     1ierr)
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
         call MPI_SEND(scs(1,2*nyzp3+1,1,m),nxvzs,mreal,kr-1,noff+5,lgrp
     1,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         if (ngc.eq.1) then
            if (jl.eq.0) then
               do 90 k = 1, nyzp3
               do 80 j = 1, nx3
               scs(j,k,2,m) = scs(j,k,2,m) - scs(j,k,1,m)
   80          continue
   90          continue
            else
               do 110 k = 1, nyzp3
               do 100 j = 1, nx3
               scs(j,k,2,m) = scs(j,k,2,m) + scs(j,k,1,m)
  100          continue
  110          continue
            endif
         else if (nyzp(1,m).eq.1) then
            do 130 k = 1, nyzp3
            do 120 j = 1, nx3
            scs(j,k+nyzp3,2,m) = 0.
  120       continue
  130       continue
         endif
      endif
c add up the guard cells
  180 do 200 k = 1, nyzp3
      do 190 j = 1, nx3
      f(j,2,k,m) = f(j,2,k,m) + scs(j,k,2,m)
      f(j,3,k,m) = f(j,3,k,m) + scs(j,k+nyzp3,2,m)
      f(j,nyzp(1,m)+1,k,m) = f(j,nyzp(1,m)+1,k,m) + scs(j,k+2*nyzp3,2,m)
      f(j,1,k,m) = 0.
      f(j,nyzp(1,m)+3,k,m) = 0.
  190 continue
  200 continue
      if (jr.lt.nvpy) then
         do 220 k = 1, nyzp3
         do 210 j = 1, nx3
         f(j,nyzp(1,m)+2,k,m) = 0.
  210    continue
  220    continue
      endif
  230 continue
  240 continue
c zero out the left edge
      do 280 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 270 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(2,m) + 3
      kl = my + js
      if (kl.eq.0) then
         do 260 k = 1, nyzp3
         do 250 j = 1, nx3
         f(j,1,k,m) = 0.
  250    continue
  260    continue
      endif
  270 continue
  280 continue
c special case for one processor in z
      if (nvpz.eq.1) then
         do 310 m = 1, mnblok
         nyzp3 = nyzp(1,m) + 3
         do 300 k = 1, nyzp3
         do 290 j = 1, nx3
         f(j,k,2,m) = f(j,k,2,m) + f(j,k,nyzp(2,m)+2,m)
         f(j,k,3,m) = f(j,k,3,m) + f(j,k,nyzp(2,m)+3,m)
         f(j,k,nyzp(2,m)+1,m) = f(j,k,nyzp(2,m)+1,m) + f(j,k,1,m)
         f(j,k,1,m) = 0.
         f(j,k,nyzp(2,m)+2,m) = 0.
         f(j,k,nyzp(2,m)+3,m) = 0.
  290    continue
  300    continue
  310    continue
         return
      endif
c buffer data in z
      do 340 m = 1, mnblok
      nyzp3 = nyzp(1,m) + 3
      do 330 k = 1, nyzp3
      do 320 j = 1, nxv
      scr(j,k,1,m) = f(j,k,nyzp(2,m)+2,m)
      scr(j,k+nyzp3,1,m) = f(j,k,nyzp(2,m)+3,m)
      scr(j,k+2*nyzp3,1,m) = f(j,k,1,m)
  320 continue
  330 continue
  340 continue
c add guard cells in z
      do 450 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 440 my = 1, mblok
      m = my + moff
      nyzp3 = nyzp(1,m) + 3
      nxvys = nxv*nyzp3
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      krr = kr + 1
      if (krr.ge.nvpz) krr = krr - nvpz
      krr = ky + nvpy*krr
      kr = ky + nvpy*kr
      kl = kz - 1
      if (kl.lt.0) kl = kl + nvpz
      kll = kl - 1
      if (kll.lt.0) kll = kll + nvpz
      kll = ky + nvpy*kll
      kl = ky + nvpy*kl
      ngc = 0
c special case of only one grid per processor
      if (nyzp(2,m).eq.1) ngc = 1
c this segment is used for shared memory computers
c     if (nyzp(2,kl).eq.1) then
c        do 360 k = 1, nyzp3
c        do 350 j = 1, nxv
c        scr(j,k,2,m) = scr(j,k,1,kl) + scr(j,k+nyzp3,1,kll)
c        scr(j,k+nyzp3,2,m) = scr(j,k+nyzp3,1,kl)
c        scr(j,k+2*nyzp3,2,m) = scr(j,k+2*nyzp3,1,kr)
c 350    continue
c 360    continue
c     else
c        do 380 k = 1, nyzp3
c        do 370 j = 1, nxv
c        scr(j,k,2,m) = scr(j,k,1,kl)
c        scr(j,k+nyzp3,2,m) = scr(j,k+nyzp3,1,kl)
c        scr(j,k+2*nyzp3,2,m) = scr(j,k+2*nyzp3,1,kr)
c 370    continue
c 380    continue
c     endif
c this segment is used for mpi computers
      call MPI_IRECV(scr(1,1,2,m),2*nxvy,mreal,kl-1,noff+5,lgrp,msid,ier
     1r)
      call MPI_SEND(scr(1,1,1,m),2*nxvys,mreal,kr-1,noff+5,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_IRECV(scr(1,2*nyzp3+1,2,m),nxvy,mreal,kr-1,noff+6,lgrp,ms
     1id,ierr)
      call MPI_SEND(scr(1,2*nyzp3+1,1,m),nxvys,mreal,kl-1,noff+6,lgrp,ie
     1rr)
      call MPI_WAIT(msid,istatus,ierr)
c special case of only one grid per processor in z
      if (nter.ge.1) go to 410
      call MPI_IRECV(ngc,1,mint,kl-1,noff+7,lgrp,msid,ierr)
      call MPI_IRECV(scr(1,1,1,m),nxvy,mreal,kll-1,noff+8,lgrp,nsid,ierr
     1)
      call MPI_SEND(nyzp(2,m),1,mint,kr-1,noff+7,lgrp,ierr)
      call MPI_SEND(scr(1,nyzp3+1,1,m),nxvys,mreal,krr-1,noff+8,lgrp,ier
     1r)
      call MPI_WAIT(msid,istatus,ierr)
      call MPI_WAIT(nsid,istatus,ierr)
      if (ngc.eq.1) then
         do 400 k = 1, nyzp3
         do 390 j = 1, nx3
         scr(j,k,2,m) = scr(j,k,2,m) + scr(j,k,1,m)
  390    continue
  400    continue
      endif
c add up the guard cells
  410 do 430 k = 1, nyzp3
      do 420 j = 1, nx3
      f(j,k,2,m) = f(j,k,2,m) + scr(j,k,2,m)
      f(j,k,3,m) = f(j,k,3,m) + scr(j,k+nyzp3,2,m)
      f(j,k,nyzp(2,m)+1,m) = f(j,k,nyzp(2,m)+1,m) + scr(j,k+2*nyzp3,2,m)
      f(j,k,1,m) = 0.
      f(j,k,nyzp(2,m)+2,m) = 0.
      f(j,k,nyzp(2,m)+3,m) = 0.
  420 continue
  430 continue
  440 continue
  450 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMACGUARDS32(f,scs,kstrt,nvpy,nx,nxv,nypmx,nzpmx,mblok,
     1nblok,kyp,kzp,ngds)
c this subroutine corrects the charge density data for particle boundary
c conditions which keep particles one grid away from the edges
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y direction, and periodic in the z direction.
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes three extra
c guard cells.
c kstrt = starting data block number
c nvpy = number of real or virtual processors in y
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      real f, scs
      integer kstrt, nvpy, nx, nxv, nypmx, nzpmx, mblok, nblok, kyp, kzp
      integer ngds
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx,ngds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, mnblok
      integer kzp1, nxvz, m, my, mz, j, k, n
      integer jr, jrr, jl, jll
      dimension istatus(lstat)
      kzp1 = kzp + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
c buffer data in y
      do 40 m = 1, mnblok
      do 30 k = 1, nzpmx
      do 20 j = 1, nxv
      do 10 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,kyp+2,k,m)
      scs(n,j,k,2,m) = f(n,j,2,k,m)
   10 continue
   20 continue
   30 continue
   40 continue
c add guard cells in y
      do 360 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 350 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr
      krr = jr + kz
      jl = ky - 1
      jll = jl
      kll = jl + kz
c special case of only one grid per processor
      if (kyp.eq.1) then
         jrr = jr + 1
         krr = jrr + kz
         jll = jl - 1
         kll = jll + kz
      endif
      kr = jr + kz
      kl = jl + kz
c fix edges if all points are on the same processor
      if (jl.eq.(-1)) then
         if (kyp.gt.2) then
            do 70 k = 1, kzp1
            do 60 j = 2, nx
            do 50 n = 1, 3
            f(n,j+1,3,k+1,m) = f(n,j+1,3,k+1,m) + 2.*f(n,j+1,2,k+1,m)
            f(n,j+1,4,k+1,m) = f(n,j+1,4,k+1,m) - f(n,j+1,2,k+1,m)
            f(n,j+1,2,k+1,m) = 0.
   50       continue
   60       continue
   70       continue
         else if (kyp.eq.2) then
            do 100 k = 1, kzp1
            do 90 j = 2, nx
            do 80 n = 1, 3
            f(n,j+1,3,k+1,m) = f(n,j+1,3,k+1,m) + 2.*f(n,j+1,2,k+1,m)
   80       continue
   90       continue
  100       continue
         endif
      endif
      if (jr.eq.nvpy) then
         if (kyp.gt.1) then
            do 130 k = 1, kzp1
            do 120 j = 2, nx
            do 110 n = 1, 3
            f(n,j+1,kyp,k+1,m) = f(n,j+1,kyp,k+1,m) - f(n,j+1,kyp+2,k+1,
     1m)
            f(n,j+1,kyp+1,k+1,m) = f(n,j+1,kyp+1,k+1,m) + 2.*f(n,j+1,kyp
     1+2,k+1,m)
            f(n,j+1,kyp+2,k+1,m) = 0.
  110       continue
  120       continue
  130       continue
         else if (kyp.eq.1) then
            do 160 k = 1, kzp1
            do 150 j = 2, nx
            do 140 n = 1, 3
            f(n,j+1,kyp+1,k+1,m) = f(n,j+1,kyp+1,k+1,m) + 2.*f(n,j+1,kyp
     1+2,k+1,m)
            f(n,j+1,kyp+2,k+1,m) = 0.
  140       continue
  150       continue
  160       continue
         endif
      endif
c this segment is used for shared memory computers
c     if (kyp.eq.2) then
c        if (jl.eq.0) then
c           do 190 k = 1, kzp1
c           do 180 j = 2, nx
c           do 170 n = 1, 3
c           f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) - scs(n,j+1,k+1,2,kl)
c 170       continue
c 180       continue
c 190       continue
c        endif
c     else if (kyp.eq.1) then
c        if (jl.eq.0) then
c           do 220 k = 1, kzp1
c           do 210 j = 2, nx
c           do 200 n = 1, 3
c           f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) + 2.*scs(n,j+1,k+1,2,kl)
c 200       continue
c 210       continue
c 220       continue
c        endif
c        if (jll.eq.0) then
c           do 250 k = 1, kzp1
c           do 240 j = 2, nx
c           do 230 n = 1, 3
c           f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) - scs(n,j+1,k+1,2,kll)
c 230       continue
c 240       continue
c 250       continue
c        endif
c        if (jr.eq.(nvpy-1)) then
c           do 280 k = 1, kzp1
c           do 270 j = 2, nx
c           do 260 n = 1, 3
c           f(n,j+1,kyp+1,k+1,m) = f(n,j+1,kyp+1,k+1,m) - scs(n,j+1,k+1,
c    11,kr)
c 260       continue
c 270       continue
c 280       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (kyp.eq.2) then
         if (jl.eq.0) then
            call MPI_IRECV(scs(1,1,1,3,m),3*nxvz,mreal,kl-1,noff+1,lgrp,
     1msid,ierr)
         endif
         if (jl.eq.(-1)) then
            call MPI_SEND(scs(1,1,1,2,m),3*nxvz,mreal,kr-1,noff+1,lgrp,i
     1err)
            do 190 k = 1, kzp1
            do 180 j = 2, nx
            do 170 n = 1, 3
            f(n,j+1,2,k+1,m) = 0.
  170       continue
  180       continue
  190       continue
         endif
         if (jl.eq.0) then
            call MPI_WAIT(msid,istatus,ierr)
            do 220 k = 1, kzp1
            do 210 j = 2, nx
            do 200 n = 1, 3
            f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) - scs(n,j+1,k+1,3,m)
            f(n,j+1,1,k+1,m) = 0.
  200       continue
  210       continue
  220       continue
         endif
      else if (kyp.eq.1) then
         if (jl.eq.0) then
            call MPI_IRECV(scs(1,1,1,3,m),3*nxvz,mreal,kl-1,noff+1,lgrp,
     1msid,ierr)
         endif
         if (jll.eq.0) then
            call MPI_IRECV(scs(1,1,1,3,m),3*nxvz,mreal,kll-1,noff+1,lgrp
     1,nsid,ierr)
         endif
         if (jl.eq.(-1)) then
            call MPI_SEND(scs(1,1,1,2,m),3*nxvz,mreal,kr-1,noff+1,lgrp,i
     1err)
            call MPI_SEND(scs(1,1,1,2,m),3*nxvz,mreal,krr-1,noff+1,lgrp,
     1ierr)
            do 250 k = 1, kzp1
            do 240 j = 2, nx
            do 230 n = 1, 3
            f(n,j+1,2,k+1,m) = 0.
  230       continue
  240       continue
  250       continue
         endif
         if (jl.eq.0) then
            call MPI_WAIT(msid,istatus,ierr)
            do 280 k = 1, kzp1
            do 270 j = 2, nx
            do 260 n = 1, 3
            f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) + 2.*scs(n,j+1,k+1,3,m)
            f(n,j+1,1,k+1,m) = 0.
  260       continue
  270       continue
  280       continue
         endif
         if (jll.eq.0) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 310 k = 1, kzp1
            do 300 j = 2, nx
            do 290 n = 1, 3
            f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) - scs(n,j+1,k+1,3,m)
            f(n,j+1,1,k+1,m) = 0.
  290       continue
  300       continue
  310       continue
         endif
         if (jr.eq.(nvpy-1)) then
            call MPI_IRECV(scs(1,1,1,3,m),3*nxvz,mreal,kr-1,noff+2,lgrp,
     1msid,ierr)
         endif
         if (jr.eq.nvpy) then
            call MPI_SEND(scs(1,1,1,1,m),3*nxvz,mreal,kl-1,noff+2,lgrp,i
     1err)
         endif
         if (jr.eq.(nvpy-1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 340 k = 1, kzp1
            do 330 j = 2, nx
            do 320 n = 1, 3
            f(n,j+1,kyp+1,k+1,m) = f(n,j+1,kyp+1,k+1,m) - scs(n,j+1,k+1,
     13,m)
            f(n,j+1,1,k+1,m) = 0.
  320       continue
  330       continue
  340       continue
         endif
      endif
  350 continue
  360 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMAGUARDS32(f,scs,kstrt,nvpy,nx,nxv,nypmx,nzpmx,mblok,n
     1blok,kyp,kzp,ngds)
c this subroutine corrects the charge density data for particle boundary
c conditions which keep particles one grid away from the edges
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y direction, and periodic in the z direction.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes three extra
c guard cells.
c kstrt = starting data block number
c nvpy = number of real or virtual processors in y
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      real f, scs
      integer kstrt, nvpy, nx, nxv, nypmx, nzpmx, mblok, nblok, kyp, kzp
      integer ngds
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,ngds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, mnblok
      integer kzp1, nxvz, m, my, mz, j, k
      integer jr, jrr, jl, jll
      dimension istatus(lstat)
      kzp1 = kzp + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
c buffer data in y
      do 30 m = 1, mnblok
      do 20 k = 1, nzpmx
      do 10 j = 1, nxv
      scs(j,k,1,m) = f(j,kyp+2,k,m)
      scs(j,k,2,m) = f(j,2,k,m)
   10 continue
   20 continue
   30 continue
c add guard cells in y
      do 250 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 240 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr
      krr = jr + kz
      jl = ky - 1
      jll = jl
      kll = jl + kz
c special case of only one grid per processor
      if (kyp.eq.1) then
         jrr = jr + 1
         krr = jrr + kz
         jll = jl - 1
         kll = jll + kz
      endif
      kr = jr + kz
      kl = jl + kz
c fix edges if all points are on the same processor
      if (jl.eq.(-1)) then
         if (kyp.gt.2) then
            do 50 k = 1, kzp1
            do 40 j = 2, nx
            f(j+1,3,k+1,m) = f(j+1,3,k+1,m) + 2.*f(j+1,2,k+1,m)
            f(j+1,4,k+1,m) = f(j+1,4,k+1,m) - f(j+1,2,k+1,m)
            f(j+1,2,k+1,m) = 0.
   40       continue
   50       continue
         else if (kyp.eq.2) then
            do 70 k = 1, kzp1
            do 60 j = 2, nx
            f(j+1,3,k+1,m) = f(j+1,3,k+1,m) + 2.*f(j+1,2,k+1,m)
   60       continue
   70       continue
         endif
      endif
      if (jr.eq.nvpy) then
         if (kyp.gt.1) then
            do 90 k = 1, kzp1
            do 80 j = 2, nx
            f(j+1,kyp,k+1,m) = f(j+1,kyp,k+1,m) - f(j+1,kyp+2,k+1,m)
            f(j+1,kyp+1,k+1,m) = f(j+1,kyp+1,k+1,m) + 2.*f(j+1,kyp+2,k+1
     1,m)
            f(j+1,kyp+2,k+1,m) = 0.
   80       continue
   90       continue
         else if (kyp.eq.1) then
            do 110 k = 1, kzp1
            do 100 j = 2, nx
            f(j+1,kyp+1,k+1,m) = f(j+1,kyp+1,k+1,m) + 2.*f(j+1,kyp+2,k+1
     1,m)
            f(j+1,kyp+2,k+1,m) = 0.
  100       continue
  110       continue
         endif
      endif
c this segment is used for shared memory computers
c     if (kyp.eq.2) then
c        if (jl.eq.0) then
c           do 130 k = 1, kzp1
c           do 120 j = 2, nx
c           f(j+1,2,k+1,m) = f(j+1,2,k+1,m) - scs(j+1,k+1,2,kl)
c 120       continue
c 130       continue
c        endif
c     else if (kyp.eq.1) then
c        if (jl.eq.0) then
c           do 150 k = 1, kzp1
c           do 140 j = 2, nx
c           f(j+1,2,k+1,m) = f(j+1,2,k+1,m) + 2.*scs(j+1,k+1,2,kl)
c 140       continue
c 150       continue
c        endif
c        if (jll.eq.0) then
c           do 170 k = 1, kzp1
c           do 160 j = 2, nx
c           f(j+1,2,k+1,m) = f(j+1,2,k+1,m) - scs(j+1,k+1,2,kll)
c 160       continue
c 170       continue
c        endif
c        if (jr.eq.(nvpy-1)) then
c           do 190 k = 1, kzp1
c           do 180 j = 2, nx
c           f(j+1,kyp+1,k+1,m) = f(j+1,kyp+1,k+1,m) - scs(j+1,k+1,1,kr)
c 180       continue
c 190       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (kyp.eq.2) then
         if (jl.eq.0) then
            call MPI_IRECV(scs(1,1,3,m),nxvz,mreal,kl-1,noff+1,lgrp,msid
     1,ierr)
         endif
         if (jl.eq.(-1)) then
            call MPI_SEND(scs(1,1,2,m),nxvz,mreal,kr-1,noff+1,lgrp,ierr)
            do 130 k = 1, kzp1
            do 120 j = 2, nx
            f(j+1,2,k+1,m) = 0.
  120       continue
  130       continue
         endif
         if (jl.eq.0) then
            call MPI_WAIT(msid,istatus,ierr)
            do 150 k = 1, kzp1
            do 140 j = 2, nx
            f(j+1,2,k+1,m) = f(j+1,2,k+1,m) - scs(j+1,k+1,3,m)
            f(j+1,1,k+1,m) = 0.
  140       continue
  150       continue
         endif
      else if (kyp.eq.1) then
         if (jl.eq.0) then
            call MPI_IRECV(scs(1,1,3,m),nxvz,mreal,kl-1,noff+1,lgrp,msid
     1,ierr)
         endif
         if (jll.eq.0) then
            call MPI_IRECV(scs(1,1,3,m),nxvz,mreal,kll-1,noff+1,lgrp,nsi
     1d,ierr)
         endif
         if (jl.eq.(-1)) then
            call MPI_SEND(scs(1,1,2,m),nxvz,mreal,kr-1,noff+1,lgrp,ierr)
            call MPI_SEND(scs(1,1,2,m),nxvz,mreal,krr-1,noff+1,lgrp,ierr
     1)
            do 170 k = 1, kzp1
            do 160 j = 2, nx
            f(j+1,2,k+1,m) = 0.
  160       continue
  170       continue
         endif
         if (jl.eq.0) then
            call MPI_WAIT(msid,istatus,ierr)
            do 190 k = 1, kzp1
            do 180 j = 2, nx
            f(j+1,2,k+1,m) = f(j+1,2,k+1,m) + 2.*scs(j+1,k+1,3,m)
            f(j+1,1,k+1,m) = 0.
  180       continue
  190       continue
         endif
         if (jll.eq.0) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 210 k = 1, kzp1
            do 200 j = 2, nx
            f(j+1,2,k+1,m) = f(j+1,2,k+1,m) - scs(j+1,k+1,3,m)
            f(j+1,1,k+1,m) = 0.
  200       continue
  210       continue
         endif
         if (jr.eq.(nvpy-1)) then
            call MPI_IRECV(scs(1,1,3,m),nxvz,mreal,kr-1,noff+2,lgrp,msid
     1,ierr)
         endif
         if (jr.eq.nvpy) then
            call MPI_SEND(scs(1,1,1,m),nxvz,mreal,kl-1,noff+2,lgrp,ierr)
         endif
         if (jr.eq.(nvpy-1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 230 k = 1, kzp1
            do 220 j = 2, nx
            f(j+1,kyp+1,k+1,m) = f(j+1,kyp+1,k+1,m) - scs(j+1,k+1,3,m)
            f(j+1,1,k+1,m) = 0.
  220       continue
  230       continue
         endif
      endif
  240 continue
  250 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNMACGUARDS32(f,scs,nyzp,kstrt,nvpy,nx,nxv,nypmx,nzpmx,
     1mblok,nblok,ngds,idds,mter)
c this subroutine corrects current density data for particle boundary
c conditions which keep particles one grid away from the edges
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y direction.
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes three extra guard cells.
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c it is assumed the nyzp(n,m) > 0.
c kstrt = starting data block number
c nvpy = number of real or virtual processors in y
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c mter = (0,1) = (no,yes) pass data to next processor only in y
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      real f, scs
      integer kstrt, nvpy, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds, mter
      integer nyzp
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx,ngds,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, nps
      integer mnblok, kzp1, nxvz, nxvzs, nxvy
      integer m, my, mz, j, k, n, jr, jrr, jl, jll
      dimension istatus(lstat)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c buffer data in y
      do 40 m = 1, mnblok
      kzp1 = nyzp(2,m) + 1
      do 30 k = 1, kzp1 + 1
      do 20 j = 1, nxv
      do 10 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,nyzp(1,m)+2,k,m)
      scs(n,j,k,2,m) = f(n,j,2,k,m)
   10 continue
   20 continue
   30 continue
   40 continue
c add guard cells in y
      do 400 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 390 my = 1, mblok
      m = my + moff
      kzp1 = nyzp(2,m) + 1
      nxvzs = nxv*(kzp1 + 1)
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr + 1
      krr = jrr + kz
      jl = ky - 1
      jll = jl - 1
      kll = jll + kz
      kr = jr + kz
      kl = jl + kz
c fix edges if all points are on the same processor
      if (jl.eq.(-1)) then
         if (nyzp(1,m).gt.2) then
            do 70 k = 1, kzp1
            do 60 j = 2, nx
            do 50 n = 1, 3
            f(n,j+1,3,k+1,m) = f(n,j+1,3,k+1,m) + 2.*f(n,j+1,2,k+1,m)
            f(n,j+1,4,k+1,m) = f(n,j+1,4,k+1,m) - f(n,j+1,2,k+1,m)
            f(n,j+1,2,k+1,m) = 0.
   50       continue
   60       continue
   70       continue
         else if (nyzp(1,m).eq.2) then
            do 100 k = 1, kzp1
            do 90 j = 2, nx
            do 80 n = 1, 3
            f(n,j+1,3,k+1,m) = f(n,j+1,3,k+1,m) + 2.*f(n,j+1,2,k+1,m)
   80       continue
   90       continue
  100       continue
         endif
      endif
      if (jr.eq.nvpy) then
         if (nyzp(1,m).gt.1) then
            do 130 k = 1, kzp1
            do 120 j = 2, nx
            do 110 n = 1, 3
            f(n,j+1,nyzp(1,m),k+1,m) = f(n,j+1,nyzp(1,m),k+1,m) - f(n,j+
     11,nyzp(1,m)+2,k+1,m)
            f(n,j+1,nyzp(1,m)+1,k+1,m) = f(n,j+1,nyzp(1,m)+1,k+1,m) + 2.
     1*f(n,j+1,nyzp(1,m)+2,k+1,m)
            f(n,j+1,nyzp(1,m)+2,k+1,m) = 0.
  110       continue
  120       continue
  130       continue
         else if (nyzp(1,m).eq.1) then
            do 160 k = 1, kzp1
            do 150 j = 2, nx
            do 140 n = 1, 3
            f(n,j+1,nyzp(1,m)+1,k+1,m) = f(n,j+1,nyzp(1,m)+1,k+1,m) + 2.
     1*f(n,j+1,nyzp(1,m)+2,k+1,m)
            f(n,j+1,nyzp(1,m)+2,k+1,m) = 0.
  140       continue
  150       continue
  160       continue
         endif
      endif
c this segment is used for shared memory computers
c     if (mter.ge.2) go to 390
c     if (mter.eq.1) go to 260
c     if (jll.eq.0) then
c        if (nyzp(1,kll).eq.1).and.(nyzp(1,kl).eq.1)) then
c           do 190 k = 1, kzp1
c           do 180 j = 2, nx
c           do 170 n = 1, 3
c           f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) - scs(n,j+1,k+1,2,kll)
c 170       continue
c 180       continue
c 190       continue
c        endif
c     endif
c     if (jr.eq.(nvpy-1)) then
c        if (nyzp(1,kr).eq.1) then
c           do 250 k = 1, kzp1
c           do 240 j = 2, nx
c           do 230 n = 1, 3
c           f(n,j+1,nyzp(1,m)+1,k+1,m) = f(n,j+1,nyzp(1,m)+1,k+1,m) - sc
c    1s(n,j+1,k+1,1,kr)
c 230       continue
c 240       continue
c 250       continue
c        endif
c     endif
c 260 if (jl.eq.0) then
c        if (nyzp(1,kl).le.2) then
c           if (nyzp(1,kl).eq.1) then
c              if (nyzp(1,m).gt.1) then
c                 do 320 k = 1, kzp1
c                 do 310 j = 2, nx
c                 do 300 n = 1, 3
c                 f(n,j+1,3,k+1,m) = f(n,j+1,3,k+1,m) - scs(n,j+1,k+1,2,
c    1kl)
c 300             continue
c 310             continue
c 320             continue
c              endif
c              do 350 k = 1, kzp1
c              do 340 j = 2, nx
c              do 330 n = 1, 3
c              f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) + 2.*scs(n,j+1,k+1,2,
c    1kl)
c 330          continue
c 340          continue
c 350          continue
c           else
c              do 380 k = 1, kzp1
c              do 370 j = 2, nx
c              do 360 n = 1, 3
c              f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) - scs(n,j+1,k+1,2,kl)
c 360          continue
c 370          continue
c 380          continue
c           endif
c        endif
c     endif
c this segment is used for mpi computers
      if (mter.ge.2) go to 390
      if (mter.eq.1) go to 260
      if (jll.eq.0) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+4,lgrp,msid,ierr)
         call MPI_IRECV(scs(1,1,1,3,m),3*nxvz,mreal,kll-1,moff+1,lgrp,ns
     1id,ierr)
      endif
      if (jl.eq.0) then
         call MPI_SEND(nyzp(1,m),1,mint,kr-1,moff+4,lgrp,ierr)
      endif
      if (jl.eq.(-1)) then
         if (nyzp(1,m).eq.1) then
            call MPI_SEND(scs(1,1,1,2,m),3*nxvzs,mreal,krr-1,moff+1,lgrp
     1,ierr)
         else
            nps = 0
            call MPI_SEND(scs(1,1,1,2,m),nps,mreal,krr-1,moff+1,lgrp,ier
     1r)
         endif
      endif
      if (jll.eq.0) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            if (ngc.eq.1) then
               do 190 k = 1, kzp1
               do 180 j = 2, nx
               do 170 n = 1, 3
               f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) - scs(n,j+1,k+1,3,m)
  170          continue
  180          continue
  190          continue
            endif
            do 220 k = 1, kzp1
            do 210 j = 2, nx
            do 200 n = 1, 3
            f(n,j+1,1,k+1,m) = 0.
  200       continue
  210       continue
  220       continue
         endif
      endif
      if (jr.eq.(nvpy-1)) then
         call MPI_IRECV(scs(1,1,1,3,m),3*nxvz,mreal,kr-1,moff+2,lgrp,msi
     1d,ierr)
      endif
      if (jr.eq.nvpy) then
         if (nyzp(1,m).eq.1) then
            call MPI_SEND(scs(1,1,1,1,m),3*nxvzs,mreal,kl-1,moff+2,lgrp,
     1ierr)
         else
            nps = 0
            call MPI_SEND(scs(1,1,1,1,m),nps,mreal,kl-1,moff+2,lgrp,ierr
     1)
         endif
      endif
      if (jr.eq.(nvpy-1)) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            do 250 k = 1, kzp1
            do 240 j = 2, nx
            do 230 n = 1, 3
            f(n,j+1,nyzp(1,m)+1,k+1,m) = f(n,j+1,nyzp(1,m)+1,k+1,m) - sc
     1s(n,j+1,k+1,3,m)
            f(n,j+1,1,k+1,m) = 0.
  230       continue
  240       continue
  250       continue
         endif
      endif
  260 if (jl.eq.0) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+5,lgrp,msid,ierr)
         call MPI_IRECV(scs(1,1,1,3,m),3*nxvz,mreal,kl-1,moff+1,lgrp,nsi
     1d,ierr)
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
         call MPI_SEND(nyzp(1,m),1,mint,kr-1,moff+5,lgrp,ierr)
         if (nyzp(1,m).le.2) then
            call MPI_SEND(scs(1,1,1,2,m),3*nxvzs,mreal,kr-1,moff+1,lgrp,
     1ierr)
            do 290 k = 1, kzp1
            do 280 j = 2, nx
            do 270 n = 1, 3
            f(n,j+1,2,k+1,m) = 0.
  270       continue
  280       continue
  290       continue
         else
            nps = 0
            call MPI_SEND(scs(1,1,1,2,m),nps,mreal,kr-1,moff+1,lgrp,ierr
     1)
         endif
      endif
      if (jl.eq.0) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            if (ngc.eq.1) then
               if (nyzp(1,m).gt.1) then
                  do 320 k = 1, kzp1
                  do 310 j = 2, nx
                  do 300 n = 1, 3
                  f(n,j+1,3,k+1,m) = f(n,j+1,3,k+1,m) - scs(n,j+1,k+1,3,
     1m)
  300             continue
  310             continue
  320             continue
               endif
               do 350 k = 1, kzp1
               do 340 j = 2, nx
               do 330 n = 1, 3
               f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) + 2.*scs(n,j+1,k+1,3,
     1m)
               f(n,j+1,1,k+1,m) = 0.
  330          continue
  340          continue
  350          continue
            else
               do 380 k = 1, kzp1
               do 370 j = 2, nx
               do 360 n = 1, 3
               f(n,j+1,2,k+1,m) = f(n,j+1,2,k+1,m) - scs(n,j+1,k+1,3,m)
               f(n,j+1,1,k+1,m) = 0.
  360          continue
  370          continue
  380          continue
            endif
         endif
      endif
  390 continue
  400 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNMAGUARDS32(f,scs,nyzp,kstrt,nvpy,nx,nxv,nypmx,nzpmx,m
     1blok,nblok,ngds,idds,mter)
c this subroutine corrects current density data for particle boundary
c conditions which keep particles one grid away from the edges
c the field is added up so as to disable quadratic interpolation
c within half a cell of the edges, and reduce it to linear interpolation
c in the y direction.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes three extra guard cells.
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c it is assumed the nyzp(n,m) > 0.
c kstrt = starting data block number
c nvpy = number of real or virtual processors in y
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c mter = (0,1) = (no,yes) pass data to next processor only in y
c quadratic interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      real f, scs
      integer kstrt, nvpy, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds, mter
      integer nyzp
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,ngds,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer ky, kz, js, ks, moff, noff, kr, krr, kl, kll, ngc, nps
      integer mnblok, kzp1, nxvz, nxvzs, nxvy, m, my, mz
      integer j, k, jr, jrr, jl, jll
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
      kzp1 = nyzp(2,m) + 1
      do 20 k = 1, kzp1 + 1
      do 10 j = 1, nxv
      scs(j,k,1,m) = f(j,nyzp(1,m)+2,k,m)
      scs(j,k,2,m) = f(j,2,k,m)
   10 continue
   20 continue
   30 continue
c add guard cells in y
      do 280 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 270 my = 1, mblok
      m = my + moff
      kzp1 = nyzp(2,m) + 1
      nxvzs = nxv*(kzp1 + 1)
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jrr = jr + 1
      krr = jrr + kz
      jl = ky - 1
      jll = jl - 1
      kll = jll + kz
      kr = jr + kz
      kl = jl + kz
c fix edges if all points are on the same processor
      if (jl.eq.(-1)) then
         if (nyzp(1,m).gt.2) then
            do 50 k = 1, kzp1
            do 40 j = 2, nx
            f(j+1,3,k+1,m) = f(j+1,3,k+1,m) + 2.*f(j+1,2,k+1,m)
            f(j+1,4,k+1,m) = f(j+1,4,k+1,m) - f(j+1,2,k+1,m)
            f(j+1,2,k+1,m) = 0.
   40       continue
   50       continue
         else if (nyzp(1,m).eq.2) then
            do 70 k = 1, kzp1
            do 60 j = 2, nx
            f(j+1,3,k+1,m) = f(j+1,3,k+1,m) + 2.*f(j+1,2,k+1,m)
   60       continue
   70       continue
         endif
      endif
      if (jr.eq.nvpy) then
         if (nyzp(1,m).gt.1) then
            do 90 k = 1, kzp1
            do 80 j = 2, nx
            f(j+1,nyzp(1,m),k+1,m) = f(j+1,nyzp(1,m),k+1,m) - f(j+1,nyzp
     1(1,m)+2,k+1,m)
            f(j+1,nyzp(1,m)+1,k+1,m) = f(j+1,nyzp(1,m)+1,k+1,m) + 2.*f(j
     1+1,nyzp(1,m)+2,k+1,m)
            f(j+1,nyzp(1,m)+2,k+1,m) = 0.
   80       continue
   90       continue
         else if (nyzp(1,m).eq.1) then
            do 110 k = 1, kzp1
            do 100 j = 2, nx
            f(j+1,nyzp(1,m)+1,k+1,m) = f(j+1,nyzp(1,m)+1,k+1,m) + 2.*f(j
     1+1,nyzp(1,m)+2,k+1,m)
            f(j+1,nyzp(1,m)+2,k+1,m) = 0.
  100       continue
  110       continue
         endif
      endif
c this segment is used for shared memory computers
c     if (mter.ge.2) go to 270
c     if (mter.eq.1) go to 180
c     if (jll.eq.0) then
c        if (nyzp(1,kll).eq.1).and.(nyzp(1,kl).eq.1)) then
c           do 130 k = 1, kzp1
c           do 120 j = 2, nx
c           f(j+1,2,k+1,m) = f(j+1,2,k+1,m) - scs(j+1,k+1,2,kll)
c 120       continue
c 130       continue
c        endif
c     endif
c     if (jr.eq.(nvpy-1)) then
c        if (nyzp(1,kr).eq.1) then
c           do 170 k = 1, kzp1
c           do 160 j = 2, nx
c           f(j+1,nyzp(1,m)+1,k+1,m) = f(j+1,nyzp(1,m)+1,k+1,m) - scs(j+
c    11,k+1,1,kr)
c 160       continue
c 170       continue
c        endif
c     endif
c 180 if (jl.eq.0) then
c        if (nyzp(1,kl).le.2) then
c           if (nyzp(1,kl).eq.1) then
c              if (nyzp(1,m).gt.1) then
c                 do 220 k = 1, kzp1
c                 do 210 j = 2, nx
c                 f(j+1,3,k+1,m) = f(j+1,3,k+1,m) - scs(j+1,k+1,2,kl)
c 210             continue
c 220             continue
c              endif
c              do 240 k = 1, kzp1
c              do 230 j = 2, nx
c              f(j+1,2,k+1,m) = f(j+1,2,k+1,m) + 2.*scs(j+1,k+1,2,kl)
c 230          continue
c 240          continue
c           else
c              do 260 k = 1, kzp1
c              do 250 j = 2, nx
c              f(j+1,2,k+1,m) = f(j+1,2,k+1,m) - scs(j+1,k+1,2,kl)
c 250          continue
c 260          continue
c           endif
c        endif
c     endif
c this segment is used for mpi computers
      if (mter.ge.2) go to 270
      if (mter.eq.1) go to 180
      if (jll.eq.0) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+4,lgrp,msid,ierr)
         call MPI_IRECV(scs(1,1,3,m),nxvz,mreal,kll-1,moff+1,lgrp,nsid,i
     1err)
      endif
      if (jl.eq.0) then
         call MPI_SEND(nyzp(1,m),1,mint,kr-1,moff+4,lgrp,ierr)
      endif
      if (jl.eq.(-1)) then
         if (nyzp(1,m).eq.1) then
            call MPI_SEND(scs(1,1,2,m),nxvzs,mreal,krr-1,moff+1,lgrp,ier
     1r)
         else
            nps = 0
            call MPI_SEND(scs(1,1,2,m),nps,mreal,krr-1,moff+1,lgrp,ierr)
         endif
      endif
      if (jll.eq.0) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            if (ngc.eq.1) then
               do 130 k = 1, kzp1
               do 120 j = 2, nx
               f(j+1,2,k+1,m) = f(j+1,2,k+1,m) - scs(j+1,k+1,3,m)
  120          continue
  130          continue
            endif
            do 150 k = 1, kzp1
            do 140 j = 2, nx
            f(j+1,1,k+1,m) = 0.
  140       continue
  150       continue
         endif
      endif
      if (jr.eq.(nvpy-1)) then
         call MPI_IRECV(scs(1,1,3,m),nxvz,mreal,kr-1,moff+2,lgrp,msid,ie
     1rr)
      endif
      if (jr.eq.nvpy) then
         if (nyzp(1,m).eq.1) then
            call MPI_SEND(scs(1,1,1,m),nxvzs,mreal,kl-1,moff+2,lgrp,ierr
     1)
         else
            nps = 0
            call MPI_SEND(scs(1,1,1,m),nps,mreal,kl-1,moff+2,lgrp,ierr)
         endif
      endif
      if (jr.eq.(nvpy-1)) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            do 170 k = 1, kzp1
            do 160 j = 2, nx
            f(j+1,nyzp(1,m)+1,k+1,m) = f(j+1,nyzp(1,m)+1,k+1,m) - scs(j+
     11,k+1,3,m)
            f(j+1,1,k+1,m) = 0.
  160       continue
  170       continue
         endif
      endif
  180 if (jl.eq.0) then
         call MPI_IRECV(ngc,1,mint,kl-1,moff+5,lgrp,msid,ierr)
         call MPI_IRECV(scs(1,1,3,m),nxvz,mreal,kl-1,moff+1,lgrp,nsid,ie
     1rr)
      endif
      if ((jl.eq.(-1)).and.(jr.lt.nvpy)) then
         call MPI_SEND(nyzp(1,m),1,mint,kr-1,moff+5,lgrp,ierr)
         if (nyzp(1,m).le.2) then
            call MPI_SEND(scs(1,1,2,m),nxvzs,mreal,kr-1,moff+1,lgrp,ierr
     1)
            do 200 k = 1, kzp1
            do 190 j = 2, nx
            f(j+1,2,k+1,m) = 0.
  190       continue
  200       continue
         else
            nps = 0
            call MPI_SEND(scs(1,1,2,m),nps,mreal,kr-1,moff+1,lgrp,ierr)
         endif
      endif
      if (jl.eq.0) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_WAIT(nsid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         if (nps.gt.0) then
            if (ngc.eq.1) then
               if (nyzp(1,m).gt.1) then
                  do 220 k = 1, kzp1
                  do 210 j = 2, nx
                  f(j+1,3,k+1,m) = f(j+1,3,k+1,m) - scs(j+1,k+1,3,m)
  210             continue
  220             continue
               endif
               do 240 k = 1, kzp1
               do 230 j = 2, nx
               f(j+1,2,k+1,m) = f(j+1,2,k+1,m) + 2.*scs(j+1,k+1,3,m)
               f(j+1,1,k+1,m) = 0.
  230          continue
  240          continue
            else
               do 260 k = 1, kzp1
               do 250 j = 2, nx
               f(j+1,2,k+1,m) = f(j+1,2,k+1,m) - scs(j+1,k+1,3,m)
               f(j+1,1,k+1,m) = 0.
  250          continue
  260          continue
            endif
         endif
      endif
  270 continue
  280 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMACGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzp
     1mx,mblok,nblok,kyp,kzp,ngds)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c no copying is done at the boundary edges.
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes one extra guard
c cell.
c scs/scr = scratch array for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c linear interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, kyp, kzp
      real f, scs, scr
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(3,nxv,nypmx,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer ky, kz, js, ks, moff, noff, jr, jl, kr, kl, mnblok
      integer nx1, nyp1, nzp1, nxvz, nxvy, m, my, mz, j, k, n
      dimension istatus(lstat)
      nx1 = nx + 1
      nyp1 = kyp + 1
      nzp1 = kzp + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c buffer data in y
      do 40 m = 1, mnblok
      do 30 k = 1, nzpmx
      do 20 j = 1, nxv
      do 10 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,kyp+1,k,m)
   10 continue
   20 continue
   30 continue
   40 continue
c add guard cells in y
      do 150 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 140 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jl = ky - 1
      kr = jr + kz
      kl = jl + kz
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 70 k = 1, nzpmx
c        do 60 j = 1, nxv
c        do 50 n = 1, 3
c        scs(n,j,k,2,m) = scs(n,j,k,1,kl)
c  50    continue
c  60    continue
c  70    continue
c     else
c        do 100 k = 1, nzp1
c        do 90 j = 1, nx1
c        do 80 n = 1, 3
c        scs(n,j,k,2,m) = 0.
c  80    continue
c  90    continue
c 100    continue
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,1,2,m),3*nxvz,mreal,kl-1,noff+1,lgrp,msi
     1d,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,1,m),3*nxvz,mreal,kr-1,noff+1,lgrp,ierr
     1)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 70 k = 1, nzp1
         do 60 j = 1, nx1
         do 50 n = 1, 3
         scs(n,j,k,2,m) = 0.
   50    continue
   60    continue
   70    continue
      endif
c add up the guard cells
      do 130 k = 1, nzp1
      do 120 j = 1, nx1
      do 110 n = 1, 3
      f(n,j,1,k,m) = f(n,j,1,k,m) + scs(n,j,k,2,m)
  110 continue
  120 continue
  130 continue
  140 continue
  150 continue
c add guard cells in z
      do 230 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 220 my = 1, mblok
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
c     do 180 k = 1, nypmx
c     do 170 j = 1, nxv
c     do 160 n = 1, 3
c     scr(n,j,k,m) = f(n,j,k,kzp+1,kl)
c 160 continue
c 170 continue
c 180 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,3*nxvy,mreal,kl-1,noff+2,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,1,kzp+1,m),3*nxvy,mreal,kr-1,noff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
      do 210 k = 1, nyp1
      do 200 j = 1, nx1
      do 190 n = 1, 3
      f(n,j,k,1,m) = f(n,j,k,1,m) + scr(n,j,k,m)
  190 continue
  200 continue
  210 continue
  220 continue
  230 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMAGUARD32L(f,scs,scr,kstrt,nvpy,nvpz,nx,nxv,nypmx,nzpm
     1x,mblok,nblok,kyp,kzp,ngds)
c this subroutine copies data from particle to field partitions, adding
c data from guard cells, where the field and particle partitions are 
c assumed to be the same.
c no copying is done at the boundary edges.
c f(j,k,l,m) = real data for grid j,k,l in particle partition m.  the
c number of grids per partition is uniform and includes one extra guard
c cell.
c scs/scr = scratch array for particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c linear interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, kyp, kzp
      real f, scs, scr
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(nxv,nypmx,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer ky, kz, js, ks, moff, noff, jr, jl, kr, kl, mnblok
      integer nx1, nyp1, nzp1, nxvz, nxvy, m, my, mz, j, k
      dimension istatus(lstat)
      nx1 = nx + 1
      nyp1 = kyp + 1
      nzp1 = kzp + 1
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
      do 110 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 100 my = 1, mblok
      m = my + moff
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jl = ky - 1
      kr = jr + kz
      kl = jl + kz
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 50 k = 1, nzpmx
c        do 40 j = 1, nxv
c        scs(j,k,2,m) = scs(j,k,1,kl)
c  40    continue
c  50    continue
c     else
c        do 70 k = 1, nzp1
c        do 60 j = 1, nx1
c        scs(j,k,2,m) = 0.
c  60    continue
c  70    continue
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,2,m),nxvz,mreal,kl-1,noff+1,lgrp,msid,ie
     1rr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,m),nxvz,mreal,kr-1,noff+1,lgrp,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 50 k = 1, nzp1
         do 40 j = 1, nx1
         scs(j,k,2,m) = 0.
   40    continue
   50    continue
      endif
c add up the guard cells
      do 90 k = 1, nzp1
      do 80 j = 1, nx1
      f(j,1,k,m) = f(j,1,k,m) + scs(j,k,2,m)
   80 continue
   90 continue
  100 continue
  110 continue
c add guard cells in z
      do 170 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 160 my = 1, mblok
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
c     do 130 k = 1, nypmx
c     do 120 j = 1, nxv
c     scr(j,k,m) = f(j,k,kzp+1,kl)
c 120 continue
c 130 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,nxvy,mreal,kl-1,noff+2,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,kzp+1,m),nxvy,mreal,kr-1,noff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
      do 150 k = 1, nyp1
      do 140 j = 1, nx1
      f(j,k,1,m) = f(j,k,1,m) + scr(j,k,m)
  140 continue
  150 continue
  160 continue
  170 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNMACGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nyp
     1mx,nzpmx,mblok,nblok,ngds,idds)
c this subroutine adds data from guard cells in non-uniform partitions
c for vector data.  no copying is done at the boundary edges in y.
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes one extra guard cell.
c scs/scr = scratch arrays for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c linear interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds
      integer nyzp
      real f, scs, scr
      dimension f(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(3,nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(3,nxv,nypmx,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer ky, kz, js, ks, moff, noff, jr, jl, kr, kl, mnblok
      integer nx1, nxvz, nxvzs, nyzp1, nxvy, nxvys, m, my, mz, j, k, n
      dimension istatus(lstat)
      nx1 = nx + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c special case for one processor in y
      if (nvpy.eq.1) go to 190
c buffer data in y
      do 40 m = 1, mnblok
      do 30 k = 1, nyzp(2,m)+1
      do 20 j = 1, nxv
      do 10 n = 1, 3
      scs(n,j,k,1,m) = f(n,j,nyzp(1,m)+1,k,m)
   10 continue
   20 continue
   30 continue
   40 continue
c add guard cells in y
      do 180 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 170 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(2,m) + 1
      nxvzs = nxv*nyzp1
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jl = ky - 1
      kr = jr + kz
      kl = jl + kz
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 70 k = 1, nyzp1
c        do 60 j = 1, nxv
c        do 50 n = 1, 3
c        scs(n,j,k,2,m) = scs(n,j,k,1,kl)
c  50    continue
c  60    continue
c  70    continue
c     else
c        do 100 k = 1, nyzp1
c        do 90 j = 1, nxv
c        do 80 n = 1, 3
c        scs(n,j,k,2,m) = 0.
c  80    continue
c  90    continue
c 100    continue
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,1,2,m),3*nxvz,mreal,kl-1,noff+1,lgrp,msi
     1d,ierr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,1,m),3*nxvzs,mreal,kr-1,noff+1,lgrp,ier
     1r)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 100 k = 1, nyzp1
         do 90 j = 1, nxv
         do 80 n = 1, 3
         scs(n,j,k,2,m) = 0.
   80    continue
   90    continue
  100    continue
      endif
c add up the guard cells
      do 130 k = 1, nyzp1
      do 120 j = 1, nx1
      do 110 n = 1, 3
      f(n,j,1,k,m) = f(n,j,1,k,m) + scs(n,j,k,2,m)
  110 continue
  120 continue
  130 continue
      if (jr.lt.nvpy) then
         do 160 k = 1, nyzp1
         do 150 j = 1, nx1
         do 140 n = 1, 3
         f(n,j,nyzp(1,m)+1,k,m) = 0.
  140    continue
  150    continue
  160    continue
      endif
  170 continue
  180 continue
c special case for one processor in z
  190 if (nvpz.eq.1) then
         do 230 m = 1, mnblok
         do 220 k = 1, nyzp(1,m)+1
         do 210 j = 1, nx1
         do 200 n = 1, 3
         f(n,j,k,1,m) = f(n,j,k,1,m) + f(n,j,k,nyzp(2,m)+1,m)
         f(n,j,k,nyzp(2,m)+1,m) = 0.
  200    continue
  210    continue
  220    continue
  230    continue
         return
      endif
c add guard cells in z
      do 310 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 300 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(1,m) + 1
      nxvys = nxv*nyzp1
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      kl = kz - 1
      if (kl.lt.0) kl = kl + nvpz
      kr = ky + nvpy*kr
      kl = ky + nvpy*kl
c this segment is used for shared memory computers
c     do 260 k = 1, nyzp1
c     do 250 j = 1, nxv
c     do 240 n = 1, 3
c     scr(n,j,k,m) = f(n,j,k,nyzp(2,m)+1,kl)
c 240 continue
c 250 continue
c 240 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,3*nxvy,mreal,kl-1,noff+2,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,1,nyzp(2,m)+1,m),3*nxvys,mreal,kr-1,noff+2,lgr
     1p,ierr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
      do 290 k = 1, nyzp1
      do 280 j = 1, nx1
      do 270 n = 1, 3
      f(n,j,k,1,m) = f(n,j,k,1,m) + scr(n,j,k,m)
      f(n,j,k,nyzp(2,m)+1,m) = 0.
  270 continue
  280 continue
  290 continue
  300 continue
  310 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNMAGUARD32L(f,scs,scr,nyzp,kstrt,nvpy,nvpz,nx,nxv,nypm
     1x,nzpmx,mblok,nblok,ngds,idds)
c this subroutine adds data from guard cells in non-uniform partitions
c for vector data.  no copying is done at the boundary edges in y.
c f(3,j,k,l,m) = real data for grid j,k,l in particle partition m.
c the grid is non-uniform and includes one extra guard cell.
c scs/scr = scratch arrays for particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nx = system length in x direction
c nxv = first dimension of f, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c mblok/nblok = number of particle partitions in y/z
c kyp/kzp = number of complex grids in y/z for each field partition.
c ngds = number of guard cells
c idds = dimensionality of domain decomposition
c linear interpolation, for distributed data
c with 2D spatial decomposition
      implicit none
      integer kstrt, nvpy, nvpz, nx, nxv, nypmx, nzpmx, mblok, nblok
      integer ngds, idds
      integer nyzp
      real f, scs, scr
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension scs(nxv,nzpmx,2*ngds,mblok*nblok)
      dimension scr(nxv,nypmx,mblok*nblok)
      dimension nyzp(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierr, msid, istatus
      integer ky, kz, js, ks, moff, noff, jr, jl, kr, kl, mnblok
      integer nx1, nxvz, nxvzs, nyzp1, nxvy, nxvys, m, my, mz, j, k
      dimension istatus(lstat)
      nx1 = nx + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      noff = nypmx*nzpmx
      mnblok = mblok*nblok
      nxvz = nxv*nzpmx
      nxvy = nxv*nypmx
c special case for one processor in y
      if (nvpy.eq.1) go to 140
c buffer data in y
      do 30 m = 1, mnblok
      do 20 k = 1, nyzp(2,m)+1
      do 10 j = 1, nxv
      scs(j,k,1,m) = f(j,nyzp(1,m)+1,k,m)
   10 continue
   20 continue
   30 continue
c add guard cells in y
      do 130 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 120 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(2,m) + 1
      nxvzs = nxv*nyzp1
      ky = my + js
      kz = nvpy*(mz + ks) + 1
      jr = ky + 1
      jl = ky - 1
      kr = jr + kz
      kl = jl + kz
c this segment is used for shared memory computers
c     if (jl.ge.0) then
c        do 50 k = 1, nyzp1
c        do 40 j = 1, nxv
c        scs(j,k,2,m) = scs(j,k,1,kl)
c  40    continue
c  50    continue
c     else
c        do 70 k = 1, nyzp1
c        do 60 j = 1, nxv
c        scs(j,k,2,m) = 0.
c  60    continue
c  70    continue
c     endif
c this segment is used for mpi computers
      if (jl.ge.0) then
         call MPI_IRECV(scs(1,1,2,m),nxvz,mreal,kl-1,noff+1,lgrp,msid,ie
     1rr)
      endif
      if (jr.lt.nvpy) then
         call MPI_SEND(scs(1,1,1,m),nxvzs,mreal,kr-1,noff+1,lgrp,ierr)
      endif
      if (jl.ge.0) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         do 70 k = 1, nyzp1
         do 60 j = 1, nxv
         scs(j,k,2,m) = 0.
   60    continue
   70    continue
      endif
c add up the guard cells
      do 90 k = 1, nyzp1
      do 80 j = 1, nx1
      f(j,1,k,m) = f(j,1,k,m) + scs(j,k,2,m)
   80 continue
   90 continue
      if (jr.lt.nvpy) then
         do 110 k = 1, nyzp1
         do 100 j = 1, nx1
         f(j,nyzp(1,m)+1,k,m) = 0.
  100    continue
  110    continue
      endif
  120 continue
  130 continue
c special case for one processor in z
  140 if (nvpz.eq.1) then
         do 170 m = 1, mnblok
         do 160 k = 1, nyzp(1,m)+1
         do 150 j = 1, nx1
         f(j,k,1,m) = f(j,k,1,m) + f(j,k,nyzp(2,m)+1,m)
         f(j,k,nyzp(2,m)+1,m) = 0.
  150    continue
  160    continue
  170    continue
         return
      endif
c add guard cells in z
      do 230 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 220 my = 1, mblok
      m = my + moff
      nyzp1 = nyzp(1,m) + 1
      nxvys = nxv*nyzp1
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      kl = kz - 1
      if (kl.lt.0) kl = kl + nvpz
      kr = ky + nvpy*kr
      kl = ky + nvpy*kl
c this segment is used for shared memory computers
c     do 190 k = 1, nyzp1
c     do 180 j = 1, nxv
c     scr(j,k,m) = f(j,k,nyzp(2,m)+1,kl)
c 180 continue
c 190 continue
c this segment is used for mpi computers
      call MPI_IRECV(scr,nxvy,mreal,kl-1,noff+2,lgrp,msid,ierr)
      call MPI_SEND(f(1,1,nyzp(2,m)+1,m),nxvys,mreal,kr-1,noff+2,lgrp,ie
     1rr)
      call MPI_WAIT(msid,istatus,ierr)
c add up the guard cells
      do 210 k = 1, nyzp1
      do 200 j = 1, nx1
      f(j,k,1,m) = f(j,k,1,m) + scr(j,k,m)
      f(j,k,nyzp(2,m)+1,m) = 0.
  200 continue
  210 continue
  220 continue
  230 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDBLSIN32C(cu,cu2,scb,scd,nx,ny,kstrt,nvpy,nxv,kyp,kzp,
     1kypd,kzpd,kyp2,kblok,lblok,k2blok)
c this subroutine creates a doubled vector array cu2 from a vector array
c cu, so that various 2d sine/cosine transforms can be performed with a
c 3d real to complex fft.  the x component is an odd function in y,
c and y component is an odd function in x.
c Asummes vector cu vanishes at end points
c linear interpolation for distributed data with 2D domain decomposition
c cu2 array may be modified
c scb/scd = scratch arrays
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nvpy = number of real or virtual processors in y
c nxv = second dimension of input array cu, must be >= nx
c kyp/kzp = number of data values per block in cu in y/z
c kypd = third dimension of input array cu, must be >= kyp
c kzpd = fourth dimension of input array cu, must be >= kzp
c kyp2 = number of data values per block in cu2 in y
c kblok/lblok = number of data blocks in y/z
c k2blok = number of data blocks in y for doubled data
      implicit none
      real cu, cu2, scb, scd
      integer nx, ny, kstrt, nvpy, nxv, kyp, kzp, kypd, kzpd, kyp2
      integer kblok, lblok, k2blok
      dimension cu(3,nxv,kypd,kzpd,kblok*lblok)
      dimension cu2(3,2*nxv,2*kypd,kzp,k2blok*lblok)
      dimension scb(3,nxv,kypd,kzpd,kblok*lblok)
      dimension scd(3,nxv,kzpd,2,kblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, lsid, msid, nsid, ierr
      integer i, j, k, l, m, my, mz, nxs, nys, ny2, kyb, kyb2, ks, js
      integer nxvy, noff, moff, koff, joff, kz, ll, lm, mm, ml, kk
      integer l1, l2, k0, k1, k2
      dimension istatus(lstat)
      nxs = nx - 1
      nys = ny - 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      kyb = ny/kyp
      ny2 = ny + ny
      kyb2 = ny2/kyp2
      nxvy = nxv*kypd
      noff = kypd*kzpd + kyb
c copy to double array in y direction
      do 140 mz = 1, lblok
      moff = k2blok*(mz - 1)
      do 130 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      ll = koff/kyp + 1
      koff = kyp*(my + js)
      lm = koff/kyp2 + 1
      kz = nvpy*(mz + ks)
      mm = ll + kz
      ml = lm + kz
c special case for one processor in y direction
      if (kyb2.eq.1) then
         do 40 l = 1, kzp
         do 20 k = 1, nys
         do 10 j = 1, nxs
         cu2(1,j+1,k+1,l,m) = cu(1,j+1,k+1,l,m)
         cu2(2,j+1,k+1,l,m) = cu(2,j+1,k+1,l,m)
         cu2(3,j+1,k+1,l,m) = cu(3,j+1,k+1,l,m)
         cu2(1,nx+j+1,k+1,l,m) = cu(1,nx-j+1,k+1,l,m)
         cu2(2,nx+j+1,k+1,l,m) = -cu(2,nx-j+1,k+1,l,m)
         cu2(3,nx+j+1,k+1,l,m) = -cu(3,nx-j+1,k+1,l,m)
         cu2(1,j+1,ny+k+1,l,m) = -cu(1,j+1,ny-k+1,l,m)
         cu2(2,j+1,ny+k+1,l,m) = cu(2,j+1,ny-k+1,l,m)
         cu2(3,j+1,ny+k+1,l,m) = -cu(3,j+1,ny-k+1,l,m)
         cu2(1,nx+j+1,ny+k+1,l,m) = -cu(1,nx-j+1,ny-k+1,l,m)
         cu2(2,nx+j+1,ny+k+1,l,m) = -cu(2,nx-j+1,ny-k+1,l,m)
         cu2(3,nx+j+1,ny+k+1,l,m) = cu(3,nx-j+1,ny-k+1,l,m)
   10    continue
         cu2(1,1,k+1,l,m) = 0.
         cu2(2,1,k+1,l,m) = 0.
         cu2(3,1,k+1,l,m) = 0.
         cu2(1,nx+1,k+1,l,m) = 0.
         cu2(2,nx+1,k+1,l,m) = 0.
         cu2(3,nx+1,k+1,l,m) = 0.
         cu2(1,1,k+ny+1,l,m) = 0.
         cu2(2,1,k+ny+1,l,m) = 0.
         cu2(3,1,k+ny+1,l,m) = 0.
         cu2(1,nx+1,k+ny+1,l,m) = 0.
         cu2(2,nx+1,k+ny+1,l,m) = 0.
         cu2(3,nx+1,k+ny+1,l,m) = 0.
   20    continue
         do 30 j = 1, nx
         cu2(1,j,1,l,m) = 0.
         cu2(2,j,1,l,m) = 0.
         cu2(3,j,1,l,m) = 0.
         cu2(1,j+nx,1,l,m) = 0.
         cu2(2,j+nx,1,l,m) = 0.
         cu2(3,j+nx,1,l,m) = 0.
         cu2(1,j,ny+1,l,m) = 0.
         cu2(2,j,ny+1,l,m) = 0.
         cu2(3,j,ny+1,l,m) = 0.
         cu2(1,j+nx,ny+1,l,m) = 0.
         cu2(2,j+nx,ny+1,l,m) = 0.
         cu2(3,j+nx,ny+1,l,m) = 0.
   30    continue
   40    continue
         return
      endif
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        do 80 l = 1, kzp
c        do 70 k = 1, kyp
c        do 60 j = 1, nx
c        do 50 i = 1, 3
c        cu2(i,j,k,l,m) = cu(i,j,k,l,mm)
c  50    continue
c  60    continue
c  70    continue
c  80   continue
c        if (kyp.lt.kyp2) then
c           do 120 l = 1, kzp
c           do 110 k = 1, kyp
c           do 100 j = 1, nx
c           do 90 i = 1, 3
c           cu2(i,j,k+kyp,l,m) = cu(i,j,k,l,mm+1)
c  90       continue
c 100       continue
c 110       continue
c 120       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         call MPI_IRECV(cu2(1,1,1,1,m),3*kzp*nxvy,mreal,mm-1,noff+1,lgrp
     1,msid,ierr)
         if (kyp.lt.kyp2) then
            call MPI_IRECV(scb(1,1,1,1,m),3*kzp*nxvy,mreal,mm,noff+1,lgr
     1p,nsid,ierr)
         endif
      endif
      if (lm.le.(kyb2/2)) then
         call MPI_SEND(cu(1,1,1,1,m),3*kzp*nxvy,mreal,ml-1,noff+1,lgrp,i
     1err)
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         call MPI_WAIT(msid,istatus,ierr)
         do 80 l = 1, kzp
         l1 = kzp - l + 1
         l2 = (l1 - 1)/4 + 1
         koff = kypd*(l1 - 4*(l2 - 1) - 1)
         do 70 k = 1, kypd
         k1 = kypd - k + 1
         k0 = k1 - 1 + koff
         k2 = k0/2 + 1
         joff = nxv*(k0 - 2*(k2 - 1))
         do 60 j = 1, nxv
         do 50 i = 1, 3
         cu2(i,j,k1,l1,m) = cu2(i,j+joff,k2,l2,m)
   50    continue
   60    continue
   70    continue
   80    continue
         if (kyp.lt.kyp2) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 120 l = 1, kzp
            do 110 k = 1, kypd
            do 100 j = 1, nxv
            do 90 i = 1, 3
            cu2(i,j,k+kyp,l,m) = scb(i,j,k,l,m)
   90       continue
  100       continue
  110       continue
  120       continue
         endif
      endif
  130 continue
  140 continue
c copy to double array in y direction
      do 360 mz = 1, lblok
      moff = k2blok*(mz - 1)
      do 350 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      ll = (ny2 - koff - 1)/kyp + 1
      koff = kyp*(my + js)
      lm = (ny2 - koff - 1)/kyp2 + 1
      koff = koff + kyp2*lm - ny2
      kz = nvpy*(mz + ks)
      mm = ll + kz
      ml = lm + kz
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        if ((ll+1).le.kyb) then
c           do 170 l = 1, kzp
c           do 160 j = 1, nx
c           do 150 i = 1, 3
c           cu2(i,j,1,l,m) = cu(i,j,1,l,mm+1)
c 150       continue
c 160       continue
c 170       continue
c        endif
c        if (kyp.lt.kyp2) then
c           do 210 l = 1, kzp
c           do 200 k = 1, kyp
c           do 190 j = 1, nx
c           do 180 i = 1, 3
c           cu2(i,j,k+kyp,l,m) = cu(i,j,k,l,mm)
c 180       continue
c 190       continue
c 200       continue
c 210       continue
c        endif
c        if (kyp.gt.1) then
c           do 250 l = 1, kzp
c           do 240 k = 2, kyp
c           do 230 j = 1, nx
c           do 220 i = 1, 3
c           cu2(i,j,k,l,m) = cu(i,j,k,l,mm-1)
c 220       continue
c 230       continue
c 240       continue
c 250       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         if ((ll+1).le.kyb) then
            call MPI_IRECV(scd(1,1,1,2,m),3*kzp*nxv,mreal,mm,noff+2,lgrp
     1,lsid,ierr)
         endif
         if (kyp.lt.kyp2) then
            call MPI_IRECV(scb(1,1,1,1,m),3*kzp*nxvy,mreal,mm-1,noff+2,l
     1grp,msid,ierr)
         endif
         if (kyp.gt.1) then
            call MPI_IRECV(cu2(1,1,1,1,m),3*kzp*nxvy,mreal,mm-2,noff+2,l
     1grp,nsid,ierr)
         endif
      endif
      if ((lm.gt.(kyb2/2)).and.(lm.le.kyb2)) then
         if (koff.eq.0) then
            if ((lm+1).le.kyb2) then
               do 170 l = 1, kzp
               do 160 j = 1, nxv
               do 150 i = 1, 3
               scd(i,j,l,1,m) = cu(i,j,1,l,m)
  150          continue
  160          continue
  170          continue
               call MPI_SEND(scd(1,1,1,1,m),3*kzp*nxv,mreal,ml,noff+2,lg
     1rp,ierr)
            endif
            if (kyp.gt.1) then
               call MPI_SEND(cu(1,1,1,1,m),3*kzp*nxvy,mreal,ml-1,noff+2,
     1lgrp,ierr)
            endif
         else
            call MPI_SEND(cu(1,1,1,1,m),3*kzp*nxvy,mreal,ml-1,noff+2,lgr
     1p,ierr)
         endif
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         if (kyp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 210 l = 1, kzp
            l1 = kzp - l + 1
            l2 = (l1 - 1)/4 + 1
            koff = kypd*(l1 - 4*(l2 - 1) - 1)
            do 200 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1 + koff
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 190 j = 1, nxv
            do 180 i = 1, 3
            cu2(i,j,k1,l1,m) = cu2(i,j+joff,k2,l2,m)
  180       continue
  190       continue
  200       continue
  210       continue
         endif
         if ((ll+1).le.kyb) then
            call MPI_WAIT(lsid,istatus,ierr)
            do 240 l = 1, kzp
            do 230 j = 1, nxv
            do 220 i = 1, 3
            cu2(i,j,1,l,m) = scd(i,j,l,2,m)
  220       continue
  230       continue
  240       continue
         endif
         if (kyp.lt.kyp2) then
            call MPI_WAIT(msid,istatus,ierr)
            do 340 l = 1, kzp
            do 330 k = 1, kypd
            do 320 j = 1, nxv
            do 310 i = 1, 3
            cu2(i,j,k+kyp,l,m) = scb(i,j,k,l,m)
  310       continue
  320       continue
  330       continue
  340       continue
         endif
      endif
  350 continue
  360 continue
c create odd array
      do 470 mz = 1, lblok
      moff = k2blok*(mz - 1)
      do 460 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      do 450 l = 1, kzp
      do 440 k = 1, kyp2
      kk = k + koff
      if ((kk.eq.1).or.(kk.eq.(ny+1))) then
         do 380 j = 1, nx
         do 370 i = 1, 3
         cu2(i,j,k,l,m) = 0.
         cu2(i,j+nx,k,l,m) = 0.
  370    continue
  380    continue
      else if (kk.le.ny) then
         do 390 j = 1, nxs
         cu2(1,nx+j+1,k,l,m) = cu2(1,nx-j+1,k,l,m)
         cu2(2,nx+j+1,k,l,m) = -cu2(2,nx-j+1,k,l,m)
         cu2(3,nx+j+1,k,l,m) = -cu2(3,nx-j+1,k,l,m)
  390    continue
         do 400 i = 1, 3
         cu2(i,1,k,l,m) = 0.
         cu2(i,nx+1,k,l,m) = 0.
  400    continue
      else if (kk.gt.(ny+1)) then
         if (k.eq.1) then
            do 410 j = 1, nxs
            cu2(1,nx+j+1,k,l,m) = -cu2(1,nx-j+1,k,l,m)
            cu2(2,nx+j+1,k,l,m) = -cu2(2,nx-j+1,k,l,m)
            cu2(3,nx+j+1,k,l,m) = cu2(3,nx-j+1,k,l,m)
  410       continue
         else
            do 420 j = 1, nxs
            cu2(1,nx+j+1,kyp2-k+2,l,m) = -cu2(1,nx-j+1,k,l,m)
            cu2(2,nx+j+1,kyp2-k+2,l,m) = -cu2(2,nx-j+1,k,l,m)
            cu2(3,nx+j+1,kyp2-k+2,l,m) = cu2(3,nx-j+1,k,l,m)
  420       continue
         endif
         do 430 i = 1, 3
         cu2(i,1,k,l,m) = 0.
         cu2(i,nx+1,k,l,m) = 0.
  430    continue
      endif
  440 continue
  450 continue
  460 continue
  470 continue
c finish odd array
      do 520 mz = 1, lblok
      moff = k2blok*(mz - 1)
      do 510 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      do 500 l = 1, kzp
      do 490 k = 1, kyp2
      kk = k + koff
      if (kk.gt.(ny+1)) then
         do 480 j = 1, nxs
         cu2(1,nx-j+1,k,l,m) = cu2(1,nx+j+1,k,l,m)
         cu2(2,nx-j+1,k,l,m) = -cu2(2,nx+j+1,k,l,m)
         cu2(3,nx-j+1,k,l,m) = -cu2(3,nx+j+1,k,l,m)
  480    continue
         cu2(1,nx+1,k,l,m) = cu2(1,nx+1,k,l,m)
         cu2(2,nx+1,k,l,m) = -cu2(2,nx+1,k,l,m)
         cu2(3,nx+1,k,l,m) = -cu2(3,nx+1,k,l,m)
      endif
  490 continue
  500 continue
  510 continue
  520 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDBLSIN32D(q,q2,scb,scd,nx,ny,kstrt,nvpy,nxv,kyp,kzp,ky
     1pd,kzpd,kyp2,kblok,lblok,k2blok)
c this subroutine creates an odd array q2 from an array q, so that
c a 2d sine transform can be performed with a 3d real to complex fft.
c linear interpolation for distributed data with 2D domain decomposition
c Asummes q vanishes at end points
c q2 array may be modified
c scb/scd = scratch arrays
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nvpy = number of real or virtual processors in y
c nxv = first dimension of input array q, must be >= nx
c kyp/kzp = number of data values per block in cu in y/z
c kypd = second dimension of input array q, must be >= kyp
c kzpd = third dimension of input array q, must be >= kzp
c kyp2 = number of data values per block in q2 in y
c kblok/lblok = number of data blocks in y/z
c k2blok = number of data blocks in y for doubled data
      implicit none
      real q, q2, scb, scd
      integer nx, ny, kstrt, nvpy, nxv, kyp, kzp, kypd, kzpd, kyp2
      integer kblok, lblok, k2blok
      dimension q(nxv,kypd,kzpd,kblok*lblok)
      dimension q2(2*nxv,2*kypd,kzp,k2blok*lblok)
      dimension scb(nxv,kypd,kzpd,kblok*lblok)
      dimension scd(nxv,kzpd,2,kblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, lsid, msid, nsid, ierr
      integer j, k, l, m, my, mz, nxs, nys, ny2, kyb, kyb2, ks, js
      integer nxvy, noff, moff, koff, joff, kz, ll, lm, mm, ml, kk
      integer l1, l2, k0, k1, k2
      dimension istatus(lstat)
      nxs = nx - 1
      nys = ny - 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      kyb = ny/kyp
      ny2 = ny + ny
      kyb2 = ny2/kyp2
      nxvy = nxv*kypd
      noff = kypd*kzpd + kyb
c copy to double array in y direction
      do 120 mz = 1, lblok
      moff = k2blok*(mz - 1)
      do 110 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      ll = koff/kyp + 1
      koff = kyp*(my + js)
      lm = koff/kyp2 + 1
      kz = nvpy*(mz + ks)
      mm = ll + kz
      ml = lm + kz
c special case for one processor in y direction
      if (kyb2.eq.1) then
         do 40 l = 1, kzp
         do 20 k = 1, nys
         do 10 j = 1, nxs
         q2(j+1,k+1,l,m) = q(j+1,k+1,l,m)
         q2(nx+j+1,k+1,l,m) = -q(nx-j+1,k+1,l,m)
         q2(j+1,ny+k+1,l,m) = -q(j+1,ny-k+1,l,m)
         q2(nx+j+1,ny+k+1,l,m) = q(nx-j+1,ny-k+1,l,m)
   10    continue
         q2(1,k+1,l,m) = 0.
         q2(nx+1,k+1,l,m) = 0.
         q2(1,k+ny+1,l,m) = 0.
         q2(nx+1,k+ny+1,l,m) = 0.
   20    continue
         do 30 j = 1, nx
         q2(j,1,l,m) = 0.
         q2(j+nx,1,l,m) = 0.
         q2(j,ny+1,l,m) = 0.
         q2(j+nx,ny+1,l,m) = 0.
   30    continue
   40    continue
         return
      endif
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        do 70 l = 1, kzp
c        do 60 k = 1, kyp
c        do 50 j = 1, nx
c        q2(j,k,l,m) = q(j,k,l,mm)
c  50    continue
c  60    continue
c  70   continue
c        if (kyp.lt.kyp2) then
c           do 100 l = 1, kzp
c           do 90 k = 1, kyp
c           do 80 j = 1, nx
c           q2(j,k+kyp,l,m) = q(j,k,l,mm+1)
c  80       continue
c  90       continue
c 100       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         call MPI_IRECV(q2(1,1,1,m),kzp*nxvy,mreal,mm-1,noff+1,lgrp,msid
     1,ierr)
         if (kyp.lt.kyp2) then
            call MPI_IRECV(scb(1,1,1,m),kzp*nxvy,mreal,mm,noff+1,lgrp,ns
     1id,ierr)
         endif
      endif
      if (lm.le.(kyb2/2)) then
         call MPI_SEND(q(1,1,1,m),kzp*nxvy,mreal,ml-1,noff+1,lgrp,ierr)
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         call MPI_WAIT(msid,istatus,ierr)
         do 70 l = 1, kzp
         l1 = kzp - l + 1
         l2 = (l1 - 1)/4 + 1
         koff = kypd*(l1 - 4*(l2 - 1) - 1)
         do 60 k = 1, kypd
         k1 = kypd - k + 1
         k0 = k1 - 1 + koff
         k2 = k0/2 + 1
         joff = nxv*(k0 - 2*(k2 - 1))
         do 50 j = 1, nxv
         q2(j,k1,l1,m) = q2(j+joff,k2,l2,m)
   50    continue
   60    continue
   70    continue
         if (kyp.lt.kyp2) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 100 l = 1, kzp
            do 90 k = 1, kypd
            do 80 j = 1, nxv
            q2(j,k+kyp,l,m) = scb(j,k,l,m)
   80       continue
   90       continue
  100       continue
         endif
      endif
  110 continue
  120 continue
c copy to double array in y direction
      do 240 mz = 1, lblok
      moff = k2blok*(mz - 1)
      do 230 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      ll = (ny2 - koff - 1)/kyp + 1
      koff = kyp*(my + js)
      lm = (ny2 - koff - 1)/kyp2 + 1
      koff = koff + kyp2*lm - ny2
      kz = nvpy*(mz + ks)
      mm = ll + kz
      ml = lm + kz
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        if ((ll+1).le.kyb) then
c           do 140 l = 1, kzp
c           do 130 j = 1, nx
c           q2(j,1,l,m) = q(j,1,l,mm+1)
c 130       continue
c 140       continue
c        endif
c        if (kyp.lt.kyp2) then
c           do 170 l = 1, kzp
c           do 160 k = 1, kyp
c           do 150 j = 1, nx
c           q2(j,k+kyp,l,m) = q(j,k,l,mm)
c 150       continue
c 160       continue
c 170       continue
c        endif
c        if (kyp.gt.1) then
c           do 200 l = 1, kzp
c           do 190 k = 2, kyp
c           do 180 j = 1, nx
c           q2(j,k,l,m) = q(j,k,l,mm-1)
c 180       continue
c 190       continue
c 200       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         if ((ll+1).le.kyb) then
            call MPI_IRECV(scd(1,1,2,m),kzp*nxv,mreal,mm,noff+2,lgrp,lsi
     1d,ierr)
         endif
         if (kyp.lt.kyp2) then
            call MPI_IRECV(scb(1,1,1,m),kzp*nxvy,mreal,mm-1,noff+2,lgrp,
     1msid,ierr)
         endif
         if (kyp.gt.1) then
            call MPI_IRECV(q2(1,1,1,m),kzp*nxvy,mreal,mm-2,noff+2,lgrp,n
     1sid,ierr)
         endif
      endif
      if ((lm.gt.(kyb2/2)).and.(lm.le.kyb2)) then
         if (koff.eq.0) then
            if ((lm+1).le.kyb2) then
               do 140 l = 1, kzp
               do 130 j = 1, nxv
               scd(j,l,1,m) = q(j,1,l,m)
  130          continue
  140          continue
               call MPI_SEND(scd(1,1,1,m),kzp*nxv,mreal,ml,noff+2,lgrp,i
     1err)
            endif
            if (kyp.gt.1) then
               call MPI_SEND(q(1,1,1,m),kzp*nxvy,mreal,ml-1,noff+2,lgrp,
     1ierr)
            endif
         else
            call MPI_SEND(q(1,1,1,m),kzp*nxvy,mreal,ml-1,noff+2,lgrp,ier
     1r)
         endif
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         if (kyp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 170 l = 1, kzp
            l1 = kzp - l + 1
            l2 = (l1 - 1)/4 + 1
            koff = kypd*(l1 - 4*(l2 - 1) - 1)
            do 160 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1 + koff
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 150 j = 1, nxv
            q2(j,k1,l1,m) = q2(j+joff,k2,l2,m)
  150       continue
  160       continue
  170       continue
         endif
         if ((ll+1).le.kyb) then
            call MPI_WAIT(lsid,istatus,ierr)
            do 190 l = 1, kzp
            do 180 j = 1, nxv
            q2(j,1,l,m) = scd(j,l,2,m)
  180       continue
  190       continue
         endif
         if (kyp.lt.kyp2) then
            call MPI_WAIT(msid,istatus,ierr)
            do 220 l = 1, kzp
            do 210 k = 1, kypd
            do 200 j = 1, nxv
            q2(j,k+kyp,l,m) = scb(j,k,l,m)
  200       continue
  210       continue
  220       continue
         endif
      endif
  230 continue
  240 continue
c create odd array
      do 320 mz = 1, lblok
      moff = k2blok*(mz - 1)
      do 310 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      do 300 l = 1, kzp
      do 290 k = 1, kyp2
      kk = k + koff
      if ((kk.eq.1).or.(kk.eq.(ny+1))) then
         do 250 j = 1, nx
         q2(j,k,l,m) = 0.
         q2(j+nx,k,l,m) = 0.
  250    continue
      else if (kk.le.ny) then
         do 260 j = 1, nxs
         q2(nx+j+1,k,l,m) = -q2(nx-j+1,k,l,m)
  260    continue
         q2(1,k,l,m) = 0.
         q2(nx+1,k,l,m) = 0.
      else if (kk.gt.(ny+1)) then
         if (k.eq.1) then
            do 270 j = 1, nxs
            q2(nx+j+1,k,l,m) = q2(nx-j+1,k,l,m)
  270       continue
         else
            do 280 j = 1, nxs
            q2(nx+j+1,kyp2-k+2,l,m) = q2(nx-j+1,k,l,m)
  280       continue
         endif
         q2(1,k,l,m) = 0.
         q2(nx+1,k,l,m) = 0.
      endif
  290 continue
  300 continue
  310 continue
  320 continue
c finish odd array
      do 370 mz = 1, lblok
      moff = k2blok*(mz - 1)
      do 360 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      do 350 l = 1, kzp
      do 340 k = 1, kyp2
      kk = k + koff
      if (kk.gt.(ny+1)) then
         do 330 j = 1, nxs
         q2(nx-j+1,k,l,m) = -q2(nx+j+1,k,l,m)
  330    continue
         q2(nx+1,k,l,m) = -q2(nx+1,k,l,m)
      endif
  340 continue
  350 continue
  360 continue
  370 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PHAFDBL32C(fxyz,fxyz2,scb,scd,nx,ny,kstrt,nvpy,nvpz,nxv
     1,kyp,kzp,kypd,kzpd,kyp2,kblok,lblok,k2blok)
c this subroutine copies data from a double array to regular array
c with guard cells for vector field and linear interpolation
c for distributed data with 2D domain decomposition
c fxyz array may be modified
c scb/scd = scratch arrays
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = second dimension of input array fxyz, must be >= nx
c kyp/kzp = number of data values per block in fxyz in y/z
c kypd = third dimension of input array fxyz, must be >= kyp
c kzpd = fourth dimension of input array fxyz, must be >= kzp
c kyp2 = number of data values per block in fxyz2 in y
c kblok/lblok = number of data blocks in y/z
c k2blok = number of data blocks in y for doubled data
      implicit none
      real fxyz, fxyz2, scb, scd
      integer nx, ny, kstrt, nvpy, nvpz, nxv, kyp, kzp, kypd, kzpd, kyp2
      integer kblok, lblok, k2blok
      dimension fxyz(3,nxv,kypd,kzpd,kblok*lblok)
      dimension fxyz2(3,2*nxv,2*kypd,kzp,k2blok*lblok)
      dimension scb(3,nxv,kypd,kzpd,kblok*lblok)
      dimension scd(3,nxv,kzpd,2,kblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer i, j, k, l, m, my, mz, nx1, ny1, kyb, kyb2, ks, js, kr, kl
      integer nxvy, noff, moff, koff, ky, kz, ll, lm, mm, ml
      dimension istatus(lstat)
      nx1 = nx + 1
      ny1 = ny + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      kyb = ny/kyp
      kyb2 = (ny + ny)/kyp2
      nxvy = nxv*kypd
      noff = kypd*kzpd + kyb
c copy from double array in x direction
      do 270 mz = 1, lblok
      moff = k2blok*(mz - 1)
      do 260 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      lm = koff/kyp + 1
      koff = kyp*(my + js)
      ll = koff/kyp2 + 1
      koff = koff - kyp2*(ll - 1)
      kz = nvpy*(mz + ks)
      mm = ll + kz
      ml = lm + kz
c special case for one processor in y direction
      if (kyb2.eq.1) then
         do 40 l = 1, kzp
         do 30 k = 1, ny1
         do 20 j = 1, nx1
         do 10 i = 1, 3
         fxyz(i,j,k,l,m) = fxyz2(i,j,k,l,m)
   10    continue
   20    continue
   30    continue
   40    continue
         go to 260
      endif
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        if ((koff.eq.0).and.(kyp.lt.kyp2)) then
c           do 80 l = 1, kzp
c           do 70 k = 1, kyp1
c           do 60 j = 1, nx1
c           do 50 i = 1, 3
c           fxyz(i,j,k,l,m) = fxyz2(i,j,k,l,mm)
c  50       continue
c  60       continue
c  70       continue
c  80       continue
c        else
c           do 120 l = 1, kzp
c           do 110 k = 1, kyp
c           do 100 j = 1, nx1
c           do 90 i = 1, 3
c           fxyz(i,j,k,l,m) = fxyz2(i,j,k+kyp,l,mm)
c  90       continue
c 100       continue
c 110       continue
c 120       continue
c           do 150 l = 1, kzp
c           do 140 j = 1, nx1
c           do 130 i = 1, 3
c           fxyz(i,j,kyp+1,l,m) = fxyz2(i,j,1,l,mm+1)
c 130       continue
c 140       continue
c 150       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_IRECV(fxyz(1,1,1,1,m),3*kzp*nxvy,mreal,mm-1,noff+1,
     1lgrp,msid,ierr)
         else
            call MPI_IRECV(fxyz(1,1,1,1,m),3*kzp*nxvy,mreal,mm-1,noff+1,
     1lgrp,msid,ierr)
            call MPI_IRECV(scd(1,1,1,2,m),3*kzp*nxv,mreal,mm,noff+1,lgrp
     1,nsid,ierr)
         endif
      endif
c pack data and send it
      if (lm.le.kyb) then
         if (kyp.lt.kyp2) then
            do 80 l = 1, kzp
            do 70 k = 1, kypd
            do 60 j = 1, nxv
            do 50 i = 1, 3
            scb(i,j,k,l,m) = fxyz2(i,j,k,l,m)
   50       continue
   60       continue
   70       continue
   80       continue
            call MPI_SEND(scb(1,1,1,1,m),3*kzp*nxvy,mreal,ml-1,noff+1,lg
     1rp,ierr)
            do 120 l = 1, kzp
            do 110 k = 1, kypd
            do 100 j = 1, nxv
            do 90 i = 1, 3
            scb(i,j,k,l,m) = fxyz2(i,j,k+kyp,l,m)
   90       continue
  100       continue
  110       continue
  120       continue
            call MPI_SEND(scb(1,1,1,1,m),3*kzp*nxvy,mreal,ml,noff+1,lgrp
     1,ierr)
         else
            do 160 l = 1, kzp
            do 150 k = 1, kypd
            do 140 j = 1, nxv
            do 130 i = 1, 3
            scb(i,j,k,l,m) = fxyz2(i,j,k,l,m)
  130       continue
  140       continue
  150       continue
  160       continue
            call MPI_SEND(scb(1,1,1,1,m),3*kzp*nxvy,mreal,ml-1,noff+1,lg
     1rp,ierr)
         endif
         if (lm.gt.1) then
            do 190 l = 1, kzp
            do 180 j = 1, nxv
            do 170 i = 1, 3
            scd(i,j,l,1,m) = fxyz2(i,j,1,l,m)
  170       continue
  180       continue
  190       continue
            call MPI_SEND(scd(1,1,1,1,m),3*kzp*nxv,mreal,ml-2,noff+1,lgr
     1p,ierr)
         endif
      else if (lm.eq.(kyb+1)) then
         do 220 l = 1, kzp
         do 210 j = 1, nxv
         do 200 i = 1, 3
         scd(i,j,l,1,m) = fxyz2(i,j,1,l,m)
  200    continue
  210    continue
  220    continue
         call MPI_SEND(scd(1,1,1,1,m),3*kzp*nxv,mreal,ml-2,noff+1,lgrp,i
     1err)
      endif
c wait for data
      if (ll.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_WAIT(msid,istatus,ierr)
         else
            call MPI_WAIT(msid,istatus,ierr)
            call MPI_WAIT(nsid,istatus,ierr)
            do 250 l = 1, kzp
            do 240 j = 1, nxv
            do 230 i = 1, 3
            fxyz(i,j,kyp+1,l,m) = scd(i,j,l,2,m)
  230       continue
  240       continue
  250       continue
         endif
      endif
  260 continue
  270 continue
c copy to guard cells in z
      do 320 mz = 1, lblok
      moff = kblok*(mz - 1)
      do 310 my = 1, kblok
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
c     do 300 k = 1, nypmx
c     do 290 j = 1, nxv
c     do 280 i = 1, 3
c     fxyz(i,j,k,kzp+1,m) = fxyz(i,j,k,1,kr)
c 280 continue
c 290 continue
c 300 continue
c this segment is used for mpi computers
      call MPI_IRECV(fxyz(1,1,1,kzp+1,m),3*nxvy,mreal,kr-1,noff+2,lgrp,m
     1sid,ierr)
      call MPI_SEND(fxyz(1,1,1,1,m),3*nxvy,mreal,kl-1,noff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
  310 continue
  320 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PHAFDBL32D(q,q2,scb,scd,nx,ny,kstrt,nvpy,nvpz,nxv,kyp,k
     1zp,kypd,kzpd,kyp2,kblok,lblok,k2blok)
c this subroutine copies data from a double array to regular array
c with guard cells for scalar field and linear interpolation
c for distributed data with 2D domain decomposition
c q array may be modified
c scb/scd = scratch arrays
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = second dimension of input array q, must be >= nx
c kyp/kzp = number of data values per block in q in y/z
c kypd = third dimension of input array q, must be >= kyp
c kzpd = fourth dimension of input array q, must be >= kzp
c kyp2 = number of data values per block in q in y
c kblok/lblok = number of data blocks in y/z
c k2blok = number of data blocks in y for doubled data
      implicit none
      real q, q2, scb, scd
      integer nx, ny, kstrt, nvpy, nvpz, nxv, kyp, kzp, kypd, kzpd, kyp2
      integer kblok, lblok, k2blok
      dimension q(nxv,kypd,kzpd,kblok*lblok)
      dimension q2(2*nxv,2*kypd,kzp,k2blok*lblok)
      dimension scb(nxv,kypd,kzpd,kblok*lblok)
      dimension scd(nxv,kzpd,2,kblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, msid, nsid, ierr
      integer j, k, l, m, my, mz, nx1, ny1, kyb, kyb2, ks, js, kr, kl
      integer nxvy, noff, moff, koff, ky, kz, ll, lm, mm, ml
      dimension istatus(lstat)
      nx1 = nx + 1
      ny1 = ny + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      kyb = ny/kyp
      kyb2 = (ny + ny)/kyp2
      nxvy = nxv*kypd
      noff = kypd*kzpd + kyb
c copy from double array in x direction
      do 200 mz = 1, lblok
      moff = k2blok*(mz - 1)
      do 190 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      lm = koff/kyp + 1
      koff = kyp*(my + js)
      ll = koff/kyp2 + 1
      koff = koff - kyp2*(ll - 1)
      kz = nvpy*(mz + ks)
      mm = ll + kz
      ml = lm + kz
c special case for one processor in y direction
      if (kyb2.eq.1) then
         do 30 l = 1, kzp
         do 20 k = 1, ny1
         do 10 j = 1, nx1
         q(j,k,l,m) = q2(j,k,l,m)
   10    continue
   20    continue
   30    continue
         go to 190
      endif
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        if ((koff.eq.0).and.(kyp.lt.kyp2)) then
c           do 60 l = 1, kzp
c           do 50 k = 1, kyp1
c           do 40 j = 1, nx1
c           q(j,k,l,m) = q2(j,k,l,mm)
c  40       continue
c  50       continue
c  60       continue
c        else
c           do 90 l = 1, kzp
c           do 80 k = 1, kyp
c           do 70 j = 1, nx1
c           q(j,k,l,m) = q2(j,k+kyp,l,mm)
c  70       continue
c  80       continue
c  90       continue
c           do 110 l = 1, kzp
c           do 100 j = 1, nx1
c           q(j,kyp+1,l,m) = q2(j,1,l,mm+1)
c 100       continue
c 110       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_IRECV(q(1,1,1,m),kzp*nxvy,mreal,mm-1,noff+1,lgrp,ms
     1id,ierr)
         else
            call MPI_IRECV(q(1,1,1,m),kzp*nxvy,mreal,mm-1,noff+1,lgrp,ms
     1id,ierr)
            call MPI_IRECV(scd(1,1,2,m),kzp*nxv,mreal,mm,noff+1,lgrp,nsi
     1d,ierr)
         endif
      endif
c pack data and send it
      if (lm.le.kyb) then
         if (kyp.lt.kyp2) then
            do 60 l = 1, kzp
            do 50 k = 1, kypd
            do 40 j = 1, nxv
            scb(j,k,l,m) = q2(j,k,l,m)
   40       continue
   50       continue
   60       continue
            call MPI_SEND(scb(1,1,1,m),kzp*nxvy,mreal,ml-1,noff+1,lgrp,i
     1err)
            do 90 l = 1, kzp
            do 80 k = 1, kypd
            do 70 j = 1, nxv
            scb(j,k,l,m) = q2(j,k+kyp,l,m)
   70       continue
   80       continue
   90       continue
            call MPI_SEND(scb(1,1,1,m),kzp*nxvy,mreal,ml,noff+1,lgrp,ier
     1r)
         else
            do 120 l = 1, kzp
            do 110 k = 1, kypd
            do 100 j = 1, nxv
            scb(j,k,l,m) = q2(j,k,l,m)
  100       continue
  110       continue
  120       continue
            call MPI_SEND(scb(1,1,1,m),kzp*nxvy,mreal,ml-1,noff+1,lgrp,i
     1err)
         endif
         if (lm.gt.1) then
            do 140 l = 1, kzp
            do 130 j = 1, nxv
            scd(j,l,1,m) = q2(j,1,l,m)
  130       continue
  140       continue
            call MPI_SEND(scd(1,1,1,m),kzp*nxv,mreal,ml-2,noff+1,lgrp,ie
     1rr)
         endif
      else if (lm.eq.(kyb+1)) then
         do 160 l = 1, kzp
         do 150 j = 1, nxv
         scd(j,l,1,m) = q2(j,1,l,m)
  150    continue
  160    continue
         call MPI_SEND(scd(1,1,1,m),kzp*nxv,mreal,ml-2,noff+1,lgrp,ierr)
      endif
c wait for data
      if (ll.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_WAIT(msid,istatus,ierr)
         else
            call MPI_WAIT(msid,istatus,ierr)
            call MPI_WAIT(nsid,istatus,ierr)
            do 180 l = 1, kzp
            do 170 j = 1, nxv
            q(j,kyp+1,l,m) = scd(j,l,2,m)
  170       continue
  180       continue
         endif
      endif
  190 continue
  200 continue
c copy to guard cells in z
      do 240 mz = 1, lblok
      moff = kblok*(mz - 1)
      do 230 my = 1, kblok
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
c     do 220 k = 1, nypmx
c     do 210 j = 1, nxv
c     q(j,k,kzp+1,m) = q(j,k,1,kr)
c 210 continue
c 220 continue
c this segment is used for mpi computers
      call MPI_IRECV(q(1,1,kzp+1,m),nxvy,mreal,kr-1,noff+2,lgrp,msid,ier
     1r)
      call MPI_SEND(q(1,1,1,m),nxvy,mreal,kl-1,noff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
  230 continue
  240 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PZTRP32D(q,q3,scb,nx,ny,nz,kstrt,nvpy,nxv,kyp,kzp,kypd,
     1kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
c this subroutine creates an zeroed array q3 from an array q, so that
c a 3d real to complex transform will perform a correct convolution.
c linear interpolation for distributed data with 2D domain decomposition
c q3 array may be modified
c scb = scratch array
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c nvpy = number of real or virtual processors in y
c nxv = first dimension of input array q, must be >= nx
c kyp/kzp = number of data values per block in q in y/z
c kypd = second dimension of input array q, must be >= kyp
c kzpd = third dimension of input array q, must be >= kzp
c kyp2/kzp2 = number of data values per block in q3 in y/z
c kblok/lblok = number of data blocks in y/z
c k2blok/l2blok = number of data blocks in y/z for tripled data
      implicit none
      real q, q3, scb
      integer nx, ny, nz, kstrt, nvpy, nxv, kyp, kzp, kypd, kzpd
      integer kyp2, kzp2, kblok, lblok, k2blok, l2blok
      dimension q(nxv,kypd,kzpd,kblok*lblok)
      dimension q3(2*nxv,2*kypd,kzp2,k2blok*l2blok)
      dimension scb(nxv,kypd,kzpd,kblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      logical lt1
      integer istatus, msid, nsid, ierr
      integer i, j, k, l, m, my, mz, ks, js, ny2, nz2, kyb, kyb2
      integer kzb, kzb2, nxvy, noff, moff, loff, koff, joff
      integer jm, km, jl, kl, kz, mm, lm, l1, l2, k0, k1, k2
      dimension istatus(lstat)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      kyb = ny/kyp
      ny2 = ny + ny
      kyb2 = ny2/kyp2
      kzb = nz/kzp
      nz2 = nz + nz
      kzb2 = nz2/kzp2
      nxvy = nxv*kypd
      noff = kypd*kzpd + kyb*kzb
c copy to triple array in y direction
      do 120 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 110 my = 1, k2blok
      m = my + moff
      loff = kyp2*(my + js)
      jm = loff/kyp + 1
      loff = kzp2*(mz + ks)
      km = loff/kzp + 1
      loff = kyp*(my + js)
      jl = loff/kyp2 + 1
      loff = kzp*(mz + ks)
      kl = loff - kzp2*(loff/kzp2)
      kz = nvpy*(mz + ks)
c special case for one processor
      if ((kyb2*kzb2).eq.1) then
         do 30 l = 1, nz
         do 20 k = 1, ny
         do 10 j = 1, nx
         q3(j,k,l,m) = q(j,k,l,m)
         q3(nx+j,k,l,m) = 0.
         q3(j,ny+k,l,m) = 0.
         q3(nx+j,ny+k,l,m) = 0.
         q3(j,k,nz+l,m) = 0.
         q3(nx+j,k,nz+l,m) = 0.
         q3(j,ny+k,nz+l,m) = 0.
         q3(nx+j,ny+k,nz+l,m) = 0.
   10    continue
   20    continue
   30    continue
         return
      endif
c copy main data in y direction
      do 100 i = 1, 2
c for odd rows in z
      if (i.eq.1) then
         mm = jm + 2*kz
         lm = jl + kz/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in z
      else if (i.eq.2) then
         if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 100
         mm = mm + nvpy
         lm = lm - nvpy/2
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kyb).and.(km.le.kzb)) then
c        do 60 l = 1, kzp
c        do 50 k = 1, kyp
c        do 40 j = 1, nx
c        q3(j,k,l+loff,m) = q(j,k,l,mm)
c  40    continue
c  50    continue
c  60    continue
c        if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
c           do 90 l = 1, kzp
c           do 80 k = 1, kyp
c           do 70 j = 1, nx
c           q3(j,k+kyp,l+loff,m) = q(j,k,l,mm+1)
c  70       continue
c  80       continue
c  90       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kyb).and.(km.le.kzb)) then
         call MPI_IRECV(q3(1,1,loff+1,m),kzp*nxvy,mreal,mm-1,noff+i,lgrp
     1,msid,ierr)
         if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
            call MPI_IRECV(scb(1,1,1,m),kzp*nxvy,mreal,mm,noff+i,lgrp,ns
     1id,ierr)
         endif
      endif
      if ((jl.le.((kyb2-1)/2+1)).and.lt1) then
         call MPI_SEND(q(1,1,1,m),kzp*nxvy,mreal,lm-1,noff+i,lgrp,ierr)
      endif
c wait for data and unpack it
      if ((jm.le.kyb).and.(km.le.kzb)) then
         call MPI_WAIT(msid,istatus,ierr)
         do 60 l = 1, kzp
         l1 = kzp - l + 1
         l2 = (l1 - 1)/4 + 1
         koff = kypd*(l1 - 4*(l2 - 1) - 1)
         do 50 k = 1, kypd
         k1 = kypd - k + 1
         k0 = k1 - 1 + koff
         k2 = k0/2 + 1
         joff = nxv*(k0 - 2*(k2 - 1))
         do 40 j = 1, nxv
         q3(j,k1,l1+loff,m) = q3(j+joff,k2,l2+loff,m)
   40    continue
   50    continue
   60    continue
         if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 90 l = 1, kzp
            do 80 k = 1, kypd
            do 70 j = 1, nxv
            q3(j,k+kyp,l+loff,m) = scb(j,k,l,m)
   70       continue
   80       continue
   90       continue
         endif
      endif
  100 continue
  110 continue
  120 continue
c zero out remainder of array
      do 200 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      loff = kzp2*(mz + ks)
      do 190 my = 1, k2blok
      koff = kyp2*(my + js)
      m = my + moff
      do 180 l = 1, kzp2
      l1 = l + loff
      if (l1.le.nz) then
         do 150 k = 1, kyp2
         k1 = k + koff
         if (k1.le.ny) then
            do 130 j = 1, nx
            q3(j+nx,k,l,m) = 0.
  130       continue
         else
            do 140 j = 1, nx
            q3(j,k,l,m) = 0.
            q3(j+nx,k,l,m) = 0.
  140       continue
         endif
  150    continue
      else
         do 170 k = 1, kyp2
         do 160 j = 1, nx
         q3(j,k,l,m) = 0.
         q3(j+nx,k,l,m) = 0.
  160    continue
  170    continue
      endif
  180 continue
  190 continue
  200 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMOVE32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,js
     1r,jsl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbmax
     2,idds,ntmax,info)
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
c info = status information
c info(1) = ierr = (0,N) = (no,yes) error condition exists
c info(2) = maximum number of particles per processor
c info(3) = minimum number of particles per processor
c info(4) = maximum number of buffer overflows in y
c info(5) = maximum number of buffer overflows in z
c info(6) = maximum number of particle passes required in y
c info(7) = maximum number of particle passes required in z
c info(8) = total number of particles on entry
c info(9) = difference of total number of particles on exit
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl
      integer npp, ihole, jsr, jsl, jss, info
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
      dimension info(9)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer iy, iz
      parameter(iy=2,iz=3)
      integer ierr, ic, js, ks, mnblok, i, n, m, my, mz, moff, nvp, iter
      integer npr, nps, npt, kb, kl, kr, j, j1, j2, nbsize, nter, mter
      integer itermax
      integer msid, istatus
      integer ibflg, iwork
      double precision bflg, work
      real an, xt
      dimension msid(4), istatus(lstat)
      dimension ibflg(4), iwork(4)
      dimension bflg(2), work(2)
      dimension kb(2)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      mnblok = mblok*nblok
      nbsize = idimp*nbmax
      do 5 j = 1, 9
      info(j) = 0
    5 continue
      itermax = 2000
c debugging section: count total number of particles before move
      npr = 0
      do 10 m = 1, mnblok
      npr = npr + npp(m)
   10 continue
c buffer outgoing particles, first in y then in z direction
      do 300 n = 1, 2
      if (n.eq.1) then
         ic = iy
         nvp = nvpy
         an = float(ny)
      elseif (n.eq.2) then
         ic = iz
         nvp = nvpz
         an = float(nz)
      endif
      iter = 2
      nter = 0
   20 mter = 0
      do 60 mz = 1, nblok
      moff = mblok*(mz - 1)
      kb(2) = mz + ks
      do 50 my = 1, mblok
      m = my + moff
      kb(1) = my + js
      jsl(1,m) = 0
      jsr(1,m) = 0
      jss(2,m) = 0
      do 30 j = 1, npp(m)
      xt = part(ic,j,m)
c particles going down or backward
      if (xt.lt.edges(2*n-1,m)) then
         if (jsl(1,m).lt.nbmax) then
            jsl(1,m) = jsl(1,m) + 1
            if (kb(n).eq.0) xt = xt + an
            do 23 i = 1, idimp
            sbufl(i,jsl(1,m),m) = part(i,j,m)
   23       continue
            sbufl(ic,jsl(1,m),m) = xt
            ihole(jsl(1,m)+jsr(1,m),m) = j
         else
            jss(2,m) = 1
            go to 40
         endif
c particles going up or forward
      else if (xt.ge.edges(2*n,m)) then
         if (jsr(1,m).lt.nbmax) then
            jsr(1,m) = jsr(1,m) + 1
            if ((kb(n)+1).eq.nvp) xt = xt - an
            do 27 i = 1, idimp
            sbufr(i,jsr(1,m),m) = part(i,j,m)
   27       continue
            sbufr(ic,jsr(1,m),m) = xt
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
      nps = max0(nps,jss(2,m))
  100 continue
      ibflg(3) = nps
c copy particle buffers
  110 iter = iter + 2
      mter = mter + 1
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
c     do 115 i = 1, idimp
c     rbufl(i,j,m) = sbufr(i,j,kl)
c 115 continue
c 120 continue
c     jsr(2,m) = jsl(1,kr)
c     do 130 j = 1, jsr(2,m)
c     do 125 i = 1, idimp
c     rbufr(i,j,m) = sbufl(i,j,kr)
c 125 continue
c 130 continue
c this segment is used for mpi computers
c post receive
      call MPI_IRECV(rbufl,nbsize,mreal,kl-1,iter-1,lgrp,msid(1),ierr)
      call MPI_IRECV(rbufr,nbsize,mreal,kr-1,iter,lgrp,msid(2),ierr)
c send particles
      call MPI_ISEND(sbufr,idimp*jsr(1,m),mreal,kr-1,iter-1,lgrp,msid(3)
     1,ierr)
      call MPI_ISEND(sbufl,idimp*jsl(1,m),mreal,kl-1,iter,lgrp,msid(4),i
     1err)
c wait for particles to arrive
      call MPI_WAIT(msid(1),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsl(2,m) = nps/idimp
      call MPI_WAIT(msid(2),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsr(2,m) = nps/idimp
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
      if (rbufr(ic,j,m).lt.edges(2*n-1,m)) jsl(1,m) = jsl(1,m) + 1
      if (rbufr(ic,j,m).ge.edges(2*n,m)) jsr(1,m) = jsr(1,m) + 1
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
      if (rbufl(ic,j,m).ge.edges(2*n,m)) jsr(1,m) = jsr(1,m) + 1
      if (rbufl(ic,j,m).lt.edges(2*n-1,m)) jss(2,m) = jss(2,m) + 1
  170 continue
      if (jss(2,m).ne.0) then
         if (n.eq.1) then
            write (2,*) 'Info:',jss(2,m),' particles returning below'
         elseif (n.eq.2) then
            write (2,*) 'Info:',jss(2,m),' particles returning back'
         endif
      endif
      jsl(1,m) = jsl(1,m) + jss(2,m)
      nps = max0(nps,jsl(1,m)+jsr(1,m))
  180 continue
      ibflg(2) = nps
c make sure sbufr and sbufl have been sent
      call MPI_WAIT(msid(3),istatus,ierr)
      call MPI_WAIT(msid(4),istatus,ierr)
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
      xt = rbufr(ic,j,m)
c particles going down or back
      if (xt.lt.edges(2*n-1,m)) then
         jsl(1,m) = jsl(1,m) + 1
         if (kb(n).eq.0) xt = xt + an
         rbufr(ic,j,m) = xt
         do 183 i = 1, idimp
         sbufl(i,jsl(1,m),m) = rbufr(i,j,m)
  183    continue
c particles going up or front, should not happen
      elseif (xt.ge.edges(2*n,m)) then
         jsr(1,m) = jsr(1,m) + 1
         if ((kb(n)+1).eq.nvp) xt = xt - an
         rbufr(ic,j,m) = xt
         do 185 i = 1, idimp
         sbufr(i,jsr(1,m),m) = rbufr(i,j,m)
  185    continue
c particles staying here
      else
         jss(2,m) = jss(2,m) + 1
         do 187 i = 1, idimp
         rbufr(i,jss(2,m),m) = rbufr(i,j,m)
  187    continue
      endif
  190 continue
      jsr(2,m) = jss(2,m)
c next check particles coming from below or back
      jss(2,m) = 0
      do 200 j = 1, jsl(2,m)
      xt = rbufl(ic,j,m)
c particles going up or front
      if (xt.ge.edges(2*n,m)) then
         if (jsr(1,m).lt.nbmax) then
            jsr(1,m) = jsr(1,m) + 1
            if ((kb(n)+1).eq.nvp) xt = xt - an
            rbufl(ic,j,m) = xt
            do 193 i = 1, idimp
            sbufr(i,jsr(1,m),m) = rbufl(i,j,m)
  193       continue 
         else
            jss(2,m) = 2*npmax
            go to 210
         endif
c particles going down back, should not happen
      elseif (xt.lt.edges(2*n-1,m)) then
         if (jsl(1,m).lt.nbmax) then
            jsl(1,m) = jsl(1,m) + 1
            if (kb(n).eq.0) xt = xt + an
            rbufl(ic,j,m) = xt
            do 195 i = 1, idimp
            sbufl(i,jsl(1,m),m) = rbufl(i,j,m)
  195       continue
         else
            jss(2,m) = 2*npmax
            go to 210
         endif
c particles staying here
      else
         jss(2,m) = jss(2,m) + 1
         do 197 i = 1, idimp
         rbufl(i,jss(2,m),m) = rbufl(i,j,m)
  197    continue
      endif
  200 continue
  210 jsl(2,m) = jss(2,m)
  220 continue
  230 continue
c check if move would overflow particle array
  240 nps = 0
      npt = npmax
      do 250 m = 1, mnblok
      jss(2,m) = npp(m) + jsl(2,m) + jsr(2,m) - jss(1,m)
      nps = max0(nps,jss(2,m))
      npt = min0(npt,jss(2,m))
  250 continue
      ibflg(1) = nps
      ibflg(4) = -npt
      call PIMAX(ibflg,iwork,4,1)
      info(2) = ibflg(1)
      info(3) = -ibflg(4)
      ierr = ibflg(1) - npmax
      if (ierr.gt.0) then
         write (2,*) 'particle overflow error, ierr = ', ierr
         info(1) = ierr
         return
      endif
c distribute incoming particles from buffers
      do 290 m = 1, mnblok
c distribute particles coming from below or back into holes
      jss(2,m) = min0(jss(1,m),jsl(2,m))
      do 260 j = 1, jss(2,m)
      do 255 i = 1, idimp
      part(i,ihole(j,m),m) = rbufl(i,j,m)
  255 continue
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
         do 263 i = 1, idimp
         part(i,ihole(j+jsl(2,m),m),m) = rbufr(i,j,m)
  263    continue
      else
c no more holes
c distribute remaining particles from below or back into bottom
         do 267 i = 1, idimp
         part(i,j+npp(m),m) = rbufl(i,j+jss(1,m),m)
  267    continue
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
            do 273 i = 1, idimp
            part(i,ihole(j2,m),m) = part(i,j1,m)
  273       continue
         endif
      else
c no more holes
c distribute remaining particles from above or front into bottom
         do 277 i = 1, idimp
         part(i,j+npp(m),m) = rbufr(i,j+jss(1,m),m)
  277    continue
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
      info(5+n) = max0(info(5+n),mter)
      if (ibflg(2).gt.0) then
         write (2,*) 'Info: particles being passed further = ', ibflg(2)
         ! flush(2)
         if (ibflg(3).gt.0) ibflg(3) = 1
         if (iter.lt.itermax) go to 110
         ierr = -((iter-2)/2)
         write (2,*) 'Iteration overflow, iter = ', ierr
         info(1) = ierr
         go to 320
      endif
c check if buffer overflowed and more particles remain to be checked
      if (ibflg(3).gt.0) then
         nter = nter + 1
         info(3+n) = nter
         go to 20
      endif
  300 continue
c debugging section: count total number of particles after move
      nps = 0
      do 310 m = 1, mnblok
      nps = nps + npp(m)
  310 continue
      bflg(2) = dble(nps)
      bflg(1) = dble(npr)
      call PDSUM(bflg,work,2,1)
      info(8) = bflg(1)
      info(9) = bflg(2) - bflg(1)
      if (bflg(1).ne.bflg(2)) then
         write (2,*) 'particle number error, old/new=',bflg(1),bflg(2)
         info(1) = 1
      endif
c information
  320 if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine DMOVE32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,js
     1r,jsl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbmax
     2,idds,ntmax,info,tinfo)
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
c info = status information
c info(1) = ierr = (0,N) = (no,yes) error condition exists
c info(2) = maximum number of particles per processor
c info(3) = minimum number of particles per processor
c info(4) = maximum number of buffer overflows in y
c info(5) = maximum number of buffer overflows in z
c info(6) = maximum number of particle passes required in y
c info(7) = maximum number of particle passes required in z
c info(8) = total number of particles on entry
c info(9) = difference of total number of particles on exit
c tinfo = timing information
c tinfo(1-2) = check particles in y/z
c tinfo(3-4) = send/receive data in y/z
c tinfo(5-6) = check if particles belong in y/z
c tinfo(7-8) = remove particles which do not belong in y/z
c tinfo(9-10) = synchronize processors
c tinfo(11-2) = distribute incoming particles
c tinfo(13) = check if particles are lost
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl
      integer npp, ihole, jsr, jsl, jss, info
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
      dimension info(9)
      real tinfo
      dimension tinfo(13)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer iy, iz
      parameter(iy=2,iz=3)
      integer ierr, ic, js, ks, mnblok, i, n, m, my, mz, moff, nvp, iter
      integer npr, nps, npt, kb, kl, kr, j, j1, j2, nbsize, nter, mter
      integer itermax
      integer msid, istatus
      integer ibflg, iwork
      double precision bflg, work
      real an, xt
      double precision dtime
      double precision MPI_WTIME
      external MPI_WTIME
      dimension msid(4), istatus(lstat)
      dimension ibflg(4), iwork(4)
      dimension bflg(2), work(2)
      dimension kb(2)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      mnblok = mblok*nblok
      nbsize = idimp*nbmax
      do 5 j = 1, 9
      info(j) = 0
    5 continue
      do 6 j = 1, 13
      tinfo(j) = 0.0
    6 continue
      itermax = 2000
c debugging section: count total number of particles before move
      npr = 0
      do 10 m = 1, mnblok
      npr = npr + npp(m)
   10 continue
c buffer outgoing particles, first in y then in z direction
      do 300 n = 1, 2
      if (n.eq.1) then
         ic = iy
         nvp = nvpy
         an = float(ny)
      elseif (n.eq.2) then
         ic = iz
         nvp = nvpz
         an = float(nz)
      endif
      iter = 2
      nter = 0
   20 mter = 0
      dtime = MPI_WTIME()
      do 60 mz = 1, nblok
      moff = mblok*(mz - 1)
      kb(2) = mz + ks
      do 50 my = 1, mblok
      m = my + moff
      kb(1) = my + js
      jsl(1,m) = 0
      jsr(1,m) = 0
      jss(2,m) = 0
      do 30 j = 1, npp(m)
      xt = part(ic,j,m)
c particles going down or backward
      if (xt.lt.edges(2*n-1,m)) then
         if (jsl(1,m).lt.nbmax) then
            jsl(1,m) = jsl(1,m) + 1
            if (kb(n).eq.0) xt = xt + an
            do 23 i = 1, idimp
            sbufl(i,jsl(1,m),m) = part(i,j,m)
   23       continue
            sbufl(ic,jsl(1,m),m) = xt
            ihole(jsl(1,m)+jsr(1,m),m) = j
         else
            jss(2,m) = 1
            go to 40
         endif
c particles going up or forward
      else if (xt.ge.edges(2*n,m)) then
         if (jsr(1,m).lt.nbmax) then
            jsr(1,m) = jsr(1,m) + 1
            if ((kb(n)+1).eq.nvp) xt = xt - an
            do 27 i = 1, idimp
            sbufr(i,jsr(1,m),m) = part(i,j,m)
   27       continue
            sbufr(ic,jsr(1,m),m) = xt
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
      tinfo(n) = tinfo(n) + real(MPI_WTIME() - dtime)
c check for full buffer condition
      nps = 0
      do 100 m = 1, mnblok
      nps = max0(nps,jss(2,m))
  100 continue
      ibflg(3) = nps
c copy particle buffers
  110 iter = iter + 2
      mter = mter + 1
      dtime = MPI_WTIME()
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
c     do 115 i = 1, idimp
c     rbufl(i,j,m) = sbufr(i,j,kl)
c 115 continue
c 120 continue
c     jsr(2,m) = jsl(1,kr)
c     do 130 j = 1, jsr(2,m)
c     do 125 i = 1, idimp
c     rbufr(i,j,m) = sbufl(i,j,kr)
c 125 continue
c 130 continue
c this segment is used for mpi computers
c post receive
      call MPI_IRECV(rbufl,nbsize,mreal,kl-1,iter-1,lgrp,msid(1),ierr)
      call MPI_IRECV(rbufr,nbsize,mreal,kr-1,iter,lgrp,msid(2),ierr)
c send particles
      call MPI_ISEND(sbufr,idimp*jsr(1,m),mreal,kr-1,iter-1,lgrp,msid(3)
     1,ierr)
      call MPI_ISEND(sbufl,idimp*jsl(1,m),mreal,kl-1,iter,lgrp,msid(4),i
     1err)
c wait for particles to arrive
      call MPI_WAIT(msid(1),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsl(2,m) = nps/idimp
      call MPI_WAIT(msid(2),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsr(2,m) = nps/idimp
  140 continue
      tinfo(2+n) = tinfo(2+n) + real(MPI_WTIME() - dtime)
  150 continue
c check if particles must be passed further
      nps = 0
      dtime = MPI_WTIME()
      do 180 m = 1, mnblok
c check if any particles coming from above or front belong here
      jsl(1,m) = 0
      jsr(1,m) = 0
      jss(2,m) = 0
      do 160 j = 1, jsr(2,m)
      if (rbufr(ic,j,m).lt.edges(2*n-1,m)) jsl(1,m) = jsl(1,m) + 1
      if (rbufr(ic,j,m).ge.edges(2*n,m)) jsr(1,m) = jsr(1,m) + 1
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
      if (rbufl(ic,j,m).ge.edges(2*n,m)) jsr(1,m) = jsr(1,m) + 1
      if (rbufl(ic,j,m).lt.edges(2*n-1,m)) jss(2,m) = jss(2,m) + 1
  170 continue
      if (jss(2,m).ne.0) then
         if (n.eq.1) then
            write (2,*) 'Info:',jss(2,m),' particles returning below'
         elseif (n.eq.2) then
            write (2,*) 'Info:',jss(2,m),' particles returning back'
         endif
      endif
      jsl(1,m) = jsl(1,m) + jss(2,m)
      nps = max0(nps,jsl(1,m)+jsr(1,m))
  180 continue
      tinfo(4+n) = tinfo(4+n) + real(MPI_WTIME() - dtime)
      ibflg(2) = nps
c make sure sbufr and sbufl have been sent
      call MPI_WAIT(msid(3),istatus,ierr)
      call MPI_WAIT(msid(4),istatus,ierr)
      if (nps.eq.0) go to 240
      dtime = MPI_WTIME()
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
      xt = rbufr(ic,j,m)
c particles going down or back
      if (xt.lt.edges(2*n-1,m)) then
         jsl(1,m) = jsl(1,m) + 1
         if (kb(n).eq.0) xt = xt + an
         rbufr(ic,j,m) = xt
         do 183 i = 1, idimp
         sbufl(i,jsl(1,m),m) = rbufr(i,j,m)
  183    continue
c particles going up or front, should not happen
      elseif (xt.ge.edges(2*n,m)) then
         jsr(1,m) = jsr(1,m) + 1
         if ((kb(n)+1).eq.nvp) xt = xt - an
         rbufr(ic,j,m) = xt
         do 185 i = 1, idimp
         sbufr(i,jsr(1,m),m) = rbufr(i,j,m)
  185    continue
c particles staying here
      else
         jss(2,m) = jss(2,m) + 1
         do 187 i = 1, idimp
         rbufr(i,jss(2,m),m) = rbufr(i,j,m)
  187    continue
      endif
  190 continue
      jsr(2,m) = jss(2,m)
c next check particles coming from below or back
      jss(2,m) = 0
      do 200 j = 1, jsl(2,m)
      xt = rbufl(ic,j,m)
c particles going up or front
      if (xt.ge.edges(2*n,m)) then
         if (jsr(1,m).lt.nbmax) then
            jsr(1,m) = jsr(1,m) + 1
            if ((kb(n)+1).eq.nvp) xt = xt - an
            rbufl(ic,j,m) = xt
            do 193 i = 1, idimp
            sbufr(i,jsr(1,m),m) = rbufl(i,j,m)
  193       continue 
         else
            jss(2,m) = 2*npmax
            go to 210
         endif
c particles going down back, should not happen
      elseif (xt.lt.edges(2*n-1,m)) then
         if (jsl(1,m).lt.nbmax) then
            jsl(1,m) = jsl(1,m) + 1
            if (kb(n).eq.0) xt = xt + an
            rbufl(ic,j,m) = xt
            do 195 i = 1, idimp
            sbufl(i,jsl(1,m),m) = rbufl(i,j,m)
  195       continue
         else
            jss(2,m) = 2*npmax
            go to 210
         endif
c particles staying here
      else
         jss(2,m) = jss(2,m) + 1
         do 197 i = 1, idimp
         rbufl(i,jss(2,m),m) = rbufl(i,j,m)
  197    continue
      endif
  200 continue
  210 jsl(2,m) = jss(2,m)
  220 continue
  230 continue
      tinfo(6+n) = tinfo(6+n) + real(MPI_WTIME() - dtime)
c check if move would overflow particle array
  240 nps = 0
      npt = npmax
      dtime = MPI_WTIME()
      do 250 m = 1, mnblok
      jss(2,m) = npp(m) + jsl(2,m) + jsr(2,m) - jss(1,m)
      nps = max0(nps,jss(2,m))
      npt = min0(npt,jss(2,m))
  250 continue
      ibflg(1) = nps
      ibflg(4) = -npt
      call PIMAX(ibflg,iwork,4,1)
      info(2) = ibflg(1)
      info(3) = -ibflg(4)
      ierr = ibflg(1) - npmax
      if (ierr.gt.0) then
         write (2,*) 'particle overflow error, ierr = ', ierr
         info(1) = ierr
         return
      endif
      tinfo(8+n) = tinfo(8+n) + real(MPI_WTIME() - dtime)
      dtime = MPI_WTIME()
c distribute incoming particles from buffers
      do 290 m = 1, mnblok
c distribute particles coming from below or back into holes
      jss(2,m) = min0(jss(1,m),jsl(2,m))
      do 260 j = 1, jss(2,m)
      do 255 i = 1, idimp
      part(i,ihole(j,m),m) = rbufl(i,j,m)
  255 continue
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
         do 263 i = 1, idimp
         part(i,ihole(j+jsl(2,m),m),m) = rbufr(i,j,m)
  263    continue
      else
c no more holes
c distribute remaining particles from below or back into bottom
         do 267 i = 1, idimp
         part(i,j+npp(m),m) = rbufl(i,j+jss(1,m),m)
  267    continue
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
            do 273 i = 1, idimp
            part(i,ihole(j2,m),m) = part(i,j1,m)
  273       continue
         endif
      else
c no more holes
c distribute remaining particles from above or front into bottom
         do 277 i = 1, idimp
         part(i,j+npp(m),m) = rbufr(i,j+jss(1,m),m)
  277    continue
      endif
  280 continue
      if (jss(2,m).gt.0) then
         npp(m) = npp(m) - jsr(2,m)
      else
         npp(m) = npp(m) + jsr(2,m)
      endif
      jss(1,m) = 0
  290 continue
      tinfo(10+n) = tinfo(10+n) + real(MPI_WTIME() - dtime)
c check if any particles have to be passed further
      info(5+n) = max0(info(5+n),mter)
      if (ibflg(2).gt.0) then
         write (2,*) 'Info: particles being passed further = ', ibflg(2)
         ! flush(2)
         if (ibflg(3).gt.0) ibflg(3) = 1
         if (iter.lt.itermax) go to 110
         ierr = -((iter-2)/2)
         write (2,*) 'Iteration overflow, iter = ', ierr
         info(1) = ierr
         go to 320
      endif
c check if buffer overflowed and more particles remain to be checked
      if (ibflg(3).gt.0) then
         nter = nter + 1
         info(3+n) = nter
         go to 20
      endif
  300 continue
c debugging section: count total number of particles after move
      nps = 0
      do 310 m = 1, mnblok
      nps = nps + npp(m)
  310 continue
      dtime = MPI_WTIME()
      bflg(2) = dble(nps)
      bflg(1) = dble(npr)
      call PDSUM(bflg,work,2,1)
      info(8) = bflg(1)
      info(9) = bflg(2) - bflg(1)
      if (bflg(1).ne.bflg(2)) then
         write (2,*) 'particle number error, old/new=',bflg(1),bflg(2)
         info(1) = 1
      endif
      tinfo(13) = tinfo(13) + real(MPI_WTIME() - dtime)
c information
  320 if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PXMOV32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,js
     1r,jsl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbmax
     2,idds,ntmax,maskp,info)
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
c info = status information
c info(1) = ierr = (0,N) = (no,yes) error condition exists
c info(2) = maximum number of particles per processor
c info(3) = minimum number of particles per processor
c info(4) = maximum number of buffer overflows in y
c info(5) = maximum number of buffer overflows in z
c info(6) = maximum number of particle passes required in y
c info(7) = maximum number of particle passes required in z
c info(8) = total number of particles on entry
c info(9) = difference of total number of particles on exit
c optimized for vector processor
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl
      integer npp, ihole, jsr, jsl, jss, maskp, info
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
      dimension info(9)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer iy, iz
      parameter(iy=2,iz=3)
      integer ierr, ic, js, ks, mnblok, i, n, m, my, mz, moff, nvp, iter
      integer npr, nps, npt, kb, kl, kr, j, j1, j2, nbsize, nter, mter
      integer itermax
      integer msid, istatus
      integer ibflg, iwork
      double precision bflg, work
      real an, xt
      dimension msid(4), istatus(lstat)
      dimension ibflg(4), iwork(4)
      dimension bflg(2), work(2)
      dimension kb(2)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      mnblok = mblok*nblok
      nbsize = idimp*nbmax
      do 5 j = 1, 9
      info(j) = 0
    5 continue
      itermax = 2000
c debugging section: count total number of particles before move
      npr = 0
      do 10 m = 1, mnblok
      npr = npr + npp(m)
   10 continue
c buffer outgoing particles, first in y then in z direction
      do 300 n = 1, 2
      if (n.eq.1) then
         ic = iy
         nvp = nvpy
         an = float(ny)
      elseif (n.eq.2) then
         ic = iz
         nvp = nvpz
         an = float(nz)
      endif
      iter = 2
      nter = 0
   20 mter = 0
      do 90 mz = 1, nblok
      moff = mblok*(mz - 1)
      kb(2) = mz + ks
      do 80 my = 1, mblok
      m = my + moff
      kb(1) = my + js
      jss(1,m) = 0
      jss(2,m) = 0
c find mask function for particles out of bounds
      do 30 j = 1, npp(m)
      xt = part(ic,j,m)
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
      xt = part(ic,ihole(j,m),m)
c particles going backward
      if (xt.lt.edges(2*n-1,m)) then
         if (kb(n).eq.0) xt = xt + an
         if (jsl(1,m).lt.nbmax) then
            jsl(1,m) = jsl(1,m) + 1
            part(ic,ihole(j,m),m) = xt
            do 53 i = 1, idimp
            sbufl(i,jsl(1,m),m) = part(i,ihole(j,m),m)
   53       continue
            ihole(jsl(1,m)+jsr(1,m),m) = ihole(j,m)
         else
            jss(2,m) = 1
c           go to 70
         endif
c particles going forward
      else
         if ((kb(n)+1).eq.nvp) xt = xt - an
         if (jsr(1,m).lt.nbmax) then
            jsr(1,m) = jsr(1,m) + 1
            do 57 i = 1, idimp
            sbufr(i,jsr(1,m),m) = part(i,ihole(j,m),m)
   57       continue
            part(ic,ihole(j,m),m) = xt
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
      nps = max0(nps,jss(2,m))
  100 continue
      ibflg(3) = nps
c copy particle buffers
  110 iter = iter + 2
      mter = mter + 1
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
c     do 115 i = 1, idimp
c     rbufl(i,j,m) = sbufr(i,j,kl)
c 115 continue
c 120 continue
c     jsr(2,m) = jsl(1,kr)
c     do 130 j = 1, jsr(2,m)
c     do 125 i = 1, idimp
c     rbufr(i,j,m) = sbufl(i,j,kr)
c 125 continue
c 130 continue
c this segment is used for mpi computers
c post receive
      call MPI_IRECV(rbufl,nbsize,mreal,kl-1,iter-1,lgrp,msid(1),ierr)
      call MPI_IRECV(rbufr,nbsize,mreal,kr-1,iter,lgrp,msid(2),ierr)
c send particles
      call MPI_ISEND(sbufr,idimp*jsr(1,m),mreal,kr-1,iter-1,lgrp,msid(3)
     1,ierr)
      call MPI_ISEND(sbufl,idimp*jsl(1,m),mreal,kl-1,iter,lgrp,msid(4),i
     1err)
c wait for particles to arrive
      call MPI_WAIT(msid(1),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsl(2,m) = nps/idimp
      call MPI_WAIT(msid(2),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsr(2,m) = nps/idimp
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
      if (rbufr(ic,j,m).lt.edges(2*n-1,m)) jsl(1,m) = jsl(1,m) + 1
      if (rbufr(ic,j,m).ge.edges(2*n,m)) jsr(1,m) = jsr(1,m) + 1
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
      if (rbufl(ic,j,m).ge.edges(2*n,m)) jsr(1,m) = jsr(1,m) + 1
      if (rbufl(ic,j,m).lt.edges(2*n-1,m)) jss(2,m) = jss(2,m) + 1
  170 continue
      if (jss(2,m).ne.0) then
         if (n.eq.1) then
            write (2,*) 'Info:',jss(2,m),' particles returning below'
         elseif (n.eq.2) then
            write (2,*) 'Info:',jss(2,m),' particles returning back'
         endif
      endif
      jsl(1,m) = jsl(1,m) + jss(2,m)
      nps = max0(nps,jsl(1,m)+jsr(1,m))
  180 continue
      ibflg(2) = nps
c make sure sbufr and sbufl have been sent
      call MPI_WAIT(msid(3),istatus,ierr)
      call MPI_WAIT(msid(4),istatus,ierr)
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
      xt = rbufr(ic,j,m)
c particles going down or back
      if (xt.lt.edges(2*n-1,m)) then
         jsl(1,m) = jsl(1,m) + 1
         if (kb(n).eq.0) xt = xt + an
         rbufr(ic,j,m) = xt
         do 183 i = 1, idimp
         sbufl(i,jsl(1,m),m) = rbufr(i,j,m)
  183    continue
c particles going up or front, should not happen
      elseif (xt.ge.edges(2*n,m)) then
         jsr(1,m) = jsr(1,m) + 1
         if ((kb(n)+1).eq.nvp) xt = xt - an
         rbufr(ic,j,m) = xt
         do 185 i = 1, idimp
         sbufr(i,jsr(1,m),m) = rbufr(i,j,m)
  185    continue
c particles staying here
      else
         jss(2,m) = jss(2,m) + 1
         do 187 i = 1, idimp
         rbufr(i,jss(2,m),m) = rbufr(i,j,m)
  187    continue
      endif
  190 continue
      jsr(2,m) = jss(2,m)
c next check particles coming from below or back
      jss(2,m) = 0
      do 200 j = 1, jsl(2,m)
      xt = rbufl(ic,j,m)
c particles going up or front
      if (xt.ge.edges(2*n,m)) then
         if (jsr(1,m).lt.nbmax) then
            jsr(1,m) = jsr(1,m) + 1
            if ((kb(n)+1).eq.nvp) xt = xt - an
            rbufl(ic,j,m) = xt
            do 193 i = 1, idimp
            sbufr(i,jsr(1,m),m) = rbufl(i,j,m)
  193       continue
         else
            jss(2,m) = 2*npmax
            go to 210
         endif
c particles going down back, should not happen
      elseif (xt.lt.edges(2*n-1,m)) then
         if (jsl(1,m).lt.nbmax) then
            jsl(1,m) = jsl(1,m) + 1
            if (kb(n).eq.0) xt = xt + an
            rbufl(ic,j,m) = xt
            do 195 i = 1, idimp
            sbufl(i,jsl(1,m),m) = rbufl(i,j,m)
  195       continue
         else
            jss(2,m) = 2*npmax
            go to 210
         endif
c particles staying here
      else
         jss(2,m) = jss(2,m) + 1
         do 197 i = 1, idimp
         rbufl(i,jss(2,m),m) = rbufl(i,j,m)
  197    continue
      endif
  200 continue
  210 jsl(2,m) = jss(2,m)
  220 continue
  230 continue
c check if move would overflow particle array
  240 nps = 0
      npt = npmax
      do 250 m = 1, mnblok
      jss(2,m) = npp(m) + jsl(2,m) + jsr(2,m) - jss(1,m)
      nps = max0(nps,jss(2,m))
      npt = min0(npt,jss(2,m))
  250 continue
      ibflg(1) = nps
      ibflg(4) = -npt
      call PIMAX(ibflg,iwork,4,1)
      info(2) = ibflg(1)
      info(3) = -ibflg(4)
      ierr = ibflg(1) - npmax
      if (ierr.gt.0) then
         write (2,*) 'particle overflow error, ierr = ', ierr
         info(1) = ierr
         return
      endif
c distribute incoming particles from buffers
      do 290 m = 1, mnblok
c distribute particles coming from below or back into holes
      jss(2,m) = min0(jss(1,m),jsl(2,m))
      do 260 j = 1, jss(2,m)
      do 255 i = 1, idimp
      part(i,ihole(j,m),m) = rbufl(i,j,m)
  255 continue
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
         do 263 i = 1, idimp
         part(i,ihole(j+jsl(2,m),m),m) = rbufr(i,j,m)
  263    continue
      else
c no more holes
c distribute remaining particles from below or back into bottom
         do 267 i = 1, idimp
         part(i,j+npp(m),m) = rbufl(i,j+jss(1,m),m)
  267    continue
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
            do 273 i = 1, idimp
            part(i,ihole(j2,m),m) = part(i,j1,m)
  273       continue
         endif
      else
c no more holes
c distribute remaining particles from above or front into bottom
         do 277 i = 1, idimp
         part(i,j+npp(m),m) = rbufr(i,j+jss(1,m),m)
  277    continue
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
      info(5+n) = max0(info(5+n),mter)
      if (ibflg(2).gt.0) then
         write (2,*) 'Info: particles being passed further = ', ibflg(2)
         ! flush(2)
         if (ibflg(3).gt.0) ibflg(3) = 1
         if (iter.lt.itermax) go to 110
         ierr = -((iter-2)/2)
         write (2,*) 'Iteration overflow, iter = ', ierr
         info(1) = ierr
         go to 320
      endif
c check if buffer overflowed and more particles remain to be checked
      if (ibflg(3).gt.0) then
         nter = nter + 1
         info(3+n) = nter
         go to 20
      endif
  300 continue
c debugging section: count total number of particles after move
      nps = 0
      do 310 m = 1, mnblok
      nps = nps + npp(m)
  310 continue
      bflg(2) = dble(nps)
      bflg(1) = dble(npr)
      call PDSUM(bflg,work,2,1)
      info(8) = bflg(1)
      info(9) = bflg(2) - bflg(1)
      if (bflg(1).ne.bflg(2)) then
         write (2,*) 'particle number error, old/new=',bflg(1),bflg(2)
         info(1) = 1
      endif
c information
  320 if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WPMOVE32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,j
     1sr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,n
     2bmax,idds,ntmax,info)
c wrapper function for particle manager
c info(8) = total number of particles on entry
c info(9) = difference of total number of particles on exit
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl, th
      integer npp, ihole, jsr, jsl, jss, info
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
      dimension info(9)
c local data
      integer mnblok, j, m, n, npr, nps, nter
      real tf
      double precision dtime
      double precision bflg, work
      dimension bflg(2), work(2)
      mnblok = mblok*nblok
      do 10 j = 1, 9
      info(j) = 0
   10 continue
c debugging section: count total number of particles before move
      npr = 0
      do 20 m = 1, mnblok
      npr = npr + npp(m)
   20 continue
c find outgoing particles, first in y then in z direction
      do 40 n = 1, 2
   30 nter = info(n+3)
      call PWTIMERA(-1,tf,dtime)
      call PMOVEH32(part,edges,npp,ihole,jss,idimp,npmax,mblok,nblok,idp
     1s,idds,ntmax,n)
      call PWTIMERA(1,tf,dtime)
      th = th + tf
c send outgoing particles
      call PMOVES32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,jsl
     1,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbmax,idds
     2,ntmax,info,n)
c particle overflow error
      if (info(1).gt.0) return
c buffer overflowed and more particles remain to be checked
      if (info(n+3).gt.nter) go to 30
c iteration overflow
      if (info(1).lt.0) go to 60
   40 continue
c debugging section: count total number of particles after move
      nps = 0
      do 50 m = 1, mnblok
      nps = nps + npp(m)
   50 continue
      bflg(2) = dble(nps)
      bflg(1) = dble(npr)
      call PDSUM(bflg,work,2,1)
      info(8) = bflg(1)
      info(9) = bflg(2) - bflg(1)
      if (bflg(1).ne.bflg(2)) then
         if (kstrt.eq.1) then
           write (2,*) 'particle number error, old/new=',bflg(1),bflg(2)
         endif
         info(1) = 1
      endif
      return
c information
   60 nter = info(n+3)
      if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, n, nbmax=', n, 
     1nbmax
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WPXMOV32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,j
     1sr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,n
     2bmax,idds,ntmax,maskp,info)
c wrapper function for particle manager
c info(8) = total number of particles on entry
c info(9) = difference of total number of particles on exit
c optimized for vector processor
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl, th
      integer npp, ihole, jsr, jsl, jss, maskp, info
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
      dimension info(9)
c local data
      integer mnblok, j, m, n, npr, nps, nter
      real tf
      double precision dtime
      double precision bflg, work
      dimension bflg(2), work(2)
      mnblok = mblok*nblok
      do 10 j = 1, 9
      info(j) = 0
   10 continue
c debugging section: count total number of particles before move
      npr = 0
      do 20 m = 1, mnblok
      npr = npr + npp(m)
   20 continue
c find outgoing particles, first in y then in z direction
      do 40 n = 1, 2
   30 nter = info(n+3)
      call PWTIMERA(-1,tf,dtime)
      call PMOVEHX32(part,edges,npp,ihole,jss,idimp,npmax,mblok,nblok,id
     1ps,idds,ntmax,maskp,n)
      call PWTIMERA(1,tf,dtime)
      th = th + tf
c send outgoing particles
      call PMOVES32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,jsl
     1,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbmax,idds
     2,ntmax,info,n)
c particle overflow error
      if (info(1).gt.0) return
c buffer overflowed and more particles remain to be checked
      if (info(n+3).gt.nter) go to 30
c iteration overflow
      if (info(1).lt.0) go to 60
   40 continue
c debugging section: count total number of particles after move
      nps = 0
      do 50 m = 1, mnblok
      nps = nps + npp(m)
   50 continue
      bflg(2) = dble(nps)
      bflg(1) = dble(npr)
      call PDSUM(bflg,work,2,1)
      info(8) = bflg(1)
      info(9) = bflg(2) - bflg(1)
      if (bflg(1).ne.bflg(2)) then
         if (kstrt.eq.1) then
           write (2,*) 'particle number error, old/new=',bflg(1),bflg(2)
         endif
         info(1) = 1
      endif
      return
c information
   60 nter = info(n+3)
      if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, n, nbmax=', n, 
     1nbmax
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine WPMOVES32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,
     1jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,
     2nbmax,idds,ntmax,info)
c wrapper function for particle manager
c info(6) = maximum number of particle passes required in y
c info(6) must be set on entry
c info(7) = maximum number of particle passes required in z
c info(7) must be set on entry
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl, th
      integer npp, ihole, jsr, jsl, jss, info
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
      dimension info(9)
c local data
      integer mnblok, j, n
      real tf
      double precision dtime
      mnblok = mblok*nblok
      do 10 j = 1, 5
      info(j) = 0
   10 continue
c find outgoing particles, first in y then in z direction
      do 40 n = 1, 2
      call PWTIMERA(-1,tf,dtime)
      call PMOVEH32(part,edges,npp,ihole,jss,idimp,npmax,mblok,nblok,idp
     1s,idds,ntmax,n)
      call PWTIMERA(1,tf,dtime)
      th = th + tf
c send outgoing particles
      call PMOVESS32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,js
     1l,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbmax,idd
     2s,ntmax,info,n)
c particle or iteration overflow error
      if (info(1).ne.0) return
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WPXMOVS32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,
     1jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,
     2nbmax,idds,ntmax,maskp,info)
c wrapper function for particle manager
c info(6) = maximum number of particle passes required in y
c info(6) must be set on entry
c info(7) = maximum number of particle passes required in z
c info(7) must be set on entry
c optimized for vector processor
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl, th
      integer npp, ihole, jsr, jsl, jss, maskp, info
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
      dimension info(9)
c local data
      integer mnblok, j, n
      real tf
      double precision dtime
      mnblok = mblok*nblok
      do 10 j = 1, 5
      info(j) = 0
   10 continue
c find outgoing particles, first in y then in z direction
      do 40 n = 1, 2
      call PWTIMERA(-1,tf,dtime)
      call PMOVEHX32(part,edges,npp,ihole,jss,idimp,npmax,mblok,nblok,id
     1ps,idds,ntmax,maskp,n)
      call PWTIMERA(1,tf,dtime)
      th = th + tf
c send outgoing particles
      call PMOVESS32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,js
     1l,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbmax,idd
     2s,ntmax,info,n)
c particle or iteration overflow error
      if (info(1).ne.0) return
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMOVEH32(part,edges,npp,ihole,jss,idimp,npmax,mblok,nbl
     1ok,idps,idds,ntmax,n)
c this subroutine determines list of particles which are leaving this
c processor, with 2D spatial decomposition
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c edges(1,m) = lower boundary in y of particle partition m
c edges(2,m) = upper boundary in y of particle partition m
c edges(3,m) = back boundary in z of particle partition m
c edges(4,m) = front boundary in z of particle partition m
c npp(m) = number of particles in partition m
c ihole = location of holes left in particle arrays
c jss(1,m) = number of holes found in partition m
c jss(2,m) = (0,1) = buffer ihole (OK, overflowed) in partition m
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mblok/nblok = number of particle partitions in y/z
c idps = number of particle partition boundaries
c idds = dimensionality of domain decomposition
c ntmax =  size of hole array for particles leaving processors
c n = partition index (n=1,2) for (y,z)
      implicit none
      real part, edges
      integer npp, ihole, jss
      integer idimp, npmax, mblok, nblok, idps, idds, ntmax, n
      dimension part(idimp,npmax,mblok*nblok)
      dimension edges(idps,mblok*nblok), npp(mblok*nblok)
      dimension jss(idds,mblok*nblok)
      dimension ihole(ntmax,mblok*nblok)
c local data
c iy, iz = partitioned co-ordinates
      integer iy, iz
      parameter(iy=2,iz=3)
      integer ic, mnblok, j, m, my, mz, moff
      real xt
      mnblok = mblok*nblok
      if (n.eq.1) then
         ic = iy
      elseif (n.eq.2) then
         ic = iz
      endif
c find particles out of bounds
      do 30 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 20 my = 1, mblok
      m = my + moff
      jss(1,m) = 0
      jss(2,m) = 0
      do 10 j = 1, npp(m)
      xt = part(ic,j,m)
      if ((xt.ge.edges(2*n,m)).or.(xt.lt.edges(2*n-1,m))) then
         if (jss(1,m).lt.ntmax) then
            jss(1,m) = jss(1,m) + 1
            ihole(jss(1,m),m) = j
         else
            jss(2,m) = 1
            go to 20
         endif
      endif
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMOVEHX32(part,edges,npp,ihole,jss,idimp,npmax,mblok,nb
     1lok,idps,idds,ntmax,maskp,n)
c this subroutine determines list of particles which are leaving this
c processor, with 2D spatial decomposition
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c edges(1,m) = lower boundary in y of particle partition m
c edges(2,m) = upper boundary in y of particle partition m
c edges(3,m) = back boundary in z of particle partition m
c edges(4,m) = front boundary in z of particle partition m
c npp(m) = number of particles in partition m
c ihole = location of holes left in particle arrays
c jss(1,m) = number of holes found in partition m
c jss(2,m) = (0,1) = buffer ihole (OK, overflowed) in partition m
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mblok/nblok = number of particle partitions in y/z
c idps = number of particle partition boundaries
c idds = dimensionality of domain decomposition
c ntmax =  size of hole array for particles leaving processors
c maskp = scratch array for particle addresses
c n = partition index (n=1,2) for (y,z)
c optimized for vector processor
      implicit none
      real part, edges
      integer npp, ihole, jss, maskp
      integer idimp, npmax, mblok, nblok, idps, idds, ntmax, n
      dimension part(idimp,npmax,mblok*nblok), maskp(npmax,mblok*nblok)
      dimension edges(idps,mblok*nblok), npp(mblok*nblok)
      dimension jss(idds,mblok*nblok)
      dimension ihole(ntmax,mblok*nblok)
c local data
c iy, iz = partitioned co-ordinates
      integer iy, iz
      parameter(iy=2,iz=3)
      integer ic, mnblok, j, m, my, mz, moff
      real xt
      mnblok = mblok*nblok
      if (n.eq.1) then
         ic = iy
      elseif (n.eq.2) then
         ic = iz
      endif
c find particles out of bounds
      do 50 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 40 my = 1, mblok
      m = my + moff
      jss(1,m) = 0
      jss(2,m) = 0
c find mask function for particles out of bounds
      do 10 j = 1, npp(m)
      xt = part(ic,j,m)
      if ((xt.ge.edges(2*n,m)).or.(xt.lt.edges(2*n-1,m))) then
         jss(1,m) = jss(1,m) + 1
         maskp(j,m) = 1
      else
         maskp(j,m) = 0
      endif
   10 continue
c set flag if hole buffer would overflow
      if (jss(1,m).gt.ntmax) then
         jss(1,m) = ntmax
         jss(2,m) = 1
      endif
c accumulate location of holes
      do 20 j = 2, npp(m)
      maskp(j,m) = maskp(j,m) + maskp(j-1,m)
   20 continue
c store addresses of particles out of bounds
      do 30 j = 2, npp(m)
      if ((maskp(j,m).gt.maskp(j-1,m)).and.(maskp(j,m).le.ntmax)) then
         ihole(maskp(j,m),m) = j
      endif
   30 continue
      if (maskp(1,m).gt.0) ihole(1,m) = 1
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PMOVES32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,j
     1sr,jsl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbma
     2x,idds,ntmax,info,n)
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
c on input, jss(1,m) = number of holes found in partition m
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
c info = status information
c info(1) = ierr = (0,N) = (no,yes) error condition exists
c info(2) = maximum number of particles per processor
c info(3) = minimum number of particles per processor
c info(4) = maximum number of buffer overflows in y
c info(5) = maximum number of buffer overflows in z
c info(6) = maximum number of particle passes required in y
c info(7) = maximum number of particle passes required in z
c n = partition index (n=1,2) for (y,z)
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl
      integer npp, ihole, jsr, jsl, jss, info
      integer ny, nz, kstrt, nvpy, nvpz, idimp, npmax, mblok, nblok
      integer idps, nbmax, idds, ntmax, n
      dimension part(idimp,npmax,mblok*nblok)
      dimension edges(idps,mblok*nblok), npp(mblok*nblok)
      dimension sbufl(idimp,nbmax,mblok*nblok)
      dimension sbufr(idimp,nbmax,mblok*nblok)
      dimension rbufl(idimp,nbmax,mblok*nblok)
      dimension rbufr(idimp,nbmax,mblok*nblok)
      dimension jsl(idds,mblok*nblok), jsr(idds,mblok*nblok)
      dimension jss(idds,mblok*nblok)
      dimension ihole(ntmax,mblok*nblok)
      dimension info(9)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
c iy, iz = partitioned co-ordinates
      integer iy, iz
      parameter(iy=2,iz=3)
      integer ierr, ic, js, ks, mnblok, i, m, my, mz, moff, nvp, iter
      integer nps, npt, kb, kl, kr, j, j1, j2, nbsize, nter, mter
      integer itermax
      integer msid, istatus
      integer ibflg, iwork
      real an, xt
      dimension msid(4), istatus(lstat)
      dimension ibflg(4), iwork(4)
      dimension kb(2)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      mnblok = mblok*nblok
      nbsize = idimp*nbmax
      if (n.eq.1) then
         ic = iy
         nvp = nvpy
         an = float(ny)
      elseif (n.eq.2) then
         ic = iz
         nvp = nvpz
         an = float(nz)
      endif
      iter = 2
      nter = info(3+n)
      itermax = 2000
      mter = 0
c buffer outgoing particles
      do 60 mz = 1, nblok
      moff = mblok*(mz - 1)
      kb(2) = mz + ks
      do 50 my = 1, mblok
      m = my + moff
      kb(1) = my + js
      jsl(1,m) = 0
      jsr(1,m) = 0
c     jss(2,m) = 0
      do 30 j = 1, jss(1,m)
      xt = part(ic,ihole(j,m),m)
c particles going down or backward
      if (xt.lt.edges(2*n-1,m)) then
         if (kb(n).eq.0) xt = xt + an
         if (jsl(1,m).lt.nbmax) then
            jsl(1,m) = jsl(1,m) + 1
            do 10 i = 1, idimp
            sbufl(i,jsl(1,m),m) = part(i,ihole(j,m),m)
   10       continue
            sbufl(ic,jsl(1,m),m) = xt
            ihole(jsl(1,m)+jsr(1,m),m) = ihole(j,m)
         else
            jss(2,m) = 1
            go to 40
         endif
c particles going up or forward
      else if (xt.ge.edges(2*n,m)) then
         if ((kb(n)+1).eq.nvp) xt = xt - an
         if (jsr(1,m).lt.nbmax) then
            jsr(1,m) = jsr(1,m) + 1
            do 20 i = 1, idimp
            sbufr(i,jsr(1,m),m) = part(i,ihole(j,m),m)
   20       continue
            sbufr(ic,jsr(1,m),m) = xt
            ihole(jsl(1,m)+jsr(1,m),m) = ihole(j,m)
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
      nps = max0(nps,jss(2,m))
  100 continue
      ibflg(3) = nps
c copy particle buffers
  110 iter = iter + 2
      mter = mter + 1
      do 170 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 160 my = 1, mblok
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
c     do 130 j = 1, jsl(2,m)
c     do 120 i = 1, idimp
c     rbufl(i,j,m) = sbufr(i,j,kl)
c 120 continue
c 130 continue
c     jsr(2,m) = jsl(1,kr)
c     do 150 j = 1, jsr(2,m)
c     do 140 i = 1, idimp
c     rbufr(i,j,m) = sbufl(i,j,kr)
c 140 continue
c 150 continue
c this segment is used for mpi computers
c post receive
      call MPI_IRECV(rbufl,nbsize,mreal,kl-1,iter-1,lgrp,msid(1),ierr)
      call MPI_IRECV(rbufr,nbsize,mreal,kr-1,iter,lgrp,msid(2),ierr)
c send particles
      call MPI_ISEND(sbufr,idimp*jsr(1,m),mreal,kr-1,iter-1,lgrp,msid(3)
     1,ierr)
      call MPI_ISEND(sbufl,idimp*jsl(1,m),mreal,kl-1,iter,lgrp,msid(4),i
     1err)
c wait for particles to arrive
      call MPI_WAIT(msid(1),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsl(2,m) = nps/idimp
      call MPI_WAIT(msid(2),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsr(2,m) = nps/idimp
  160 continue
  170 continue
c check if particles must be passed further
      nps = 0
      do 200 m = 1, mnblok
c check if any particles coming from above or front belong here
      jsl(1,m) = 0
      jsr(1,m) = 0
      jss(2,m) = 0
      do 180 j = 1, jsr(2,m)
      if (rbufr(ic,j,m).lt.edges(2*n-1,m)) jsl(1,m) = jsl(1,m) + 1
      if (rbufr(ic,j,m).ge.edges(2*n,m)) jsr(1,m) = jsr(1,m) + 1
  180 continue
      if (jsr(1,m).ne.0) then
         if (n.eq.1) then
            write (2,*) 'Info:',jsr(1,m),' particles returning above'
         elseif (n.eq.2) then
            write (2,*) 'Info:',jsr(1,m),' particles returning front'
         endif
      endif
c check if any particles coming from below or back belong here
      do 190 j = 1, jsl(2,m)
      if (rbufl(ic,j,m).ge.edges(2*n,m)) jsr(1,m) = jsr(1,m) + 1
      if (rbufl(ic,j,m).lt.edges(2*n-1,m)) jss(2,m) = jss(2,m) + 1
  190 continue
      if (jss(2,m).ne.0) then
         if (n.eq.1) then
            write (2,*) 'Info:',jss(2,m),' particles returning below'
         elseif (n.eq.2) then
            write (2,*) 'Info:',jss(2,m),' particles returning back'
         endif
      endif
      jsl(1,m) = jsl(1,m) + jss(2,m)
      nps = max0(nps,jsl(1,m)+jsr(1,m))
  200 continue
      ibflg(2) = nps
c make sure sbufr and sbufl have been sent
      call MPI_WAIT(msid(3),istatus,ierr)
      call MPI_WAIT(msid(4),istatus,ierr)
      if (nps.eq.0) go to 320
c remove particles which do not belong here
      do 310 mz = 1, nblok
      moff = mblok*(mz - 1)
      kb(2) = mz + ks
      do 300 my = 1, mblok
      m = my + moff
      kb(1) = my + js
c first check particles coming from above or front
      jsl(1,m) = 0
      jsr(1,m) = 0
      jss(2,m) = 0
      do 240 j = 1, jsr(2,m)
      xt = rbufr(ic,j,m)
c particles going down or back
      if (xt.lt.edges(2*n-1,m)) then
         jsl(1,m) = jsl(1,m) + 1
         if (kb(n).eq.0) xt = xt + an
         rbufr(ic,j,m) = xt
         do 210 i = 1, idimp
         sbufl(i,jsl(1,m),m) = rbufr(i,j,m)
  210    continue
c particles going up or front, should not happen
      elseif (xt.ge.edges(2*n,m)) then
         jsr(1,m) = jsr(1,m) + 1
         if ((kb(n)+1).eq.nvp) xt = xt - an
         rbufr(ic,j,m) = xt
         do 220 i = 1, idimp
         sbufr(i,jsr(1,m),m) = rbufr(i,j,m)
  220    continue
c particles staying here
      else
         jss(2,m) = jss(2,m) + 1
         do 230 i = 1, idimp
         rbufr(i,jss(2,m),m) = rbufr(i,j,m)
  230    continue
      endif
  240 continue
      jsr(2,m) = jss(2,m)
c next check particles coming from below or back
      jss(2,m) = 0
      do 280 j = 1, jsl(2,m)
      xt = rbufl(ic,j,m)
c particles going up or front
      if (xt.ge.edges(2*n,m)) then
         if (jsr(1,m).lt.nbmax) then
            jsr(1,m) = jsr(1,m) + 1
            if ((kb(n)+1).eq.nvp) xt = xt - an
            rbufl(ic,j,m) = xt
            do 250 i = 1, idimp
            sbufr(i,jsr(1,m),m) = rbufl(i,j,m)
  250       continue 
         else
            jss(2,m) = 2*npmax
            go to 290
         endif
c particles going down back, should not happen
      elseif (xt.lt.edges(2*n-1,m)) then
         if (jsl(1,m).lt.nbmax) then
            jsl(1,m) = jsl(1,m) + 1
            if (kb(n).eq.0) xt = xt + an
            rbufl(ic,j,m) = xt
            do 260 i = 1, idimp
            sbufl(i,jsl(1,m),m) = rbufl(i,j,m)
  260       continue
         else
            jss(2,m) = 2*npmax
            go to 290
         endif
c particles staying here
      else
         jss(2,m) = jss(2,m) + 1
         do 270 i = 1, idimp
         rbufl(i,jss(2,m),m) = rbufl(i,j,m)
  270    continue
      endif
  280 continue
  290 jsl(2,m) = jss(2,m)
  300 continue
  310 continue
c check if move would overflow particle array
  320 nps = 0
      npt = npmax
      do 330 m = 1, mnblok
      jss(2,m) = npp(m) + jsl(2,m) + jsr(2,m) - jss(1,m)
      nps = max0(nps,jss(2,m))
      npt = min0(npt,jss(2,m))
  330 continue
      ibflg(1) = nps
      ibflg(4) = -npt
      call PIMAX(ibflg,iwork,4,1)
      info(2) = ibflg(1)
      info(3) = -ibflg(4)
      ierr = ibflg(1) - npmax
      if (ierr.gt.0) then
         write (2,*) 'particle overflow error, ierr = ', ierr
         info(1) = ierr
         return
      endif
c distribute incoming particles from buffers
      do 420 m = 1, mnblok
c distribute particles coming from below or back into holes
      jss(2,m) = min0(jss(1,m),jsl(2,m))
      do 350 j = 1, jss(2,m)
      do 340 i = 1, idimp
      part(i,ihole(j,m),m) = rbufl(i,j,m)
  340 continue
  350 continue
      if (jss(1,m).gt.jsl(2,m)) then
         jss(2,m) = min0(jss(1,m)-jsl(2,m),jsr(2,m))
      else
         jss(2,m) = jsl(2,m) - jss(1,m)
      endif
      do 380 j = 1, jss(2,m)
c no more particles coming from below or back
c distribute particles coming from above or front into holes
      if (jss(1,m).gt.jsl(2,m)) then
         do 360 i = 1, idimp
         part(i,ihole(j+jsl(2,m),m),m) = rbufr(i,j,m)
  360    continue
      else
c no more holes
c distribute remaining particles from below or back into bottom
         do 370 i = 1, idimp
         part(i,j+npp(m),m) = rbufl(i,j+jss(1,m),m)
  370    continue
      endif
  380 continue
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
      do 410 j = 1, jsr(2,m)
c holes left over
c fill up remaining holes in particle array with particles from bottom
      if (jss(2,m).gt.0) then
         j1 = npp(m) - j + 1
         j2 = jss(1,m) + jss(2,m) - j + 1
         if (j1.gt.ihole(j2,m)) then
c move particle only if it is below current hole
            do 390 i = 1, idimp
            part(i,ihole(j2,m),m) = part(i,j1,m)
  390       continue
         endif
      else
c no more holes
c distribute remaining particles from above or front into bottom
         do 400 i = 1, idimp
         part(i,j+npp(m),m) = rbufr(i,j+jss(1,m),m)
  400    continue
      endif
  410 continue
      if (jss(2,m).gt.0) then
         npp(m) = npp(m) - jsr(2,m)
      else
         npp(m) = npp(m) + jsr(2,m)
      endif
      jss(1,m) = 0
  420 continue
c check if any particles have to be passed further
      info(5+n) = max0(info(5+n),mter)
      if (ibflg(2).gt.0) then
         write (2,*) 'Info: particles being passed further = ', ibflg(2)
         ! flush(2)
         if (ibflg(3).gt.0) ibflg(3) = 1
         if (iter.lt.itermax) go to 110
         ierr = -((iter-2)/2)
         write (2,*) 'Iteration overflow, iter = ', ierr
         info(1) = ierr
         return
      endif
c check if buffer overflowed and more particles remain to be checked
      if (ibflg(3).gt.0) then
         nter = nter + 1
         info(3+n) = nter
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PMOVESS32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,
     1jsr,jsl,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbm
     2ax,idds,ntmax,info,n)
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
c on input, jss(1,m) = number of holes found in partition m
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
c info = status information
c info(1) = ierr = (0,N) = (no,yes) error condition exists
c info(2) = maximum number of particles per processor
c info(3) = minimum number of particles per processor
c info(4) = maximum number of buffer overflows in y
c info(5) = maximum number of buffer overflows in z
c info(6) = maximum number of particle passes required in y
c info(6) must be set on entry
c info(7) = maximum number of particle passes required in z
c info(7) must be set on entry
c n = partition index (n=1,2) for (y,z)
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl
      integer npp, ihole, jsr, jsl, jss, info
      integer ny, nz, kstrt, nvpy, nvpz, idimp, npmax, mblok, nblok
      integer idps, nbmax, idds, ntmax, n
      dimension part(idimp,npmax,mblok*nblok)
      dimension edges(idps,mblok*nblok), npp(mblok*nblok)
      dimension sbufl(idimp,nbmax,mblok*nblok)
      dimension sbufr(idimp,nbmax,mblok*nblok)
      dimension rbufl(idimp,nbmax,mblok*nblok)
      dimension rbufr(idimp,nbmax,mblok*nblok)
      dimension jsl(idds,mblok*nblok), jsr(idds,mblok*nblok)
      dimension jss(idds,mblok*nblok)
      dimension ihole(ntmax,mblok*nblok)
      dimension info(9)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
c iy, iz = partitioned co-ordinates
      integer iy, iz
      parameter(iy=2,iz=3)
      integer ierr, ic, js, ks, mnblok, i, m, my, mz, moff, nvp, iter
      integer nps, npt, kb, kl, kr, j, j1, j2, nbsize, nter, mter
      integer itermax
      integer msid, istatus
      integer ibflg
      real an, xt
      dimension msid(4), istatus(lstat)
      dimension ibflg(4)
      dimension kb(2)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      mnblok = mblok*nblok
      nbsize = idimp*nbmax
      if (n.eq.1) then
         ic = iy
         nvp = nvpy
         an = float(ny)
      elseif (n.eq.2) then
         ic = iz
         nvp = nvpz
         an = float(nz)
      endif
      iter = 2
c     nter = info(3+n)
      nter = 0
      itermax = 2000
      mter = 0
c buffer outgoing particles
      do 60 mz = 1, nblok
      moff = mblok*(mz - 1)
      kb(2) = mz + ks
      do 50 my = 1, mblok
      m = my + moff
      kb(1) = my + js
      jsl(1,m) = 0
      jsr(1,m) = 0
c     jss(2,m) = 0
      do 30 j = 1, jss(1,m)
      xt = part(ic,ihole(j,m),m)
c particles going down or backward
      if (xt.lt.edges(2*n-1,m)) then
         if (kb(n).eq.0) xt = xt + an
         if (jsl(1,m).lt.nbmax) then
            jsl(1,m) = jsl(1,m) + 1
            do 10 i = 1, idimp
            sbufl(i,jsl(1,m),m) = part(i,ihole(j,m),m)
   10       continue
            sbufl(ic,jsl(1,m),m) = xt
            ihole(jsl(1,m)+jsr(1,m),m) = ihole(j,m)
         else
            jss(2,m) = 1
            go to 40
         endif
c particles going up or forward
      else if (xt.ge.edges(2*n,m)) then
         if ((kb(n)+1).eq.nvp) xt = xt - an
         if (jsr(1,m).lt.nbmax) then
            jsr(1,m) = jsr(1,m) + 1
            do 20 i = 1, idimp
            sbufr(i,jsr(1,m),m) = part(i,ihole(j,m),m)
   20       continue
            sbufr(ic,jsr(1,m),m) = xt
            ihole(jsl(1,m)+jsr(1,m),m) = ihole(j,m)
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
      nps = max0(nps,jss(2,m))
  100 continue
      ibflg(3) = nps
c copy particle buffers
  110 iter = iter + 2
      mter = mter + 1
      do 170 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 160 my = 1, mblok
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
c     do 130 j = 1, jsl(2,m)
c     do 120 i = 1, idimp
c     rbufl(i,j,m) = sbufr(i,j,kl)
c 120 continue
c 130 continue
c     jsr(2,m) = jsl(1,kr)
c     do 150 j = 1, jsr(2,m)
c     do 140 i = 1, idimp
c     rbufr(i,j,m) = sbufl(i,j,kr)
c 140 continue
c 150 continue
c this segment is used for mpi computers
c post receive
      call MPI_IRECV(rbufl,nbsize,mreal,kl-1,iter-1,lgrp,msid(1),ierr)
      call MPI_IRECV(rbufr,nbsize,mreal,kr-1,iter,lgrp,msid(2),ierr)
c send particles
      call MPI_ISEND(sbufr,idimp*jsr(1,m),mreal,kr-1,iter-1,lgrp,msid(3)
     1,ierr)
      call MPI_ISEND(sbufl,idimp*jsl(1,m),mreal,kl-1,iter,lgrp,msid(4),i
     1err)
c wait for particles to arrive
      call MPI_WAIT(msid(1),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsl(2,m) = nps/idimp
      call MPI_WAIT(msid(2),istatus,ierr)
      call MPI_GET_COUNT(istatus,mreal,nps,ierr)
      jsr(2,m) = nps/idimp
  160 continue
  170 continue
c check if particles must be passed further
      nps = 0
      do 200 m = 1, mnblok
c check if any particles coming from above or front belong here
      jsl(1,m) = 0
      jsr(1,m) = 0
      jss(2,m) = 0
      do 180 j = 1, jsr(2,m)
      if (rbufr(ic,j,m).lt.edges(2*n-1,m)) jsl(1,m) = jsl(1,m) + 1
      if (rbufr(ic,j,m).ge.edges(2*n,m)) jsr(1,m) = jsr(1,m) + 1
  180 continue
      if (jsr(1,m).ne.0) then
         if (n.eq.1) then
            write (2,*) 'Info:',jsr(1,m),' particles returning above'
         elseif (n.eq.2) then
            write (2,*) 'Info:',jsr(1,m),' particles returning front'
         endif
      endif
c check if any particles coming from below or back belong here
      do 190 j = 1, jsl(2,m)
      if (rbufl(ic,j,m).ge.edges(2*n,m)) jsr(1,m) = jsr(1,m) + 1
      if (rbufl(ic,j,m).lt.edges(2*n-1,m)) jss(2,m) = jss(2,m) + 1
  190 continue
      if (jss(2,m).ne.0) then
         if (n.eq.1) then
            write (2,*) 'Info:',jss(2,m),' particles returning below'
         elseif (n.eq.2) then
            write (2,*) 'Info:',jss(2,m),' particles returning back'
         endif
      endif
      jsl(1,m) = jsl(1,m) + jss(2,m)
      nps = max0(nps,jsl(1,m)+jsr(1,m))
  200 continue
      ibflg(2) = nps
c make sure sbufr and sbufl have been sent
      call MPI_WAIT(msid(3),istatus,ierr)
      call MPI_WAIT(msid(4),istatus,ierr)
      if ((nps.eq.0).or.(mter.eq.info(5+n))) go to 320
c remove particles which do not belong here
      do 310 mz = 1, nblok
      moff = mblok*(mz - 1)
      kb(2) = mz + ks
      do 300 my = 1, mblok
      m = my + moff
      kb(1) = my + js
c first check particles coming from above or front
      jsl(1,m) = 0
      jsr(1,m) = 0
      jss(2,m) = 0
      do 240 j = 1, jsr(2,m)
      xt = rbufr(ic,j,m)
c particles going down or back
      if (xt.lt.edges(2*n-1,m)) then
         jsl(1,m) = jsl(1,m) + 1
         if (kb(n).eq.0) xt = xt + an
         rbufr(ic,j,m) = xt
         do 210 i = 1, idimp
         sbufl(i,jsl(1,m),m) = rbufr(i,j,m)
  210    continue
c particles going up or front, should not happen
      elseif (xt.ge.edges(2*n,m)) then
         jsr(1,m) = jsr(1,m) + 1
         if ((kb(n)+1).eq.nvp) xt = xt - an
         rbufr(ic,j,m) = xt
         do 220 i = 1, idimp
         sbufr(i,jsr(1,m),m) = rbufr(i,j,m)
  220    continue
c particles staying here
      else
         jss(2,m) = jss(2,m) + 1
         do 230 i = 1, idimp
         rbufr(i,jss(2,m),m) = rbufr(i,j,m)
  230    continue
      endif
  240 continue
      jsr(2,m) = jss(2,m)
c next check particles coming from below or back
      jss(2,m) = 0
      do 280 j = 1, jsl(2,m)
      xt = rbufl(ic,j,m)
c particles going up or front
      if (xt.ge.edges(2*n,m)) then
         if (jsr(1,m).lt.nbmax) then
            jsr(1,m) = jsr(1,m) + 1
            if ((kb(n)+1).eq.nvp) xt = xt - an
            rbufl(ic,j,m) = xt
            do 250 i = 1, idimp
            sbufr(i,jsr(1,m),m) = rbufl(i,j,m)
  250       continue 
         else
            jss(2,m) = 2*npmax
            go to 290
         endif
c particles going down back, should not happen
      elseif (xt.lt.edges(2*n-1,m)) then
         if (jsl(1,m).lt.nbmax) then
            jsl(1,m) = jsl(1,m) + 1
            if (kb(n).eq.0) xt = xt + an
            rbufl(ic,j,m) = xt
            do 260 i = 1, idimp
            sbufl(i,jsl(1,m),m) = rbufl(i,j,m)
  260       continue
         else
            jss(2,m) = 2*npmax
            go to 290
         endif
c particles staying here
      else
         jss(2,m) = jss(2,m) + 1
         do 270 i = 1, idimp
         rbufl(i,jss(2,m),m) = rbufl(i,j,m)
  270    continue
      endif
  280 continue
  290 jsl(2,m) = jss(2,m)
  300 continue
  310 continue
c check if move would overflow particle array
  320 nps = 0
      npt = npmax
      do 330 m = 1, mnblok
      jss(2,m) = npp(m) + jsl(2,m) + jsr(2,m) - jss(1,m)
      nps = max0(nps,jss(2,m))
      npt = min0(npt,jss(2,m))
  330 continue
      ibflg(1) = nps
      ibflg(4) = -npt
c     call PIMAX(ibflg,iwork,4,1)
      info(2) = ibflg(1)
      info(3) = -ibflg(4)
      ierr = ibflg(1) - npmax
      if (ierr.gt.0) then
         write (2,*) 'particle overflow error, ierr = ', ierr
         info(1) = ierr
         return
      endif
c distribute incoming particles from buffers
      do 420 m = 1, mnblok
c distribute particles coming from below or back into holes
      jss(2,m) = min0(jss(1,m),jsl(2,m))
      do 350 j = 1, jss(2,m)
      do 340 i = 1, idimp
      part(i,ihole(j,m),m) = rbufl(i,j,m)
  340 continue
  350 continue
      if (jss(1,m).gt.jsl(2,m)) then
         jss(2,m) = min0(jss(1,m)-jsl(2,m),jsr(2,m))
      else
         jss(2,m) = jsl(2,m) - jss(1,m)
      endif
      do 380 j = 1, jss(2,m)
c no more particles coming from below or back
c distribute particles coming from above or front into holes
      if (jss(1,m).gt.jsl(2,m)) then
         do 360 i = 1, idimp
         part(i,ihole(j+jsl(2,m),m),m) = rbufr(i,j,m)
  360    continue
      else
c no more holes
c distribute remaining particles from below or back into bottom
         do 370 i = 1, idimp
         part(i,j+npp(m),m) = rbufl(i,j+jss(1,m),m)
  370    continue
      endif
  380 continue
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
      do 410 j = 1, jsr(2,m)
c holes left over
c fill up remaining holes in particle array with particles from bottom
      if (jss(2,m).gt.0) then
         j1 = npp(m) - j + 1
         j2 = jss(1,m) + jss(2,m) - j + 1
         if (j1.gt.ihole(j2,m)) then
c move particle only if it is below current hole
            do 390 i = 1, idimp
            part(i,ihole(j2,m),m) = part(i,j1,m)
  390       continue
         endif
      else
c no more holes
c distribute remaining particles from above or front into bottom
         do 400 i = 1, idimp
         part(i,j+npp(m),m) = rbufr(i,j+jss(1,m),m)
  400    continue
      endif
  410 continue
      if (jss(2,m).gt.0) then
         npp(m) = npp(m) - jsr(2,m)
      else
         npp(m) = npp(m) + jsr(2,m)
      endif
      jss(1,m) = 0
  420 continue
c check if any particles have to be passed further
c     info(5+n) = max0(info(5+n),mter)
      if (mter.eq.info(5+n)) then
         if (ibflg(2).gt.0) then
            write (2,*) 'Exceeded maximum number of passes = ', n, mter
            info(1) = -1
            return
         endif
c        write (2,*) 'Info: particles being passed further = ', ibflg(2)
      else
         if (ibflg(3).gt.0) ibflg(3) = 1
         if (iter.lt.itermax) go to 110
         ierr = -((iter-2)/2)
         write (2,*) 'Iteration overflow, iter = ', ierr
         info(1) = ierr
         return
      endif
c check if buffer overflowed and more particles remain to be checked
      if (ibflg(3).gt.0) then
         nter = nter + 1
         info(3+n) = nter
         write (2,*) 'Buffer overflow = ', n, nter
         info(1) = -2
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PFMOVE32(f,g,h,noff,nyzp,noffs,nyzps,noffd,nyzpd,jsr,js
     1l,isign,kyp,kzp,kstrt,nvpy,nvpz,nxv,nypmx,nzpmx,mblok,nblok,idds,m
     2ter,nter,ierr)
c this subroutine moves fields into appropriate spatial regions,
c between non-uniform and uniform partitions
c f(j,k,l,m) = real data for grid j,k,l in field partition m.
c the grid is non-uniform and includes extra guard cells.
c g(j,k,l,m) = scratch data for grid j,k,l in field partition m.
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c noffs(m)/nzyps(m) = source or scratch arrays for field partition m
c noffd(m)/nyzpd(m) = destination or scratch arrays for field partition m
c jsl(idds,m) = number of particles going down in field partition m
c jsr(idds,m) = number of particles going up in field partition m
c isign = -1, move from non-uniform (noff/nyzp) to uniform (kyp/kzp)
c    fields
c isign = 1, move from uniform (kyp/kzp) to non-uniform (noff/nyzp)
c    fields
c if isign = 0, the noffs/nyzps contains the source partition, 
c    noffd/nyzpd contains the destination partition, and  kyp, kzp are
c    not used.  the source partitions noffs/nyzps and noff/nyzp are
c    modified.  isign = 2 is a special optimized version of isign = 0,
c    which can be used when the partitions in z are the same for each y,
c    where noff/nyzp are not modified.
c kyp/kzp = number of complex grids in y/z for each field partition.
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = first dimension of f, must be >= nx
c nypmx = maximum size of particle partition in y, must be >= kyp + 1
c nzpmx = maximum size of particle partition in z, must be >= kzp + 1
c mblok/nblok = number of particle partitions in y/z
c idds = dimensionality of domain decomposition
c mter/nter = number of shifts required in y/z
c if mter/nter = 0, then number of shifts is determined and returned
c ierr = (0,1) = (no,yes) error condition exists
      implicit none
      real f, g, h
      integer noff, nyzp, noffs, nyzps, noffd, nyzpd, jsr, jsl
      integer isign, kyp, kzp, kstrt, nvpy, nvpz, nxv, nypmx, nzpmx
      integer mblok, nblok, idds, mter, nter, ierr
      dimension f(nxv,nypmx,nzpmx,mblok*nblok)
      dimension g(nxv,nypmx*nzpmx,mblok*nblok)
      dimension h(nxv,nypmx*nzpmx,mblok*nblok)
      dimension noff(idds,mblok*nblok), nyzp(idds,mblok*nblok)
      dimension noffs(idds,mblok*nblok), nyzps(idds,mblok*nblok)
      dimension noffd(idds,mblok*nblok), nyzpd(idds,mblok*nblok)
      dimension jsl(idds,mblok*nblok), jsr(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer j, k, l, my, mz, m, n, js, ks, moff, koff, loff, mnblok
      integer mnter, nxvyz, nbsize, iter, npr, nps, nnter, kl, kr, kk
      integer ll, nn, ne, nypm, nzpm, nyzpmn
      integer msid, istatus
      integer ibflg, iwork
      dimension istatus(lstat)
      dimension ibflg(2), iwork(2)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      mnblok = mblok*nblok
      nbsize = nxv*nypmx*nzpmx
      ne = 2
      ierr = 0
      do 20 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 10 my = 1, mblok
      m = my + moff
c move from non-uniform to uniform fields
      if (isign.lt.0) then
         noffs(1,m) = noff(1,m)
         noffs(2,m) = noff(2,m)
         nyzps(1,m) = nyzp(1,m)
         nyzps(2,m) = nyzp(2,m)
         noffd(1,m) = kyp*(my + js)
         noffd(2,m) = kzp*(mz + ks)
         nyzpd(1,m) = kyp
         nyzpd(2,m) = kzp
c move from uniform to non-uniform fields
      else if (isign.eq.1) then
         noffs(2,m) = kzp*(mz + ks)
         noffs(1,m) = kyp*(my + js)
         nyzps(1,m) = kyp
         nyzps(2,m) = kzp
         noffd(1,m) = noff(1,m)
         noffd(2,m) = noff(2,m)
         nyzpd(1,m) = nyzp(1,m)
         nyzpd(2,m) = nyzp(2,m)
c move from non-uniform to non-uniform fields
      else if (isign.eq.0) then
         noff(2,m) = noffd(2,m)
         nyzp(2,m) = nyzpd(2,m)
         noffd(2,m) = kzp*(mz + ks)
         nyzpd(2,m) = kzp
         ne = ne + 1
      endif
c extend partitions to include (ny+1,nz+1) grids
      if ((my+js).eq.(nvpy-1)) then
         nyzps(1,m) = nyzps(1,m) + 1
         nyzpd(1,m) = nyzpd(1,m) + 1
      endif
      if ((mz+ks).eq.(nvpz-1)) then
         nyzps(2,m) = nyzps(2,m) + 1
         nyzpd(2,m) = nyzpd(2,m) + 1
      endif
   10 continue
   20 continue
c main loop over decompositions
      do 750 nn = 1, ne
      if (isign.le.0) then
         n = 3 - nn
      else
         n = nn
      endif
c restore previous parameters
      if (nn.eq.3) then
         do 40 mz = 1, nblok
         moff = mblok*(mz - 1)
         do 30 my = 1, mblok
         m = my + moff
         noffd(2,m) = noff(2,m)
         nyzpd(2,m) = nyzp(2,m)
c extend partition to include (nz+1) grids
         if ((mz+ks).eq.(nvpz-1)) then
            nyzpd(2,m) = nyzpd(2,m) + 1
         endif
   30    continue
   40    continue
         n = 2
      endif
      iter = 2
      if (n.eq.1) then
         nyzpmn = nypmx
         mnter = mter
      else if (n.eq.2) then
         nyzpmn = nzpmx
         mnter = nter
      endif
      if (nn.eq.3) mnter = 0
c exit if certain flags are set
      if (mnter.lt.0) go to 750
c determine number of outgoing grids
   50 do 70 m = 1, mnblok
      kl = noffd(n,m)
      kr = kl + nyzpd(n,m)
      jsl(1,m) = 0
      jsr(1,m) = 0
      do 60 l = 1, nyzps(n,m)
      kk = l + noffs(n,m)
c fields going up
      if (kk.gt.kr) then
         jsr(1,m) = jsr(1,m) + 1
c fields going down
      else if (kk.le.kl) then
         jsl(1,m) = jsl(1,m) + 1
      endif
   60 continue
   70 continue
c copy fields
      iter = iter + 2
      npr = 0
      nnter = 0
c get fields from below
      do 270 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 260 my = 1, mblok
      m = my + moff
      kl = my + js + nvpy*(mz + ks) + 1
      if (n.eq.1) then
         kr = kl + 1
         kl = kl - 1
         nxvyz = nxv*nyzps(2,m)
      else if (n.eq.2) then
         kr = kl + nvpy
         kl = kl - nvpy
         nxvyz = nxv*nyzps(1,m)
      endif
      jsl(2,m) = 0
      jsr(2,m) = 0
c this segment is used for shared memory computers
c     if (noffs(n,m).gt.noffd(n,m)) then     
c        jsl(2,m) = jsr(1,kl)
c        if (n.eq.1) then
c           do 100 l = 1, nyzps(2,m)
c           koff = nypmx*(l - 1)
c           do 90 k = 1, jsl(2,m)
c           kk = k + koff
c           do 80 j = 1, nxv
c           g(j,kk,m) = f(j,k+nyzps(n,kl)-jsr(1,kl),l,kl)
c  80       continue
c  90       continue
c 100       continue
c        else if (n.eq.2) then
c           do 130 l = 1, jsl(2,m)
c           koff = nypmx*(l - 1)
c           do 120 k = 1, nyzps(1,m)
c           kk = k + koff
c           do 110 j = 1, nxv
c           g(j,kk,m) = f(j,k,l+nyzps(n,kl)-jsr(1,kl),kl)
c 110       continue
c 120       continue
c 130       continue
c        endif
c     endif
c this segment is used for mpi computers
c post receive from left
      if (noffs(n,m).gt.noffd(n,m)) then  
         call MPI_IRECV(h,nbsize,mreal,kl-1,iter-1,lgrp,msid,ierr)
      endif
c send fields to right
      if (jsr(1,m).gt.0) then
         if (n.eq.1) then
            do 160 l = 1, nyzps(2,m)
            koff = jsr(1,m)*(l - 1)
            do 150 k = 1, jsr(1,m)
            kk = k + koff
            do 140 j = 1, nxv
            g(j,kk,m) = f(j,k+nyzps(n,m)-jsr(1,m),l,m)
  140       continue
  150       continue
  160       continue
         else if (n.eq.2) then
            do 190 l = 1, jsr(1,m)
            koff = nyzps(1,m)*(l - 1)
            do 180 k = 1, nyzps(1,m)
            kk = k + koff
            do 170 j = 1, nxv
            g(j,kk,m) = f(j,k,l+nyzps(n,m)-jsr(1,m),m)
  170       continue
  180       continue
  190       continue
         endif
         call MPI_SEND(g,nxvyz*jsr(1,m),mreal,kr-1,iter-1,lgrp,ierr)
      endif
c wait for fields to arrive
      if (noffs(n,m).gt.noffd(n,m)) then 
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         jsl(2,m) = nps/nxvyz
c shift received data
         if (n.eq.1) then
            do 220 l = 1, nyzps(2,m)
            koff = nypmx*(l - 1)
            loff = jsl(2,m)*(l - 1)
            do 210 k = 1, jsl(2,m)
            kk = k + koff
            ll = k + loff
            do 200 j = 1, nxv
            g(j,kk,m) = h(j,ll,m)
  200       continue
  210       continue
  220       continue
         else if (n.eq.2) then
            do 250 l = 1, jsl(2,m)
            koff = nypmx*(l - 1)
            loff = nyzps(1,m)*(l - 1)
            do 240 k = 1, nyzps(1,m)
            kk = k + koff
            ll = k + loff
            do 230 j = 1, nxv
            g(j,kk,m) = h(j,ll,m)
  230       continue
  240       continue
  250       continue
         endif
      endif
  260 continue
  270 continue
c adjust field
      do 430 m = 1, mnblok
c adjust field size
      nyzps(n,m) = nyzps(n,m) - jsr(1,m)
c do not allow move to overflow field array
      jsr(1,m) = max0((nyzps(n,m)+jsl(2,m)-nyzpmn),0)
      nyzps(n,m) = nyzps(n,m) - jsr(1,m)
      if (jsr(1,m).gt.0) then
         npr = max0(npr,jsr(1,m))
c save whatever is possible into end of g
         kk = min0(jsr(1,m),nyzpmn-jsl(2,m))
         if (n.eq.1) then
            do 300 l = 1, nyzps(2,m)
            koff = nyzpmn - kk + nypmx*(l - 1)
            do 290 k = 1, kk
            do 280 j = 1, nxv
            g(j,k+koff,m) = f(j,nyzps(n,m)+k,l,m)
  280       continue
  290       continue
  300       continue
         else if (n.eq.2) then
            do 330 l = 1, kk
            koff = nypmx*(nyzpmn - kk + l - 1)
            do 320 k = 1, nyzps(1,m)
            do 310 j = 1, nxv
            g(j,k+koff,m) = f(j,k,nyzps(n,m)+l,m)
  310       continue
  320       continue
  330       continue
         endif
      endif
c shift data which is staying, if necessary
      if ((nyzps(n,m).gt.0).and.(jsl(2,m).gt.0)) then
         if (n.eq.1) then
            do 360 l = 1, nyzps(2,m)
            do 350 k = 1, nyzps(n,m)
            kk = nyzps(n,m) - k + 1
            do 340 j = 1, nxv
            f(j,kk+jsl(2,m),l,m) = f(j,kk,l,m)
  340       continue
  350       continue
  360       continue
         else if (n.eq.2) then
            do 390 l = 1, nyzps(n,m)
            ll = nyzps(n,m) - l + 1
            do 380 k = 1, nyzps(1,m)
            do 370 j = 1, nxv
            f(j,k,ll+jsl(2,m),m) = f(j,k,ll,m)
  370       continue
  380       continue
  390       continue
         endif
      endif
c insert data coming from left
      if (n.eq.1) then
         nypm = jsl(2,m)
         nzpm = nyzps(2,m)
      else if (n.eq.2) then
         nypm = nyzps(1,m)
         nzpm = jsl(2,m)
      endif
      do 420 l = 1, nzpm
      koff = nypmx*(l - 1)
      do 410 k = 1, nypm
      kk = k + koff
      do 400 j = 1, nxv
      f(j,k,l,m) = g(j,kk,m)
  400 continue
  410 continue
  420 continue
c adjust field size and offset
      nyzps(n,m) = nyzps(n,m) + jsl(2,m)
      noffs(n,m) = noffs(n,m) - jsl(2,m)
  430 continue
c get fields from above
      do 600 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 590 my = 1, mblok
      m = my + moff
      kl = my + js + nvpy*(mz + ks) + 1
      if (n.eq.1) then
         kr = kl + 1
         kl = kl - 1
         nxvyz = nxv*nyzps(2,m)
      else if (n.eq.2) then
         kr = kl + nvpy
         kl = kl - nvpy
         nxvyz = nxv*nyzps(1,m)
      endif
c this segment is used for shared memory computers
c     if ((noffs(n,m)+nyzps(n,m)).lt.(noffd(n,m)+nyzpd(n,m))) then 
c        jsr(2,m) = jsl(1,kr)
c        if (n.eq.1) then
c           nypm = jsr(2,m)
c           nzpm = nyzps(2,m)
c        else if (n.eq.2) then
c           nypm = nyzps(1,m)
c           nzpm = jsr(2,m)
c        endif
c        do 460 l = 1, nzpm
c        koff = nypmx*(l - 1)
c        do 450 k = 1, nypm
c        kk = k + koff
c        do 440 j = 1, nxv
c        g(j,kk,m) =  f(j,k,l,kr)
c 440    continue
c 450    continue
c 460    continue
c     endif
c this segment is used for mpi computers
c post receive from right
      if ((noffs(n,m)+nyzps(n,m)).lt.(noffd(n,m)+nyzpd(n,m))) then   
         call MPI_IRECV(h,nbsize,mreal,kr-1,iter,lgrp,msid,ierr)
      endif
c send fields to left
      if (jsl(1,m).gt.0) then
         if (n.eq.1) then
            do 490 l = 1, nyzps(2,m)
            koff = jsl(1,m)*(l - 1)
            do 480 k = 1, jsl(1,m)
            kk = k + koff
            do 470 j = 1, nxv
            g(j,kk,m) = f(j,k,l,m)
  470       continue
  480       continue
  490       continue
         else if (n.eq.2) then
            do 520 l = 1, jsl(1,m)
            koff = nyzps(1,m)*(l - 1)
            do 510 k = 1, nyzps(1,m)
            kk = k + koff
            do 500 j = 1, nxv
            g(j,kk,m) = f(j,k,l,m)
  500       continue
  510       continue
  520       continue
         endif
         call MPI_SEND(g,nxvyz*jsl(1,m),mreal,kl-1,iter,lgrp,ierr)
      endif
c wait for fields to arrive
      if ((noffs(n,m)+nyzps(n,m)).lt.(noffd(n,m)+nyzpd(n,m))) then  
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nps,ierr)
         jsr(2,m) = nps/nxvyz
c shift received data
         if (n.eq.1) then
            do 550 l = 1, nyzps(2,m)
            koff = nypmx*(l - 1)
            loff = jsr(2,m)*(l - 1)
            do 540 k = 1, jsr(2,m)
            kk = k + koff
            ll = k + loff 
            do 530 j = 1, nxv
            g(j,kk,m) = h(j,ll,m)
  530       continue
  540       continue
  550       continue
         else if (n.eq.2) then
            do 580 l = 1, jsr(2,m)
            koff = nypmx*(l - 1)
            loff = nyzps(1,m)*(l - 1)
            do 570 k = 1, nyzps(1,m)
            kk = k + koff
            ll = k + loff 
            do 560 j = 1, nxv
            g(j,kk,m) = h(j,ll,m)
  560       continue
  570       continue
  580       continue
         endif
      endif
  590 continue
  600 continue
c adjust field
      do 740 m = 1, mnblok
c adjust field size
      nyzps(n,m) = nyzps(n,m) - jsl(1,m)
      noffs(n,m) = noffs(n,m) + jsl(1,m)
c shift data which is staying, if necessary
      if ((nyzps(n,m).gt.0).and.(jsl(1,m).gt.0)) then
        if (n.eq.1) then
            do 630 l = 1, nyzps(2,m)
            do 620 k = 1, nyzps(n,m)
            do 610 j = 1, nxv
            f(j,k,l,m) = f(j,k+jsl(1,m),l,m)
  610       continue
  620       continue
  630       continue
         else if (n.eq.2) then
            do 660 l = 1, nyzps(n,m)
            do 650 k = 1, nyzps(1,m)
            do 640 j = 1, nxv
            f(j,k,l,m) = f(j,k,l+jsl(1,m),m)
  640       continue
  650       continue
  660       continue
         endif
      endif
c do not allow move to overflow field array
      jsl(1,m) = max0((nyzps(n,m)+jsr(2,m)-nyzpmn),0)
      if (jsl(1,m).gt.0) then
         npr = max0(npr,jsl(1,m))
         jsr(2,m) = jsr(2,m) - jsl(1,m)
c do not process if prior error
      else if (jsr(1,m).gt.0) then
         go to 730
      endif
c insert data coming from right
      if (n.eq.1) then
         do 690 l = 1, nyzps(2,m)
         koff = nypmx*(l - 1)
         do 680 k = 1, jsr(2,m)
         kk = k + koff
         do 670 j = 1, nxv
         f(j,k+nyzps(n,m),l,m) = g(j,kk,m)
  670    continue
  680    continue
  690    continue
      else if (n.eq.2) then
         do 720 l = 1, jsr(2,m)
         koff = nypmx*(l - 1)
         do 710 k = 1, nyzps(1,m)
         kk = k + koff
         do 700 j = 1, nxv
         f(j,k,l+nyzps(n,m),m) = g(j,kk,m)
  700    continue
  710    continue
  720    continue
      endif
c adjust field size and offset
      nyzps(n,m) = nyzps(n,m) + jsr(2,m)
c check if new partition is uniform
  730 nnter = nnter + abs(nyzps(n,m)-nyzpd(n,m)) + abs(noffs(n,m)-noffd(
     1n,m))
  740 continue
c calculate number of iterations
      nps = iter/2 - 1
      if (nps.le.mnter) then
c process errors
         if (npr.ne.0) then
            ierr = npr
            write (2,*) 'local field overflow error, ierr = ', ierr
            go to 760
         endif
         if (nps.lt.mnter) go to 50
         go to 750
      endif
c process errors
      ibflg(1) = npr
      ibflg(2) = nnter
      call PIMAX(ibflg,iwork,2,1)
c field overflow error
      if (ibflg(1).ne.0) then
         ierr = ibflg(1)
         write (2,*) 'global field overflow error, ierr = ', ierr
         go to 760
      endif
c check if any fields have to be passed further
      if (ibflg(2).gt.0) then
         write (2,*) 'Info: fields being passed further = ', ibflg(2)
         go to 50
      endif
      mnter = nps
      if (n.eq.1) then
         mter = mnter
      else if (n.eq.2) then
         if (nn.lt.3) nter = mnter
      endif
  750 continue
c restore partitions to normal
  760 do 780 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 770 my = 1, mblok
      m = my + moff
      if ((my+js).eq.(nvpy-1)) then
         nyzps(1,m) = nyzps(1,m) - 1
         nyzpd(1,m) = nyzpd(1,m) - 1
      endif
      if ((mz+ks).eq.(nvpz-1)) then
         nyzps(2,m) = nyzps(2,m) - 1
         nyzpd(2,m) = nyzpd(2,m) - 1
      endif
  770 continue
  780 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine REPARTD32(edges,edg,eds,eg,es,et2,npicyz,noff,nyzp,anpa
     1v,nypmin,nypmax,nzpmin,nzpmax,kstrt,nvpy,nvpz,mblok,nblok,idps,idd
     2s,myzpm1)
c this subroutines finds new partitions boundaries (edges,noff,nyzp)
c from old partition information (npic,nyzp).
c edges(1,m) = lower boundary in y of particle partition m
c edges(2,m) = upper boundary in y of particle partition m
c edges(3,m) = back boundary in z of particle partition m
c edges(4,m) = front boundary in z of particle partition m
c edg/eds/eg/es/et2 = scratch arrays
c npicyz(m) = number of particles per grid in y and z in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c anpav = average number of particles per partition desired
c nypmin/nypmax = minimum/maximum value of nyzp(1,m) in new partition
c nzpmin/nzpmax = minimum/maximum value of nyzp(2,m) in new partition
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c mblok/nblok = number of field partitions in y/z.
c idps = number of partition boundaries
c idds = dimensionality of domain decompositions
c myzpm1 = maximum size of particle partition in either direction
      implicit none
      real edges, edg, eds, eg, es, et2
      integer npicyz, noff, nyzp
      real anpav
      integer nypmin, nypmax, nzpmin, nzpmax, kstrt, nvpy, nvpz
      integer mblok, nblok, idps, idds, myzpm1
      dimension edges(idps,mblok*nblok)
      dimension edg(myzpm1,mblok*nblok), eds(myzpm1,mblok*nblok)
      dimension eg(idds,mblok*nblok), es(idds,mblok*nblok)
      dimension et2(2*idds,mblok*nblok)
      dimension npicyz(myzpm1,idds,mblok*nblok)
      dimension nyzp(idds,mblok*nblok), noff(idds,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, iter, nter, nyzp1, npav, k1, kb, kl, kr, k, my, mz
      integer m, n, nvp, moff, mnblok, ierr
      real sum1, at1, at2, anpv, apav, anpl, anpr
      integer msid, istatus
      integer ibflg, iwork
      dimension istatus(lstat)
      dimension ibflg(4), iwork(4)
      dimension kb(2)
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      mnblok = mblok*nblok
      iter = 2
c main loop over decompositions
      do 260 n = 1, 2
      if (n.eq.1) then
         nvp = nvpy
         anpv = anpav*real(nvpz)
      elseif (n.eq.2) then
         nvp = nvpz
         anpv = anpav*real(nvpy)
      endif
c copy number of particles and grid in current partition
      do 20 m = 1, mnblok
      sum1 = 0.
      do 10 k = 1, nyzp(n,m)
      at1 = npicyz(k,n,m)
      sum1 = sum1 + at1
      eds(k,m) = at1
   10 continue
      edges(2*n-1,m) = sum1
      edges(2*n,m) = nyzp(n,m)
      et2(1,m) = edges(2*n-1,m)
      et2(2,m) = edges(2*n,m)
   20 continue
c perform running sum
      call PSCAN2(edges(2*n-1,1),eg,es,2,kstrt,nvpy,nvpz,n,mblok,nblok)
      do 30 m = 1, mnblok
      es(1,m) = et2(1,m)
      es(2,m) = et2(2,m)
      et2(1,m) = edges(2*n-1,m)
      et2(2,m) = edges(2*n,m)
      et2(3,m) = et2(1,m)
      et2(4,m) = et2(2,m)
      eg(1,m) = 0.
      eg(2,m) = 0.
      edges(2*n,m) = 1.0
   30 continue
c move partitions
   40 iter = iter + 2
c get partition from left
      do 70 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 60 my = 1, mblok
      m = my + moff
      kb(1) = my + js
      kb(2) = mz + ks
      kl = kb(1) + nvpy*kb(2) + 1
      if (n.eq.1) then
         kr = kl + 1
         kl = kl - 1
      else if (n.eq.2) then
         kr = kl + nvpy
         kl = kl - nvpy
      endif
c apav = desired number of particles on processor to left
      npav = real(kb(n))*anpv + .5
      apav = real(npav)
c anpl = deficit of particles on processor to left
      anpl = apav - et2(1,m) + es(1,m)
c anpr = excess of particles on current processor
      npav = real(kb(n)+1)*anpv + .5
      anpr = et2(1,m) - real(npav)
c this segment is used for shared memory computers
c     if (anpl.lt.0.) then
c        nyzp1 = es(2,kl)
c        do 50 k = 1, nyzp1
c        edg(k,m) = eds(k,kl)
c  50    continue
c        eg(1,m) = es(1,kl)
c        eg(2,m) = es(2,kl)
c     endif
c this segment is used for mpi computers
c post receive from left
      if (anpl.lt.0.) then
         call MPI_IRECV(edg,myzpm1,mreal,kl-1,iter-1,lgrp,msid,ierr)
      endif
c send partition to right
      if (anpr.gt.0.) then
         nyzp1 = es(2,m)
         call MPI_SEND(eds,nyzp1,mreal,kr-1,iter-1,lgrp,ierr)
      endif
c wait for partition to arrive
      if (anpl.lt.0.) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nyzp1,ierr)
         eg(2,m) = nyzp1
         sum1 = 0.
         do 50 k = 1, nyzp1
         sum1 = sum1 + edg(k,m)
   50    continue
         eg(1,m) = sum1
      endif
   60 continue
   70 continue
c find new partitions
      nter = 0
      do 120 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 110 my = 1, mblok
      m = my + moff
      kb(1) = my + js
      kb(2) = mz + ks
      kl = kb(1) + nvpy*kb(2) + 1
      if (n.eq.1) then
         kl = kl - 1
      else if (n.eq.2) then
         kl = kl - nvpy
      endif
      npav = real(kb(n))*anpv + .5
      apav = real(npav)
      anpl = apav - et2(1,m) + es(1,m)
      npav = real(kb(n)+1)*anpv + .5
      anpr = et2(1,m) - real(npav)
c left boundary is on the left
      if (anpl.lt.0.) then
         if ((anpl+eg(1,m)).ge.0.) then
            nyzp1 = eg(2,m)
            k1 = nyzp1
            sum1 = 0.
   80       at1 = sum1
            sum1 = sum1 - edg(k1,m)
            k1 = k1 - 1
            if ((sum1.gt.anpl).and.(k1.gt.0)) go to 80
            at1 = real(nyzp1 - k1 - 1) + (anpl - at1)/(sum1 - at1)
            edges(2*n-1,m) = (et2(2,m) - es(2,m)) - at1
c left boundary is even further to left
         else
            nter = nter + 1
         endif
c left boundary is inside
      else if (et2(1,m).ge.apav) then
         nyzp1 = es(2,m)
         k1 = 1
         sum1 = 0.
   90    at1 = sum1
         sum1 = sum1 + eds(k1,m)
         k1 = k1 + 1
         if ((sum1.lt.anpl).and.(k1.le.nyzp1)) go to 90
         at2 = real(k1 - 2)
         if (sum1.gt.at1) at2 = at2 + (anpl - at1)/(sum1 - at1)
         edges(2*n-1,m) = (et2(2,m) - es(2,m)) + at2
      endif
c right going data will need to be sent
      if (anpr.gt.es(1,m)) nter = nter + 1
      if (kb(n).gt.0) then
         nyzp1 = eg(2,m)
         do 100 k = 1, nyzp1
         eds(k,m) = edg(k,m)
  100    continue
         et2(1,m) = et2(1,m) - es(1,m)
         et2(2,m) = et2(2,m) - es(2,m)
         es(1,m) = eg(1,m)
         es(2,m) = eg(2,m)
      endif
  110 continue
  120 continue
c get more data from left
      if (nter.gt.0) go to 40
      iter = nvp + 2
c restore partition data
      do 140 m = 1, mnblok
      sum1 = 0.
      do 130 k = 1, nyzp(n,m)
      at1 = npicyz(k,n,m)
      sum1 = sum1 + at1
      eds(k,m) = at1
  130 continue
      et2(1,m) = et2(3,m)
      et2(2,m) = et2(4,m)
      es(1,m) = sum1
      es(2,m) = nyzp(n,m)
      eg(1,m) = 0.
      eg(2,m) = 0.
  140 continue
c continue moving partitions
  150 iter = iter + 2
c get partition from right
      do 180 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 170 my = 1, mblok
      m = my + moff
      kb(1) = my + js
      kb(2) = mz + ks
      kl = kb(1) + nvpy*kb(2) + 1
      if (n.eq.1) then
         kr = kl + 1
         kl = kl - 1
      else if (n.eq.2) then
         kr = kl + nvpy
         kl = kl - nvpy
      endif
      npav = real(kb(n))*anpv + .5
      apav = real(npav)
      npav = real(kb(n)-1)*anpv + .5
      anpl = real(npav) - et2(1,m) + es(1,m)
c this segment is used for shared memory computers
c     if (et2(1,m).lt.apav) then
c        nyzp1 = es(2,kr)
c        do 160 k = 1, nyzp1
c        edg(k,m) = eds(k,kr)
c 160    continue
c        eg(1,m) = es(1,kr)
c        eg(2,m) = es(2,kr)
c     endif
c this segment is used for mpi computers
c post receive from right
      if (et2(1,m).lt.apav) then
         call MPI_IRECV(edg,myzpm1,mreal,kr-1,iter,lgrp,msid,ierr)
      endif
c send partition to left
      if (anpl.gt.0.) then
         nyzp1 = es(2,m)
         call MPI_SEND(eds,nyzp1,mreal,kl-1,iter,lgrp,ierr)
      endif
c wait for partition to arrive
      if (et2(1,m).lt.apav) then
         call MPI_WAIT(msid,istatus,ierr)
         call MPI_GET_COUNT(istatus,mreal,nyzp1,ierr)
         eg(2,m) = nyzp1
         sum1 = 0.
         do 160 k = 1, nyzp1
         sum1 = sum1 + edg(k,m)
  160    continue
         eg(1,m) = sum1
      endif
  170 continue
  180 continue
c find new partitions
      nter = 0
      do 220 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 210 my = 1, mblok
      m = my + moff
      kb(1) = my + js
      kb(2) = mz + ks
      kl = kb(1) + nvpy*kb(2) + 1
      if (n.eq.1) then
         kr = kl + 1
         kl = kl - 1
      else if (n.eq.2) then
         kr = kl + nvpy
         kl = kl - nvpy
      endif
      npav = real(kb(n))*anpv + .5
      apav = real(npav)
      npav = real(kb(n)-1)*anpv + .5
      anpl = real(npav) - et2(1,m) + es(1,m)
      anpr = et2(1,m) - apav
c left boundary is on the right
      if (et2(1,m).lt.apav) then
         if ((et2(1,m)+eg(1,m)).ge.apav) then
            nyzp1 = eg(2,m)
            k1 = 1
            sum1 = 0.
            at2 = -anpr
  190       at1 = sum1
            sum1 = sum1 + edg(k1,m)
            k1 = k1 + 1
            if ((sum1.lt.at2).and.(k1.le.nyzp1)) go to 190
            at1 = real(k1 - 2) + (at2 - at1)/(sum1 - at1)
            edges(2*n-1,m) = et2(2,m) + at1
c left boundary is even further to right
         else
            nter = nter + 1
         endif
      endif
c left going data will need to be sent
      if (anpl.gt.es(1,m)) nter = nter + 1
      if ((kb(n)+2).le.nvp) then
         nyzp1 = eg(2,m)
         do 200 k = 1, nyzp1
         eds(k,m) = edg(k,m)
  200    continue
         et2(1,m) = et2(1,m) + eg(1,m)
         et2(2,m) = et2(2,m) + eg(2,m)
         es(1,m) = eg(1,m)
         es(2,m) = eg(2,m)
      endif
  210 continue
  220 continue
c get more data from right
      if (nter.gt.0) go to 150
c send left edge to processor on right
      iter = 2
      do 240 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 230 my = 1, mblok
      m = my + moff
      kb(1) = my + js
      kb(2) = mz + ks
      kl = kb(1) + nvpy*kb(2) + 1
      if (n.eq.1) then
         kr = kl + 1
         kl = kl - 1
      else if (n.eq.2) then
         kr = kl + nvpy
         kl = kl - nvpy
      endif
c this segment is used for shared memory computers
c     if ((kb(n)+2).le.nvp) then
c        edges(2*n,m) = edges(2*n-1,kr)
c     else
c        edges(2*n,m) = et2(4,m)
c     endif
c this segment is used for mpi computers
c post receive from right
      if ((kb(n)+2).le.nvp) then
         call MPI_IRECV(edges(2*n,m),1,mreal,kr-1,iter,lgrp,msid,ierr)
      endif
c send left edge to left
      if (kb(n).gt.0) then
         call MPI_SEND(edges(2*n-1,m),1,mreal,kl-1,iter,lgrp,ierr)
      endif
c wait for edge to arrive
      if ((kb(n)+2).le.nvp) then
         call MPI_WAIT(msid,istatus,ierr)
      else
         edges(2*n,m) = et2(4,m)
      endif
  230 continue
  240 continue
c calculate number of grids and offsets in new partitions
      do 250 m = 1, mnblok
      kl = edges(2*n-1,m) + .5
      noff(n,m) = kl
      kr = edges(2*n,m) + .5
      nyzp(n,m) = kr - kl
      edges(2*n-1,m) = real(kl)
      edges(2*n,m) = real(kr)
  250 continue
  260 continue
c find minimum and maximum partition size
      nypmin = nyzp(1,1)
      nypmax = nyzp(1,1)
      nzpmin = nyzp(2,1)
      nzpmax = nyzp(2,1)
      do 270 m = 1, mnblok
      nypmin = min0(nypmin,nyzp(1,m))
      nypmax = max0(nypmax,nyzp(1,m))
      nzpmin = min0(nzpmin,nyzp(2,m))
      nzpmax = max0(nzpmax,nyzp(2,m))
  270 continue
      ibflg(1) = -nypmin
      ibflg(2) = nypmax
      ibflg(3) = -nzpmin
      ibflg(4) = nzpmax
      call PIMAX(ibflg,iwork,4,1)
      nypmin = -ibflg(1)
      nypmax = ibflg(2)
      nzpmin = -ibflg(3)
      nzpmax = ibflg(4)
      return
      end
c-----------------------------------------------------------------------
      subroutine FNOFF32(edges,noff,nyzp,nypmin,nypmax,nzpmin,nzpmax,mnb
     1lok,idps,idds)
c this subroutines finds new partitions arrays (noff,nzyp) from edges
c edges(1,m) = lower boundary in y of particle partition m
c edges(2,m) = upper boundary in y of particle partition m
c edges(3,m) = back boundary in z of particle partition m
c edges(4,m) = front boundary in z of particle partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c nypmin/nypmax = minimum/maximum value of nyzp(1,m) in new partition
c nzpmin/nzpmax = minimum/maximum value of nyzp(2,m) in new partition
c mnblok = number of field partitions.
c idps = number of partition boundaries
c idds = dimensionality of domain decompositions
      implicit none
      real edges
      integer noff, nyzp, nypmin, nypmax, nzpmin, nzpmax, mnblok, idps
      integer idds
      dimension edges(idps,mnblok)
      dimension nyzp(idds,mnblok), noff(idds,mnblok)
c local data
      integer kl, kr, m, n
      integer ibflg, iwork
      dimension ibflg(4), iwork(4)
c calculate number of grids and offsets in new partitions
      do 20 n = 1, 2
      do 10 m = 1, mnblok
      kl = edges(2*n-1,m) + .5
      noff(n,m) = kl
      kr = edges(2*n,m) + .5
      nyzp(n,m) = kr - kl
   10 continue
   20 continue
c find minimum and maximum partition size
      nypmin = nyzp(1,1)
      nypmax = nyzp(1,1)
      nzpmin = nyzp(2,1)
      nzpmax = nyzp(2,1)
      do 30 m = 1, mnblok
      nypmin = min0(nypmin,nyzp(1,m))
      nypmax = max0(nypmax,nyzp(1,m))
      nzpmin = min0(nzpmin,nyzp(2,m))
      nzpmax = max0(nzpmax,nyzp(2,m))
   30 continue
      ibflg(1) = -nypmin
      ibflg(2) = nypmax
      ibflg(3) = -nzpmin
      ibflg(4) = nzpmax
      call PIMAX(ibflg,iwork,4,1)
      nypmin = -ibflg(1)
      nypmax = ibflg(2)
      nzpmin = -ibflg(3)
      nzpmax = ibflg(4)
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
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, ls, kxb, kyb, kzb
      integer jkblok, kxym, mtr, ntr, mntr
      integer m, mx, mz, l, i, moff, ioff, joff, koff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
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
         call MPI_IRECV(t(1,1,1,m),kxyp*kyp*kzp,mcplx,ir-1,ir+kxym+1,lgr
     1p,msid,ierr)
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
         call MPI_SEND(s(1,1,1,m),kxyp*kyp*kzp,mcplx,is-1,m+kstrt+kxym,l
     1grp,ierr)
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
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, kxb, kyzb, kzb
      integer mlblok, kyzm, mtr, ntr, mntr
      integer m, mx, my, l, i, moff, ioff, koff, loff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
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
         call MPI_IRECV(t(1,1,1,m),kxyp*kyzp*kzp,mcplx,ir-1,ir+kyzm+1,lg
     1rp,msid,ierr)
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
         call MPI_SEND(s(1,1,1,m),kxyp*kyzp*kzp,mcplx,is-1,m+kstrt+kyzm,
     1lgrp,ierr)
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
c g(1:3,k+kyp*(m-1),j,l,mx,mz) = f(1:3,j+kxyp*(l-1),k,l,my,mz), where
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
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, ls, kxb, kyb, kzb
      integer jkblok, kxym, mtr, ntr, mntr
      integer m, mx, mz, l, i, moff, ioff, joff, koff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
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
         call MPI_IRECV(t(1,1,1,1,m),3*kxyp*kyp*kzp,mcplx,ir-1,ir+kxym+1
     1,lgrp,msid,ierr)
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
         call MPI_SEND(s(1,1,1,1,m),3*kxyp*kyp*kzp,mcplx,is-1,m+kstrt+kx
     1ym,lgrp,ierr)
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
c h(1:3,j+kzp*(l-1),k,l,my,mz) = g(1:3,k+kyzp*(m-1),j,l,mx,mz), where
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
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, kxb, kyzb, kzb
      integer mlblok, kyzm, mtr, ntr, mntr
      integer m, mx, my, l, i, moff, ioff, koff, loff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
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
         call MPI_IRECV(t(1,1,1,1,m),3*kxyp*kyzp*kzp,mcplx,ir-1,ir+kyzm+
     11,lgrp,msid,ierr)
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
         call MPI_SEND(s(1,1,1,1,m),3*kxyp*kyzp*kzp,mcplx,is-1,m+kstrt+k
     1yzm,lgrp,ierr)
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
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, ls, kxb, kyb, kzb
      integer jkblok, kxym, mtr, ntr, mntr
      integer m, mx, mz, l, i, moff, ioff, joff, koff, loff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
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
         call MPI_IRECV(f(1+kxyp*kyp*kzp*is),kxyp*kyp*kzp,mcplx,ir-1,ir+
     1kxym+1,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kyb*kzb)).and.(ii.le.ntr)) then
         ir = is0 + kxym*(ii - 1)
         is = ir + kxb*(mz + ks) + 1
         call MPI_SEND(g(1+kxyp*kyp*kzp*ir),kxyp*kyp*kzp,mcplx,is-1,m+ks
     1trt+kxym,lgrp,ierr)
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
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, kxb, kyzb, kzb
      integer mlblok, kyzm, mtr, ntr, mntr, nzyv
      integer m, mx, my, l, i, moff, ioff, joff, koff, loff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
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
         call MPI_IRECV(g(1+kxyp*kyzp*kzp*is),kxyp*kyzp*kzp,mcplx,ir-1,i
     1r+kyzm+1,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.ntr)) then
         ir = is0 + kyzm*(ii - 1)
         is = mx + js + kxb*ir + 1
         call MPI_SEND(h(1+kxyp*kyzp*kzp*ir),kxyp*kyzp*kzp,mcplx,is-1,m+
     1kstrt+kyzm,lgrp,ierr)
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
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, ls, kxb, kyb, kzb
      integer jkblok, kxym, mtr, ntr, mntr
      integer m, mx, mz, l, i, moff, ioff, joff, koff, loff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
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
         call MPI_IRECV(f(1+3*kxyp*kyp*kzp*is),3*kxyp*kyp*kzp,mcplx,ir-1
     1,ir+kxym+1,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kyb*kzb)).and.(ii.le.ntr)) then
         ir = is0 + kxym*(ii - 1)
         is = ir + kxb*(mz + ks) + 1
         call MPI_SEND(g(1+3*kxyp*kyp*kzp*ir),3*kxyp*kyp*kzp,mcplx,is-1,
     1m+kstrt+kxym,lgrp,ierr)
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
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, kxb, kyzb, kzb
      integer mlblok, kyzm, mtr, ntr, mntr, nzyv
      integer m, mx, my, l, i, moff, ioff, joff, koff, loff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
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
         call MPI_IRECV(g(1+3*kxyp*kyzp*kzp*is),3*kxyp*kyzp*kzp,mcplx,ir
     1-1,ir+kyzm+1,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.ntr)) then
         ir = is0 + kyzm*(ii - 1)
         is = mx + js + kxb*ir + 1
         call MPI_SEND(h(1+3*kxyp*kyzp*kzp*ir),3*kxyp*kyzp*kzp,mcplx,is-
     11,m+kstrt+kyzm,lgrp,ierr)
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
      subroutine PNTPOS3A(f,g,s,t,nx,ny,nz,kstrt,nxv,nyv,kxyp,kyp,kzp,kx
     1ypd,kypd,kzpd,jblok,kblok,lblok,ndim)
c this subroutine performs a transpose of a matrix f, distributed in y
c and z to a matrix g, distributed in x and z, that is,
c g(1:ndim,k+kyp*(m-1),j,l,mx,mz) = f(1:ndim,j+kxyp*(l-1),k,l,my,mz),
c where, 1 <= j <= kxyp, 1 <= k <= kyp, 1 <= l <= kzp, and
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
c ndim = leading dimension of arrays f and g
      implicit none
      integer nx, ny, nz, kstrt, nxv, nyv, kxyp, kyp, kzp
      integer kxypd, kypd, kzpd, jblok, kblok, lblok, ndim
      complex f, g, s, t
      dimension f(ndim,nxv,kypd,kzpd,kblok*lblok)
      dimension g(ndim,nyv,kxypd,kzpd,jblok*lblok)
      dimension s(ndim,kxyp,kyp,kzp,kblok*lblok)
      dimension t(ndim,kxyp,kyp,kzp,jblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, ls, kxb, kyb, kzb
      integer jkblok, kxym, mtr, ntr, mntr
      integer m, mx, mz, l, i, moff, ioff, joff, koff, k, j, n
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyb = ny/kyp
      kzb = nz/kzp
      ks = (kstrt - 1)/kyb
      js = kstrt - kyb*ks - 2
      ks = ks - 1
      ls = (kstrt - 1)/kxb - 1
c this segment is used for shared memory computers
c     do 70 mz = 1, lblok
c     moff = jblok*(mz - 1)
c     ioff = kblok*(mz - 1)
c     do 60 mx = 1, jblok
c     joff = kxyp*(mx + js)
c     m = mx + moff
c     do 50 i = 1, kyb
c     koff = kyp*(i - 1)
c     do 40 l = 1, kzp
c     do 30 k = 1, kyp
c     do 20 j = 1, kxyp
c     do 10 n = 1, ndim
c     g(n,k+koff,j,l,m) = f(n,j+joff,k,l,i+ioff)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c  70 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      do 120 mz = 1, lblok
      moff = jkblok*(mz - 1)
      do 110 mx = 1, jkblok
      m = mx + moff
      do 100 i = 1, kxym
      ir0 = iand(kxym-1,ieor(mx+js,i-1))
      is0 = ir0
      do 90 ii = 1, mntr
c post receive
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         koff = kyp*ir
         ir = ir + kyb*(mz + ls) + 1
         call MPI_IRECV(t(1,1,1,1,m),ndim*kxyp*kyp*kzp,mcplx,ir-1,ir+kxy
     1m+1,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kyb*kzb)).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         joff = kxyp*is
         is = is + kxb*(mz + ks) + 1
         do 40 l = 1, kzp
         do 30 k = 1, kyp
         do 20 j = 1, kxyp
         do 10 n = 1, ndim
         s(n,j,k,l,m) = f(n,j+joff,k,l,m)
   10    continue
   20    continue
   30    continue
   40    continue
         call MPI_SEND(s(1,1,1,1,m),ndim*kxyp*kyp*kzp,mcplx,is-1,m+kstrt
     1+kxym,lgrp,ierr)
      endif
c receive data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
         do 80 l = 1, kzp
         do 70 k = 1, kyp
         do 60 j = 1, kxyp
         do 50 n = 1, ndim
         g(n,k+koff,j,l,m) = t(n,j,k,l,m)
   50    continue
   60    continue
   70    continue
   80    continue
      endif
   90 continue
  100 continue
  110 continue
  120 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNTPOS3B(g,h,s,t,nx,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,k
     1xypd,kyzpd,kzpd,jblok,mblok,lblok,ndim)
c this subroutine performs a transpose of a matrix g, distributed in x
c and z to a matrix h, distributed in x and y, that is,
c h(1:ndim,j+kzp*(l-1),k,l,my,mz) = g(1:ndim,k+kyzp*(m-1),j,l,mx,mz),
c where, 1 <= j <= kxyp, 1 <= k <= kyzp, 1 <= l <= kzp, and
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
c ndim = leading dimension of arrays g and h
      implicit none
      integer nx, ny, nz, kstrt, nyv, nzv, kxyp, kyzp, kzp
      integer kxypd, kyzpd, kzpd, jblok, mblok, lblok, ndim
      complex g, h, s, t
      dimension g(ndim,nyv,kxypd,kzpd,jblok*lblok)
      dimension h(ndim,nzv,kxypd,kyzpd,jblok*mblok)
      dimension s(ndim,kyzp,kxyp,kzp,jblok*lblok)
      dimension t(ndim,kyzp,kxyp,kzp,jblok*mblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, kxb, kyzb, kzb
      integer mlblok, kyzm, mtr, ntr, mntr
      integer m, mx, my, l, i, moff, ioff, koff, loff, k, j, n
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyzb = ny/kyzp
      kzb = nz/kzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
c this segment is used for shared memory computers
c     do 70 my = 1, mblok
c     moff = jblok*(my - 1)
c     koff = kyzp*(my + ks)
c     do 60 mx = 1, jblok
c     m = mx + moff
c     do 50 i = 1, kzb
c     loff = kzp*(i - 1)
c     ioff = jblok*(i - 1)
c     do 40 l = 1, kzp
c     do 30 j = 1, kxyp
c     do 20 k = 1, kyzp
c     do 10 n = 1, ndim
c     h(n,l+loff,j,k,m) = g(n,k+koff,j,l,mx+ioff)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c  70 continue
c this segment is used for mpi computers
      mlblok = max0(mblok,lblok)
      kyzm = min0(kyzb,kzb)
      mtr = kzb/kyzm
      ntr = kyzb/kyzm
      mntr = max0(mtr,ntr)
      do 120 my = 1, mlblok
      moff = jblok*(my - 1)
      do 110 mx = 1, jblok
      m = mx + moff
      do 100 i = 1, kyzm
      ir0 = iand(kyzm-1,ieor(my+ks,i-1))
      is0 = ir0
      do 90 ii = 1, mntr
c post receive
      if ((kstrt.le.(kxb*kyzb)).and.(ii.le.mtr)) then
         ir = ir0 + kyzm*(ii - 1)
         loff = kzp*ir
         ir = mx + js + kxb*ir + 1
         call MPI_IRECV(t(1,1,1,1,m),ndim*kxyp*kyzp*kzp,mcplx,ir-1,ir+ky
     1zm+1,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.ntr)) then
         is = is0 + kyzm*(ii - 1)
         koff = kyzp*is
         is = mx + js + kxb*is + 1
         do 40 l = 1, kzp
         do 30 j = 1, kxyp
         do 20 k = 1, kyzp
         do 10 n = 1, ndim
         s(n,k,j,l,m) = g(n,k+koff,j,l,m)
   10    continue
   20    continue
   30    continue
   40    continue
         call MPI_SEND(s(1,1,1,1,m),ndim*kxyp*kyzp*kzp,mcplx,is-1,m+kstr
     1t+kyzm,lgrp,ierr)
      endif
c receive data
      if ((kstrt.le.(kxb*kyzb)).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
         do 80 l = 1, kzp
         do 70 j = 1, kxyp
         do 60 k = 1, kyzp
         do 50 n = 1, ndim
         h(n,l+loff,j,k,m) = t(n,k,j,l,m)
   50    continue
   60    continue
   70    continue
   80    continue
      endif
   90 continue
  100 continue
  110 continue
  120 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNTPOS3AX(f,g,nx,ny,nz,kstrt,nxv,nyv,kxyp,kyp,kzp,kxypd
     1,kypd,kzpd,jblok,kblok,lblok,ndim)
c this subroutine performs a transpose of a matrix f, distributed in y
c and z to a matrix g, distributed in x and z, that is,
c g(1:ndim,k+kyp*(m-1),j,l,mx,mz) = f(1:ndim,j+kxyp*(l-1),k,l,my,mz),
c where, 1 <= j <= kxyp, 1 <= k <= kyp, 1 <= l <= kzp, and
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
c ndim = leading dimension of arrays f and g
c optimized version
      implicit none
      integer nx, ny, nz, kstrt, nxv, nyv, kxyp, kyp, kzp
      integer kxypd, kypd, kzpd, jblok, kblok, lblok, ndim
      complex f, g
      dimension f(ndim*nxv*kypd*kzpd*kblok*lblok)
      dimension g(ndim*nyv*kxypd*kzpd*jblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, ls, kxb, kyb, kzb
      integer jkblok, kxym, mtr, ntr, mntr
      integer m, mx, mz, l, i, moff, ioff, joff, koff, loff, k, j, n
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyb = ny/kyp
      kzb = nz/kzp
      ks = (kstrt - 1)/kyb
      js = kstrt - kyb*ks - 2
      ks = ks - 1
      ls = (kstrt - 1)/kxb - 1
c this segment is used for shared memory computers
c     do 70 mz = 1, lblok
c     moff = jblok*(mz - 1)
c     ioff = kblok*(mz - 1) - 1
c     do 60 mx = 1, jblok
c     joff = kxyp*(mx + js) - 1
c     m = mx + moff
c     do 50 i = 1, kyb
c     koff = kyp*(i - 1) - 1
c     do 40 l = 1, kzp
c     do 30 k = 1, kyp
c     do 20 j = 1, kxyp
c     do 10 n = 1, ndim
c     g(n+ndim*(k+koff+nyv*(j-1+kxypd*(l-1+kzpd*(m-1))))) = f(n+ndim*(j+
c    1joff+nxv*(k-1+kypd*(l-1+kzpd*(i+ioff)))))
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c  70 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
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
         joff = ndim*kxyp*is
         is = kzp*(is + ioff) - 1
         do 30 l = 1, kzp
         ir = kyp*(l + is) - 1
         koff = kypd*(l + loff) - 1
         do 20 k = 1, kyp
         do 10 j = 1, ndim*kxyp
         g(j+ndim*kxyp*(k+ir)) = f(j+joff+ndim*nxv*(k+koff))
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
         call MPI_IRECV(f(1+ndim*kxyp*kyp*kzp*is),ndim*kxyp*kyp*kzp,mcpl
     1x,ir-1,ir+kxym+1,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kyb*kzb)).and.(ii.le.ntr)) then
         ir = is0 + kxym*(ii - 1)
         is = ir + kxb*(mz + ks) + 1
         call MPI_SEND(g(1+ndim*kxyp*kyp*kzp*ir),ndim*kxyp*kyp*kzp,mcplx
     1,is-1,m+kstrt+kxym,lgrp,ierr)
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
      do 190 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 180 mx = 1, jblok
      m = mx + moff
      ioff = kyb*(m - 1) 
      loff = kzpd*(m - 1) - 1
      do 170 i = 1, kxym
      ir0 = iand(kxym-1,ieor(mx+js,i-1))
      do 160 ii = 1, mtr
      if (kstrt.le.(kxb*kzb)) then
         ir = ir0 + kxym*(ii - 1)
         koff = kyp*ir - 1
         ir = kzp*(ir + ioff) - 1
         do 150 l = 1, kzp
         is = kyp*(l + ir) - 1
         joff = kxypd*(l + loff) - 1
         do 140 k = 1, kyp
         do 130 j = 1, kxyp
         do 120 n = 1, ndim
         g(ndim*(k+koff+nyv*(j+joff))+n) = f(ndim*(j+kxyp*(k+is)-1)+n)
  120    continue
  130    continue
  140    continue
  150    continue
      endif
  160 continue
  170 continue
  180 continue
  190 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNTPOS3BX(g,h,nx,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxyp
     1d,kyzpd,kzpd,jblok,mblok,lblok,ndim)
c this subroutine performs a transpose of a matrix g, distributed in x
c and z to a matrix h, distributed in x and y, that is,
c h(1:ndim,j+kzp*(l-1),k,l,my,mz) = g(1:ndim,k+kyzp*(m-1),j,l,mx,mz),
c where, 1 <= j <= kxyp, 1 <= k <= kyzp, 1 <= l <= kzp, and
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
c ndim = leading dimension of arrays g and h
c optimized version
      implicit none
      integer nx, ny, nz, kstrt, nyv, nzv, kxyp, kyzp, kzp
      integer kxypd, kyzpd, kzpd, jblok, mblok, lblok, ndim
      complex g, h
      dimension g(ndim*nyv*kxypd*kzpd*jblok*lblok)
      dimension h(ndim*nzv*kxypd*kyzpd*jblok*mblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, kxb, kyzb, kzb
      integer mlblok, kyzm, mtr, ntr, mntr, nzyv
      integer m, mx, my, l, i, moff, ioff, joff, koff, loff, k, j, n
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyzb = ny/kyzp
      kzb = nz/kzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
c this segment is used for shared memory computers
c     do 70 my = 1, mblok
c     moff = jblok*(my - 1)
c     koff = kyzp*(my + ks) - 1
c     do 60 mx = 1, jblok
c     m = mx + moff
c     do 50 i = 1, kzb
c     loff = kzp*(i - 1) - 1
c     ioff = jblok*(i - 1) - 1
c     do 40 l = 1, kzp
c     do 30 j = 1, kxyp
c     do 20 k = 1, kyzp
c     do 10 n = 1, ndim
c     h(n+ndim*(l+loff+nzv*(j-1+kxypd*(k-1+kyzpd*(m-1))))) = g(n+ndim*(k
c    1+koff+nyv*(j-1+kxypd*(l-1+kzpd*(mx+ioff)))))
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c  70 continue
c this segment is used for mpi computers
      mlblok = max0(mblok,lblok)
      kyzm = min0(kyzb,kzb)
      mtr = kzb/kyzm
      ntr = kyzb/kyzm
      mntr = max0(mtr,ntr)
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
         koff = ndim*kyzp*is
         is = kzp*(is + ioff) - 1
         do 30 l = 1, kzp
         ir = kxyp*(l + is) - 1
         joff = kxypd*(l + loff) - 1
         do 20 j = 1, kxyp
         do 10 k = 1, ndim*kyzp
         h(k+ndim*kyzp*(j+ir)) = g(k+koff+ndim*nyv*(j+joff))
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
         call MPI_IRECV(g(1+ndim*kxyp*kyzp*kzp*is),ndim*kxyp*kyzp*kzp,mc
     1plx,ir-1,ir+kyzm+1,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.ntr)) then
         ir = is0 + kyzm*(ii - 1)
         is = mx + js + kxb*ir + 1
         call MPI_SEND(h(1+ndim*kxyp*kyzp*kzp*ir),ndim*kxyp*kyzp*kzp,mcp
     1lx,is-1,m+kstrt+kyzm,lgrp,ierr)
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
      do 190 my = 1, mblok
      moff = jblok*(my - 1)
      do 180 mx = 1, jblok
      m = mx + moff
      ioff = kzb*(m - 1)
      koff = kyzpd*(m - 1) - 1
      do 170 i = 1, kyzm
      ir0 = iand(kyzm-1,ieor(my+ks,i-1))
      do 160 ii = 1, mtr
      if (kstrt.le.(kxb*kyzb)) then
         ir = ir0 + kyzm*(ii - 1)
         joff = kzp*ir - 1
         ir = kzp*(ir + ioff) - 1
         do 150 l = 1, kzp
         is = kxyp*(l + ir) - 1
         do 140 j = 1, kxyp
         loff = nzv*(j - 1) + joff
         do 130 k = 1, kyzp
         do 120 n = 1, ndim
         h(ndim*(l+loff+nzyv*(k+koff))+n) = g(ndim*(k+kyzp*(j+is)-1)+n)
  120    continue
  130    continue
  140    continue
  150    continue
      endif
  160 continue
  170 continue
  180 continue
  190 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PN2TPOS3A(f1,f2,g1,g2,s,t,nx,ny,nz,kstrt,nxv,nyv,kxyp,k
     1yp,kzp,kxypd,kypd,kzpd,jblok,kblok,lblok,ndim1,ndim2)
c this subroutine performs a transpose of two matrices f1 and f2,
c distributed in y and z to two matrices g1 and g2,
c distributed in x and z, that is,
c g1(1:ndim1,k+kyp*(m-1),j,l,mx,mz) = f1(1:ndim1,j+kxyp*(l-1),k,l,my,mz)
c g2(1:ndim2,k+kyp*(m-1),j,l,mx,mz) = f2(1:ndim2,j+kxyp*(l-1),k,l,my,mz)
c where, 1 <= j <= kxyp, 1 <= k <= kyp, 1 <= l <= kzp, and
c 1 <= mx <= nx/kxyp, 1 <= my <= ny/kyp, 1 <= mz <= nz/kzp
c and where indices mx, my, and mz can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f1, f2 = complex input arrays
c g1, g2 = complex output arrays
c s, t = complex scratch arrays
c nx/ny/nz = number of points in x/y/z
c kstrt = starting data block number
c nxv/nyv = first dimension of f/g
c kypd/kxypd = second dimension of f/g
c kzpd = third dimension of f and g
c kxyp/kyp/kzp = number of data values per block in x/y/z
c jblok/kblok/lblok = number of data blocks in x/y/z
c ndim1 = leading dimension of arrays f1 and g1
c ndim2 = leading dimension of arrays f2 and g2
      implicit none
      integer nx, ny, nz, kstrt, nxv, nyv, kxyp, kyp, kzp
      integer kxypd, kypd, kzpd, jblok, kblok, lblok, ndim1, ndim2
      complex f1, f2, g1, g2, s, t
      dimension f1(ndim1,nxv,kypd,kzpd,kblok*lblok)
      dimension f2(ndim2,nxv,kypd,kzpd,kblok*lblok)
      dimension g1(ndim1,nyv,kxypd,kzpd,jblok*lblok)
      dimension g2(ndim2,nyv,kxypd,kzpd,jblok*lblok)
      dimension s(ndim1+ndim2,kxyp,kyp,kzp,kblok*lblok)
      dimension t(ndim1+ndim2,kxyp,kyp,kzp,jblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, ls, kxb, kyb, kzb, ndim
      integer jkblok, kxym, mtr, ntr, mntr
      integer m, mx, mz, l, i, moff, ioff, joff, koff, k, j, n
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyb = ny/kyp
      kzb = nz/kzp
      ks = (kstrt - 1)/kyb
      js = kstrt - kyb*ks - 2
      ks = ks - 1
      ls = (kstrt - 1)/kxb - 1
      ndim = ndim1 + ndim2
c this segment is used for shared memory computers
c     do 120 mz = 1, lblok
c     moff = jblok*(mz - 1)
c     ioff = kblok*(mz - 1)
c     do 100 mx = 1, jblok
c     joff = kxyp*(mx + js)
c     m = mx + moff
c     do 50 i = 1, kyb
c     koff = kyp*(i - 1)
c     do 40 l = 1, kzp
c     do 30 k = 1, kyp
c     do 20 j = 1, kxyp
c     do 10 n = 1, ndim1
c     g1(n,k+koff,j,l,m) = f1(n,j+joff,k,l,i+ioff)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c     do 100 i = 1, kyb
c     koff = kyp*(i - 1)
c     do 90 l = 1, kzp
c     do 80 k = 1, kyp
c     do 70 j = 1, kxyp
c     do 60 n = 1, ndim2
c     g2(n,k+koff,j,l,m) = f2(n,j+joff,k,l,i+ioff)
c  60 continue
c  70 continue
c  80 continue
c  90 continue
c 100 continue
c 110 continue
c 120 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      do 170 mz = 1, lblok
      moff = jkblok*(mz - 1)
      do 160 mx = 1, jkblok
      m = mx + moff
      do 150 i = 1, kxym
      ir0 = iand(kxym-1,ieor(mx+js,i-1))
      is0 = ir0
      do 140 ii = 1, mntr
c post receive
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         koff = kyp*ir
         ir = ir + kyb*(mz + ls) + 1
         call MPI_IRECV(t(1,1,1,1,m),ndim*kxyp*kyp*kzp,mcplx,ir-1,ir+kxy
     1m+1,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kyb*kzb)).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         joff = kxyp*is
         is = is + kxb*(mz + ks) + 1
         do 40 l = 1, kzp
         do 30 k = 1, kyp
         do 20 j = 1, kxyp
         do 10 n = 1, ndim1
         s(n,j,k,l,m) = f1(n,j+joff,k,l,m)
   10    continue
   20    continue
   30    continue
   40    continue
         do 80 l = 1, kzp
         do 70 k = 1, kyp
         do 60 j = 1, kxyp
         do 50 n = 1, ndim2
         s(n+ndim1,j,k,l,m) = f2(n,j+joff,k,l,m)
   50    continue
   60    continue
   70    continue
   80    continue
         call MPI_SEND(s(1,1,1,1,m),ndim*kxyp*kyp*kzp,mcplx,is-1,m+kstrt
     1+kxym,lgrp,ierr)
      endif
c receive data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
         do 130 l = 1, kzp
         do 120 k = 1, kyp
         do 110 j = 1, kxyp
         do 90 n = 1, ndim1
         g1(n,k+koff,j,l,m) = t(n,j,k,l,m)
   90    continue
         do 100 n = 1, ndim2
         g2(n,k+koff,j,l,m) = t(n+ndim1,j,k,l,m)
  100    continue
  110    continue
  120    continue
  130    continue
      endif
  140 continue
  150 continue
  160 continue
  170 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PN2TPOS3B(g1,g2,h1,h2,s,t,nx,ny,nz,kstrt,nyv,nzv,kxyp,k
     1yzp,kzp,kxypd,kyzpd,kzpd,jblok,mblok,lblok,ndim1,ndim2)
c this subroutine performs a transpose of two matrices g1 and g2,
c distributed in x and z to two matrices h1 and h2,
c distributed in x and y, that is,
c h1(1:ndim1,j+kzp*(l-1),k,l,my,mz) = g1(1:ndim1,k+kyzp*(m-1),j,l,mx,mz)
c h2(1:ndim2,j+kzp*(l-1),k,l,my,mz) = g2(1:ndim2,k+kyzp*(m-1),j,l,mx,mz)
c where, 1 <= j <= kxyp, 1 <= k <= kyzp, 1 <= l <= kzp, and
c 1 <= mx <= nx/kxyp, 1 <= my <= ny/kyzp, 1 <= mz <= nz/kzp
c and where indices mx, my, and mz can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c g1, g2 = complex input arrays
c h1, h2 = complex output arrays
c s, t = complex scratch arrays
c nx/ny/nz = number of points in x/y/z
c kstrt = starting data block number
c nyv/nzv = first dimension of g/h
c kxypd = second dimension of f and g
c kzpd/kyzpd = third dimension of g/h
c kxyp/kyzp/kzp = number of data values per block in x/y/z
c jblok/mblok/lblok = number of data blocks in x/y/z
c ndim1 = leading dimension of arrays g1 and h1
c ndim2 = leading dimension of arrays g2 and h2
      implicit none
      integer nx, ny, nz, kstrt, nyv, nzv, kxyp, kyzp, kzp
      integer kxypd, kyzpd, kzpd, jblok, mblok, lblok, ndim1, ndim2
      complex g1, g2, h1, h2, s, t
      dimension g1(ndim1,nyv,kxypd,kzpd,jblok*lblok)
      dimension g2(ndim2,nyv,kxypd,kzpd,jblok*lblok)
      dimension h1(ndim1,nzv,kxypd,kyzpd,jblok*mblok)
      dimension h2(ndim2,nzv,kxypd,kyzpd,jblok*mblok)
      dimension s(ndim1+ndim2,kyzp,kxyp,kzp,jblok*lblok)
      dimension t(ndim1+ndim2,kyzp,kxyp,kzp,jblok*mblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, kxb, kyzb, kzb, ndim
      integer mlblok, kyzm, mtr, ntr, mntr
      integer m, mx, my, l, i, moff, ioff, koff, loff, k, j, n
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyzb = ny/kyzp
      kzb = nz/kzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      ndim = ndim1 + ndim2
c this segment is used for shared memory computers
c     do 120 my = 1, mblok
c     moff = jblok*(my - 1)
c     koff = kyzp*(my + ks)
c     do 110 mx = 1, jblok
c     m = mx + moff
c     do 50 i = 1, kzb
c     loff = kzp*(i - 1)
c     ioff = jblok*(i - 1)
c     do 40 l = 1, kzp
c     do 30 j = 1, kxyp
c     do 20 k = 1, kyzp
c     do 10 n = 1, ndim1
c     h1(n,l+loff,j,k,m) = g1(n,k+koff,j,l,mx+ioff)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c     do 100 i = 1, kzb
c     loff = kzp*(i - 1)
c     ioff = jblok*(i - 1)
c     do 90 l = 1, kzp
c     do 80 j = 1, kxyp
c     do 70 k = 1, kyzp
c     do 60 n = 1, ndim2
c     h2(n,l+loff,j,k,m) = g2(n,k+koff,j,l,mx+ioff)
c  60 continue
c  70 continue
c  80 continue
c  90 continue
c 100 continue
c 110 continue
c 120 continue
c this segment is used for mpi computers
      mlblok = max0(mblok,lblok)
      kyzm = min0(kyzb,kzb)
      mtr = kzb/kyzm
      ntr = kyzb/kyzm
      mntr = max0(mtr,ntr)
      do 170 my = 1, mlblok
      moff = jblok*(my - 1)
      do 160 mx = 1, jblok
      m = mx + moff
      do 150 i = 1, kyzm
      ir0 = iand(kyzm-1,ieor(my+ks,i-1))
      is0 = ir0
      do 140 ii = 1, mntr
c post receive
      if ((kstrt.le.(kxb*kyzb)).and.(ii.le.mtr)) then
         ir = ir0 + kyzm*(ii - 1)
         loff = kzp*ir
         ir = mx + js + kxb*ir + 1
         call MPI_IRECV(t(1,1,1,1,m),ndim*kxyp*kyzp*kzp,mcplx,ir-1,ir+ky
     1zm+1,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.ntr)) then
         is = is0 + kyzm*(ii - 1)
         koff = kyzp*is
         is = mx + js + kxb*is + 1
         do 40 l = 1, kzp
         do 30 j = 1, kxyp
         do 20 k = 1, kyzp
         do 10 n = 1, ndim1
         s(n,k,j,l,m) = g1(n,k+koff,j,l,m)
   10    continue
   20    continue
   30    continue
   40    continue
         do 80 l = 1, kzp
         do 70 j = 1, kxyp
         do 60 k = 1, kyzp
         do 50 n = 1, ndim2
         s(n+ndim1,k,j,l,m) = g2(n,k+koff,j,l,m)
   50    continue
   60    continue
   70    continue
   80    continue
         call MPI_SEND(s(1,1,1,1,m),ndim*kxyp*kyzp*kzp,mcplx,is-1,m+kstr
     1t+kyzm,lgrp,ierr)
      endif
c receive data
      if ((kstrt.le.(kxb*kyzb)).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
         do 130 l = 1, kzp
         do 120 j = 1, kxyp
         do 110 k = 1, kyzp
         do 90 n = 1, ndim1
         h1(n,l+loff,j,k,m) = t(n,k,j,l,m)
   90    continue
         do 100 n = 1, ndim2
         h2(n,l+loff,j,k,m) = t(n+ndim1,k,j,l,m)
  100    continue
  110    continue
  120    continue
  130    continue
      endif
  140 continue
  150 continue
  160 continue
  170 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSUM2(f,g,nxp,kstrt,nvpy,nvpz,nd,mblok,nblok)
c this subroutine performs multiple parallel sums of a vector, that is:
c f(j,k) = sum over k of f(j,k) in y or z direction
c assumes the number of processors nvpy, nvpz are a power of two.
c the algorithm performs partial sums in binary pairs, as follows:
c first, adjacent processors exchange vectors and sum them.  next,
c processors separated by 2 exchange the new vectors and sum them, then
c those separated by 4, up to processors separated by nvpy/2, or nvpz/2.
c at the end, all processors contain the same summation.
c f = input and output data
c g = scratch array
c nxp = number of data values in vector
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nd = (0,1) = integrate in (y,z) direction.
c mblok/nblok = number of field partitions in y/z.
c written by viktor k. decyk, ucla
      implicit none
      real f, g
      integer nxp, kstrt, nvpy, nvpz, nd, mblok, nblok
      dimension f(nxp,mblok*nblok), g(nxp,mblok*nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, ierr, msid
      integer js, ks, nvp, l, kxs, kys, ky, kz, k, kb, lb, j, mnblok
      dimension istatus(lstat)
      if (kstrt.gt.(nvpy*nvpz)) return
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      if (nd.eq.1) then
         nvp = nvpy
         kys = 1
      else if (nd.eq.2) then
         nvp = nvpz
         kys = nvpy
      else
         return
      endif
      mnblok = mblok*nblok
      l = 1
      kxs = 1
c main iteration loop
   10 if (kxs.ge.nvp) go to 70
c shift data
      do 40 kz = 1, nblok
      do 30 ky = 1, mblok
      kb = ky + js + nvpy*(kz + ks)
      lb = kb/kys
      kb = kb + 1
      lb = lb - 2*(lb/2)
c this loop is used for shared memory computers
c     do 20 j = 1, nxp
c     if (lb.eq.0) then
c        g(j,k) = f(j,kb+kys)
c     else
c        g(j,k) = f(j,kb-kys)
c     endif
c  20 continue
c this segment is used for mpi computers
      if (lb.eq.0) then
         call MPI_IRECV(g,nxp,mreal,kb+kys-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(f,nxp,mreal,kb+kys-1,l+nxp,lgrp,ierr)
      else
         call MPI_IRECV(g,nxp,mreal,kb-kys-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(f,nxp,mreal,kb-kys-1,l+nxp,lgrp,ierr)
      endif
      call MPI_WAIT(msid,istatus,ierr)
   30 continue
   40 continue
c perform sum
      do 60 k = 1, mnblok
      do 50 j = 1, nxp
      f(j,k) = f(j,k) + g(j,k)
   50 continue
   60 continue
      l = l + 1
      kxs = kxs + kxs
      kys = kys + kys
      go to 10
   70 return
      end
c-----------------------------------------------------------------------
      subroutine PISUM2(if,ig,nxp,kstrt,nvpy,nvpz,nd,mblok,nblok)
c this subroutine performs multiple parallel sums of a vector, that is:
c if(j,k) = sum over k of if(j,k) in y or z direction
c assumes the number of processors nvpy, nvpz are a power of two.
c the algorithm performs partial sums in binary pairs, as follows:
c first, adjacent processors exchange vectors and sum them.  next,
c processors separated by 2 exchange the new vectors and sum them, then
c those separated by 4, up to processors separated by nvpy/2, or nvpz/2.
c at the end, all processors contain the same summation.
c if = input and output integer data
c ig = scratch integer array
c nxp = number of data values in vector
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nd = (0,1) = integrate in (y,z) direction.
c mblok/nblok = number of field partitions in y/z.
c written by viktor k. decyk, ucla
      implicit none
      integer if, ig
      integer nxp, kstrt, nvpy, nvpz, nd, mblok, nblok
      dimension if(nxp,nblok), ig(nxp,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mint = default datatype for integers
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, ierr, msid
      integer js, ks, nvp, l, kxs, kys, ky, kz, k, kb, lb, j, mnblok
      dimension istatus(lstat)
      if (kstrt.gt.(nvpy*nvpz)) return
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      if (nd.eq.1) then
         nvp = nvpy
         kys = 1
      else if (nd.eq.2) then
         nvp = nvpz
         kys = nvpy
      else
         return
      endif
      mnblok = mblok*nblok
      l = 1
      kxs = 1
c main iteration loop
   10 if (kxs.ge.nvp) go to 70
c shift data
      do 40 kz = 1, nblok
      do 30 ky = 1, mblok
      kb = ky + js + nvpy*(kz + ks)
      lb = kb/kys
      kb = kb + 1
      lb = lb - 2*(lb/2)
c this loop is used for shared memory computers
c     do 20 j = 1, nxp
c     if (lb.eq.0) then
c        ig(j,k) = if(j,kb+kys)
c     else
c        ig(j,k) = if(j,kb-kys)
c     endif
c  20 continue
c this segment is used for mpi computers
      if (lb.eq.0) then
         call MPI_IRECV(ig,nxp,mint,kb+kys-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(if,nxp,mint,kb+kys-1,l+nxp,lgrp,ierr)
      else
         call MPI_IRECV(ig,nxp,mint,kb-kys-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(if,nxp,mint,kb-kys-1,l+nxp,lgrp,ierr)
      endif
      call MPI_WAIT(msid,istatus,ierr)
   30 continue
   40 continue
c perform sum
      do 60 k = 1, mnblok
      do 50 j = 1, nxp
      if(j,k) = if(j,k) + ig(j,k)
   50 continue
   60 continue
      l = l + 1
      kxs = kxs + kxs
      kys = kys + kys
      go to 10
   70 return
      end
c-----------------------------------------------------------------------
      subroutine PSCAN2(f,g,s,nxp,kstrt,nvpy,nvpz,nd,mblok,nblok)
c this subroutine performs multiple parallel prefix reductions of a
c vector, that is: f(j,k) = sum over k of f(j,k), where the sum is over
c k values less than the processor id in the appropriate direction.
c assumes the number of processors nvpy, nvpz are a power of two.
c the algorithm performs partial sums in binary pairs, as follows:
c first, adjacent processors exchange vectors and sum them.  next,
c processors separated by 2 exchange the new vectors and sum them, then
c those separated by 4, up to processors separated by separated by
c nvpy/2, or nvpz/2.
c at the end, all processors contain the same summation.
c f = input and output data
c g, s = scratch array
c nxp = number of data values in vector
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nd = (0,1) = integrate in (y,z) direction.
c mblok/nblok = number of field partitions in y/z.
c written by viktor k. decyk, ucla
      implicit none
      real f, g, s
      integer nxp, kstrt, nvpy, nvpz, nd, mblok, nblok
      dimension f(nxp,mblok*nblok), g(nxp,mblok*nblok)
      dimension s(nxp,mblok*nblok)
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
      integer js, ks, nvp, l, kxs, kys, ky, kz, koff, k, kb, lb, j
      integer mnblok
      dimension istatus(lstat)
      if (kstrt.gt.(nvpy*nvpz)) return
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      if (nd.eq.1) then
         nvp = nvpy
         kys = 1
      else if (nd.eq.2) then
         nvp = nvpz
         kys = nvpy
      else
         return
      endif
      mnblok = mblok*nblok
      l = 1
      kxs = 1
c initialize global sum
      do 20 k = 1, mnblok
      do 10 j = 1, nxp
      s(j,k) = f(j,k)
   10 continue
   20 continue
c main iteration loop
   30 if (kxs.ge.nvp) go to 100
c shift data
      do 70 kz = 1, nblok
      koff = mblok*(kz - 1)
      do 60 ky = 1, mblok
      k = ky + koff
      kb = ky + js + nvpy*(kz + ks)
      lb = kb/kys
      kb = kb + 1
      lb = lb - 2*(lb/2)
c this loop is used for shared memory computers
c     do 40 j = 1, nxp
c     if (lb.eq.0) then
c        g(j,k) = s(j,kb+kys)
c     else
c        g(j,k) = s(j,kb-kys)
c     endif
c  40 continue
c this segment is used for mpi computers
      if (lb.eq.0) then
         call MPI_IRECV(g,nxp,mreal,kb+kys-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(s,nxp,mreal,kb+kys-1,l+nxp,lgrp,ierr)
      else
         call MPI_IRECV(g,nxp,mreal,kb-kys-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(s,nxp,mreal,kb-kys-1,l+nxp,lgrp,ierr)
      endif
      call MPI_WAIT(msid,istatus,ierr)
c perform prefix scan
      if (lb.ne.0) then
         do 50 j = 1, nxp
         f(j,k) = f(j,k) + g(j,k)
   50    continue
      endif
   60 continue
   70 continue
c perform sum
      do 90 k = 1, mnblok
      do 80 j = 1, nxp
      s(j,k) = s(j,k) + g(j,k)
   80 continue
   90 continue
      l = l + 1
      kxs = kxs + kxs
      kys = kys + kys
      go to 30
  100 return
      end
c-----------------------------------------------------------------------
      subroutine P0COPY2(f,g,nxp,kstrt,nvpy,nvpz,nd)
c this subroutine copies real data f on node 0 in appropriate direction
c to g on other nodes
c f/g = data to be sent/received
c nxp = size of data f
c input: f, nxp
c output: g
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nd = (0,1) = integrate in (y,z) direction.
      implicit none
      real f, g
      integer nxp, kstrt, nvpy, nvpz, nd
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
      integer istatus, js, ks, idproc, idr, nvp, nxs, i, id, ierr
      dimension istatus(lstat)
      if (kstrt.gt.(nvpy*nvpz)) return
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 1
      if (nd.eq.1) then
         idproc = js
         idr = nvpy*ks
         nvp = nvpy
         nxs = 1
      else if (nd.eq.2) then
         idproc = ks
         idr = js
         nvp = nvpz
         nxs = nvpy
      else
         return
      endif
c this segment is used for mpi computers
c node 0 sends messages to other nodes
      if (idproc.eq.0) then
         do 10 i = 2, nvp
            id = idr + nxs*(i - 1)
            call MPI_SEND(f,nxp,mreal,id,95,lgrp,ierr)
   10    continue
c other nodes receive data from node 0
      else
         call MPI_RECV(g,nxp,mreal,idr,95,lgrp,istatus,ierr)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PWRITE32(f,nx,kyp,kzp,nxv,kypmx,kzpmx,mnblok,iunit,nrec
     1,lrec,name)
c this subroutine collects distributed real 3d data f and writes to a
c direct access binary file with 2D spatial decomposition
c f = input data to be written, modified on node 0
c nx/kyp/kzp = length of data f in x/y/z on each processor to write
c nxv = first dimension of data array f, must be >= nx
c kypmx = second dimension of data array f, must be >= kyp
c kzpmx = third dimension of data array f, must be >= kzp
c mnblok = number of parallel partitions
c iunit = fortran unit number
c nrec = current record number for write, if nrec > 0
c if nrec < 0, open new file and write first record
c if nrec = 0, open old file, do not write
c lrec = record length (used only if nrec <= 0)
c name = file name (used only if nrec <= 0)
c input: f,nx,kyp,kzp,nxv,kypmx,kzpmx,mnblok,iunit,nrec,lrec,name
c output: nrec
      implicit none
      integer nx, kyp, kzp, nxv, kypmx, kzpmx, mnblok, iunit, nrec, lrec
      real f
      character*(*) name
      dimension f(nxv,kypmx,kzpmx,mnblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, ierr
      integer nvp, idproc, np, ioff, id, nrec0, i, j, k, l, m, nxyv
      dimension istatus(lstat)
      nxyv = nxv*kypmx
c this segment is used for shared memory computers
c     if (nrec.lt.0) then
c        open(unit=iunit,file=name,form='unformatted',access='direct',re
c    1cl=lrec,status='replace')
c        nrec = 1
c open old file
c     else if (nrec.eq.0) then
c        open(unit=iunit,file=name,form='unformatted',access='direct',re
c    1cl=lrec,status='old')
c     endif
c     write (unit=iunit,rec=nrec) ((((f(j,k,l,m),j=1,nx),k=1,kyp),l=1,kz
c    1p),m=1,mnblok)
c     nrec = nrec + 1
c this segment is used for mpi computers
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)
c node 0 receives messages from other nodes
      if (idproc.eq.0) then
         if (nrec.lt.0) then
            open(unit=iunit,file=name,form='unformatted',access='direct'
     1,recl=lrec,status='replace')
            nrec = 1
c open old file
         else if (nrec.eq.0) then
            open(unit=iunit,file=name,form='unformatted',access='direct'
     1,recl=lrec,status='old')
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
            call MPI_RECV(f,nxyv*kzp,mreal,id,99,lworld,istatus,ierr)
         endif
c first write data for node 0
         nrec0 = nrec
         write (unit=iunit,rec=nrec) ((((f(j,k,l,m),j=1,nx),k=1,kyp),l=1
     1,kzp),m=1,mnblok)
         nrec = nrec + 1
c then write data from remaining nodes
         do 10 i = 2, np
         id = i - ioff
         call MPI_RECV(f,nxyv*kzp,mreal,id,99,lworld,istatus,ierr)
         write (unit=iunit,rec=nrec) ((((f(j,k,l,m),j=1,nx),k=1,kyp),l=1
     1,kzp),m=1,mnblok)
         nrec = nrec + 1
   10    continue
c read data back for node 0
         read (unit=iunit,rec=nrec0) ((((f(j,k,l,m),j=1,nx),k=1,kyp),l=1
     1,kzp),m=1,mnblok)
c other nodes send data to node 0
      elseif (idproc.le.(nproc+1)) then
         call MPI_SEND(f,nxyv*kzp,mreal,0,99,lworld,ierr)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PREAD32(f,nx,kyp,kzp,nxv,kypmx,kzpmx,mnblok,iunit,nrec,
     1lrec,name,ierror)
c this subroutine reads real 3d data f from a direct access binary file
c and distributes it with 2D spatial decomposition
c f = output data to be read
c nx/kyp/kzp = length of data f in x/y/z on each processor to read
c nxv = first dimension of data array f, must be >= nx
c kzpmx = second dimension of data array f, must be >= kyp
c kzpmx = third dimension of data array f, must be >= kzp
c mnblok = number of parallel partitions
c iunit = fortran unit number
c nrec = current record number for read, if nrec > 0
c if nrec < 0, open new file and read first record
c if nrec = 0, open old file, do not read
c lrec = record length (used only if nrec <= 0)
c name = file name (used only if nrec <= 0)
c ierror = error indicator
c input: f,nx,kyp,kzp,nxv,kypmx,kzpmx,mnblok,iunit,nrec,lrec,name
c output: f, nrec, ierror
      implicit none
      integer nx, kyp, kzp, nxv, kypmx, kzpmx, mnblok, iunit, nrec, lrec
      integer ierror
      real f
      character*(*) name
      dimension f(nxv,kypmx,kzpmx,mnblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, ierr
      integer nvp, idproc, np, ioff, id, nrec0, i, j, k, l, m, nxyv
      dimension istatus(lstat)
      nxyv = nxv*kypmx
c this segment is used for shared memory computers
c     if (nrec.lt.0) then
c        open(unit=iunit,file=name,form='unformatted',access='direct',re
c    1cl=lrec,status='old')
c        if (nrec.eq.0) return
c        nrec = 1
c     endif
c     read (unit=iunit,rec=nrec,err=10) ((((f(j,k,l,m),j=1,nx),k=1,kyp),
c    1l=1,kzp),m=1,mnblok)
c     nrec = nrec + 1
c     return
c 10 ierror = 1
c     return
c this segment is used for mpi computers
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)
c node 0 receives messages from other nodes
      if (idproc.eq.0) then
         if (nrec.lt.0) then
            open(unit=iunit,file=name,form='unformatted',access='direct'
     1,recl=lrec,status='old')
            if (nrec.eq.0) return
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
         do 30 i = 2, np
            read (unit=iunit,rec=nrec,err=10) ((((f(j,k,l,m),j=1,nx),k=1
     1,kyp),l=1,kzp),m=1,mnblok)
            go to 20
   10       ierror = ierror + 1
   20       nrec = nrec + 1
            id = i - ioff
            call MPI_SEND(f,nxyv*kzp,mreal,id,98,lworld,ierr)
   30    continue
c then read data from node 0
         read (unit=iunit,rec=nrec0,err=40) ((((f(j,k,l,m),j=1,nx),k=1,k
     1yp),l=1,kzp),m=1,mnblok)
         go to 50
   40    ierror = ierror + 1
   50    if (ioff.eq.0) then
            id = 1
            call MPI_SEND(f,nxyv*kzp,mreal,id,98,lworld,ierr)
         endif
c other nodes receive data from node 0
      else if (idproc.le.(nproc+1)) then
         if (nrec.eq.0) return
         call MPI_RECV(f,nxyv*kzp,mreal,0,98,lworld,istatus,ierr)
      endif
c check for error condition
      call PIMAX(ierror,ierr,1,1)
      return
      end
c-----------------------------------------------------------------------
      subroutine PCWRITE32(f,g,nx,ny,nz,kxyp,kyzp,nzv,kxypd,kyzpd,jblok,
     1mblok,iunit,nrec,lrec,name)
c this subroutine collects distributed complex 3d data f and writes to a
c direct access binary file with 2D spatial decomposition
c f = input data to be written
c g = scratch data
c nx/ny/nz = total length of data f in x/y/z to write
c kxyp/kyzp = maximum length of data f in x/y on each processor to write
c nzv = first dimension of data array f, must be >= nz
c kxypd = second dimension of data array f, must be >= kxyp
c kyzpd = third dimension of data array f, must be >= kyzp
c jblok/mblok = number of field partitions in x/y
c iunit = fortran unit number
c nrec = current record number for write, if nrec > 0
c if nrec < 0, open new file and write first record
c if nrec = 0, open old file, do not write
c lrec = record length (used only if nrec <= 0)
c name = file name (used only if nrec <= 0)
c input: f, nx, ny, nz, kxyp, kyzp, nzv, kxypd, kyzpd, jblok, mblok,
c        iunit ,nrec, lrec, name
c output: nrec
      implicit none
      integer nx, ny, nz, kxyp, kyzp, nzv, kxypd, kyzpd, jblok, mblok
      integer iunit, nrec, lrec
      complex f, g
      character*(*) name
      dimension f(nzv,kxypd,kyzpd,jblok*mblok)
      dimension g(nzv,kxypd*kyzpd,jblok*mblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus
      integer js, ks, kxb, kxpp, kxypp, kyzpp, nvp, idproc, np, ioff, id
      integer i, j, k, l, mx, my, m, joff, moff, ierr
      dimension istatus(lstat)
      do 80 my = 1, mblok
      moff = jblok*(my - 1)
      do 70 mx = 1, jblok
      m = mx + moff
c this segment is used for shared memory computers
c     if (nrec.lt.0) then
c        open(unit=iunit,file=name,form='unformatted',access='direct',re
c    1cl=lrec,status='replace')
c        nrec = 1
c open old file
c     else if (nrec.eq.0) then
c        open(unit=iunit,file=name,form='unformatted',access='direct',re
c    1cl=lrec,status='old')
c     endif
c     do 20 k = 1, ny
c     do 10 j = 1, nx
c     write (unit=iunit,rec=nrec) (f(l,j,k,m),l=1,nz)
c     nrec = nrec + 1
c  10 continue
c  20 continue
c this segment is used for mpi computers
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)
c node 0 receives messages from other nodes
      if (idproc.eq.0) then
         if (nrec.lt.0) then
            open(unit=iunit,file=name,form='unformatted',access='direct'
     1,recl=lrec,status='replace')
            nrec = 1
c open old file
         else if (nrec.eq.0) then
            open(unit=iunit,file=name,form='unformatted',access='direct'
     1,recl=lrec,status='old')
         return
         endif
c no special diagnostic node
         if (nvp.eq.nproc) then
            np = nvp
            ioff = 1
            kxypp = min(nx,kxyp)
            kyzpp = min(ny,kyzp)
            kxpp = kxypp*kyzpp
c special diagnostic node present
         else
            np = nvp - nproc
            ioff = 0
            id = 1
            call MPI_RECV(g,nzv*kxyp*kyzp,mcplx,id,99,lworld,istatus,ier
     1r)
            call MPI_GET_COUNT(istatus,mcplx,kxpp,ierr)
            kxpp = kxpp/nzv
         endif
c first write data for node 0
         do 10 j = 1, kxpp
            write (unit=iunit,rec=nrec) (g(l,j,m),l=1,nz)
            nrec = nrec + 1
   10    continue
c then write data from remaining nodes
         do 30 i = 2, np
            id = i - ioff
            call MPI_RECV(g,nzv*kxyp*kyzp,mcplx,id,99,lworld,istatus,ier
     1r)
            call MPI_GET_COUNT(istatus,mcplx,kxpp,ierr)
            kxpp = kxpp/nzv
            do 20 j = 1, kxpp
               write (unit=iunit,rec=nrec) (g(l,j,m),l=1,nz)
               nrec = nrec + 1
   20       continue
   30    continue
c other nodes send data to node 0
      elseif (idproc.le.(nproc+1)) then
         if (nrec.eq.0) return
c find amount of data to write
         kxb = (nx/2)/kxyp
         js = idproc
         if (nvp.eq.nproc) js = js + 1
         ks = js/kxb
         js = js - kxb*ks - 1
         ks = ks - 1
         kxypp = nx - kxyp*(mx + js)
         kyzpp = ny - kyzp*(my + ks)
         if (kxypp.gt.kxyp) then
            kxypp = kxyp
         else if (kxypp.le.0) then
            kxypp = 0
         endif
         if (kyzpp.gt.kyzp) then
            kyzpp = kyzp
         else if (kyzpp.le.0) then
            kyzpp = 0
         endif
         kxpp = kxypp*kyzpp
c copy the data
         do 60 k = 1, kyzpp
         joff = kxypp*(k - 1)
         do 50 j = 1, kxypp
         do 40 l = 1, nzv
         g(l,j+joff,m) = f(l,j,k,m)
   40    continue
   50    continue
   60    continue
c send the data
         call MPI_SEND(g,nzv*kxpp,mcplx,0,99,lworld,ierr)
      endif
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCREAD32(f,g,nx,ny,nz,kxyp,kyzp,nzv,kxypd,kyzpd,jblok,m
     1blok,iunit,nrec,lrec,name,ierror)
c this subroutine reads complex 3d data f from a direct access binary
c file and distributes it with 2D spatial decomposition
c f = output data to be read
c g = scratch data
c nx/ny/nz = total length of data f in x/y/z to read
c kxyp/kyzp = maximum length of data f in x/y on each processor to write
c nzv = first dimension of data array f, must be >= nz
c kxypd = second dimension of data array f, must be >= kxyp
c kyzpd = third dimension of data array f, must be >= kyzp
c jblok/mblok = number of field partitions in x/y
c iunit = fortran unit number
c nrec = current record number for write, if nrec > 0
c if nrec < 0, open old file and read first record
c if nrec = 0, open old file, do not read
c lrec = record length (used only if nrec <= 0)
c name = file name (used only if nrec <= 0)
c input: nx, ny, nz, kxyp, kyzp, nzv, kxypd, kyzpd, jblok, mblok,
c        iunit ,nrec, lrec, name
c output: f, nrec, ierror
      implicit none
      integer nx, ny, nz, kxyp, kyzp, nzv, kxypd, kyzpd, jblok, mblok
      integer iunit, nrec, lrec, ierror
      complex f, g
      character*(*) name
      dimension f(nzv,kxypd,kyzpd,jblok*mblok)
      dimension g(nzv,kxypd*kyzpd,jblok*mblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c lworld = MPI_COMM_WORLD communicator
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus
      integer js, ks, kxb, kxpp, kxypp, kyzpp, nvp, idproc, np, ioff, id
      integer nrec0, i, j, k, l, mx, my, m, joff, moff, ierr
      dimension istatus(lstat)
      ierror = 0
      do 150 my = 1, mblok
      moff = jblok*(my - 1)
      do 140 mx = 1, jblok
      m = mx + moff
c this segment is used for shared memory computers
c     if (nrec.le.0) then
c        open(unit=iunit,file=name,form='unformatted',access='direct',re
c    1cl=lrec,status='old')
c        if (nrec.eq.0) return
c        nrec = 1
c     endif
c     do 20 k = 1, ny
c     do 10 j = 1, nx
c     read (unit=iunit,rec=nrec,err=30) (f(l,j,k,m),l=1,nz)
c     nrec = nrec + 1
c  10 continue
c  20 continue
c     return
c  30 ierror = 1
c     return
c this segment is used for mpi computers
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lworld,idproc,ierr)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lworld,nvp,ierr)
c node 0 sends messages to other nodes
      if (idproc.eq.0) then
         if (nrec.le.0) then
            open(unit=iunit,file=name,form='unformatted',access='direct'
     1,recl=lrec,status='old')
            if (nrec.eq.0) return
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
         nrec = nrec + min(nx,kxyp)*min(ny,kyzp)
         do 40 i = 2, np
c find amount of data to read
            kxb = (nx/2)/kxyp
            js = i - 2
            ks = js/kxb
            js = js - kxb*ks - 1
            ks = ks - 1
            kxypp = nx - kxyp*(mx + js)
            kyzpp = ny - kyzp*(my + ks)
            if (kxypp.gt.kxyp) then
               kxypp = kxyp
            else if (kxypp.le.0) then
               kxypp = 0
            endif
            if (kyzpp.gt.kyzp) then
               kyzpp = kyzp
            else if (kyzpp.le.0) then
               kyzpp = 0
            endif
            kxpp = kxypp*kyzpp
            do 10 j = 1, kxpp
               read (unit=iunit,rec=nrec,err=20) (g(l,j,m),l=1,nz)
               nrec = nrec + 1
   10       continue
            go to 30
   20       ierror = ierror + 1
            nrec = nrec - 1
   30       id = i - ioff
            call MPI_SEND(g,nzv*kxpp,mcplx,id,98,lworld,ierr)
   40    continue
c then read data from node 0
         kxpp = min(nx,kxyp)*min(ny,kyzp)
         do 50 j = 1, kxpp
            read (unit=iunit,rec=nrec0,err=60) (g(l,j,m),l=1,nz)
            nrec0 = nrec0 + 1
   50    continue
         go to 70
   60    ierror = ierror + 1
   70    if (ioff.eq.0) then
            id = 1
            call MPI_SEND(g,nzv*kxpp,mcplx,id,98,lworld,ierr)
         else
            do 100 k = 1, kyzpp
            joff = kxypp*(k - 1)
            do 90 j = 1, kxypp
            do 80 l = 1, nzv
            f(l,j,k,m) = g(l,j+joff,m)
   80       continue
   90       continue
  100       continue
         endif
c other nodes receive data from node 0
      else if (idproc.le.(nproc+1)) then
         if (nrec.eq.0) return
         call MPI_RECV(g,nzv*kxyp*kyzp,mcplx,0,98,lworld,istatus,ierr)
c copy the data
         do 130 k = 1, kyzpp
         joff = kxypp*(k - 1)
         do 120 j = 1, kxypp
         do 110 l = 1, nzv
         f(l,j,k,m) = g(l,j+joff,m)
  110    continue
  120    continue
  130    continue
      endif
  140 continue
  150 continue
c check for error condition
      call PIMAX(ierror,ierr,1,1)
      return
      end
