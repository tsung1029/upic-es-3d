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
c info(9) = total number of particles on exit
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
      real an, xt
      dimension msid(4), istatus(lstat)
      dimension ibflg(4), iwork(4)
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
      ibflg(2) = nps
      ibflg(1) = npr
      call PISUM(ibflg,iwork,2,1)
      info(8) = ibflg(1)
      info(9) = ibflg(2)
      if (ibflg(1).ne.ibflg(2)) then
         write (2,*) 'particle number error, old/new=',ibflg(1),ibflg(2)
         info(1) = 1
      endif
c information
  320 if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, nbmax=', nbmax
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
c nypmx = second dimension of f, must be >= ny
c nzpmx = maximum size of field partition, including guard cells.
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
         nyzpmn = nypmx - 1
         mnter = mter
      else if (n.eq.2) then
         nyzpmn = nzpmx - 1
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
      subroutine REPARTD32(edges,edg,eds,eg,es,et2,npicyz,noff,nyzp,npav
     1,nypmin,nypmax,nzpmin,nzpmax,kstrt,nvpy,nvpz,mblok,nblok,idps,idds
     2,myzpm1)
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
c npav = average number of particles per partition desired
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
      integer npav, nypmin, nypmax, nzpmin, nzpmax, kstrt, nvpy, nvpz
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
      integer js, ks, iter, nter, nyzp1, k1, kb, kl, kr, k, my, mz, m, n
      integer nvp, moff, mnblok, ierr
      real sum1, at1, at2, anpav, apav, anpl, anpr
      integer msid, istatus
      integer ibflg, iwork
      dimension istatus(lstat)
      dimension ibflg(4), iwork(4)
      dimension kb(2)
c exit if flag is set
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      mnblok = mblok*nblok
      iter = 2
c main loop over decompositions
      do 260 n = 1, 2
      if (n.eq.1) then
         nvp = nvpy
         anpav = real(npav*nvpz)
      elseif (n.eq.2) then
         nvp = nvpz
         anpav = real(npav*nvpy)
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
      apav = real(kb(n))*anpav
c anpl = deficit of particles on processor to left
      anpl = apav - et2(1,m) + es(1,m)
c anpr = excess of particles on current processor
      anpr = et2(1,m) - apav - anpav
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
      apav = real(kb(n))*anpav
      anpl = apav - et2(1,m) + es(1,m)
      anpr = et2(1,m) - apav - anpav
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
      apav = real(kb(n))*anpav
      anpl = apav - et2(1,m) + es(1,m)
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
      if (anpl.gt.anpav) then
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
      apav = real(kb(n))*anpav
      anpl = apav - et2(1,m) + es(1,m)
      anpr = et2(1,m) - apav - anpav
c left boundary is on the right
      if (et2(1,m).lt.apav) then
         if ((et2(1,m)+eg(1,m)).ge.apav) then
            nyzp1 = eg(2,m)
            k1 = 1
            sum1 = 0.
            at2 = - (anpr + anpav)
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
      if ((anpl-es(1,m)).gt.anpav) nter = nter + 1
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
