c-----------------------------------------------------------------------
      program pfft32xtime
      implicit none
      integer indx, indy, indz, indnvpy, indnvpz, mshare
      integer nx, ny, nz, nxh, nvpy, nvpz, nvp, kyp, kzp
      integer kxyp, kyzp, kyzmx, kzyp, kyb, kzb
      integer kblok, jblok, lblok, mblok, jkbmx, mlbmx, jkblok, mlblok
      integer nmyz, nyz, nmx, nxyz, nmxh, nxyzh, nxhyz
      integer nxv, nxvh, nyv, nzv, ndv, nvrp, nblok, ndim, nloop
c indnvpy/indnvpz = exponent determining number of real or virtual
c processors in y/z, indnvpy must be <= indy, indnvpz must be <= indz
c mshare = (0,1) = (no,yes) architecture is shared memory
c     parameter( indx =   7, indy =   8, indz =   6)
c     parameter( indx =   6, indy =   5, indz =   7)
      parameter( indx =   5, indy =   5, indz =   6)
c     parameter( indx =  11, indy =  11, indz =  11)
      parameter( indnvpy =   1, indnvpz =   1, mshare =   0)
c vector dimension, number of time iterations
      parameter(ndim=3,nloop=10)
      parameter(nx=2**indx,ny=2**indy,nz=2**indz,nxh=nx/2)
      parameter(nvpy=2**indnvpy,nvpz=2**indnvpz,nvp=nvpy*nvpz)
      parameter(kyp=(ny-1)/nvpy+1,kzp=(nz-1)/nvpz+1)
      parameter(kxyp=(nxh-1)/nvpy+1,kyzp=(ny-1)/nvpz+1)
c kzyp = maximum(kyzp,kyp)
      parameter(kyzmx=kyzp*(kyp/kyzp)+kyp*(kyzp/kyp))
      parameter(kzyp=kyzmx/(2-kyzp/kyzmx-kyp/kyzmx))
      parameter(kyb=ny/kyp,kzb=nz/kzp)
      parameter(kblok=1+mshare*(ny/kyp-1))
      parameter(jblok=1+mshare*(nxh/kxyp-1))
      parameter(lblok=1+mshare*(nz/kzp-1))
      parameter(mblok=1+mshare*(ny/kyzp-1))
      parameter(nblok=kblok*lblok)
c jkblok = maximum(jblok,kblok)
      parameter(jkbmx=jblok*(kblok/jblok)+kblok*(jblok/kblok))
      parameter(jkblok=jkbmx/(2-jblok/jkbmx-kblok/jkbmx))
c mlblok = maximum(mblok,lblok)
      parameter(mlbmx=mblok*(lblok/mblok)+lblok*(mblok/lblok))
      parameter(mlblok=mlbmx/(2-mblok/mlbmx-lblok/mlbmx))
      parameter(nmyz=ny*(nz/ny)+nz*(ny/nz))
      parameter(nyz=nmyz/(2-ny/nmyz-nz/nmyz))
      parameter(nmx=nx*(nyz/nx)+nyz*(nx/nyz))
      parameter(nxyz=nmx/(2-nx/nmx-nyz/nmx))
      parameter(nmxh=nxh*(nyz/nxh)+nyz*(nxh/nyz))
      parameter(nxhyz=nmxh/(2-nxh/nmxh-nyz/nmxh))
      parameter(nxv=nx+2,nxvh=nxv/2,nyv=ny+2,nzv=nz+2,nxyzh=nxyz/2)
      parameter(ndv=256,nvrp=(ndv-1)/nvp+1)
c
      integer ntpose, idproc, kstrt, ks, js, isign
      integer i, j, k, l, m,  my, mz, moff, nd
      real time, wtime, timee, timep, ttp, tsp, epsmax, eps, hh
      real sum1, sum2
      double precision dtime, etime
      dimension time(2), wtime(2), timee(2*nloop), timep(2*nloop)
c
      integer mixup
      complex sct
      double precision vran
      dimension mixup(nxhyz), sct(nxyzh)
      dimension vran(nvrp,nblok)
c
      real f, e
      complex g, h, bs, br
      dimension f(nxv,kyp,kzp,kblok*lblok), e(nxv,kyp,kzp,kblok*lblok)
      dimension g(nyv,kxyp,kzp,jblok*lblok)
      dimension h(nzv,kxyp,kyzp,jblok*mblok)
      dimension bs(kxyp,kzyp,kzp,jkblok*lblok)
      dimension br(kxyp,kzyp,kzp,jblok*mlblok)
c
      real fn, en
      complex gn, hn, bsn, brn
      dimension fn(ndim,nxv,kyp,kzp,kblok*lblok)
      dimension en(ndim,nxv,kyp,kzp,kblok*lblok)
      dimension gn(ndim,nyv,kxyp,kzp,jblok*lblok)
      dimension hn(ndim,nzv,kxyp,kyzp,jblok*mblok)
      dimension bsn(ndim,kxyp,kzyp,kzp,jkblok*lblok)
      dimension brn(ndim,kxyp,kzyp,kzp,jblok*mlblok)
c
c ntpose = (0,1) = (no,yes) input, output data are transposed in PFFT32R
      data ntpose /1/
c initialize for parallel processing
      call PPINIT0(idproc,nvp)
      kstrt = idproc + 1
      ks = (kstrt - 1)/kyb
      js = kstrt - kyb*ks - 2
      ks = ks - 1
c prepare fft tables
      call WPFFT32RINIT(mixup,sct,indx,indy,indz,nxhyz,nxyzh)
c generate test data
      call PTIMERA(-1,time,dtime)
      if (kstrt.gt.(kyb*kzb)) go to 60
      do 50 mz = 1, lblok
      moff = kblok*(mz - 1)
      do 40 my = 1, kblok
      m = my + moff
      do 30 l = 1, kzp
      do 20 k = 1, kyp
      do 10 j = 1, nx
      call pranorm(vran,kstrt,nvp,nvrp,nd,nvrp,nblok)
      hh = vran(1,m)
      f(j,k,l,m) = hh
      e(j,k,l,m) = hh
   10 continue
   20 continue
   30 continue
   40 continue
   50 continue
   60 continue
      call PTIMERA(1,time,dtime)
      if (idproc==0) then
         write (71,*) 'indx,indy,indz=',indx,indy,indz
         write (71,*) 'nvp,npvy,nvpz=',nvp,nvpy,nvpz
         write (71,*) 'init max/min time=',time(1),time(2), 'sec'
      endif
c timings
      tsp = 0.0
      call PTIMERA(-1,time,dtime)
      do 70 i = 1, nloop
      call PWTIMERA(-1,timee(2*i-1),etime)
c scalar real to complex fourier space transform
      isign = -1
      call WPFFT32R(f,g,h,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,ind
     1z,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxyp,kyp,kyzp,kzp,kzyp,jblo
     2k,kblok,jkblok,lblok,mblok,mlblok,nxhyz,nxyzh)
      call PWTIMERA(1,timee(2*i-1),etime)
      tsp = tsp + ttp
      timep(2*i-1) = ttp
c fourier space to real transform
      isign = 1
      call WPFFT32R(f,g,h,bs,br,isign,ntpose,mixup,sct,ttp,indx,indy,ind
     1z,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxyp,kyp,kyzp,kzp,kzyp,jblo
     2k,kblok,jkblok,lblok,mblok,mlblok,nxhyz,nxyzh)
      call PWTIMERA(1,timee(2*i),etime)
      tsp = tsp + ttp
      timep(2*i) = ttp
   70 continue
      call PTIMERA(1,time,dtime)
      if (idproc==0) then
         write (71,*) 'scalar max/min time=',time(1),time(2), 'sec'
      endif
      time(1) = tsp
      time(2) = -tsp
      call PMAX(time,wtime,2,1)
      time(2) = -time(2)
      if (idproc==0) then
         write (71,*) 'transpose max/min time=', time(1), time(2), 'sec'
         sum1 = 0.0
         sum2 = 0.0
         do 75 i = 1, 2*nloop
         write (71,*) 'node 0 total, transpose time=',timee(i),timep(i)
         sum1 = sum1 + timee(i)
         sum2 = sum2 + timep(i)
   75    continue
         write (71,*) 'node 0 scalar sum totals=', sum1, sum2
      endif
c verify real scalar data
      epsmax = 0.
      if (kstrt.gt.(kyb*kzb)) go to 130
      do 120 mz = 1, lblok
      moff = kblok*(mz - 1)
      do 110 my = 1, kblok
      m = my + moff
      do 100 l = 1, kzp
      do 90 k = 1, kyp
      do 80 j = 1, nx
      eps = abs(f(j,k,l,m) - e(j,k,l,m))
      if (eps.gt.epsmax) then
c        write (71,*) j,k,l,f(j,k,l,m),e(j,k,l,m),eps
         epsmax = eps
      endif
   80 continue
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
      if (idproc==0) then
         write (71,*) 'local real epsmax=',epsmax
      endif
      call PMAX(epsmax,eps,1,1)
      if (idproc==0) then
         write (71,*) 'global real epsmax=',epsmax
      endif
c generate multi-dimensional test data
      call PTIMERA(-1,time,dtime)
      if (kstrt.gt.(kyb*kzb)) go to 200
      do 190 mz = 1, lblok
      moff = kblok*(mz - 1)
      do 180 my = 1, kblok
      m = my + moff
      do 170 l = 1, kzp
      do 160 k = 1, kyp
      do 150 j = 1, nx
      do 140 i = 1, ndim
      call pranorm(vran,kstrt,nvp,nvrp,nd,nvrp,nblok)
      hh = vran(1,m)
      fn(i,j,k,l,m) = hh
      en(i,j,k,l,m) = hh
  140 continue
  150 continue
  160 continue
  170 continue
  180 continue
  190 continue
  200 continue
      call PTIMERA(1,time,dtime)
      if (idproc==0) then
         write (71,*) 'multi init time=',time(1),time(2), 'sec'
      endif
c timings
      tsp = 0.0
      call PTIMERA(-1,time,dtime)
      do 210 i = 1, nloop
      call PWTIMERA(-1,timee(2*i-1),etime)
c vector real to complex fourier space transform
      isign = -1
      call WPFFT32R3(fn,gn,hn,bsn,brn,isign,ntpose,mixup,sct,ttp,indx,in
     1dy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxyp,kyp,kyzp,kzp,kzy
     2p,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyz,nxyzh)
      call PWTIMERA(1,timee(2*i),etime)
      tsp = tsp + ttp
      timep(2*i-1) = ttp
c fourier space to real transform
      isign = 1
      call WPFFT32R3(fn,gn,hn,bsn,brn,isign,ntpose,mixup,sct,ttp,indx,in
     1dy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxyp,kyp,kyzp,kzp,kzy
     2p,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyz,nxyzh)
      call PWTIMERA(1,timee(2*i),etime)
      tsp = tsp + ttp
      timep(2*i) = ttp
  210 continue
      call PTIMERA(1,time,dtime)
      if (idproc==0) then
         write (71,*) 'vector max/min time=',time(1),time(2), 'sec'
      endif
      time(1) = tsp
      time(2) = -tsp
      call PMAX(time,wtime,2,1)
      time(2) = -time(2)
      if (idproc==0) then
         write (71,*) 'transpose max/min time=',time(1),time(2), 'sec'
         sum1 = 0.0
         sum2 = 0.0
         do 215 i = 1, 2*nloop
         write (71,*) 'node 0 total, transpose time=',timee(i), timep(i)
         sum1 = sum1 + timee(i)
         sum2 = sum2 + timep(i)
  215    continue
         write (71,*) 'node 0 vector sum totals=', sum1, sum2
      endif
c verify real vector data
      epsmax = 0.
      if (kstrt.gt.(kyb*kzb)) go to 280
      do 270 mz = 1, lblok
      moff = kblok*(mz - 1)
      do 260 my = 1, kblok
      m = my + moff
      do 250 l = 1, kzp
      do 240 k = 1, kyp
      do 230 j = 1, nx
      do 220 i = 1, ndim
      eps = abs(fn(i,j,k,l,m) - en(i,j,k,l,m))
      if (eps.gt.epsmax) then
c        write (71,*) i,j,k,l,fn(i,j,k,l,m),en(i,j,k,l,m),eps
         epsmax = eps
      endif
  220 continue
  230 continue
  240 continue
  250 continue
  260 continue
  270 continue
  280 continue
      if (idproc==0) then
         write (71,*) 'local real epsmax=',epsmax
      endif
      call PMAX(epsmax,eps,1,1)
      if (idproc==0) then
         write (71,*) 'global real epsmax=',epsmax
      endif
c
      call PPEXIT
      stop
      end
c-----------------------------------------------------------------------
      function vresult(prec)
      implicit none
      real prec, vresult
      vresult = prec
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
      subroutine WPFFT32RINIT(mixup,sct,indx,indy,indz,nxhyzd,nxyzhd)
c this subroutine calculates tables needed by a three dimensional
c real to complex fast fourier transform and its inverse.
c input: indx, indy, indz, nxhyzd, nxyzhd
c output: mixup, sct
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c written by viktor k. decyk, ucla
      implicit none
      integer indx, indy, indz, nxhyzd, nxyzhd
      integer mixup
      complex sct
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, ny, nz, nxyz, nxhyz, nxyzh
      integer j, k, lb, ll, jb, it
      real dnxyz, arg
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      ny = 2**indy
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
      do 20 j = 1, nxhyz
      lb = j - 1
      ll = 0
      do 10 k = 1, ndx1yz
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   10 continue
      mixup(j) = ll + 1
   20 continue
c sine/cosine table for the angles 2*n*pi/nxyz
      nxyzh = nxyz/2
      dnxyz = 6.28318530717959/float(nxyz)
      do 30 j = 1, nxyzh
      arg = dnxyz*float(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFFT32R(f,g,h,bs,br,isign,ntpose,mixup,sct,ttp,indx,in
     1dy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd
     2,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd)
c wrapper function for real to complex fft
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd, kzyp
      integer jblok, kblok, jkblok, lblok, mblok, mlblok, nxhyzd, nxyzhd
      integer mixup
      real ttp
      complex f, g, h, bs, br, sct
      dimension f(nxvh,kypd,kzpd,kblok*lblok)
      dimension g(nyv,kxypd,kzpd,jblok*lblok)
      dimension h(nzv,kxypd,kyzpd,jblok*mblok)
      dimension bs(kxyp,kzyp,kzp,jkblok*lblok)
      dimension br(kxyp,kzyp,kzp,jblok*mlblok)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer nxh, ny, nz, kypi, kxypi
      real tp, tf
      double precision dtime
      data kypi, kxypi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
      nz = 2**indz
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call PFFT32RXX(f,isign,mixup,sct,indx,indy,indz,kstrt,kypi,kyp,
     1nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,nxhyzd,nxyzhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PTPOS3A(f,g,bs,br,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,kx
     1ypd,kypd,kzpd,jblok,kblok,lblok)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PFFT32RXY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kxy
     1p,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,nxhyzd,nxyzhd)
c transpose g array to h
         call PWTIMERA(-1,tp,dtime)
         call PTPOS3B(g,h,bs,br,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kx
     1ypd,kyzpd,kzpd,jblok,mblok,lblok)
         call PWTIMERA(1,tp,dtime)
c perform z fft
         call PFFT32RXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kxy
     1p,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,nxhyzd,nxyzhd)
c transpose h array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PTPOS3B(h,g,br,bs,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp
     1,kxypd,kzpd,kyzpd,jblok,lblok,mblok)
            call PTPOS3A(g,f,br,bs,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp
     1,kypd,kxypd,kzpd,kblok,jblok,lblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to h
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PTPOS3A(f,g,bs,br,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp
     1,kxypd,kypd,kzpd,jblok,kblok,lblok)
            call PTPOS3B(g,h,bs,br,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp
     1,kxypd,kyzpd,kzpd,jblok,mblok,lblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform z fft
         call PFFT32RXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kxy
     1p,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,nxhyzd,nxyzhd)
c transpose h array to g
         call PWTIMERA(-1,tp,dtime)
         call PTPOS3B(h,g,br,bs,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kx
     1ypd,kzpd,kyzpd,jblok,lblok,mblok)
         call PWTIMERA(1,tp,dtime)
c perform y fft
         call PFFT32RXY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kxy
     1p,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,nxhyzd,nxyzhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PTPOS3A(g,f,br,bs,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,ky
     1pd,kxypd,kzpd,kblok,jblok,lblok)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PFFT32RXX(f,isign,mixup,sct,indx,indy,indz,kstrt,kypi,kyp,
     1nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,nxhyzd,nxyzhd)
      endif
      ttp = ttp + tp
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT32RXX(f,isign,mixup,sct,indx,indy,indz,kstrt,kypi,k
     1ypp,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,nxhyzd,nxyzhd)
c this subroutine performs the x part of a three dimensional real to
c complex fast fourier transform and its inverse, for a subset of y,
c using complex arithmetic, for data which is distributed in blocks,
c with 2D spatial decomposition
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform is performed
c f(n,k,i,id) = (1/nx*ny*nz)*sum(f(j,k,i,id)*exp(-sqrt(-1)*2pi*n*j/nx)
c if isign = 1, a forward fourier transform is performed
c f(n,k,i,id) = (sum(f(j,k,i,id)*exp(sqrt(-1)*2pi*n*j/nx)
c kstrt = starting data block number
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f
c kyp/kzp = number of data values per block in y/z
c kypd = second dimension of f
c kzpd = third dimension of f
c kblok/lblok = number of data blocks in y/z
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c the real data is stored in a complex array of length nx/2, ny, nz
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c f(j,k,l,i) = mode j-1,kk-1,ll-1, where kk = k + kyp*(ix - 1),
c ll = l + kzp*(iz - 1), and i = iy + nprocy*(iz - 1)
c 1 <= j <= nx/2, 1 <= kk <= ny, and 1 <= ll <= nz, except for
c f(1,k,l,i) = mode nx/2,kk-1,ll-1, where
c ny/2+2 <= kk <= ny, 1 <= ll <= nz, and
c f(1,1,l,i) = mode nx/2,0,ll-1
c f(1,ny/2+1,l,i) = mode nx/2,ny/2,ll-1, where nz/2+2 <= ll <= nz, and
c imaginary part of f(1,1,1,1) = real part of mode nx/2,0,0
c imaginary part of f(1,ny/2+1,1,1) = real part of mode nx/2,ny/2,0
c imaginary part of f(1,1,1,(nz/2)/kzp+1) =
c    real part of mode nx/2,0,nz/2
c imaginary part of f(1,ny/2+1,1,(nz/2)/kzp+1) =
c    real part of mode nx/2,ny/2,nz/2
c using jpl storage convention
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, kypi, kypp, nxvh, mixup
      integer kyp, kzp, kypd, kzpd, kblok, lblok, nxhyzd, nxyzhd
      complex f, sct
      dimension f(nxvh,kypd,kzpd,kblok*lblok)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nz, nxyz, nxhyz
      integer j, k, l, i, m, n, ns, ns2, km, kmr, k1, k2, j1, j2, klblok
      integer nrx, nry, kyb, kzb, kypt
      real ani
      complex s, t, t1
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kyb = ny/kyp
      kzb = nz/kzp
      kypt = kypi + kypp - 1
      if (kstrt.gt.(kyb*kzb)) return
      klblok = kblok*lblok
      if (isign.gt.0) go to 120
c inverse fourier transform
      ani = 0.5/(float(nx)*float(ny)*float(nz))
      do 110 m = 1, klblok
      do 100 n = 1, kzp
      nrx = nxhyz/nxh
c bit-reverse array elements in x
      do 20 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 20
      do 10 i = kypi, kypt
      t = f(j1,i,n,m)
      f(j1,i,n,m) = f(j,i,n,m)
      f(j,i,n,m) = t
   10 continue
   20 continue
c first transform in x
      nrx = nxyz/nxh
      do 60 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 50 i = kypi, kypt
      do 40 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 30 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t = s*f(j2,i,n,m)
      f(j2,i,n,m) = f(j1,i,n,m) - t
      f(j1,i,n,m) = f(j1,i,n,m) + t
   30 continue
   40 continue
   50 continue
   60 continue
c unscramble coefficients and normalize
      kmr = nxyz/nx
      nry = nxhyz/ny
      do 80 k = kypi, kypt
      do 70 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      t = conjg(f(nxh2-j,k,n,m))
      s = f(j,k,n,m) + t
      t = (f(j,k,n,m) - t)*t1
      f(j,k,n,m) = ani*(s + t)
      f(nxh2-j,k,n,m) = ani*conjg(s - t)
   70 continue
   80 continue
      do 90 k = kypi, kypt
      f(nxhh+1,k,n,m) = 2.*ani*conjg(f(nxhh+1,k,n,m))
      f(1,k,n,m) = 2.*ani*cmplx(real(f(1,k,n,m)) + aimag(f(1,k,n,m)),rea
     1l(f(1,k,n,m)) - aimag(f(1,k,n,m)))
   90 continue
  100 continue
  110 continue
      return
c forward fourier transform
  120 do 230 m = 1, klblok
      do 220 n = 1, kzp
c scramble coefficients
      kmr = nxyz/nx
      do 140 k = kypi, kypt
      do 130 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      t = conjg(f(nxh2-j,k,n,m))
      s = f(j,k,n,m) + t
      t = (f(j,k,n,m) - t)*t1
      f(j,k,n,m) = s + t
      f(nxh2-j,k,n,m) = conjg(s - t)
  130 continue
  140 continue
      do 150 k = kypi, kypt
      f(nxhh+1,k,n,m) = 2.*conjg(f(nxhh+1,k,n,m))
      f(1,k,n,m) = cmplx(real(f(1,k,n,m)) + aimag(f(1,k,n,m)),real(f(1,k
     1,n,m)) - aimag(f(1,k,n,m)))
  150 continue
      nrx = nxhyz/nxh
c bit-reverse array elements in x
      do 170 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 170
      do 160 i = kypi, kypt
      t = f(j1,i,n,m)
      f(j1,i,n,m) = f(j,i,n,m)
      f(j,i,n,m) = t
  160 continue
  170 continue
c finally transform in x
      nrx = nxyz/nxh
      do 210 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 200 i = kypi, kypt
      do 190 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 180 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t = s*f(j2,i,n,m)
      f(j2,i,n,m) = f(j1,i,n,m) - t
      f(j1,i,n,m) = f(j1,i,n,m) + t
  180 continue
  190 continue
  200 continue
  210 continue
  220 continue
  230 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT32RXY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,
     1kxypp,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,nxhyzd,nxyzhd)
c this subroutine performs the y part of a three dimensional real to
c complex fast fourier transform and its inverse, for a subset of x,
c using complex arithmetic, for data which is distributed in blocks,
c with 2D spatial decomposition
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform is performed
c g(m,j,i,ld) = sum(g(k,j,i,id)*exp(-sqrt(-1)*2pi*mm*kk/ny)
c where mm = m + kyp*(ldy - 1) and kk = k + kyp*(idy - 1),
c if isign = 1, a forward fourier transform is performed
c g(m,j,i,id) = (sum(g(k,j,i,id)*exp(sqrt(-1)*2pi*mm*kk/ny)
c kstrt = starting data block number
c kxypi = initial x index used
c kxypp = number of x indices used
c nyv = first dimension of g
c kxyp/kzp = number of data values per block in x/z
c kxypd = second dimension of g
c kzpd = third dimension of g
c jblok/lblok = number of data blocks in x/z
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c the real data is stored in a complex array of length nx/2, ny, nz
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c f(j,k,l,i) = mode j-1,kk-1,ll-1, where kk = k + kyp*(ix - 1),
c ll = l + kzp*(iz - 1), and i = iy + nprocy*(iz - 1)
c 1 <= j <= nx/2, 1 <= kk <= ny, and 1 <= ll <= nz, except for
c f(1,k,l,i) = mode nx/2,kk-1,ll-1, where
c ny/2+2 <= kk <= ny, 1 <= ll <= nz, and
c f(1,1,l,i) = mode nx/2,0,ll-1
c f(1,ny/2+1,l,i) = mode nx/2,ny/2,ll-1, where nz/2+2 <= ll <= nz, and
c imaginary part of f(1,1,1,1) = real part of mode nx/2,0,0
c imaginary part of f(1,ny/2+1,1,1) = real part of mode nx/2,ny/2,0
c imaginary part of f(1,1,1,(nz/2)/kzp+1) =
c    real part of mode nx/2,0,nz/2
c imaginary part of f(1,ny/2+1,1,(nz/2)/kzp+1) =
c    real part of mode nx/2,ny/2,nz/2
c using jpl storage convention
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, kxypi, kxypp, nyv, mixup
      integer kxyp, kzp, kxypd, kzpd, jblok, lblok, nxhyzd, nxyzhd
      complex g, sct
      dimension g(nyv,kxypd,kzpd,jblok*lblok)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh, ny2, nz, nxyz, nxhyz
      integer j, k, l, i, m, n, ns, ns2, km, kmr, k1, k2, j1, j2, jlblok
      integer js, ks, nry, kzb, kxb, kxypt, mx, mz, moff
      complex s, t
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kxb = nxh/kxyp
      kzb = nz/kzp
      kxypt = kxypi + kxypp - 1
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      if (kstrt.gt.(kxb*kzb)) return
      jlblok = jblok*lblok
      if (isign.gt.0) go to 130
c inverse fourier transform
      do 80 m = 1, jlblok
      do 70 n = 1, kzp
      nry = nxhyz/ny
c bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 i = kxypi, kxypt
      t = g(k1,i,n,m)
      g(k1,i,n,m) = g(k,i,n,m)
      g(k,i,n,m) = t
   10 continue
   20 continue
c then transform in y
      nry = nxyz/ny
      do 60 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 50 i = kxypi, kxypt
      do 40 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 30 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t = s*g(j2,i,n,m)
      g(j2,i,n,m) = g(j1,i,n,m) - t
      g(j1,i,n,m) = g(j1,i,n,m) + t
   30 continue
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble modes kx = 0, nx/2
      do 120 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 110 mx = 1, jblok
      if ((mx+js).gt.0) go to 110
      m = mx + moff
      if (kxypi.eq.1) then
         do 100 n = 1, kzp
         do 90 k = 2, nyh
         s = g(ny2-k,1,n,m)
         g(ny2-k,1,n,m) = .5*cmplx(aimag(g(k,1,n,m) + s),real(g(k,1,n,m)
     1 - s))
         g(k,1,n,m) = .5*cmplx(real(g(k,1,n,m) + s),aimag(g(k,1,n,m) - s
     1))
   90    continue
  100    continue
      endif
  110 continue
  120 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
  130 do 170 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 160 mx = 1, jblok
      if ((mx+js).gt.0) go to 160
      m = mx + moff
      if (kxypi.eq.1) then
         do 150 n = 1, kzp
         do 140 k = 2, nyh
         s = cmplx(aimag(g(ny2-k,1,n,m)),real(g(ny2-k,1,n,m)))
         g(ny2-k,1,n,m) = conjg(g(k,1,n,m) - s)
         g(k,1,n,m) = g(k,1,n,m) + s
  140    continue
  150    continue
      endif
  160 continue
  170 continue
      do 250 m = 1, jlblok
      do 240 n = 1, kzp
      nry = nxhyz/ny
c bit-reverse array elements in y
      do 190 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 190
      do 180 i = kxypi, kxypt
      t = g(k1,i,n,m)
      g(k1,i,n,m) = g(k,i,n,m)
      g(k,i,n,m) = t
  180 continue
  190 continue
c then transform in y
      nry = nxyz/ny
      do 230 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 220 i = kxypi, kxypt
      do 210 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 200 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t = s*g(j2,i,n,m)
      g(j2,i,n,m) = g(j1,i,n,m) - t
      g(j1,i,n,m) = g(j1,i,n,m) + t
  200 continue
  210 continue
  220 continue
  230 continue
  240 continue
  250 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT32RXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,
     1kxypp,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,nxhyzd,nxyzhd)
c this subroutine performs the z part of a three dimensional real to
c complex fast fourier transform and its inverse, for a subset of x,
c using complex arithmetic, for data which is distributed in blocks,
c with 2D spatial decomposition
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform is performed
c h(l,n,m,ld) = sum(h(i,j,k,id)*exp(sqrt(-1)*2pi*ll*ii/nz)
c where ll = l + kzp*(ldz - 1) and ii = i + kzp*(idz - 1),
c if isign = 1, a forward fourier transform is performed
c h(l,n,m,ld) = (sum(h(i,j,k,id)*exp(sqrt(-1)*2pi*ll*ii/nz))
c kstrt = starting data block number
c kxypi = initial x index used
c kxypp = number of x indices used
c nzv = first dimension of h
c kxyp/kyzp = number of data values per block in x/y
c kxypd = second dimension of h
c kyzpd = third dimension of h
c jblok/mblok = number of data blocks in x/y
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c the real data is stored in a complex array of length nx/2, ny, nz
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c f(j,k,l,i) = mode j-1,kk-1,ll-1, where kk = k + kyp*(ix - 1),
c ll = l + kzp*(iz - 1), and i = iy + nprocy*(iz - 1)
c 1 <= j <= nx/2, 1 <= kk <= ny, and 1 <= ll <= nz, except for
c f(1,k,l,i) = mode nx/2,kk-1,ll-1, where
c ny/2+2 <= kk <= ny, 1 <= ll <= nz, and
c f(1,1,l,i) = mode nx/2,0,ll-1
c f(1,ny/2+1,l,i) = mode nx/2,ny/2,ll-1, where nz/2+2 <= ll <= nz, and
c imaginary part of f(1,1,1,1) = real part of mode nx/2,0,0
c imaginary part of f(1,ny/2+1,1,1) = real part of mode nx/2,ny/2,0
c imaginary part of f(1,1,1,(nz/2)/kzp+1) =
c    real part of mode nx/2,0,nz/2
c imaginary part of f(1,ny/2+1,1,(nz/2)/kzp+1) =
c    real part of mode nx/2,ny/2,nz/2
c using jpl storage convention
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, kxypi, kxypp, nzv, mixup
      integer kxyp, kyzp, kxypd, kyzpd, jblok, mblok, nxhyzd, nxyzhd
      complex h, sct
      dimension h(nzv,kxypd,kyzpd,jblok*mblok)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh, nz, nzh, nz2, nxyz, nxhyz
      integer j, k, l, i, m, n, ns, ns2, km, kmr, k1, k2, j1, j2, jmblok
      integer l1, js, ks, nrz, kxb, kyzb, kxypt, mx, my, moff
      complex s, t
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      nz = 2**indz
      nzh = nz/2
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      kxypt = kxypi + kxypp - 1
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      if (kstrt.gt.(kxb*kyzb)) return
      jmblok = jblok*mblok
      if (isign.gt.0) go to 140
c inverse fourier transform
      do 80 m = 1, jmblok
      do 70 n = 1, kyzp
      nrz = nxhyz/nz
c bit-reverse array elements in z
      do 20 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 20
      do 10 i = kxypi, kxypt
      t = h(l1,i,n,m)
      h(l1,i,n,m) = h(l,i,n,m)
      h(l,i,n,m) = t
   10 continue
   20 continue
c finally transform in z
      nrz = nxyz/nz
      do 60 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 50 i = kxypi, kxypt
      do 40 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 30 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t = s*h(j2,i,n,m)
      h(j2,i,n,m) = h(j1,i,n,m) - t
      h(j1,i,n,m) = h(j1,i,n,m) + t
   30 continue
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble modes kx = 0, nx/2
      do 130 my = 1, mblok
      moff = jblok*(my - 1)
      do 120 mx = 1, jblok
      if ((mx+js).gt.0) go to 120
      m = mx + moff
      if ((my+ks).eq.0) then
         if (kxypi.eq.1) then
            do 90 n = 2, nzh
            s = h(nz2-n,1,1,m)
            h(nz2-n,1,1,m) = .5*cmplx(aimag(h(n,1,1,m) + s),real(h(n,1,1
     1,m) - s))
            h(n,1,1,m) = .5*cmplx(real(h(n,1,1,m) + s),aimag(h(n,1,1,m) 
     1- s))
   90       continue
         endif
      endif
      kyzb = nyh/kyzp
      if ((my+ks).eq.kyzb) then
         k1 = nyh - kyzb*kyzp + 1
         if (kxypi.eq.1) then
            do 100 n = 2, nzh
            s = h(nz2-n,1,k1,m)
            h(nz2-n,1,k1,m) = .5*cmplx(aimag(h(n,1,k1,m) + s),real(h(n,1
     1,k1,m) - s))
            h(n,1,k1,m) = .5*cmplx(real(h(n,1,k1,m) + s),aimag(h(n,1,k1,
     1m) - s))
  100       continue
        endif
      endif
  120 continue
  130 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
  140 do 180 my = 1, mblok
      moff = jblok*(my - 1)
      do 170 mx = 1, jblok
      if ((mx+js).gt.0) go to 170
      m = mx + moff
      if ((my+ks).eq.0) then
         if (kxypi.eq.1) then
            do 150 n = 2, nzh
            s = cmplx(aimag(h(nz2-n,1,1,m)),real(h(nz2-n,1,1,m)))
            h(nz2-n,1,1,m) = conjg(h(n,1,1,m) - s)
            h(n,1,1,m) = h(n,1,1,m) + s
  150       continue
         endif
      endif
      kyzb = nyh/kyzp
      if ((my+ks).eq.kyzb) then
         k1 = nyh - kyzb*kyzp + 1
         if (kxypi.eq.1) then
            do 160 n = 2, nzh
            s = cmplx(aimag(h(nz2-n,1,k1,m)),real(h(nz2-n,1,k1,m)))
            h(nz2-n,1,k1,m) = conjg(h(n,1,k1,m) - s)
            h(n,1,k1,m) = h(n,1,k1,m) + s
  160       continue
         endif
      endif
  170 continue
  180 continue
      do 260 m = 1, jmblok
      do 250 n = 1, kyzp
      nrz = nxhyz/nz
c bit-reverse array elements in z
      do 200 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 200
      do 190 i = kxypi, kxypt
      t = h(l1,i,n,m)
      h(l1,i,n,m) = h(l,i,n,m)
      h(l,i,n,m) = t
  190 continue
  200 continue
c first transform in z
      nrz = nxyz/nz
      do 240 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 230 i = kxypi, kxypt
      do 220 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 210 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t = s*h(j2,i,n,m)
      h(j2,i,n,m) = h(j1,i,n,m) - t
      h(j1,i,n,m) = h(j1,i,n,m) + t
  210 continue
  220 continue
  230 continue
  240 continue
  250 continue
  260 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT32C(f,g,h,bs,br,isign,ntpose,mixup,sct,indx,indy,in
     1dz,kstrt,nxv,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,kzyp,
     2jblok,kblok,jkblok,lblok,mblok,mlblok,nxyzd,nxyzhd)
c this subroutine performs a three dimensional complex to complex fast
c fourier transform and its inverse, using complex arithmetic,
c for data which is distributed in blocks, with 2D spatial decomposition
c for isign = 0, input: isign, indx, indy, indz, kstrt, nxyzd, nxyzhd
c output: mixup, sct
c for isign = (-1,1), input: all, output: f, g, bs, br
c approximate flop count: 5*N*log2(N)/nvp
c where N = nx*ny*nz, and nvp = number of procs
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c ntpose = (0,1) = (no,yes) input, output data are transposed
c if isign = 0, the fft tables are prepared
c if isign = -1, an inverse fourier transform is performed
c if ntpose = 0, f is the input and output array, g, h are scratch arrays
c f(n,m,l,ld) = (1/nx*ny*nz)*sum(f(j,k,i,id)*exp(-sqrt(-1)*2pi*n*j/nx)*
c       exp(-sqrt(-1)*2pi*mm*kk/ny)*exp(-sqrt(-1)*2pi*ll*ii/nz))
c where ll = l + kzp*(ldz - 1) and ii = i + kzp*(idz - 1),
c and mm = m + kyp*(ldy - 1) and kk = k + kyp*(idy - 1),
c and ld = ldy + nprocy*(ldz - 1), id = idy + nprocy*(idz - 1),
c and nprocy = number of processors in y.
c if ntpose = 1, f is the input and h is the output, and g is scratch
c h(l,n,m,nd) = (1/nx*ny*nz)*sum(f(j,k,i,id)*exp(-sqrt(-1)*2pi*nn*j/nx)*
c       exp(-sqrt(-1)*2pi*mm*kk/ny)*exp(-sqrt(-1)*2pi*l*ii/nz))
c where nn = n + kxyp*(ndy - 1) and ii = i + kzp*(idz - 1),
c and mm = m + kyzp*(mdz - 1) and kk = k + kyp*(idy - 1),
c and nd = ndy + nprocy*(ndz - 1), and id = idy + nprocy*(idz - 1)
c if isign = 1, a forward fourier transform is performed
c if ntpose = 0, f is the input and output array, g, h are scratch arrays
c f(n,m,l,ld) = (sum(f(j,k,i,id)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*mm*kk/ny)*exp(sqrt(-1)*2pi*ll*ii/nz))
c if ntpose = 1, h is the input and f is the output, and g is scratch
c f(j,k,i,id) = sum(h(l,n,m,nd)*exp(sqrt(-1)*2pi*nn*j/nx)*
c       exp(sqrt(-1)*2pi*mm*kk/ny)*exp(sqrt(-1)*2pi*l*ii/nz))
c bs, br = scratch arrays
c kstrt = starting data block number
c nxv/nyv/nzv = first dimension of f/g/h
c kypd/kxypd = second dimension of f/g,h
c kzpd/kyzpd = third dimension of f,g/h
c kxyp/kyp,kyzp/kzp = number of data values per block in x/y/z
c kzyp = maximum(kyzp,kyp)
c jblok/kblok,mblok/lblok = number of data blocks in x/y/z
c jkblok = maximum(jblok,kblok)
c mlblok = maximum(mblok,lblok)
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxyz = maximum of (nx,ny,nz)
c nxyzh = one half of maximum of (nx,ny,nz)
c written by viktor k. decyk, ucla
c parallel version
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nxv, nyv, nzv
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd, kzyp
      integer jblok, kblok, jkblok, lblok, mblok, mlblok, nxyzd, nxyzhd
      integer mixup
      complex f, g, h, bs, br, sct
      dimension f(nxv,kypd,kzpd,kblok*lblok)
      dimension g(nyv,kxypd,kzpd,jblok*lblok)
      dimension h(nzv,kxypd,kyzpd,jblok*mblok)
      dimension bs(kxyp,kzyp,kzp,jkblok*lblok)
      dimension br(kxyp,kzyp,kzp,jblok*mlblok)
      dimension mixup(nxyzd), sct(nxyzhd)
c local data
      integer indxyz, nx, nxh, ny, nyh, nz, nzh, nxyz, nxyzh, i, l, m, n
      integer j, j1, j2, k, k1, k2, l1, lb, ll, jb, it, ns, ns2, km, kmr
      integer nrx, nry, nrz, kyb, kzb, kxb, kyzb, klblok, jlblok, jmblok
      real dnxyz, arg, ani
      complex s, t
      indxyz = max0(indx,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      nz = 2**indz
      nzh = nz/2
      nxyz = 2**indxyz
      if (isign) 50, 10, 430
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
   10 do 30 j = 1, nxyz
      lb = j - 1
      ll = 0
      do 20 k = 1, indxyz
      jb = lb/2
      it = lb - 2*jb
      lb = jb
      ll = 2*ll + it
   20 continue
      mixup(j) = ll + 1
   30 continue
c sine/cosine table for the angles 2*n*pi/nxyz
      nxyzh = nxyz/2
      dnxyz = 6.28318530717959/float(nxyz)
      do 40 j = 1, nxyzh
      arg = dnxyz*float(j - 1)
      sct(j) = cmplx(cos(arg),-sin(arg))
   40 continue
      return
c inverse fourier transform
   50 kyb = ny/kyp
      kzb = nz/kzp
      if (kstrt.gt.(kyb*kzb)) go to 160
      klblok = kblok*lblok
      nrx = nxyz/nx
      nry = nxyz/ny
      do 90 m = 1, klblok
c bit-reverse array elements in x
      do 80 n = 1, kzp
      do 70 j = 1, nx
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 70
      do 60 i = 1, kyp
      t = f(j1,i,n,m)
      f(j1,i,n,m) = f(j,i,n,m)
      f(j,i,n,m) = t
   60 continue
   70 continue
   80 continue
   90 continue
c first transform in x
      do 150 l = 1, indx
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxh/ns
      kmr = km*nrx
      do 140 m = 1, klblok
      do 130 n = 1, kzp
      do 120 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 110 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 100 i = 1, kyp
      t = s*f(j2,i,n,m)
      f(j2,i,n,m) = f(j1,i,n,m) - t
      f(j1,i,n,m) = f(j1,i,n,m) + t
  100 continue
  110 continue
  120 continue
  130 continue
  140 continue
  150 continue
c transpose f array to g
  160 call PTPOS3A(f,g,bs,br,nx,ny,nz,kstrt,nxv,nyv,kxyp,kyp,kzp,kxypd,k
     1ypd,kzpd,jblok,kblok,lblok)
      kxb = nx/kxyp
      if (kstrt.gt.(kxb*kzb)) go to 270
      jlblok = jblok*lblok
      do 200 m = 1, jlblok
c bit-reverse array elements in y
      do 190 n = 1, kzp
      do 180 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 180
      do 170 i = 1, kxyp
      t = g(k1,i,n,m)
      g(k1,i,n,m) = g(k,i,n,m)
      g(k,i,n,m) = t
  170 continue
  180 continue
  190 continue
  200 continue
c then transform in y
      do 260 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 250 m = 1, jlblok
      do 240 n = 1, kzp
      do 230 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 220 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 210 i = 1, kxyp
      t = s*g(j2,i,n,m)
      g(j2,i,n,m) = g(j1,i,n,m) - t
      g(j1,i,n,m) = g(j1,i,n,m) + t
  210 continue
  220 continue
  230 continue
  240 continue
  250 continue
  260 continue
c transpose g array to h
  270 call PTPOS3B(g,h,bs,br,nx,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxypd,
     1kyzpd,kzpd,jblok,mblok,lblok)
      kyzb = ny/kyzp
      if (kstrt.gt.(kxb*kyzb)) go to 420
      jmblok = jblok*mblok
      nrz = nxyz/nz
      do 310 m = 1, jmblok
c bit-reverse array elements in z
      do 300 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 300
      do 290 n = 1, kyzp
      do 280 i = 1, kxyp
      t = h(l1,i,n,m)
      h(l1,i,n,m) = h(l,i,n,m)
      h(l,i,n,m) = t
  280 continue
  290 continue
  300 continue
  310 continue
c finally transform in z
      do 370 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 360 m = 1, jmblok
      do 350 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 340 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 330 n = 1, kyzp
      do 320 i = 1, kxyp
      t = s*h(j2,i,n,m)
      h(j2,i,n,m) = h(j1,i,n,m) - t
      h(j1,i,n,m) = h(j1,i,n,m) + t
  320 continue
  330 continue
  340 continue
  350 continue
  360 continue
  370 continue
c normalize result
      ani = 1.0/(float(nx)*float(ny)*float(nz))
      do 410 m = 1, jmblok
      do 400 k = 1, kyzp
      do 390 j = 1, kxyp
      do 380 l = 1, nz
      h(l,j,k,m) = h(l,j,k,m)*ani
  380 continue
  390 continue
  400 continue
  410 continue
c transpose h array to f
  420 if (ntpose.eq.0) then
         call PTPOS3B(h,g,br,bs,nx,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kxy
     1pd,kzpd,kyzpd,jblok,lblok,mblok)
         call PTPOS3A(g,f,br,bs,ny,nx,nz,kstrt,nyv,nxv,kyp,kxyp,kzp,kypd
     1,kxypd,kzpd,kblok,jblok,lblok)
      endif
      return
c forward fourier transform
c transpose f array to h
  430 if (ntpose.eq.0) then
         call PTPOS3A(f,g,bs,br,nx,ny,nz,kstrt,nxv,nyv,kxyp,kyp,kzp,kxyp
     1d,kypd,kzpd,jblok,kblok,lblok)
         call PTPOS3B(g,h,bs,br,nx,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxy
     1pd,kyzpd,kzpd,jblok,mblok,lblok)
      endif
      kxb = nx/kxyp
      kyzb = ny/kyzp
      if (kstrt.gt.(kxb*kyzb)) go to 540
      jmblok = jblok*mblok
      nrz = nxyz/nz
      do 470 m = 1, jmblok
c bit-reverse array elements in z
      do 460 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 460
      do 450 n = 1, kyzp
      do 440 i = 1, kxyp
      t = h(l1,i,n,m)
      h(l1,i,n,m) = h(l,i,n,m)
      h(l,i,n,m) = t
  440 continue
  450 continue
  460 continue
  470 continue
c first transform in z
      do 530 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 520 m = 1, jmblok
      do 510 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 500 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 490 n = 1, kyzp
      do 480 i = 1, kxyp
      t = s*h(j2,i,n,m)
      h(j2,i,n,m) = h(j1,i,n,m) - t
      h(j1,i,n,m) = h(j1,i,n,m) + t
  480 continue
  490 continue
  500 continue
  510 continue
  520 continue
  530 continue
c transpose h array to g
  540 call PTPOS3B(h,g,br,bs,nx,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kxypd,
     1kzpd,kyzpd,jblok,lblok,mblok)
      kzb = nz/kzp
      if (kstrt.gt.(kxb*kzb)) go to 650
      jlblok = jblok*lblok
      nrx = nxyz/nx
      nry = nxyz/ny
      do 580 m = 1, jlblok
c bit-reverse array elements in y
      do 570 n = 1, kzp
      do 560 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 560
      do 550 i = 1, kxyp
      t = g(k1,i,n,m)
      g(k1,i,n,m) = g(k,i,n,m)
      g(k,i,n,m) = t
  550 continue
  560 continue
  570 continue
  580 continue
c then transform in y
      do 640 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 630 m = 1, jlblok
      do 620 n = 1, kzp
      do 610 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 600 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 590 i = 1, kxyp
      t = s*g(j2,i,n,m)
      g(j2,i,n,m) = g(j1,i,n,m) - t
      g(j1,i,n,m) = g(j1,i,n,m) + t
  590 continue
  600 continue
  610 continue
  620 continue
  630 continue
  640 continue
c transpose g array to f
  650 call PTPOS3A(g,f,br,bs,ny,nx,nz,kstrt,nyv,nxv,kyp,kxyp,kzp,kypd,kx
     1ypd,kzpd,kblok,jblok,lblok)
      kyb = ny/kyp
      if (kstrt.gt.(kyb*kzb)) return
      klblok = kblok*lblok
      do 690 m = 1, klblok
c bit-reverse array elements in x
      do 680 n = 1, kzp
      do 670 j = 1, nx
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 670
      do 660 i = 1, kyp
      t = f(j1,i,n,m)
      f(j1,i,n,m) = f(j,i,n,m)
      f(j,i,n,m) = t
  660 continue
  670 continue
  680 continue
  690 continue
c finally transform in x
      do 750 l = 1, indx
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxh/ns
      kmr = km*nrx
      do 740 m = 1, klblok
      do 730 n = 1, kzp
      do 720 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 710 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 700 i = 1, kyp
      t = s*f(j2,i,n,m)
      f(j2,i,n,m) = f(j1,i,n,m) - t
      f(j1,i,n,m) = f(j1,i,n,m) + t
  700 continue
  710 continue
  720 continue
  730 continue
  740 continue
  750 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFFT32R3(f,g,h,bs,br,isign,ntpose,mixup,sct,ttp,indx,i
     1ndy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzp
     2d,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd)
c wrapper function for real to complex fft
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd, kzyp
      integer jblok, kblok, jkblok, lblok, mblok, mlblok, nxhyzd, nxyzhd
      integer mixup
      real ttp
      complex f, g, h, bs, br, sct
      dimension f(3,nxvh,kypd,kzpd,kblok*lblok)
      dimension g(3,nyv,kxypd,kzpd,jblok*lblok)
      dimension h(3,nzv,kxypd,kyzpd,jblok*mblok)
      dimension bs(3,kxyp,kzyp,kzp,jkblok*lblok)
      dimension br(3,kxyp,kzyp,kzp,jblok*mlblok)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer nxh, ny, nz, kypi, kxypi
      real tp, tf
      double precision dtime
      data kypi, kxypi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
      nz = 2**indz
c inverse fourier transform
      if (isign.lt.0) then
c perform x fft
         call PFFT32R3XX(f,isign,mixup,sct,indx,indy,indz,kstrt,kypi,kyp
     1,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,nxhyzd,nxyzhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call P3TPOS3A(f,g,bs,br,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,k
     1xypd,kypd,kzpd,jblok,kblok,lblok)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PFFT32R3XY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1yp,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,nxhyzd,nxyzhd)
c transpose g array to h
         call PWTIMERA(-1,tp,dtime)
         call P3TPOS3B(g,h,bs,br,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,k
     1xypd,kyzpd,kzpd,jblok,mblok,lblok)
         call PWTIMERA(1,tp,dtime)
c perform z fft
         call PFFT32R3XZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1yp,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,nxhyzd,nxyzhd)
c transpose h array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call P3TPOS3B(h,g,br,bs,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyz
     1p,kxypd,kzpd,kyzpd,jblok,lblok,mblok)
            call P3TPOS3A(g,f,br,bs,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kz
     1p,kypd,kxypd,kzpd,kblok,jblok,lblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to h
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call P3TPOS3A(f,g,bs,br,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kz
     1p,kxypd,kypd,kzpd,jblok,kblok,lblok)
            call P3TPOS3B(g,h,bs,br,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kz
     1p,kxypd,kyzpd,kzpd,jblok,mblok,lblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform z fft
         call PFFT32R3XZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1yp,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,nxhyzd,nxyzhd)
c transpose h array to g
         call PWTIMERA(-1,tp,dtime)
         call P3TPOS3B(h,g,br,bs,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,k
     1xypd,kzpd,kyzpd,jblok,lblok,mblok)
         call PWTIMERA(1,tp,dtime)
c perform y fft
         call PFFT32R3XY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1yp,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,nxhyzd,nxyzhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call P3TPOS3A(g,f,br,bs,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,k
     1ypd,kxypd,kzpd,kblok,jblok,lblok)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PFFT32R3XX(f,isign,mixup,sct,indx,indy,indz,kstrt,kypi,kyp
     1,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,nxhyzd,nxyzhd)
      endif
      ttp = ttp + tp
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT32R3XX(f,isign,mixup,sct,indx,indy,indz,kstrt,kypi,
     1kypp,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,nxhyzd,nxyzhd)
c this subroutine performs the x part of 3 three dimensional real to
c complex fast fourier transforms and their inverses, for a subset of y,
c using complex arithmetic, for data which is distributed in blocks,
c with 2D spatial decomposition
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform is performed
c f(1:3,n,k,i,id) = (1/nx*ny*nz)*sum(f(1:3,j,k,i,id)*
c       exp(-sqrt(-1)*2pi*n*j/nx)
c if isign = 1, a forward fourier transform is performed
c f(1:3,n,k,i,id) = (sum(f(1:3,j,k,i,id)*exp(sqrt(-1)*2pi*n*j/nx)
c kstrt = starting data block number
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f
c kyp/kzp = number of data values per block in y/z
c kypd = second dimension of f
c kzpd = third dimension of f
c kblok/lblok = number of data blocks in y/z
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c the real data is stored in a complex array of length nx/2, ny, nz
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c f(j,k,l,i) = mode j-1,kk-1,ll-1, where kk = k + kyp*(ix - 1),
c ll = l + kzp*(iz - 1), and i = iy + nprocy*(iz - 1)
c 1 <= j <= nx/2, 1 <= kk <= ny, and 1 <= ll <= nz, except for
c f(1,k,l,i) = mode nx/2,kk-1,ll-1, where
c ny/2+2 <= kk <= ny, 1 <= ll <= nz, and
c f(1,1,l,i) = mode nx/2,0,ll-1
c f(1,ny/2+1,l,i) = mode nx/2,ny/2,ll-1, where nz/2+2 <= ll <= nz, and
c imaginary part of f(1,1,1,1) = real part of mode nx/2,0,0
c imaginary part of f(1,ny/2+1,1,1) = real part of mode nx/2,ny/2,0
c imaginary part of f(1,1,1,(nz/2)/kzp+1) =
c    real part of mode nx/2,0,nz/2
c imaginary part of f(1,ny/2+1,1,(nz/2)/kzp+1) =
c    real part of mode nx/2,ny/2,nz/2
c using jpl storage convention
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, kypi, kypp, nxvh, mixup
      integer kyp, kzp, kypd, kzpd, kblok, lblok, nxhyzd, nxyzhd
      complex f, sct
      dimension f(3,nxvh,kypd,kzpd,kblok*lblok)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nz, nxyz, nxhyz
      integer j, k, l, i, m, n, ns, ns2, km, kmr, k1, k2, j1, j2, jj
      integer klblok, nrx, nry, kyb, kzb, kypt
      real ani, at1, at2
      complex s, t, t1, t2, t3
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kyb = ny/kyp
      kzb = nz/kzp
      kypt = kypi + kypp - 1
      if (kstrt.gt.(kyb*kzb)) return
      klblok = kblok*lblok
      if (isign.gt.0) go to 160
c inverse fourier transform
      ani = 0.5/(float(nx)*float(ny)*float(nz))
c swap complex components
      do 150 m = 1, klblok
      do 140 n = 1, kzp
      do 20 i = kypi, kypt
      do 10 j = 1, nxh
      at1 = real(f(3,j,i,n,m))
      f(3,j,i,n,m) = cmplx(real(f(2,j,i,n,m)),aimag(f(3,j,i,n,m)))
      at2 = aimag(f(2,j,i,n,m))
      f(2,j,i,n,m) = cmplx(aimag(f(1,j,i,n,m)),at1)
      f(1,j,i,n,m) = cmplx(real(f(1,j,i,n,m)),at2)
   10 continue
   20 continue
      nrx = nxhyz/nxh
c bit-reverse array elements in x
      do 40 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 40
      do 30 i = kypi, kypt
      t1 = f(1,j1,i,n,m)
      t2 = f(2,j1,i,n,m)
      t3 = f(3,j1,i,n,m)
      f(1,j1,i,n,m) = f(1,j,i,n,m)
      f(2,j1,i,n,m) = f(2,j,i,n,m)
      f(3,j1,i,n,m) = f(3,j,i,n,m)
      f(1,j,i,n,m) = t1
      f(2,j,i,n,m) = t2
      f(3,j,i,n,m) = t3
   30 continue
   40 continue
c first transform in x
      nrx = nxyz/nxh
      do 80 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 70 i = kypi, kypt
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t1 = s*f(1,j2,i,n,m)
      t2 = s*f(2,j2,i,n,m)
      t3 = s*f(3,j2,i,n,m)
      f(1,j2,i,n,m) = f(1,j1,i,n,m) - t1
      f(2,j2,i,n,m) = f(2,j1,i,n,m) - t2
      f(3,j2,i,n,m) = f(3,j1,i,n,m) - t3
      f(1,j1,i,n,m) = f(1,j1,i,n,m) + t1
      f(2,j1,i,n,m) = f(2,j1,i,n,m) + t2
      f(3,j1,i,n,m) = f(3,j1,i,n,m) + t3
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
      kmr = nxyz/nx
      nry = nxhyz/ny
      do 110 k = kypi, kypt
      do 100 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 90 jj = 1, 3
      t = conjg(f(jj,nxh2-j,k,n,m))
      s = f(jj,j,k,n,m) + t
      t = (f(jj,j,k,n,m) - t)*t1
      f(jj,j,k,n,m) = ani*(s + t)
      f(jj,nxh2-j,k,n,m) = ani*conjg(s - t)
   90 continue
  100 continue
  110 continue
      do 130 k = kypi, kypt
      do 120 jj = 1, 3
      f(jj,nxhh+1,k,n,m) = 2.*ani*conjg(f(jj,nxhh+1,k,n,m))
      f(jj,1,k,n,m) = 2.*ani*cmplx(real(f(jj,1,k,n,m)) + aimag(f(jj,1,k,
     1n,m)),real(f(jj,1,k,n,m)) - aimag(f(jj,1,k,n,m)))
  120 continue
  130 continue
  140 continue
  150 continue
      return
c forward fourier transform
  160 do 310 m = 1, klblok
      do 300 n = 1, kzp
c scramble coefficients
      kmr = nxyz/nx
      do 190 k = kypi, kypt
      do 180 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 170 jj = 1, 3
      t = conjg(f(jj,nxh2-j,k,n,m))
      s = f(jj,j,k,n,m) + t
      t = (f(jj,j,k,n,m) - t)*t1
      f(jj,j,k,n,m) = s + t
      f(jj,nxh2-j,k,n,m) = conjg(s - t)
  170 continue
  180 continue
  190 continue
      do 210 k = kypi, kypt
      do 200 jj = 1, 3
      f(jj,nxhh+1,k,n,m) = 2.*conjg(f(jj,nxhh+1,k,n,m))
      f(jj,1,k,n,m) = cmplx(real(f(jj,1,k,n,m)) + aimag(f(jj,1,k,n,m)),r
     1eal(f(jj,1,k,n,m)) - aimag(f(jj,1,k,n,m)))
  200 continue
  210 continue
      nrx = nxhyz/nxh
c bit-reverse array elements in x
      do 230 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 230
      do 220 i = kypi, kypt
      t1 = f(1,j1,i,n,m)
      t2 = f(2,j1,i,n,m)
      t3 = f(3,j1,i,n,m)
      f(1,j1,i,n,m) = f(1,j,i,n,m)
      f(2,j1,i,n,m) = f(2,j,i,n,m)
      f(3,j1,i,n,m) = f(3,j,i,n,m)
      f(1,j,i,n,m) = t1
      f(2,j,i,n,m) = t2
      f(3,j,i,n,m) = t3
  220 continue
  230 continue
c finally transform in x
      nrx = nxyz/nxh
      do 270 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 260 i = kypi, kypt
      do 250 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 240 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t1 = s*f(1,j2,i,n,m)
      t2 = s*f(2,j2,i,n,m)
      t3 = s*f(3,j2,i,n,m)
      f(1,j2,i,n,m) = f(1,j1,i,n,m) - t1
      f(2,j2,i,n,m) = f(2,j1,i,n,m) - t2
      f(3,j2,i,n,m) = f(3,j1,i,n,m) - t3
      f(1,j1,i,n,m) = f(1,j1,i,n,m) + t1
      f(2,j1,i,n,m) = f(2,j1,i,n,m) + t2
      f(3,j1,i,n,m) = f(3,j1,i,n,m) + t3
  240 continue
  250 continue
  260 continue
  270 continue
c swap complex components
      do 290 i = kypi, kypt
      do 280 j = 1, nxh
      at1 = real(f(3,j,i,n,m))
      f(3,j,i,n,m) = cmplx(aimag(f(2,j,i,n,m)),aimag(f(3,j,i,n,m)))
      at2 = real(f(2,j,i,n,m))
      f(2,j,i,n,m) = cmplx(at1,aimag(f(1,j,i,n,m)))
      f(1,j,i,n,m) = cmplx(real(f(1,j,i,n,m)),at2)
  280 continue
  290 continue
  300 continue
  310 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT32R3XY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi
     1,kxypp,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,nxhyzd,nxyzhd)
c this subroutine performs the y part of 3 three dimensional real to
c complex fast fourier transforms and their inverses, for a subset of x,
c using complex arithmetic, for data which is distributed in blocks,
c with 2D spatial decomposition
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform is performed
c g(1:3,m,j,i,ld) = sum(g(1:3,k,j,i,id)*exp(-sqrt(-1)*2pi*mm*kk/ny)
c where mm = m + kyp*(ldy - 1) and kk = k + kyp*(idy - 1),
c if isign = 1, a forward fourier transform is performed
c g(1:3,m,j,i,id) = (sum(g(1:3,k,j,i,id)*exp(sqrt(-1)*2pi*mm*kk/ny)
c kstrt = starting data block number
c kxypi = initial x index used
c kxypp = number of x indices used
c nyv = first dimension of g
c kxyp/kzp = number of data values per block in x/z
c kxypd = second dimension of g
c kzpd = third dimension of g
c jblok/lblok = number of data blocks in x/z
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c the real data is stored in a complex array of length nx/2, ny, nz
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c f(j,k,l,i) = mode j-1,kk-1,ll-1, where kk = k + kyp*(ix - 1),
c ll = l + kzp*(iz - 1), and i = iy + nprocy*(iz - 1)
c 1 <= j <= nx/2, 1 <= kk <= ny, and 1 <= ll <= nz, except for
c f(1,k,l,i) = mode nx/2,kk-1,ll-1, where
c ny/2+2 <= kk <= ny, 1 <= ll <= nz, and
c f(1,1,l,i) = mode nx/2,0,ll-1
c f(1,ny/2+1,l,i) = mode nx/2,ny/2,ll-1, where nz/2+2 <= ll <= nz, and
c imaginary part of f(1,1,1,1) = real part of mode nx/2,0,0
c imaginary part of f(1,ny/2+1,1,1) = real part of mode nx/2,ny/2,0
c imaginary part of f(1,1,1,(nz/2)/kzp+1) =
c    real part of mode nx/2,0,nz/2
c imaginary part of f(1,ny/2+1,1,(nz/2)/kzp+1) =
c    real part of mode nx/2,ny/2,nz/2
c using jpl storage convention
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, kxypi, kxypp, nyv, mixup
      integer kxyp, kzp, kxypd, kzpd, jblok, lblok, nxhyzd, nxyzhd
      complex g, sct
      dimension g(3,nyv,kxypd,kzpd,jblok*lblok)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh, ny2, nz, nxyz, nxhyz
      integer j, k, l, i, m, n, ns, ns2, km, kmr, k1, k2, j1, j2, jj
      integer jlblok, js, ks, nry, kzb, kxb, kxypt, mx, mz, moff
      complex s, t1, t2, t3
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kxb = nxh/kxyp
      kzb = nz/kzp
      kxypt = kxypi + kxypp - 1
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      if (kstrt.gt.(kxb*kzb)) return
      jlblok = jblok*lblok
      if (isign.gt.0) go to 140
c inverse fourier transform
      do 80 m = 1, jlblok
      do 70 n = 1, kzp
      nry = nxhyz/ny
c bit-reverse array elements in y
      do 20 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 20
      do 10 i = kxypi, kxypt
      t1 = g(1,k1,i,n,m)
      t2 = g(2,k1,i,n,m)
      t3 = g(3,k1,i,n,m)
      g(1,k1,i,n,m) = g(1,k,i,n,m)
      g(2,k1,i,n,m) = g(2,k,i,n,m)
      g(3,k1,i,n,m) = g(3,k,i,n,m)
      g(1,k,i,n,m) = t1
      g(2,k,i,n,m) = t2
      g(3,k,i,n,m) = t3
   10 continue
   20 continue
c then transform in y
      nry = nxyz/ny
      do 60 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 50 i = kxypi, kxypt
      do 40 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 30 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t1 = s*g(1,j2,i,n,m)
      t2 = s*g(2,j2,i,n,m)
      t3 = s*g(3,j2,i,n,m)
      g(1,j2,i,n,m) = g(1,j1,i,n,m) - t1
      g(2,j2,i,n,m) = g(2,j1,i,n,m) - t2
      g(3,j2,i,n,m) = g(3,j1,i,n,m) - t3
      g(1,j1,i,n,m) = g(1,j1,i,n,m) + t1
      g(2,j1,i,n,m) = g(2,j1,i,n,m) + t2
      g(3,j1,i,n,m) = g(3,j1,i,n,m) + t3
   30 continue
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble modes kx = 0, nx/2
      do 130 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 120 mx = 1, jblok
      if ((mx+js).gt.0) go to 120
      m = mx + moff
      if (kxypi.eq.1) then
         do 110 n = 1, kzp
         do 100 k = 2, nyh
         do 90 jj = 1, 3
         s = g(jj,ny2-k,1,n,m)
         g(jj,ny2-k,1,n,m) = .5*cmplx(aimag(g(jj,k,1,n,m) + s),real(g(jj
     1,k,1,n,m) - s))
         g(jj,k,1,n,m) = .5*cmplx(real(g(jj,k,1,n,m) + s),aimag(g(jj,k,1
     1,n,m) - s))
   90    continue
  100    continue
  110    continue
      endif
  120 continue
  130 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
  140 do 190 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 180 mx = 1, jblok
      if ((mx+js).gt.0) go to 180
      m = mx + moff
      if (kxypi.eq.1) then
         do 170 n = 1, kzp
         do 160 k = 2, nyh
         do 150 jj = 1, 3
         s = cmplx(aimag(g(jj,ny2-k,1,n,m)),real(g(jj,ny2-k,1,n,m)))
         g(jj,ny2-k,1,n,m) = conjg(g(jj,k,1,n,m) - s)
         g(jj,k,1,n,m) = g(jj,k,1,n,m) + s
  150    continue
  160    continue
  170    continue
      endif
  180 continue
  190 continue
      do 270 m = 1, jlblok
      do 260 n = 1, kzp
      nry = nxhyz/ny
c bit-reverse array elements in y
      do 210 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 210
      do 200 i = kxypi, kxypt
      t1 = g(1,k1,i,n,m)
      t2 = g(2,k1,i,n,m)
      t3 = g(3,k1,i,n,m)
      g(1,k1,i,n,m) = g(1,k,i,n,m)
      g(2,k1,i,n,m) = g(2,k,i,n,m)
      g(3,k1,i,n,m) = g(3,k,i,n,m)
      g(1,k,i,n,m) = t1
      g(2,k,i,n,m) = t2
      g(3,k,i,n,m) = t3
  200 continue
  210 continue
c then transform in y
      nry = nxyz/ny
      do 250 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 240 i = kxypi, kxypt
      do 230 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 220 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t1 = s*g(1,j2,i,n,m)
      t2 = s*g(2,j2,i,n,m)
      t3 = s*g(3,j2,i,n,m)
      g(1,j2,i,n,m) = g(1,j1,i,n,m) - t1
      g(2,j2,i,n,m) = g(2,j1,i,n,m) - t2
      g(3,j2,i,n,m) = g(3,j1,i,n,m) - t3
      g(1,j1,i,n,m) = g(1,j1,i,n,m) + t1
      g(2,j1,i,n,m) = g(2,j1,i,n,m) + t2
      g(3,j1,i,n,m) = g(3,j1,i,n,m) + t3
  220 continue
  230 continue
  240 continue
  250 continue
  260 continue
  270 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT32R3XZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi
     1,kxypp,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,nxhyzd,nxyzhd)
c this subroutine performs the z part of 3 three dimensional real to
c complex fast fourier transforms and their inverses, for a subset of x,
c using complex arithmetic, for data which is distributed in blocks,
c with 2D spatial decomposition
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform is performed
c h(1:3,l,n,m,ld) = sum(h(1:3,i,j,k,id)*exp(sqrt(-1)*2pi*ll*ii/nz)
c where ll = l + kzp*(ldz - 1) and ii = i + kzp*(idz - 1),
c if isign = 1, a forward fourier transform is performed
c h(1:3,l,n,m,ld) = (sum(h(1:3,i,j,k,id)*exp(sqrt(-1)*2pi*ll*ii/nz))
c kstrt = starting data block number
c kxypi = initial x index used
c kxypp = number of x indices used
c nzv = first dimension of h
c kxyp/kyzp = number of data values per block in x/y
c kxypd = second dimension of h
c kyzpd = third dimension of h
c jblok/mblok = number of data blocks in x/y
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c the real data is stored in a complex array of length nx/2, ny, nz
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c f(j,k,l,i) = mode j-1,kk-1,ll-1, where kk = k + kyp*(ix - 1),
c ll = l + kzp*(iz - 1), and i = iy + nprocy*(iz - 1)
c 1 <= j <= nx/2, 1 <= kk <= ny, and 1 <= ll <= nz, except for
c f(1,k,l,i) = mode nx/2,kk-1,ll-1, where
c ny/2+2 <= kk <= ny, 1 <= ll <= nz, and
c f(1,1,l,i) = mode nx/2,0,ll-1
c f(1,ny/2+1,l,i) = mode nx/2,ny/2,ll-1, where nz/2+2 <= ll <= nz, and
c imaginary part of f(1,1,1,1) = real part of mode nx/2,0,0
c imaginary part of f(1,ny/2+1,1,1) = real part of mode nx/2,ny/2,0
c imaginary part of f(1,1,1,(nz/2)/kzp+1) =
c    real part of mode nx/2,0,nz/2
c imaginary part of f(1,ny/2+1,1,(nz/2)/kzp+1) =
c    real part of mode nx/2,ny/2,nz/2
c using jpl storage convention
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, kxypi, kxypp, nzv, mixup
      integer kxyp, kyzp, kxypd, kyzpd, jblok, mblok, nxhyzd, nxyzhd
      complex h, sct
      dimension h(3,nzv,kxypd,kyzpd,jblok*mblok)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh, nz, nzh, nz2, nxyz, nxhyz
      integer j, k, l, i, m, n, ns, ns2, km, kmr, k1, k2, j1, j2, jj
      integer jmblok, l1, js, ks, nrz, kxb, kyzb, kxypt, mx, my, moff
      complex s, t1, t2, t3
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      nz = 2**indz
      nzh = nz/2
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      kxypt = kxypi + kxypp - 1
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      if (kstrt.gt.(kxb*kyzb)) return
      jmblok = jblok*mblok
      if (isign.gt.0) go to 160
c inverse fourier transform
      do 80 m = 1, jmblok
      do 70 n = 1, kyzp
      nrz = nxhyz/nz
c bit-reverse array elements in z
      do 20 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 20
      do 10 i = kxypi, kxypt
      t1 = h(1,l1,i,n,m)
      t2 = h(2,l1,i,n,m)
      t3 = h(3,l1,i,n,m)
      h(1,l1,i,n,m) = h(1,l,i,n,m)
      h(2,l1,i,n,m) = h(2,l,i,n,m)
      h(3,l1,i,n,m) = h(3,l,i,n,m)
      h(1,l,i,n,m) = t1
      h(2,l,i,n,m) = t2
      h(3,l,i,n,m) = t3
   10 continue
   20 continue
c finally transform in z
      nrz = nxyz/nz
      do 60 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 50 i = kxypi, kxypt
      do 40 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 30 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t1 = s*h(1,j2,i,n,m)
      t2 = s*h(2,j2,i,n,m)
      t3 = s*h(3,j2,i,n,m)
      h(1,j2,i,n,m) = h(1,j1,i,n,m) - t1
      h(2,j2,i,n,m) = h(2,j1,i,n,m) - t2
      h(3,j2,i,n,m) = h(3,j1,i,n,m) - t3
      h(1,j1,i,n,m) = h(1,j1,i,n,m) + t1
      h(2,j1,i,n,m) = h(2,j1,i,n,m) + t2
      h(3,j1,i,n,m) = h(3,j1,i,n,m) + t3
   30 continue
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble modes kx = 0, nx/2
      do 140 my = 1, mblok
      moff = jblok*(my - 1)
      do 130 mx = 1, jblok
      if ((mx+js).gt.0) go to 130
      m = mx + moff
      if ((my+ks).eq.0) then
         if (kxypi.eq.1) then
            do 100 n = 2, nzh
            do 90 jj = 1, 3
            s = h(jj,nz2-n,1,1,m)
            h(jj,nz2-n,1,1,m) = .5*cmplx(aimag(h(jj,n,1,1,m) + s),real(h
     1(jj,n,1,1,m) - s))
            h(jj,n,1,1,m) = .5*cmplx(real(h(jj,n,1,1,m) + s),aimag(h(jj,
     1n,1,1,m) - s))
   90       continue
  100       continue
         endif
      endif
      kyzb = nyh/kyzp
      if ((my+ks).eq.kyzb) then
         k1 = nyh - kyzb*kyzp + 1
         if (kxypi.eq.1) then
            do 120 n = 2, nzh
            do 110 jj = 1, 3
            s = h(jj,nz2-n,1,k1,m)
            h(jj,nz2-n,1,k1,m) = .5*cmplx(aimag(h(jj,n,1,k1,m) + s),real
     1(h(jj,n,1,k1,m) - s))
            h(jj,n,1,k1,m) = .5*cmplx(real(h(jj,n,1,k1,m) + s),aimag(h(j
     1j,n,1,k1,m) - s))
  110       continue
  120       continue
        endif
      endif
  130 continue
  140 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
  160 do 220 my = 1, mblok
      moff = jblok*(my - 1)
      do 210 mx = 1, jblok
      if ((mx+js).gt.0) go to 210
      m = mx + moff
      if ((my+ks).eq.0) then
         if (kxypi.eq.1) then
            do 180 n = 2, nzh
            do 170 jj = 1, 3
            s = cmplx(aimag(h(jj,nz2-n,1,1,m)),real(h(jj,nz2-n,1,1,m)))
            h(jj,nz2-n,1,1,m) = conjg(h(jj,n,1,1,m) - s)
            h(jj,n,1,1,m) = h(jj,n,1,1,m) + s
  170       continue
  180       continue
         endif
      endif
      kyzb = nyh/kyzp
      if ((my+ks).eq.kyzb) then
         k1 = nyh - kyzb*kyzp + 1
         if (kxypi.eq.1) then
            do 200 n = 2, nzh
            do 190 jj = 1, 3
            s = cmplx(aimag(h(jj,nz2-n,1,k1,m)),real(h(jj,nz2-n,1,k1,m))
     1)
            h(jj,nz2-n,1,k1,m) = conjg(h(jj,n,1,k1,m) - s)
            h(jj,n,1,k1,m) = h(jj,n,1,k1,m) + s
  190       continue
  200       continue
         endif
      endif
  210 continue
  220 continue
      do 300 m = 1, jmblok
      do 290 n = 1, kyzp
      nrz = nxhyz/nz
c bit-reverse array elements in z
      do 240 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 240
      do 230 i = kxypi, kxypt
      t1 = h(1,l1,i,n,m)
      t2 = h(2,l1,i,n,m)
      t3 = h(3,l1,i,n,m)
      h(1,l1,i,n,m) = h(1,l,i,n,m)
      h(2,l1,i,n,m) = h(2,l,i,n,m)
      h(3,l1,i,n,m) = h(3,l,i,n,m)
      h(1,l,i,n,m) = t1
      h(2,l,i,n,m) = t2
      h(3,l,i,n,m) = t3
  230 continue
  240 continue
c first transform in z
      nrz = nxyz/nz
      do 280 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 270 i = kxypi, kxypt
      do 260 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 250 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t1 = s*h(1,j2,i,n,m)
      t2 = s*h(2,j2,i,n,m)
      t3 = s*h(3,j2,i,n,m)
      h(1,j2,i,n,m) = h(1,j1,i,n,m) - t1
      h(2,j2,i,n,m) = h(2,j1,i,n,m) - t2
      h(3,j2,i,n,m) = h(3,j1,i,n,m) - t3
      h(1,j1,i,n,m) = h(1,j1,i,n,m) + t1
      h(2,j1,i,n,m) = h(2,j1,i,n,m) + t2
      h(3,j1,i,n,m) = h(3,j1,i,n,m) + t3
  250 continue
  260 continue
  270 continue
  280 continue
  290 continue
  300 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine WP2FFT32RN(f1,f2,g1,g2,h1,h2,bs,br,ss,isign,ntpose,mixu
     1p,sct,ttp,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxyp
     2d,kypd,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,ndim1
     3,ndim2,nxhyzd,nxyzhd)
c wrapper function for two real to complex ffts
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd, kzyp
      integer jblok, kblok, jkblok, lblok, mblok, mlblok, ndim1, ndim2
      integer nxhyzd, nxyzhd
      integer mixup
      real ttp
      complex f1, f2, g1, g2, h1, h2, bs, br, ss, sct
      dimension f1(ndim1,nxvh,kypd,kzpd,kblok*lblok)
      dimension f2(ndim2,nxvh,kypd,kzpd,kblok*lblok)
      dimension g1(ndim1,nyv,kxypd,kzpd,jblok*lblok)
      dimension g2(ndim2,nyv,kxypd,kzpd,jblok*lblok)
      dimension h1(ndim1,nzv,kxypd,kyzpd,jblok*mblok)
      dimension h2(ndim2,nzv,kxypd,kyzpd,jblok*mblok)
      dimension bs(ndim1+ndim2,kxyp,kzyp,kzp,jkblok*lblok)
      dimension br(ndim1+ndim2,kxyp,kzyp,kzp,jblok*mlblok)
      dimension ss(ndim1+ndim2,nxvh)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer nxh, ny, nz, kypi, kxypi
      real tp, tf
      double precision dtime
      data kypi, kxypi /1,1/
c calculate range of indices
      nxh = 2**(indx - 1)
      ny = 2**indy
      nz = 2**indz
c inverse fourier transform
      if (isign.lt.0) then
c perform x ffts
         call PFFT32RNXX(f1,ss,isign,mixup,sct,indx,indy,indz,kstrt,kypi
     1,kyp,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,ndim1,nxhyzd,nxyzhd)
         call PFFT32RNXX(f2,ss,isign,mixup,sct,indx,indy,indz,kstrt,kypi
     1,kyp,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,ndim2,nxhyzd,nxyzhd)
c transpose f arrays to g
         call PWTIMERA(-1,ttp,dtime)
         call PN2TPOS3A(f1,f2,g1,g2,bs,br,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,
     1kyp,kzp,kxypd,kypd,kzpd,jblok,kblok,lblok,ndim1,ndim2)
         call PWTIMERA(1,ttp,dtime)
c perform y ffts
         call PFFT32RNXY(g1,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,k
     1xyp,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,ndim1,nxhyzd,nxyzhd)
         call PFFT32RNXY(g2,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,k
     1xyp,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,ndim2,nxhyzd,nxyzhd)
c transpose g arrays to h
         call PWTIMERA(-1,tp,dtime)
         call PN2TPOS3B(g1,g2,h1,h2,bs,br,nxh,ny,nz,kstrt,nyv,nzv,kxyp,k
     1yzp,kzp,kxypd,kyzpd,kzpd,jblok,mblok,lblok,ndim1,ndim2)
         call PWTIMERA(1,tp,dtime)
c perform z ffts
         call PFFT32RNXZ(h1,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,k
     1xyp,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,ndim1,nxhyzd,nxyzhd)
         call PFFT32RNXZ(h2,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,k
     1xyp,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,ndim2,nxhyzd,nxyzhd)
c transpose h arrays to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PN2TPOS3B(h1,h2,g1,g2,br,bs,nxh,nz,ny,kstrt,nzv,nyv,kxy
     1p,kzp,kyzp,kxypd,kzpd,kyzpd,jblok,lblok,mblok,ndim1,ndim2)
            call PN2TPOS3A(g1,g2,f1,f2,br,bs,ny,nxh,nz,kstrt,nyv,nxvh,ky
     1p,kxyp,kzp,kypd,kxypd,kzpd,kblok,jblok,lblok,ndim1,ndim2)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f arrays to h
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PN2TPOS3A(f1,f2,g1,g2,bs,br,nxh,ny,nz,kstrt,nxvh,nyv,kx
     1yp,kyp,kzp,kxypd,kypd,kzpd,jblok,kblok,lblok,ndim1,ndim2)
            call PN2TPOS3B(g1,g2,h1,h2,bs,br,nxh,ny,nz,kstrt,nyv,nzv,kxy
     1p,kyzp,kzp,kxypd,kyzpd,kzpd,jblok,mblok,lblok,ndim1,ndim2)
            call PWTIMERA(1,tf,dtime)
         endif
c perform z ffts
         call PFFT32RNXZ(h1,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,k
     1xyp,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,ndim1,nxhyzd,nxyzhd)
         call PFFT32RNXZ(h2,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,k
     1xyp,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,ndim2,nxhyzd,nxyzhd)
c transpose h arrays to g
         call PWTIMERA(-1,tp,dtime)
         call PN2TPOS3B(h1,h2,g1,g2,br,bs,nxh,nz,ny,kstrt,nzv,nyv,kxyp,k
     1zp,kyzp,kxypd,kzpd,kyzpd,jblok,lblok,mblok,ndim1,ndim2)
         call PWTIMERA(1,tp,dtime)
c perform y ffts
         call PFFT32RNXY(g1,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,k
     1xyp,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,ndim1,nxhyzd,nxyzhd)
         call PFFT32RNXY(g2,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,k
     1xyp,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,ndim2,nxhyzd,nxyzhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PN2TPOS3A(g1,g2,f1,f2,br,bs,ny,nxh,nz,kstrt,nyv,nxvh,kyp,k
     1xyp,kzp,kypd,kxypd,kzpd,kblok,jblok,lblok,ndim1,ndim2)
         call PWTIMERA(1,ttp,dtime)
c perform x ffts
         call PFFT32RNXX(f1,ss,isign,mixup,sct,indx,indy,indz,kstrt,kypi
     1,kyp,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,ndim1,nxhyzd,nxyzhd)
         call PFFT32RNXX(f2,ss,isign,mixup,sct,indx,indy,indz,kstrt,kypi
     1,kyp,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,ndim2,nxhyzd,nxyzhd)
      endif
      ttp = ttp + tp
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT32RNXX(f,ss,isign,mixup,sct,indx,indy,indz,kstrt,ky
     1pi,kypp,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,ndim,nxhyzd,nxyzhd)
c this subroutine performs the x part of N three dimensional real to
c complex fast fourier transforms and their inverses, for a subset of y,
c using complex arithmetic, for data which is distributed in blocks,
c with 2D spatial decomposition
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform is performed
c f(1:N,n,k,i,id) = (1/nx*ny*nz)*sum(f(1:N,j,k,i,id)*
c       exp(-sqrt(-1)*2pi*n*j/nx)
c if isign = 1, a forward fourier transform is performed
c f(1:N,n,k,i,id) = (sum(f(1:N,j,k,i,id)*exp(sqrt(-1)*2pi*n*j/nx)
c kstrt = starting data block number
c kypi = initial y index used
c kypp = number of y indices used
c nxvh = first dimension of f
c kyp/kzp = number of data values per block in y/z
c kypd = second dimension of f
c kzpd = third dimension of f
c kblok/lblok = number of data blocks in y/z
c ndim = leading dimension of array f
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c the real data is stored in a complex array of length nx/2, ny, nz
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c f(j,k,l,i) = mode j-1,kk-1,ll-1, where kk = k + kyp*(ix - 1),
c ll = l + kzp*(iz - 1), and i = iy + nprocy*(iz - 1)
c 1 <= j <= nx/2, 1 <= kk <= ny, and 1 <= ll <= nz, except for
c f(1,k,l,i) = mode nx/2,kk-1,ll-1, where
c ny/2+2 <= kk <= ny, 1 <= ll <= nz, and
c f(1,1,l,i) = mode nx/2,0,ll-1
c f(1,ny/2+1,l,i) = mode nx/2,ny/2,ll-1, where nz/2+2 <= ll <= nz, and
c imaginary part of f(1,1,1,1) = real part of mode nx/2,0,0
c imaginary part of f(1,ny/2+1,1,1) = real part of mode nx/2,ny/2,0
c imaginary part of f(1,1,1,(nz/2)/kzp+1) =
c    real part of mode nx/2,0,nz/2
c imaginary part of f(1,ny/2+1,1,(nz/2)/kzp+1) =
c    real part of mode nx/2,ny/2,nz/2
c using jpl storage convention
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, kypi, kypp, nxvh, mixup
      integer kyp, kzp, kypd, kzpd, kblok, lblok, ndim, nxhyzd, nxyzhd
      complex f, ss, sct
      dimension f(ndim,nxvh,kypd,kzpd,kblok*lblok), ss(ndim,nxvh)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nz, nxyz, nxhyz
      integer j, k, l, i, m, n, ns, ns2, km, kmr, k1, k2, j1, j2, jj
      integer klblok, nrx, nry, kyb, kzb, kypt
      real ani
      complex s, t, t1
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kyb = ny/kyp
      kzb = nz/kzp
      kypt = kypi + kypp - 1
      if (kstrt.gt.(kyb*kzb)) return
      klblok = kblok*lblok
      if (isign.gt.0) go to 160
c inverse fourier transform
      ani = 0.5/(float(nx)*float(ny)*float(nz))
c swap complex components
      call PSWAPC32N(f,ss,isign,nxh,kypi,kypt,kzp,nxvh,kypd,kzpd,kblok,l
     1blok,ndim)
      do 150 m = 1, klblok
      do 140 n = 1, kzp
      nrx = nxhyz/nxh
c bit-reverse array elements in x
      do 30 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 30
      do 20 i = kypi, kypt
      do 10 jj = 1, ndim
      t1 = f(jj,j1,i,n,m)
      f(jj,j1,i,n,m) = f(jj,j,i,n,m)
      f(jj,j,i,n,m) = t1
   10 continue
   20 continue
   30 continue
c first transform in x
      nrx = nxyz/nxh
      do 80 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 70 i = kypi, kypt
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 40 jj = 1, ndim
      t1 = s*f(jj,j2,i,n,m)
      f(jj,j2,i,n,m) = f(jj,j1,i,n,m) - t1
      f(jj,j1,i,n,m) = f(jj,j1,i,n,m) + t1
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
c unscramble coefficients and normalize
      kmr = nxyz/nx
      nry = nxhyz/ny
      do 110 k = kypi, kypt
      do 100 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 90 jj = 1, ndim
      t = conjg(f(jj,nxh2-j,k,n,m))
      s = f(jj,j,k,n,m) + t
      t = (f(jj,j,k,n,m) - t)*t1
      f(jj,j,k,n,m) = ani*(s + t)
      f(jj,nxh2-j,k,n,m) = ani*conjg(s - t)
   90 continue
  100 continue
  110 continue
      do 130 k = kypi, kypt
      do 120 jj = 1, ndim
      f(jj,nxhh+1,k,n,m) = 2.*ani*conjg(f(jj,nxhh+1,k,n,m))
      f(jj,1,k,n,m) = 2.*ani*cmplx(real(f(jj,1,k,n,m)) + aimag(f(jj,1,k,
     1n,m)),real(f(jj,1,k,n,m)) - aimag(f(jj,1,k,n,m)))
  120 continue
  130 continue
  140 continue
  150 continue
      return
c forward fourier transform
  160 do 310 m = 1, klblok
      do 300 n = 1, kzp
c scramble coefficients
      kmr = nxyz/nx
      do 190 k = kypi, kypt
      do 180 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 170 jj = 1, ndim
      t = conjg(f(jj,nxh2-j,k,n,m))
      s = f(jj,j,k,n,m) + t
      t = (f(jj,j,k,n,m) - t)*t1
      f(jj,j,k,n,m) = s + t
      f(jj,nxh2-j,k,n,m) = conjg(s - t)
  170 continue
  180 continue
  190 continue
      do 210 k = kypi, kypt
      do 200 jj = 1, ndim
      f(jj,nxhh+1,k,n,m) = 2.*conjg(f(jj,nxhh+1,k,n,m))
      f(jj,1,k,n,m) = cmplx(real(f(jj,1,k,n,m)) + aimag(f(jj,1,k,n,m)),r
     1eal(f(jj,1,k,n,m)) - aimag(f(jj,1,k,n,m)))
  200 continue
  210 continue
      nrx = nxhyz/nxh
c bit-reverse array elements in x
      do 240 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 240
      do 230 i = kypi, kypt
      do 220 jj = 1, ndim
      t1 = f(jj,j1,i,n,m)
      f(jj,j1,i,n,m) = f(jj,j,i,n,m)
      f(jj,j,i,n,m) = t1
  220 continue
  230 continue
  240 continue
c finally transform in x
      nrx = nxyz/nxh
      do 290 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 280 i = kypi, kypt
      do 270 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 260 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 250 jj = 1, ndim
      t1 = s*f(jj,j2,i,n,m)
      f(jj,j2,i,n,m) = f(jj,j1,i,n,m) - t1
      f(jj,j1,i,n,m) = f(jj,j1,i,n,m) + t1
  250 continue
  260 continue
  270 continue
  280 continue
  290 continue
  300 continue
  310 continue
c swap complex components
      call PSWAPC32N(f,ss,isign,nxh,kypi,kypt,kzp,nxvh,kypd,kzpd,kblok,l
     1blok,ndim)
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT32RNXY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi
     1,kxypp,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,ndim,nxhyzd,nxyzhd)
c this subroutine performs the y part of N three dimensional real to
c complex fast fourier transforms and their inverses, for a subset of x,
c using complex arithmetic, for data which is distributed in blocks,
c with 2D spatial decomposition
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform is performed
c g(1:N,m,j,i,ld) = sum(g(1:N,k,j,i,id)*exp(-sqrt(-1)*2pi*mm*kk/ny)
c where mm = m + kyp*(ldy - 1) and kk = k + kyp*(idy - 1),
c if isign = 1, a forward fourier transform is performed
c g(1:N,m,j,i,id) = (sum(g(1:N,k,j,i,id)*exp(sqrt(-1)*2pi*mm*kk/ny)
c kstrt = starting data block number
c kxypi = initial x index used
c kxypp = number of x indices used
c nyv = first dimension of g
c kxyp/kzp = number of data values per block in x/z
c kxypd = second dimension of g
c kzpd = third dimension of g
c jblok/lblok = number of data blocks in x/z
c ndim = leading dimension of array g
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c the real data is stored in a complex array of length nx/2, ny, nz
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c f(j,k,l,i) = mode j-1,kk-1,ll-1, where kk = k + kyp*(ix - 1),
c ll = l + kzp*(iz - 1), and i = iy + nprocy*(iz - 1)
c 1 <= j <= nx/2, 1 <= kk <= ny, and 1 <= ll <= nz, except for
c f(1,k,l,i) = mode nx/2,kk-1,ll-1, where
c ny/2+2 <= kk <= ny, 1 <= ll <= nz, and
c f(1,1,l,i) = mode nx/2,0,ll-1
c f(1,ny/2+1,l,i) = mode nx/2,ny/2,ll-1, where nz/2+2 <= ll <= nz, and
c imaginary part of f(1,1,1,1) = real part of mode nx/2,0,0
c imaginary part of f(1,ny/2+1,1,1) = real part of mode nx/2,ny/2,0
c imaginary part of f(1,1,1,(nz/2)/kzp+1) =
c    real part of mode nx/2,0,nz/2
c imaginary part of f(1,ny/2+1,1,(nz/2)/kzp+1) =
c    real part of mode nx/2,ny/2,nz/2
c using jpl storage convention
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, kxypi, kxypp, nyv, mixup
      integer kxyp, kzp, kxypd, kzpd, jblok, lblok, ndim, nxhyzd, nxyzhd
      complex g, sct
      dimension g(ndim,nyv,kxypd,kzpd,jblok*lblok)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh, ny2, nz, nxyz, nxhyz
      integer j, k, l, i, m, n, ns, ns2, km, kmr, k1, k2, j1, j2, jj
      integer jlblok, js, ks, nry, kzb, kxb, kxypt, mx, mz, moff
      complex s, t1
      if (isign.eq.0) return
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nz = 2**indz
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kxb = nxh/kxyp
      kzb = nz/kzp
      kxypt = kxypi + kxypp - 1
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      if (kstrt.gt.(kxb*kzb)) return
      jlblok = jblok*lblok
      if (isign.gt.0) go to 160
c inverse fourier transform
      do 100 m = 1, jlblok
      do 90 n = 1, kzp
      nry = nxhyz/ny
c bit-reverse array elements in y
      do 30 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 30
      do 20 i = kxypi, kxypt
      do 10 jj = 1, ndim
      t1 = g(jj,k1,i,n,m)
      g(jj,k1,i,n,m) = g(jj,k,i,n,m)
      g(jj,k,i,n,m) = t1
   10 continue
   20 continue
   30 continue
c then transform in y
      nry = nxyz/ny
      do 80 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 70 i = kxypi, kxypt
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 40 jj = 1, ndim
      t1 = s*g(jj,j2,i,n,m)
      g(jj,j2,i,n,m) = g(jj,j1,i,n,m) - t1
      g(jj,j1,i,n,m) = g(jj,j1,i,n,m) + t1
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
   90 continue
  100 continue
c unscramble modes kx = 0, nx/2
      do 150 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 140 mx = 1, jblok
      if ((mx+js).gt.0) go to 140
      m = mx + moff
      if (kxypi.eq.1) then
         do 130 n = 1, kzp
         do 120 k = 2, nyh
         do 110 jj = 1, ndim
         s = g(jj,ny2-k,1,n,m)
         g(jj,ny2-k,1,n,m) = .5*cmplx(aimag(g(jj,k,1,n,m) + s),real(g(jj
     1,k,1,n,m) - s))
         g(jj,k,1,n,m) = .5*cmplx(real(g(jj,k,1,n,m) + s),aimag(g(jj,k,1
     1,n,m) - s))
  110    continue
  120    continue
  130    continue
      endif
  140 continue
  150 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
  160 do 210 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 200 mx = 1, jblok
      if ((mx+js).gt.0) go to 200
      m = mx + moff
      if (kxypi.eq.1) then
         do 190 n = 1, kzp
         do 180 k = 2, nyh
         do 170 jj = 1, ndim
         s = cmplx(aimag(g(jj,ny2-k,1,n,m)),real(g(jj,ny2-k,1,n,m)))
         g(jj,ny2-k,1,n,m) = conjg(g(jj,k,1,n,m) - s)
         g(jj,k,1,n,m) = g(jj,k,1,n,m) + s
  170    continue
  180    continue
  190    continue
      endif
  200 continue
  210 continue
      do 310 m = 1, jlblok
      do 300 n = 1, kzp
      nry = nxhyz/ny
c bit-reverse array elements in y
      do 240 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 240
      do 230 i = kxypi, kxypt
      do 220 jj = 1, ndim
      t1 = g(jj,k1,i,n,m)
      g(jj,k1,i,n,m) = g(jj,k,i,n,m)
      g(jj,k,i,n,m) = t1
  220 continue
  230 continue
  240 continue
c then transform in y
      nry = nxyz/ny
      do 290 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 280 i = kxypi, kxypt
      do 270 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 260 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 250 jj = 1, ndim
      t1 = s*g(jj,j2,i,n,m)
      g(jj,j2,i,n,m) = g(jj,j1,i,n,m) - t1
      g(jj,j1,i,n,m) = g(jj,j1,i,n,m) + t1
  250 continue
  260 continue
  270 continue
  280 continue
  290 continue
  300 continue
  310 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT32RNXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi
     1,kxypp,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,ndim,nxhyzd,nxyzhd)
c this subroutine performs the z part of N three dimensional real to
c complex fast fourier transforms and their inverses, for a subset of x,
c using complex arithmetic, for data which is distributed in blocks,
c with 2D spatial decomposition
c for isign = (-1,1), input: all, output: f
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c if isign = -1, an inverse fourier transform is performed
c h(1:N,l,n,m,ld) = sum(h(1:N,i,j,k,id)*exp(sqrt(-1)*2pi*ll*ii/nz)
c where ll = l + kzp*(ldz - 1) and ii = i + kzp*(idz - 1),
c if isign = 1, a forward fourier transform is performed
c h(1:N,l,n,m,ld) = (sum(h(1:N,i,j,k,id)*exp(sqrt(-1)*2pi*ll*ii/nz))
c kstrt = starting data block number
c kxypi = initial x index used
c kxypp = number of x indices used
c nzv = first dimension of h
c kxyp/kyzp = number of data values per block in x/y
c kxypd = second dimension of h
c kyzpd = third dimension of h
c jblok/mblok = number of data blocks in x/y
c ndim = leading dimension of array h
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c the real data is stored in a complex array of length nx/2, ny, nz
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c f(j,k,l,i) = mode j-1,kk-1,ll-1, where kk = k + kyp*(ix - 1),
c ll = l + kzp*(iz - 1), and i = iy + nprocy*(iz - 1)
c 1 <= j <= nx/2, 1 <= kk <= ny, and 1 <= ll <= nz, except for
c f(1,k,l,i) = mode nx/2,kk-1,ll-1, where
c ny/2+2 <= kk <= ny, 1 <= ll <= nz, and
c f(1,1,l,i) = mode nx/2,0,ll-1
c f(1,ny/2+1,l,i) = mode nx/2,ny/2,ll-1, where nz/2+2 <= ll <= nz, and
c imaginary part of f(1,1,1,1) = real part of mode nx/2,0,0
c imaginary part of f(1,ny/2+1,1,1) = real part of mode nx/2,ny/2,0
c imaginary part of f(1,1,1,(nz/2)/kzp+1) =
c    real part of mode nx/2,0,nz/2
c imaginary part of f(1,ny/2+1,1,(nz/2)/kzp+1) =
c    real part of mode nx/2,ny/2,nz/2
c using jpl storage convention
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, indx, indy, indz, kstrt, kxypi, kxypp, nzv, mixup
      integer kxyp, kyzp, kxypd, kyzpd, jblok, mblok, ndim
      integer nxhyzd, nxyzhd
      complex h, sct
      dimension h(ndim,nzv,kxypd,kyzpd,jblok*mblok)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, ny, nyh, nz, nzh, nz2, nxyz, nxhyz
      integer j, k, l, i, m, n, ns, ns2, km, kmr, k1, k2, j1, j2, jj
      integer jmblok, l1, js, ks, nrz, kxb, kyzb, kxypt, mx, my, moff
      complex s, t1
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      ny = 2**indy
      nyh = ny/2
      nz = 2**indz
      nzh = nz/2
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      kxypt = kxypi + kxypp - 1
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      if (kstrt.gt.(kxb*kyzb)) return
      jmblok = jblok*mblok
      if (isign.gt.0) go to 180
c inverse fourier transform
      do 100 m = 1, jmblok
      do 90 n = 1, kyzp
      nrz = nxhyz/nz
c bit-reverse array elements in z
      do 30 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 30
      do 20 i = kxypi, kxypt
      do 10 jj = 1, ndim
      t1 = h(jj,l1,i,n,m)
      h(jj,l1,i,n,m) = h(jj,l,i,n,m)
      h(jj,l,i,n,m) = t1
   10 continue
   20 continue
   30 continue
c finally transform in z
      nrz = nxyz/nz
      do 80 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 70 i = kxypi, kxypt
      do 60 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 50 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 40 jj = 1, ndim
      t1 = s*h(jj,j2,i,n,m)
      h(jj,j2,i,n,m) = h(jj,j1,i,n,m) - t1
      h(jj,j1,i,n,m) = h(jj,j1,i,n,m) + t1
   40 continue
   50 continue
   60 continue
   70 continue
   80 continue
   90 continue
  100 continue
c unscramble modes kx = 0, nx/2
      do 160 my = 1, mblok
      moff = jblok*(my - 1)
      do 150 mx = 1, jblok
      if ((mx+js).gt.0) go to 150
      m = mx + moff
      if ((my+ks).eq.0) then
         if (kxypi.eq.1) then
            do 120 n = 2, nzh
            do 110 jj = 1, ndim
            s = h(jj,nz2-n,1,1,m)
            h(jj,nz2-n,1,1,m) = .5*cmplx(aimag(h(jj,n,1,1,m) + s),real(h
     1(jj,n,1,1,m) - s))
            h(jj,n,1,1,m) = .5*cmplx(real(h(jj,n,1,1,m) + s),aimag(h(jj,
     1n,1,1,m) - s))
  110       continue
  120       continue
         endif
      endif
      kyzb = nyh/kyzp
      if ((my+ks).eq.kyzb) then
         k1 = nyh - kyzb*kyzp + 1
         if (kxypi.eq.1) then
            do 140 n = 2, nzh
            do 130 jj = 1, ndim
            s = h(jj,nz2-n,1,k1,m)
            h(jj,nz2-n,1,k1,m) = .5*cmplx(aimag(h(jj,n,1,k1,m) + s),real
     1(h(jj,n,1,k1,m) - s))
            h(jj,n,1,k1,m) = .5*cmplx(real(h(jj,n,1,k1,m) + s),aimag(h(j
     1j,n,1,k1,m) - s))
  130       continue
  140       continue
        endif
      endif
  150 continue
  160 continue
      return
c forward fourier transform
c scramble modes kx = 0, nx/2
  180 do 240 my = 1, mblok
      moff = jblok*(my - 1)
      do 230 mx = 1, jblok
      if ((mx+js).gt.0) go to 230
      m = mx + moff
      if ((my+ks).eq.0) then
         if (kxypi.eq.1) then
            do 200 n = 2, nzh
            do 190 jj = 1, ndim
            s = cmplx(aimag(h(jj,nz2-n,1,1,m)),real(h(jj,nz2-n,1,1,m)))
            h(jj,nz2-n,1,1,m) = conjg(h(jj,n,1,1,m) - s)
            h(jj,n,1,1,m) = h(jj,n,1,1,m) + s
  190       continue
  200       continue
         endif
      endif
      kyzb = nyh/kyzp
      if ((my+ks).eq.kyzb) then
         k1 = nyh - kyzb*kyzp + 1
         if (kxypi.eq.1) then
            do 220 n = 2, nzh
            do 210 jj = 1, ndim
            s = cmplx(aimag(h(jj,nz2-n,1,k1,m)),real(h(jj,nz2-n,1,k1,m))
     1)
            h(jj,nz2-n,1,k1,m) = conjg(h(jj,n,1,k1,m) - s)
            h(jj,n,1,k1,m) = h(jj,n,1,k1,m) + s
  210       continue
  220       continue
         endif
      endif
  230 continue
  240 continue
      do 340 m = 1, jmblok
      do 330 n = 1, kyzp
      nrz = nxhyz/nz
c bit-reverse array elements in z
      do 270 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 270
      do 260 i = kxypi, kxypt
      do 250 jj = 1, ndim
      t1 = h(jj,l1,i,n,m)
      h(jj,l1,i,n,m) = h(jj,l,i,n,m)
      h(jj,l,i,n,m) = t1
  250 continue
  260 continue
  270 continue
c first transform in z
      nrz = nxyz/nz
      do 320 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 310 i = kxypi, kxypt
      do 300 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 290 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 280 jj = 1, ndim
      t1 = s*h(jj,j2,i,n,m)
      h(jj,j2,i,n,m) = h(jj,j1,i,n,m) - t1
      h(jj,j1,i,n,m) = h(jj,j1,i,n,m) + t1
  280 continue
  290 continue
  300 continue
  310 continue
  320 continue
  330 continue
  340 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSWAPC32N(f,s,isign,nxh,kypi,kypt,kzp,nxvh,kypd,kzpd,kb
     1lok,lblok,ndim)
c this subroutine swaps components for multiple ffts
c f = input  array
c s = scratch array
c isign = (-1,1) = swap (real-to-complex,complex-to-real)
c nxh = complex dimension in x direction
c kypi/kypt = initial/final y index used
c kzp = complex dimension in z direction
c nxvh = half of the second dimension of f
c kypd = third dimension of f
c kzpd = fourth dimension of f
c lblok = number of data blocks in z
c ndim = leading dimension of array f
      implicit none
      integer isign, nxh, kypi, kypt, kzp, nxvh, kypd, kzpd
      integer kblok, lblok, ndim
      real f, s
      dimension f(ndim,2*nxvh,kypd,kzpd,kblok*lblok), s(2*ndim*nxvh)
c local data
      integer i, j, k, l, m, ioff, klblok
      klblok = kblok*lblok
c swap complex components
c real to complex
      if (isign.lt.0) then
         do 80 m = 1, klblok
         do 70 l = 1, kzp
         do 60 k = kypi, kypt
         do 20 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 10 i = 1, ndim
         s(2*i+ioff-1) = f(i,2*j-1,k,l,m)
         s(2*i+ioff) = f(i,2*j,k,l,m)
   10    continue
   20    continue
         do 50 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 30 i = 1, ndim
         f(i,2*j-1,k,l,m) = s(i+ioff)
   30    continue
         ioff = ioff + ndim
         do 40 i = 1, ndim
         f(i,2*j,k,l,m) = s(i+ioff)
   40    continue
   50    continue
   60    continue
   70    continue
   80    continue
      else if (isign.gt.0) then
c swap complex components
         do 160 m = 1, klblok
         do 150 l = 1, kzp
         do 140 k = kypi, kypt
         do 110 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 90 i = 1, ndim
         s(i+ioff) = f(i,2*j-1,k,l,m)
   90    continue
         ioff = ioff + ndim
         do 100 i = 1, ndim
         s(i+ioff) = f(i,2*j,k,l,m)
  100    continue
  110    continue
         do 130 j = 1, nxh
         ioff = 2*ndim*(j - 1)
         do 120 i = 1, ndim
         f(i,2*j-1,k,l,m) = s(2*i+ioff-1)
         f(i,2*j,k,l,m) = s(2*i+ioff)
  120    continue
  130    continue
  140    continue
  150    continue
  160    continue
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
      function ranorm()
c this program calculates a random number y from a gaussian distribution
c with zero mean and unit variance, according to the method of
c mueller and box:
c    y(k) = (-2*ln(x(k)))**1/2*sin(2*pi*x(k+1))
c    y(k+1) = (-2*ln(x(k)))**1/2*cos(2*pi*x(k+1)),
c where x is a random number uniformly distributed on (0,1).
c written for the ibm by viktor k. decyk, ucla
      integer r1,r2,r4,r5
      double precision ranorm,h1l,h1u,h2l,r0,r3,asc,bsc,temp
      save iflg,r1,r2,r4,r5,h1l,h1u,h2l,r0
      data r1,r2,r4,r5 /885098780,1824280461,1396483093,55318673/
      data h1l,h1u,h2l /65531.0d0,32767.0d0,65525.0d0/
      data iflg,r0 /0,0.0d0/
      if (iflg.eq.0) go to 10
      ranorm = r0
      r0 = 0.0d0
      iflg = 0
      return
   10 isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1 - (r1/isc)*isc
      r3 = h1l*dble(r1) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2/isc
      isc = r2 - i1*isc
      r0 = h1l*dble(r2) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1 = r3 - dble(isc)*bsc
      temp = dsqrt(-2.0d0*dlog((dble(r1) + dble(r2)*asc)*asc))
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r4 - (r4/isc)*isc
      r3 = h2l*dble(r4) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r5/isc
      isc = r5 - i1*isc
      r0 = h2l*dble(r5) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r5 = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r4 = r3 - dble(isc)*bsc
      r0 = 6.28318530717959d0*((dble(r4) + dble(r5)*asc)*asc)
      ranorm = temp*dsin(r0)
      r0 = temp*dcos(r0)
      iflg = 1
      return
      end
c-----------------------------------------------------------------------
      subroutine pranorm(vran,kstrt,nvp,nvrp,ndp,nvrd,nblok)
c this program calculates nvrp random numbers for nvp processors from a 
c gaussian distribution with zero mean and unit variance, according to
c the method of mueller and box:
c    y(k) = (-2*ln(x(k)))**1/2*sin(2*pi*x(k+1))
c    y(k+1) = (-2*ln(x(k)))**1/2*cos(2*pi*x(k+1)),
c where x is a random number uniformly distributed on (0,1).
c written for the ibm by viktor k. decyk, ucla
c parallel version
c each random number is generated from a separate seed, which is just
c the normal seed calculated every 100,000,000 million numbers apart.
c thus if more than 100,000,000 arrays of numbers are requested, the
c different numbers will no longer be unique.
c each processor uses no more than (ndv-1)/nvp+1 seeds.
c if nvp > ndv, adjacent mdp processors will share the same seed, and
c return the same random numbers, where mdp = nvp/min0(nvp,ndv)
c vran = output array of random numbers
c kstrt = starting data block number
c nvp = number of real or virtual processors
c nvrp = number of random numbers requested per processor, <= ndv/nvp
c ndp = number of random numbers returned = min((ndv-1)/nvp+1,nvrp)
c nvrd = first dimension of vran array
c nblok = number of data blocks
c ndv = total maximum number of random seeds, currently 256
      parameter(ndvb=64,ndv=4*ndvb)
      integer r1,r2,r4,r5
      integer r1a,r1b,r1c,r1d,r2a,r2b,r2c,r2d
      integer r4a,r4b,r4c,r4d,r5a,r5b,r5c,r5d
      double precision vran,vran0,h1l,h1u,h2l,r0,r3,asc,bsc,temp
      dimension vran(nvrd,nblok)
      dimension r1(ndv), r2(ndv), r4(ndv), r5(ndv), vran0(ndv)
      dimension r1a(ndvb), r1b(ndvb), r1c(ndvb), r1d(ndvb)
      dimension r2a(ndvb), r2b(ndvb), r2c(ndvb), r2d(ndvb)
      dimension r4a(ndvb), r4b(ndvb), r4c(ndvb), r4d(ndvb)
      dimension r5a(ndvb), r5b(ndvb), r5c(ndvb), r5d(ndvb)
      equivalence (r1a(1),r1(1)), (r1b(1),r1(ndvb+1))
      equivalence (r1c(1),r1(2*ndvb+1)), (r1d(1),r1(3*ndvb+1))
      equivalence (r2a(1),r2(1)), (r2b(1),r2(ndvb+1))
      equivalence (r2c(1),r2(2*ndvb+1)), (r2d(1),r2(3*ndvb+1))
      equivalence (r4a(1),r4(1)), (r4b(1),r4(ndvb+1))
      equivalence (r4c(1),r4(2*ndvb+1)), (r4d(1),r4(3*ndvb+1))
      equivalence (r5a(1),r5(1)), (r5b(1),r5(ndvb+1))
      equivalence (r5c(1),r5(2*ndvb+1)), (r5d(1),r5(3*ndvb+1))
      save iflg,h1l,h1u,h2l
      save r1,r2,r4,r5,vran0
      data iflg /0/
      data h1l,h1u,h2l /65531.0d0,32767.0d0,65525.0d0/
      data r1a /359740401,579253173,631138885,301320988,616045397,197865
     1899,320616057,1483498427,94225675,293599141,119335525,1401265604,4
     238232031,258875505,1516636373,2083742725,1049217694,678455002,2039
     3395477,737755465,2088659986,2037589420,402632816,1340891228,103903
     41049,1282693042,1660116788,464751626,1886797,1266915258,880133036,
     51015492382,111053996,868360335,1180991652,204559462,240784624,1406
     692733,1701763304,1141745325,8877543,1168775120,605400900,805925055
     7,632316507,334371796,229344900,1638112616,1487940890,2056688597,14
     898710389,1659680628,968083694,320658079,486957793,1770449327,16311
     976869,1007414014,8949141,161873637,319971082,1275092926,1886651485
     a,1762615778/
      data r1b /1902075624,1151890347,1362323973,1914314098,798243597,96
     1731537,839057753,112122977,1165362756,2393148,2146479398,860350843
     2,572224344,1756173944,953937974,1643196300,1592579351,1800565521,5
     311044950,222837089,1923358284,1578824916,1765893970,1937519044,168
     42760131,1948161675,2073215223,1229841282,2118406728,1887892227,344
     5980751,1829747047,20529831,1502561544,53526120,430292544,146251102
     60,730097511,1497253836,428186839,1181683234,197448938,508787720,20
     742011641,614837975,2076223582,48029161,1020375667,163928726,200907
     8920,320376889,234970895,51708971,251857994,1779848646,1946275994,1
     9177370210,1283584601,1689751898,69259745,1677212679,875358602,1420
     a772803,356371886/
      data r1c /2075250663,281816231,864163276,1925828303,1843883212,993
     1629822,1434862656,264732686,259532163,1298178521,484684781,9696164
     241,1553247448,1430153606,747890590,1274723866,733959575,132229199,
     3459407902,1645068415,1623986828,1953151234,652920379,189619762,534
     4804100,923506282,1721538240,2106876417,7582793,738712402,212595263
     52,1213150903,1078894248,582271624,1397557297,125560607,1988864109,
     6276375255,994086582,380247559,1184293989,1137541898,860286163,1140
     7532794,321783130,389804527,283551060,52554372,828399130,610283089,
     8897047301,543356336,462348655,1994887932,1092898994,1939108492,611
     9066406,1066474171,1172793639,320744419,704402188,120414300,2056136
     a496,117017729/
      data r1d /1217742571,1637785258,1840721305,2075238194,232681617,55
     11071857,953335598,1969426114,735109832,132190333,857915067,4588492
     25,733433950,2110708122,616214284,1878839215,556665886,1087030932,2
     3038412013,2046203556,1562577107,420426295,1948814634,17606056,5508
     485707,199759695,1631428751,1010026375,1356312656,246616232,1243148
     5264,1811048205,671834288,970938638,669392386,224431878,1789222198,
     62079181645,379743655,220719678,422295726,1282466353,136102308,2055
     7455617,594945507,1594669190,1995807941,977339803,464385955,1450641
     8689,429859287,968768344,1766727865,2004299444,358421412,569085829,
     92082681841,1541809188,826594553,172674284,2135647255,632194101,155
     a019427,737107546/
      data r2a /28358029,1188633485,1412792717,1103488909,663375245,4951
     104909,1001331085,437223309,1352918413,2003585933,644395405,1972967
     2309,2096987533,1419109261,341985677,1415753613,748098957,889158541
     3,94101901,913065869,1601219981,413733773,2048227725,464904077,3613
     483309,2140318605,1909395853,71268237,1323556237,1773945741,1825089
     5933,1879641997,192771469,1462098829,1795309965,1595058061,12639963
     601,1204777869,1820055949,1365000077,242263437,1001982861,189932788
     75,1189468045,1422540173,853713805,2033125773,1068461965,509859213,
     8759970701,73965965,1001981837,1799187853,720753549,316815757,99002
     97661,995558797,736062349,614191501,1032599437,246455693,805897101,
     a966093197,1129697165/
      data r2b /1699362189,930257805,1372520845,1281320845,1059310989,11
     109144461,1833474445,1487470477,473785741,1342557069,201470349,1748
     2146061,2090270093,1630495629,771475853,2063347597,1613796749,19729
     360141,1396007309,285591437,1191849357,222466957,2075064717,7098448
     477,824427917,673983373,661164429,1188624269,511532429,1180025741,1
     5449273741,1721929613,253162893,1740594061,144425357,162277261,4931
     69309,208204685,1041586573,804634509,2047485325,877824909,199327374
     71,1501517709,1952693645,1601971085,852003213,105443213,1912427917,
     8233159565,1912742285,911378317,1926688141,1066357645,880523661,177
     91839373,1995474317,1954081677,2050314637,539342733,2118786445,7488
     a48013,1127147917,1508855693/
      data r2c / 149140877,1745623949,258507149,385410957,381504909,6494
     142189,1591875981,1463975821,668394893,1755270029,832287117,4495829
     289,1009810829,768140173,127224205,1637199757,1405752717,1983019917
     3,1624170893,731858829,1856220557,1104941965,1028159885,2028527501,
     4213730701,281389965,486674829,1232238477,773250445,1659847565,2147
     5199373,490475405,1387296141,945347469,1714766221,1950721933,205586
     67789,285373325,1336859021,1318010765,631481741,1827408781,10134777
     773,739825549,1409105293,1276486541,744622477,216166285,93771149,78
     80090253,530293133,1894516621,980446605,338219917,370489741,1479909
     9261,1921648013,2098359181,265212301,1119827853,769891725,176554074
     a9,214460813,814272397/
      data r2d /1820145037,1487248269,218235277,563242893,777440653,1263
     1481741,276535693,366739341,1936745869,1094241165,389362061,2247617
     241,1003093389,979526541,556714381,137310093,123966861,919337869,77
     38592653,104384397,1446849933,913675149,1054996877,125984653,676775
     4309,962538381,1385927053,202110861,2108710285,1065927565,177138318
     51,332763021,1447687565,1223842701,63881613,517941133,841190797,143
     66283789,558389645,757645197,289219981,1703250829,1107423629,105187
     75213,1939258765,2024743821,1710983565,1400631181,1496339853,253279
     8117,221585805,1803913101,1107946893,683824013,934197645,114237325,
     9774079885,1168894861,1701335437,626571149,494738829,1708491661,375
     a515533,1193430925/
      data r4a /1494108262,1989145426,1552962897,694486778,1754614022,25
     11970638,1939753416,1149145743,2023686104,1040210686,1683330424,114
     23900382,581089000,293341551,1273074283,260197364,1916786754,878264
     3750,253082287,2146393565,272494066,425825733,91072590,206790821,18
     437933060,642168769,604837205,1945376344,1849933443,1587706286,6836
     562288,1769955586,346238332,251055511,1157232075,1040925359,7394408
     692,2007877002,1222081873,202759530,1948470029,2014942097,131050200
     70,1652912961,1267494332,963938746,1664474419,1873307198,1529793301
     8,1507924495,909609381,741685359,1091759396,392812189,1071533011,10
     95548766,996769269,2063692145,1854856650,208609944,1401636755,58899
     a605,1364597140,412212905/
      data r4b /150549516,1353755158,1122564414,2064237085,2086423914,13
     112181247,346330836,630583039,444972475,1587991646,965618785,144897
     28269,1198444297,108893502,1901270003,1123044673,154725217,38784221
     38,1934149525,1128620729,858627279,1507043503,426565970,2092730831,
     4280923487,986857306,191175479,220038288,277186333,1949669045,17699
     551057,899167496,1866757140,956408013,2008119210,96771395,165919905
     68,1612869550,1218031183,111337293,1165651938,258749796,1753099596,
     71704672562,1077801679,1721451514,139284141,528357630,1989300518,71
     80151550,681413757,730068156,822220019,1864573946,1152221386,642609
     9786,138405954,1149057404,2103698846,2095273026,270963999,160636110
     a2,1146066183,397835297/
      data r4c /45680022,1194183101,111642895,94366755,1797407667,209651
     15221,1481300931,652265045,1354070401,707645090,515950383,147864288
     20,1642983694,1238478320,1010345311,494264691,527373412,1405573259,
     3982118495,1297417849,117730705,1742495871,1983136315,1000431070,10
     450781343,2026272727,1230984446,1688512941,1063422363,1998811297,92
     50322742,403753971,9109649,2088320231,354370925,603965117,579278221
     6,1983246006,961296177,2146679829,434827931,1353693210,573182476,20
     710574344,1517478342,572844831,589901835,1459338914,33770780,143838
     82545,258083642,320647151,61390503,1451456060,1687502245,2002543610
     9,1385685109,980554348,1982562648,1329962193,221084816,1835442284,4
     a06171742,298531967/
      data r4d /669940997,1546129994,1248642245,206062863,411560928,5267
     144080,1664211095,226482318,13562669,649464267,982394338,278739208,
     21505811704,1670976391,1429407408,1748224138,511906636,1953892688,2
     359650830,1765738781,2002985015,1335655730,1214548744,371367120,169
     41790787,1849958713,211583658,1235895726,1786480126,2005714801,1098
     5102272,1644815811,974623474,1803444709,1340350230,1809678106,15870
     676780,1309213347,1334826999,1294545857,1729467192,1376050101,83473
     78865,1884898754,298580642,70401833,1818902862,1866602282,199874943
     80,1983487753,623232475,1042295301,883403328,1772850530,1547059812,
     91452810335,403842699,2063645760,394686709,1656142916,1245756901,12
     a85160367,229190554,1743839168/
      data r5a /1454526097,1497138321,585808529,1270673553,1806902929,44
     19666193,1896583825,107858065,1928593041,1318990993,829188753,86183
     29505,1819596433,1957629073,1678590609,1385134225,1479913105,218096
     3785,149822097,1677742225,909543057,395361425,537850513,1739663505,
     4108486289,341939345,695192209,1570898065,1224226449,57830545,62184
     57185,1171445905,2109279889,1690518673,317815441,541307025,61616296
     61,945036433,1930580625,1827965073,1039842961,2116351121,1165175441
     7,736452753,1232836241,909495441,169083537,1561737361,1195142801,16
     819436689,1089788561,8851601,926762641,2098691217,1779806865,372762
     9769,427695761,199775377,91654801,505987217,1845425809,217656465,32
     a0299665,408524945/
      data r5b /884985489,4850833,318257809,80375953,1841342097,17088421
     129,85529233,1669023889,419528337,1034663057,1769597585,879501457,9
     214511505,129797265,1075495569,2006775953,1178807953,1141728401,150
     3706833,755880081,1212417681,1922972817,1142715025,1421781137,10153
     440689,326046865,1904036497,1856995473,587576977,645917841,28718760
     51,2061523089,2076610193,735102097,587135633,2035363985,1187473041,
     6593599633,656396945,1778518161,67649169,221410449,494971537,129098
     75617,864622225,1766018193,102859409,572766353,1430908561,932455569
     8,1627544209,1771344017,1766508177,2015689873,774058641,591751313,1
     9871421073,720753809,1837370001,1328955537,1745647249,1342614673,52
     a2510993,1835473041/
      data r5c /1389186705,1733788817,1124448913,2111303825,802039441,18
     194276241,1495700113,8964241,2131689105,1824076945,1636264593,19709
     205233,1083168401,1523190929,1546142353,1554675857,1951444625,99161
     38193,1225333393,907759761,441550481,229358737,673837713,30156945,8
     448453265,1383896209,2039138961,1069351057,1024669329,160263313,102
     56269841,1877858449,970198673,853427345,1930197649,308195473,685041
     6297,1315904657,455955089,655329425,169197201,1547695249,898509457,
     7771776657,1570150033,1548799121,1110377105,657537169,592932497,131
     89216273,1091558033,312610961,1532511889,858946705,842052241,188448
     91681,93920913,167990417,361859729,1078182033,572126865,1393831057,
     a1798464145,41195665/
      data r5d / 819646097,241501329,856898193,921006225,836478609,10059
     168529,1832129169,1570130065,622624401,1539749009,429189777,1988567
     2185,178083473,1842842769,943047313,28833937,1650339473,1915249809,
     31226218129,2133381265,744425105,1756970129,1278702225,1859758225,1
     4755307665,1368003729,1100499601,1355448465,388019857,748350609,691
     5610257,620451985,937528977,2045494417,52034193,1802252433,12563513
     677,964467857,1329255057,605882513,1344487057,1800238225,228305553,
     71326309521,1201936017,257838225,1044152977,1816049809,828698257,63
     82235153,1629313681,2075103377,224773777,775945361,1983787665,21034
     970225,1537646225,688968849,2107574929,1901150353,472348305,3713056
     a17,2000675473,1468143761/
      ks = kstrt - 2
      ndvp = (ndv - 1)/nvp + 1
      ndp = min0(ndvp,nvrp)
      mdp = nvp/min0(nvp,ndv)
      if (iflg.eq.0) go to 30
      do 20 k = 1, nblok
      id = (k + ks)/mdp
      do 10 j = 1, ndp
      l = ndp*id + j
      vran(j,k) = vran0(l)
      vran0(l) = 0.0d0
   10 continue
   20 continue
      iflg = 0
      return
   30 do 50 k = 1, nblok
      id = (k + ks)/mdp
      do 40 j = 1, ndp
      l = ndp*id + j
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r1(l) - (r1(l)/isc)*isc
      r3 = h1l*dble(r1(l)) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r2(l)/isc
      isc = r2(l) - i1*isc
      r0 = h1l*dble(r2(l)) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r2(l) = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r1(l) = r3 - dble(isc)*bsc
      temp = dsqrt(-2.0d0*dlog((dble(r1(l)) + dble(r2(l))*asc)*asc))
      isc = 65536
      asc = dble(isc)
      bsc = asc*asc
      i1 = r4(l) - (r4(l)/isc)*isc
      r3 = h2l*dble(r4(l)) + asc*h1u*dble(i1)
      i1 = r3/bsc
      r3 = r3 - dble(i1)*bsc
      bsc = 0.5d0*bsc
      i1 = r5(l)/isc
      isc = r5(l) - i1*isc
      r0 = h2l*dble(r5(l)) + asc*h1u*dble(isc)
      asc = 1.0d0/bsc
      isc = r0*asc
      r5(l) = r0 - dble(isc)*bsc
      r3 = r3 + (dble(isc) + 2.0d0*h1u*dble(i1))
      isc = r3*asc
      r4(l) = r3 - dble(isc)*bsc
      r0 = 6.28318530717959d0*((dble(r4(l)) + dble(r5(l))*asc)*asc)
      vran(j,k) = temp*dsin(r0)
      vran0(l) = temp*dcos(r0)
   40 continue
   50 continue
      iflg = 1
      return
      end
