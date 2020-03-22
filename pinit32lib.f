c 3d parallel PIC library for initialization
c with 2D domain decomposition
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: march 14, 2008
c-----------------------------------------------------------------------
      subroutine PISTR32(part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,
     1npy,npz,nx,ny,nz,idimp,npmax,mblok,nblok,idps,ipbc,ierr)
c for 3d code, this subroutine calculates initial particle co-ordinates
c and velocities with uniform density and maxwellian velocity with drift
c for distributed data with 2D spatial decomposition
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
c nps(m) = starting address of particles in partition m
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c npx/npy/npz = initial number of particles distributed in x/y/z
c direction
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mblok/nblok = number of particle partitions in y/z
c idps = number of partition boundaries
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
c ierr = (0,1) = (no,yes) error condition exists
c ranorm = gaussian random number with zero mean and unit variance
c     integer*8 npxy, npxyz
      double precision dnpxy, dnpxyz, dt1
      double precision ranorm
      double precision sum0, sum1, sum2, sum4, work4
      dimension part(idimp,npmax,mblok*nblok)
      dimension edges(idps,mblok*nblok)
      dimension npp(mblok*nblok), nps(mblok*nblok)
      dimension sum4(4), work4(4)
      ierr = 0
c particle distribution constant
      dnpxy = dble(npx)*dble(npy)
c set boundary values
      edgelx = 0.
      edgely = 0.
      edgelz = 0.
      at1 = float(nx)/float(npx)
      at2 = float(ny)/float(npy)
      at3 = float(nz)/float(npz)
      if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgelz = 1.
         at1 = float(nx-2)/float(npx)
         at2 = float(ny-2)/float(npy)
         at3 = float(nz-2)/float(npz)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgely = 1.
         edgelz = 0.
         at1 = float(nx-2)/float(npx)
         at2 = float(ny-2)/float(npy)
         at3 = float(nz)/float(npz)
      endif
      do 50 l = 1, npz
      zt = edgelz + at3*(float(l) - .5)
      do 40 k = 1, npy
      yt = edgely + at2*(float(k) - .5)
      do 30 j = 1, npx
c uniform density profile
      xt = edgelx + at1*(float(j) - .5)
c maxwellian velocity distribution
      vxt = vtx*ranorm()
      vyt = vty*ranorm()
      vzt = vtz*ranorm()
      do 20 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 10 my = 1, mblok
      m = my + moff
      if ((yt.ge.edges(1,m)).and.(yt.lt.edges(2,m))) then
         if ((zt.ge.edges(3,m)).and.(zt.lt.edges(4,m))) then
            npt = npp(m) + 1
            if (npt.le.npmax) then
               part(1,npt,m) = xt
               part(2,npt,m) = yt
               part(3,npt,m) = zt
               part(4,npt,m) = vxt
               part(5,npt,m) = vyt
               part(6,npt,m) = vzt
               npp(m) = npt
            else
               ierr = ierr + 1
            endif
         endif
      endif
   10 continue
   20 continue
   30 continue
   40 continue
   50 continue
      npxyzp = 0
c add correct drift
      sum4(1) = 0.0d0
      sum4(2) = 0.0d0
      sum4(3) = 0.0d0
      mnblok = mblok*nblok
      do 70 m = 1, mnblok
      sum0 = 0.0d0
      sum1 = 0.0d0
      sum2 = 0.0d0
      do 60 j = nps(m), npp(m)
      npxyzp = npxyzp + 1
      sum0 = sum0 + part(4,j,m)
      sum1 = sum1 + part(5,j,m)
      sum2 = sum2 + part(6,j,m)
   60 continue
      sum4(1) = sum4(1) + sum0
      sum4(2) = sum4(2) + sum1
      sum4(3) = sum4(3) + sum2
   70 continue
      sum4(4) = npxyzp
      call PDSUM(sum4,work4,4,1)
      dnpxyz = sum4(4)
      call PIMAX(ierr,iwork,1,1)
      dt1 = 1.0d0/dnpxyz
      sum4(1) = dt1*sum4(1) - vdx
      sum4(2) = dt1*sum4(2) - vdy
      sum4(3) = dt1*sum4(3) - vdz
      do 90 m = 1, mnblok
      do 80 j = nps(m), npp(m)
      part(4,j,m) = part(4,j,m) - sum4(1)
      part(5,j,m) = part(5,j,m) - sum4(2)
      part(6,j,m) = part(6,j,m) - sum4(3)
   80 continue
   90 continue
c process errors
      if (ierr.gt.0) then
         write (2,*) 'particle overflow error, ierr = ', ierr
      else if (dnpxyz.ne.(dnpxy*dble(npz))) then
         ierr = dnpxyz - dnpxy*dble(npz)
         write (2,*) 'particle distribution truncated, np = ', dnpxyz
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PVISTR32(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,n
     1pz,nx,ny,nz,idimp,npmax,mblok,nblok,ipbc,vranx,vrany,vranz,kstrt,n
     2vp,ndv,nvrp,ierr)
c for 3 code, this subroutine calculates initial particle co-ordinates
c and velocities with uniform density and maxwellian velocity with drift
c for distributed data with 2D spatial decomposition
c on input, the array npp contains the number of particles already
c stored in part, normally zero.
c on output, contains the number of particles stored in part.
c the total number of particles should be a multiple of ndv, if the
c number of processors is <= ndv, or else a multiple of the number of
c processors.
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c npp(m) = number of particles in partition m
c nps(m) = starting address of particles in partition m
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c npx/npy/npz = initial number of particles distributed in x/y/z
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mblok/nblok = number of particle partitions in y/z
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
c vranx,vrany,vranz = output arrays for parallel gaussian random numbers
c with zero mean and unit variance
c kstrt = starting data block number
c nvp = number of real or virtual processors
c ndv = total maximum number of random seeds, currently 256
c nvrp = number of parallel seeds per processor
c ierr = (0,1) = (no,yes) error condition exists
c npx*npy*npz should be a multiple of nvrp*nvp
c     integer*8 npxy, npxyz, ipp, nppv, loff, nng
      double precision vranx, vrany, vranz
      double precision dnpxy, dnpxyz, dnpx, dnppv, dloff, dnng, dt1
      double precision sum0, sum1, sum2, sum4, work4
      dimension part(idimp,npmax,mblok*nblok)
      dimension npp(mblok*nblok), nps(mblok*nblok)
      dimension vranx(nvrp,mblok*nblok), vrany(nvrp,mblok*nblok)
      dimension vranz(nvrp,mblok*nblok)
      dimension sum4(4), work4(4)
      ierr = 0
c particle distribution constants
      ks = kstrt - 2
      dnpx = dble(npx)
      dnpxy = dnpx*dble(npy)
      dnpxyz = dnpxy*dble(npz)
c npd = total number of seeds used
      npd = nvrp*nvp
c ipp = number of particles per seed
      ipp = dnpxyz/dble(npd)
c nppv = number of particles per processor
      nppv = nvrp*ipp
      dnppv = dble(nppv)
c adjacent mdp processors share the same seed
      mdp = nvp/min0(nvp,ndv)
c check unsupported particle number
      if (dnpxyz.ne.(dble(ipp)*dble(npd))) ierr = 1
      call PIMAX(ierr,iwork,1,1)
      if (ierr.gt.0) then
         if (kstrt.eq.1) then
           write (2,*) 'number of seeds not multiple of particle number'
     1, npd, dnpxyz
         endif
         return
      endif
c set boundary values
      edgelx = 0.
      edgely = 0.
      edgelz = 0.
      at1 = float(nx)/float(npx)
      at2 = float(ny)/float(npy)
      at3 = float(nz)/float(npz)
      if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgelz = 1.
         at1 = float(nx-2)/float(npx)
         at2 = float(ny-2)/float(npy)
         at3 = float(nz-2)/float(npz)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgely = 1.
         edgelz = 0.
         at1 = float(nx-2)/float(npx)
         at2 = float(ny-2)/float(npy)
         at3 = float(nz)/float(npz)
      endif
      mnblok = mblok*nblok
c outer loop over processor blocks which share the same seed
      do 50 ii = 1, mdp
c particles in each block get random numbers from the same seed
      do 40 n = 1, ipp
      call pranorm(vranx,kstrt,nvp,nvrp,nd,nvrp,mnblok)
      call pranorm(vrany,kstrt,nvp,nvrp,nd,nvrp,mnblok)
      call pranorm(vranz,kstrt,nvp,nvrp,nd,nvrp,mnblok)
c inner loop over independent seeds
      do 30 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 20 my = 1, mblok
      m = my + moff
      id = (m + ks)/mdp
      is = m + ks - id*mdp + 1
      dloff = dnppv*dble(m + ks)
      do 10 jj = 1, nvrp
c i = local particle number belonging to jj seed
      i = n + ipp*(jj - 1)
c nng = global particle number
      dnng = dble(i) + dloff
c j, k, l = used to determine spatial location of this particle
      l = (dnng - 1.0d0)/dnpxy
      dnng = dnng - dnpxy*dble(l)
      l = l + 1
      k = (dnng - 1.0d0)/dnpx
      j = dnng - dnpx*dble(k)
      k = k + 1
      vxt = vtx*vranx(jj,m)
      vyt = vty*vrany(jj,m)
      vzt = vtz*vranz(jj,m)
c nth group keeps nth block of random numbers
      if (ii.eq.is) then
         npt = nps(m) + i - 1
         if (npt.le.npmax) then
c uniform density profile
            part(1,npt,m) = edgelx + at1*(float(j) - .5)
            part(2,npt,m) = edgely + at2*(float(k) - .5)
            part(3,npt,m) = edgelz + at3*(float(l) - .5)
c maxwellian velocity distribution
            part(4,npt,m) = vxt
            part(5,npt,m) = vyt
            part(6,npt,m) = vzt
         else
            ierr = ierr + 1
         endif
      endif
   10 continue
c update particle number
      if (ii.eq.is) npp(m) = npp(m) + nvrp
   20 continue
   30 continue
   40 continue
   50 continue
      npxyzp = 0
c add correct drift
      sum4(1) = 0.0d0
      sum4(2) = 0.0d0
      sum4(3) = 0.0d0
      do 70 m = 1, mnblok
      sum0 = 0.0d0
      sum1 = 0.0d0
      sum2 = 0.0d0
      do 60 j = nps(m), npp(m)
      npxyzp = npxyzp + 1
      sum0 = sum0 + part(4,j,m)
      sum1 = sum1 + part(5,j,m)
      sum2 = sum2 + part(6,j,m)
   60 continue
      sum4(1) = sum4(1) + sum0
      sum4(2) = sum4(2) + sum1
      sum4(3) = sum4(3) + sum2
   70 continue
      sum4(4) = npxyzp
      call PDSUM(sum4,work4,4,1)
      dnpxyz = sum4(4)
      call PIMAX(ierr,iwork,1,1)
      dt1 = 1.0d0/dnpxyz
      sum4(1) = dt1*sum4(1) - vdx
      sum4(2) = dt1*sum4(2) - vdy
      sum4(3) = dt1*sum4(3) - vdz
      do 90 m = 1, mnblok
      do 80 j = nps(m), npp(m)
      part(4,j,m) = part(4,j,m) - sum4(1)
      part(5,j,m) = part(5,j,m) - sum4(2)
      part(6,j,m) = part(6,j,m) - sum4(3)
   80 continue
   90 continue
c process errors
      if (ierr.gt.0) then
         if (kstrt.eq.1) then
            write (2,*) 'particle overflow error, ierr = ', ierr
         endif
      else if (dnpxyz.ne.(dnpxy*dble(npz))) then
         ierr = dnpxyz - dnpxy*dble(npz)
         if (kstrt.eq.1) then
            write (2,*) 'particle distribution truncated, np = ', dnpxyz
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PLDISTR32(part,nps,anlx,anly,anlz,npx,npy,npz,nx,ny,nz,
     1idimp,npmax,mblok,nblok,kstrt,nvp,ipbc,ierr)
c for 3d code, this subroutine calculates initial particle co-ordinates
c with the following tri-linear density profile:
c n(x,y,z) = n(x)*n(y)*n(z), 
c where n(x) = n0x*(1. + anlx*(x/nx - .5))
c   and n(y) = n0y*(1. + anly*(y/ny - .5))
c   and n(z) = n0z*(1. + anlz*(z/nz - .5))
c n0x = npx/(nx - 2*edgelx), n0y = npy/(ny - 2*edgely)
c and n0z = npz/(nz - 2*edgelz)
c for distributed data with 2D spatial decomposition
c particles are not necessarily in the correct processor.
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c nps(m) = starting address of particles in partition m
c anlx/anly/anlz = initial linear density weight in x/y/z direction
c npx/npy/npz = initial number of particles distributed in x/y/z
c direction
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mblok/nblok = number of particle partitions in y/z
c kstrt = starting data block number
c nvp = number of real or virtual processors
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
c ierr = (0,1) = (no,yes) error condition exists
      implicit none
      integer npx, npy, npz, nx, ny, nz, idimp, npmax, mblok, nblok
      integer kstrt, nvp, ipbc, ierr
      real anlx, anly, anlz
      integer nps
      real part
      dimension part(idimp,npmax,mblok*nblok)
      dimension nps(mblok*nblok)
c local data
c     integer*8 npxy, npxyz, nppv, loff, nng
      double precision dnpxy, dnpxyz, dnpx, dnppv, dloff, dnng
      integer nppv, ks, j, k, l, m, my, mz, n, moff, iwork, noff
      real edgelx, edgely, edgelz, at1, at2, at3, bt1, bt2, bt3
      real antx, anty, antz, xt, yt, zt
      ierr = 0
c particle distribution constants
      ks = kstrt - 2
      dnpx = dble(npx)
      dnpxy = dnpx*dble(npy)
      nppv = dnpxy*dble(npz)/dble(nvp)
      nppv = min(nppv,npmax)
      dnppv = dble(nppv)
      dnpxyz = dnppv*dble(nvp)
c check for errors
      if (dnpxyz.ne.(dnpxy*dble(npz))) ierr = 1
      call PIMAX(ierr,iwork,1,1)
      if (ierr.gt.0) then
         ierr = dnpxyz - dnpxy*dble(npz)
         if (kstrt.eq.1) then
            write (2,*) 'particle distribution truncated, np = ', dnpxyz
         endif
         return
      endif
c set boundary values
      edgelx = 0.
      edgely = 0.
      edgelz = 0.
      at1 = float(nx)/float(npx)
      at2 = float(ny)/float(npy)
      at3 = float(nz)/float(npz)
      if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgelz = 1.
         at1 = float(nx-2)/float(npx)
         at2 = float(ny-2)/float(npy)
         at3 = float(nz-2)/float(npz)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgely = 1.
         edgelz = 0.
         at1 = float(nx-2)/float(npx)
         at2 = float(ny-2)/float(npy)
         at3 = float(nz)/float(npz)
      endif
      if (anlz.ne.0.) then
         antz = anlz/float(nz)
         at3 = 2.*antz*at3
         bt3 = 1. - .5*antz*(float(nz) - 2.*edgelz)
      endif
      if (anly.ne.0.) then
         anty = anly/float(ny)
         at2 = 2.*anty*at2
         bt2 = 1. - .5*anty*(float(ny) - 2.*edgely)
      endif
      if (anlx.ne.0.) then
         antx = anlx/float(nx)
         at1 = 2.*antx*at1
         bt1 = 1. - .5*antx*(float(nx) - 2.*edgelx)
      endif
c uniform density profile
      do 30 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 20 my = 1, mblok
      m = my + moff
      dloff = dnppv*dble(m + ks)
      noff = nps(m) - 1
      do 10 n = 1, nppv
c j, k, l = used to determine spatial location of this particle
      dnng = dble(n) + dloff
      l = (dnng - 1.0d0)/dnpxy
      dnng = dnng - dnpxy*dble(l)
      l = l + 1
      k = (dnng - 1.0d0)/dnpx
      j = dnng - dnpx*dble(k)
      k = k + 1
c linear density in x
      if (anlx.ne.0.) then
         xt = edgelx + (sqrt(bt1*bt1 + at1*(float(j) - .5)) - bt1)/antx
c uniform density in x
      else
         xt = edgelx + at1*(float(j) - .5)
      endif
c linear density in y
      if (anly.ne.0.) then
         yt = edgely + (sqrt(bt2*bt2 + at2*(float(k) - .5)) - bt2)/anty
c uniform density in y
      else
         yt = edgely + at2*(float(k) - .5)
      endif
c linear density in z
      if (anlz.ne.0.) then
         zt = edgelz + (sqrt(bt3*bt3 + at3*(float(l) - .5)) - bt3)/antz
c uniform density in z
      else
         zt = edgelz + at3*(float(l) - .5)
      endif
      part(1,n+noff,m) = xt
      part(2,n+noff,m) = yt
      part(3,n+noff,m) = zt
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFDISTR32(part,nps,fnx,argx1,argx2,argx3,fny,argy1,argy
     12,argy3,fnz,argz1,argz2,argz3,npx,npy,npz,nx,ny,nz,idimp,npmax,mbl
     2ok,nblok,kstrt,nvp,ipbc,ierr)
c for 3d code, this subroutine calculates initial particle co-ordinates
c with general density profile n(x,y) = n(x)*n(y)*n(z), 
c where density in x is given by n(x) = fnx(x,argx1,argx2,argx3,0)
c and integral of the density is given by = fnx(x,argx1,argx2,argx3,1)
c and where density in y is given by n(y) = fny(y,argy1,argy2,argy3,0)
c and integral of the density is given by = fny(y,argy1,argy2,argy3,1)
c and where density in z is given by n(z) = fnz(z,argz1,argz2,argz3,0)
c and integral of the density is given by = fnz(z,argz1,argz2,argz3,1)
c for distributed data with 2D spatial decomposition.
c particles are not necessarily in the correct processor, and the total
c number of particles should be a multiple of the number of processors.
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c nps(m) = starting address of particles in partition m
c fnx/fny/fnz = density and density integral function in x/y/z direction
c argx1,argx2,argx3 = arguments to fnx
c argy1,argy2,argy3 = arguments to fny
c argz1,argz2,argz3 = arguments to fnz
c npx/npy/npz = initial number of particles distributed in x/y/z
c direction
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mblok/nblok = number of particle partitions in y/z
c kstrt = starting data block number
c nvp = number of real or virtual processors
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
c ierr = (0,1) = (no,yes) error condition exists
      implicit none
      integer nps, npx, npy, npz, nx, ny, nz, idimp, npmax, mblok, nblok
      integer kstrt, nvp, ipbc, ierr
      real argx1, argx2, argx3, argy1, argy2, argy3, argz1, argz2, argz3
      real part
      dimension part(idimp,npmax,mblok*nblok)
      dimension nps(mblok*nblok)
      real fnx, fny, fnz
      external fnx, fny, fnz
c local data
c     integer*8 npxy, npxyz, nppv, loff, koff, nng
      double precision dnpxy, dnpxyz, dnpx, dnppv, dloff, dkoff, dnng
      integer nppv, mppv
      integer ks, lc, kc, jc, i, j, k, l, m, my, mz, moff, n, nn, mm
      integer noff, imax, iwork
      real eps, big, edgelx, edgely, edgelz, anx, any, anz
      real bnx, bny, bnz, xt, yt, zt, xt0, yt0, zt0, x0, y0, z0
      real xn, yn, zn, f, fp
      ierr = 0
c particle distribution constants
      ks = kstrt - 2
      dnpx = dble(npx)
      dnpxy = dnpx*dble(npy)
      nppv = dnpxy*dble(npz)/dble(nvp)
      nppv = min(nppv,npmax)
      dnppv = dble(nppv)
      dnpxyz = dnppv*dble(nvp)
      mppv = min(dnpxy,dnppv)
c eps = convergence criterion
      imax = 20
      eps = 0.0002
      big = 0.5
c check for errors
      if (dnpxyz.ne.(dnpxy*dble(npz))) ierr = 1
      call PIMAX(ierr,iwork,1,1)
      if (ierr.gt.0) then
         ierr = dnpxyz - dnpxy*dble(npz)
         if (kstrt.eq.1) then
            write (2,*) 'particle distribution truncated, np = ', dnpxyz
         endif
         return
      endif
c set boundary values
      edgelx = 0.
      edgely = 0.
      edgelz = 0.
      if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgelz = 1.
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgely = 1.
         edgelz = 0.
      endif
c find normalization for function
      anx = float(nx) - edgelx
      any = float(ny) - edgely
      anz = float(nz) - edgelz
      x0 = fnx(edgelx,argx1,argx2,argx3,1)
      y0 = fny(edgely,argy1,argy2,argy3,1)
      z0 = fnz(edgelz,argz1,argz2,argz3,1)
      bnx = float(npx)/(fnx(anx,argx1,argx2,argx3,1) - x0)
      bny = float(npy)/(fny(any,argy1,argy2,argy3,1) - y0)
      bnz = float(npz)/(fnz(anz,argz1,argz2,argz3,1) - z0)
      x0 = bnx*x0 - .5
      y0 = bny*y0 - .5
      z0 = bnz*z0 - .5
c calculate density profile
      do 140 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 130 my = 1, mblok
      m = my + moff
      dloff = dnppv*dble(m + ks)
      noff = nps(m) - 1
c integrate to find starting point in z
      lc = dloff/dnpxy
      zt0 = edgelz
      zt = zt0 + 0.5/(bnz*fnz(zt0,argz1,argz2,argz3,0))
      do 20 l = 1, lc
      zn = float(l) + z0
c guess next value for zt
      if (l.gt.1) zt = zt + 1.0/(bnz*fnz(zt,argz1,argz2,argz3,0))
      zt = max(edgelz,min(zt,anz))
      i = 0
   10 f = bnz*fnz(zt,argz1,argz2,argz3,1) - zn
c find improved value for zt
      if (abs(f).ge.eps) then
c newton's method
         if (abs(f).lt.big) then
            fp = bnz*fnz(zt,argz1,argz2,argz3,0)
            zt0 = zt
            zt = zt - f/fp
            zt = max(edgelz,min(zt,anz))
c bisection method
         else if (f.gt.0.) then
            fp = .5*(zt - zt0)
            zt = zt0 + fp
         else
            fp = zt - zt0
c           zt0 = zt
            zt = zt + fp
         endif
         i = i + 1
         if (i.lt.imax) go to 10
         if (ierr.eq.0) ierr = 1
      endif
      zt0 = zt
   20 continue
c quit if error
      if (ierr.ne.0) go to 130
c integrate to find starting point in y
      dkoff = dloff - dnpxy*dble(lc)
      kc = dkoff/dnpx
      yt0 = edgely
      yt = yt0 + 0.5/(bny*fny(yt0,argy1,argy2,argy3,0))
      do 40 k = 1, kc
      yn = float(k) + y0
c guess next value for yt
      if (k.gt.1) yt = yt + 1.0/(bny*fny(yt,argy1,argy2,argy3,0))
      yt = max(edgely,min(yt,any))
      i = 0
   30 f = bny*fny(yt,argy1,argy2,argy3,1) - yn
c find improved value for yt
      if (abs(f).ge.eps) then
c newton's method
         if (abs(f).lt.big) then
            fp = bny*fny(yt,argy1,argy2,argy3,0)
            yt0 = yt
            yt = yt - f/fp
            yt = max(edgely,min(yt,any))
c bisection method
         else if (f.gt.0.) then
            fp = .5*(yt - yt0)
            yt = yt0 + fp
         else
            fp = yt - yt0
c           yt0 = yt
            yt = yt + fp
         endif
         i = i + 1
         if (i.lt.imax) go to 30
         if (ierr.eq.0) ierr = 2
      endif
      yt0 = yt
   40 continue
c quit if error
      if (ierr.ne.0) go to 130
c integrate to find starting point in x
      jc = dkoff - dnpx*dble(kc)
      xt0 = edgelx
      xt = xt0 + 0.5/(bnx*fnx(xt0,argx1,argx2,argx3,0))
      do 60 j = 1, jc
      xn = float(j) + x0
c guess next value for xt
      if (j.gt.1) xt = xt + 1.0/(bnx*fnx(xt,argx1,argx2,argx3,0))
      xt = max(edgelx,min(xt,anx))
      i = 0
   50 f = bnx*fnx(xt,argx1,argx2,argx3,1) - xn
c find improved value for xt
      if (abs(f).ge.eps) then
c newton's method
         if (abs(f).lt.big) then
            fp = bnx*fnx(xt,argx1,argx2,argx3,0)
            xt0 = xt
            xt = xt - f/fp
            xt = max(edgelx,min(xt,anx))
c bisection method
         else if (f.gt.0.) then
            fp = .5*(xt - xt0)
            xt = xt0 + fp
         else
            fp = xt - xt0
c           xt0 = xt
            xt = xt + fp
         endif
         i = i + 1
         if (i.lt.imax) go to 50
         if (ierr.eq.0) ierr = 3
      endif
      xt0 = xt
   60 continue
c quit if error
      if (ierr.ne.0) go to 130
c density profile in x
      do 80 n = 1, min(npx,nppv)
c j, k, l = used to determine spatial location of this particle
      dnng = dble(n) + dloff
      l = (dnng - 1.0d0)/dnpxy
      dnng = dnng - dnpxy*dble(l)
      l = l + 1
      k = (dnng - 1.0d0)/dnpx
      j = dnng - dnpx*dble(k)
      k = k + 1
      xn = float(j) + x0
c guess next value for xt
      if (j.eq.1) then
         xt0 = edgelx
         xt = xt0 + 0.5/(bnx*fnx(xt0,argx1,argx2,argx3,0))
      else
         xt = xt + 1.0/(bnx*fnx(xt,argx1,argx2,argx3,0))
      endif
      xt = max(edgelx,min(xt,anx))
      i = 0
   70 f = bnx*fnx(xt,argx1,argx2,argx3,1) - xn
c find improved value for xt
      if (abs(f).ge.eps) then
c newton's method
         if (abs(f).lt.big) then
            fp = bnx*fnx(xt,argx1,argx2,argx3,0)
            xt0 = xt
            xt = xt - f/fp
            xt = max(edgelx,min(xt,anx))
c bisection method
         else if (f.gt.0.) then
            fp = .5*(xt - xt0)
            xt = xt0 + fp
         else
            fp = xt - xt0
c           xt0 = xt
            xt = xt + fp
         endif
         i = i + 1
         if (i.lt.imax) go to 70
         if (ierr.eq.0) ierr = 4
      endif
      part(1,n+noff,m) = xt
      xt0 = xt
   80 continue
c quit if error
      if (ierr.ne.0) go to 130
c density profile in y
      kc = kc + 1
      do 100 n = 1, mppv
c j, k, l = used to determine spatial location of this particle
      dnng = dble(n) + dloff
      l = (dnng - 1.0d0)/dnpxy
      dnng = dnng - dnpxy*dble(l)
      l = l + 1
      k = (dnng - 1.0d0)/dnpx
      k = k + 1
      if (k.eq.kc) then
         yn = float(k) + y0
c guess next value for yt
         if (k.eq.1) then
            yt = edgely
            yt = yt + 0.5/(bny*fny(yt,argy1,argy2,argy3,0))
         else
            yt = yt + 1.0/(bny*fny(yt,argy1,argy2,argy3,0))
         endif
         yt = max(edgely,min(yt,any))
         i = 0
   90    f = bny*fny(yt,argy1,argy2,argy3,1) - yn
c find improved value for yt
         if (abs(f).ge.eps) then
c newton's method
            if (abs(f).lt.big) then
               fp = bny*fny(yt,argy1,argy2,argy3,0)
               yt0 = yt
               yt = yt - f/fp
               yt = max(edgely,min(yt,any))
c bisection method
            else if (f.gt.0.) then
               fp = .5*(yt - yt0)
               yt = yt0 + fp
            else
               fp = yt - yt0
c              yt0 = yt
               yt = yt + fp
            endif
            i = i + 1
            if (i.lt.imax) go to 90
            if (ierr.eq.0) ierr = 5
         endif
c store co-ordinates
         kc = kc + 1
         if (kc.gt.npy) kc = 1
         mm = ((n + jc - 1)/npx)*npx - jc + 1
         if (mm.lt.1) mm = 1
         part(2,mm+noff,m) = yt
         yt0 = yt
      endif
  100 continue
c quit if error
      if (ierr.ne.0) go to 130
c density profile in z
      lc = lc + 1
      do 120 n = 1, nppv
c j, k, l = used to determine spatial location of this particle
      dnng = dble(n) + dloff
      l = (dnng - 1.0d0)/dnpxy
      dnng = dnng - dnpxy*dble(l)
      l = l + 1
      k = (dnng - 1.0d0)/dnpx
      j = dnng - dnpx*dble(k)
      k = k + 1
      nn = (dble(n) - 1.0d0)/dnpxy
      nn = dble(n) - dble(nn)*dnpxy
      mm = ((nn + jc - 1)/npx)*npx - jc + 1
      nn = nn - ((nn - 1)/npx)*npx
      if (mm.lt.1) mm = 1
      if (l.eq.lc) then
         zn = float(l) + z0
c guess next value for zt
         if (l.gt.1) zt = zt + 1.0/(bnz*fnz(zt,argz1,argz2,argz3,0))
         zt = max(edgelz,min(zt,anz))
         i = 0
  110    f = bnz*fnz(zt,argz1,argz2,argz3,1) - zn
c find improved value for zt
         if (abs(f).ge.eps) then
c newton's method
            if (abs(f).lt.big) then
               fp = bnz*fnz(zt,argz1,argz2,argz3,0)
               zt0 = zt
               zt = zt - f/fp
               zt = max(edgelz,min(zt,anz))
c bisection method
            else if (f.gt.0.) then
               fp = .5*(zt - zt0)
               zt = zt0 + fp
            else
               fp = zt - zt0
c              zt0 = zt
               zt = zt + fp
            endif
            i = i + 1
            if (i.lt.imax) go to 110
            if (ierr.eq.0) ierr = 6
         endif
         lc = lc + 1
         zt0 = zt
      endif
      part(1,n+noff,m) = part(1,nn+noff,m)
      part(2,n+noff,m) = part(2,mm+noff,m)
      part(3,n+noff,m) = zt
  120 continue
  130 continue
  140 continue
c handle errors
      call PIMAX(ierr,iwork,1,1)
      if (ierr.gt.0) then
         ierr = -ierr
         if (kstrt.eq.1) then
            write (2,*) 'newton iteration max exceeded, ierr = ', ierr
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PRDISTR32(part,nps,fnx,argx1,argx2,argx3,fny,argy1,argy
     12,argy3,fnz,argz1,argz2,argz3,npx,npy,npz,nx,ny,nz,idimp,npmax,mnb
     2lok,kstrt,nvp,ipbc)
c for 3d code, this subroutine randomizes initial particle co-ordinates
c which have been previously calculated using a general density profile
c n(x,y,z) = n(x)*n(y)*n(z),  
c where density in x is given by n(x) = fnx(x,argx1,argx2,argx3,0)
c and integral of the density is given by = fnx(x,argx1,argx2,argx3,1)
c and where density in y is given by n(y) = fny(y,argy1,argy2,argy3,0)
c and integral of the density is given by = fny(y,argy1,argy2,argy3,1)
c and where density in z is given by n(z) = fnz(z,argz1,argz2,argz3,0)
c and integral of the density is given by = fnz(z,argz1,argz2,argz3,1)
c for distributed data with 2D spatial decomposition.
c particles are not necessarily in the correct processor.
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c nps(m) = starting address of particles in partition m
c fnx/fny/fnz = density and density integral function in x/y/z direction
c argx1,argx2,argx3 = arguments to fnx
c argy1,argy2,argy3 = arguments to fny
c argz1,argz2,argz3 = arguments to fnz
c npx/npy/npz = initial number of particles distributed in x/y/z
c direction
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions
c kstrt = starting data block number
c nvp = number of real or virtual processors
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
c with spatial decomposition
      implicit none
      integer npx, npy, npz, nx, ny, nz, idimp, npmax, mnblok, kstrt
      integer nvp, ipbc
      integer nps
      real argx1, argx2, argx3, argy1, argy2, argy3, argz1, argz2, argz3
      real part
      dimension part(idimp,npmax,mnblok)
      dimension nps(mnblok)
      real fnx, fny, fnz
      external fnx, fny, fnz
c local data
c     integer*8 npxy, nppv, loff, lc, kc, nn
      double precision dnpxy, dnpx, dnppv, dloff, dlc, dkc, nng
      integer nppv, ks, j, k, l, m, n, noff
      real edgelx, edgely, edgelz, anx, any, anz, bnx, bny, bnz
      real xt, yt, zt, xt0, yt0, zt0, x0, y0, z0, xn, yn, zn
      double precision randum
c particle distribution constants
      ks = kstrt - 2
      dnpx = dble(npx)
      dnpxy = dnpx*dble(npy)
      nppv = dnpxy*dble(npz)/dble(nvp)
      nppv = min(nppv,npmax)
      dnppv = dble(nppv)
c set boundary values
      edgelx = 0.
      edgely = 0.
      edgelz = 0.
      if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgelz = 1.
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgely = 1.
         edgelz = 0.
      endif
c find normalization for function
      anx = float(nx) - edgelx
      any = float(ny) - edgely
      anz = float(nz) - edgelz
      x0 = fnx(edgelx,argx1,argx2,argx3,1)
      y0 = fny(edgely,argy1,argy2,argy3,1)
      z0 = fnz(edgelz,argz1,argz2,argz3,1)
      bnx = float(npx)/(fnx(anx,argx1,argx2,argx3,1) - x0)
      bny = float(npy)/(fny(any,argy1,argy2,argy3,1) - y0)
      bnz = float(npz)/(fnz(anz,argz1,argz2,argz3,1) - z0)
c randomize co-ordinates
      do 40 l = 1, npz
      dlc = dnpxy*dble(l - 1)
      do 30 k = 1, npy
      dkc = dnpx*dble(k - 1) + dlc
      do 20 j = 1, npx
      nng = dble(j) + dkc
      xn = randum()
      yn = randum()
      zn = randum()
      do 10 m = 1, mnblok
      dloff = dnppv*dble(m + ks)
      n = nng - dloff
      if ((n.gt.0).and.(n.le.nppv)) then
         noff = nps(m) - 1
         xt = part(1,n+noff,m)
         xt0 = 1.0/(bnx*fnx(xt,argx1,argx2,argx3,0))
         xn = xt + xt0*(xn - 0.5)
         xn = max(edgelx,min(xn,anx))
         part(1,n+noff,m) = xn
         yt = part(2,n+noff,m)
         yt0 = 1.0/(bny*fny(yt,argy1,argy2,argy3,0))
         yn = yt + yt0*(yn - 0.5)
         yn = max(edgely,min(yn,any))
         part(2,n+noff,m) = yn
         zt = part(3,n+noff,m)
         zt0 = 1.0/(bnz*fnz(zt,argz1,argz2,argz3,0))
         zn = zt + zt0*(zn - 0.5)
         zn = max(edgelz,min(zn,anz))
         part(3,n+noff,m) = zn
      endif
   10 continue
   20 continue
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PVRDISTR32(part,nps,fnx,argx1,argx2,argx3,fny,argy1,arg
     1y2,argy3,fnz,argz1,argz2,argz3,npx,npy,npz,nx,ny,nz,idimp,npmax,mn
     2blok,vranx,vrany,vranz,kstrt,nvp,ndv,nvrp,ipbc,ierr)
c for 3d code, this subroutine randomizes initial particle co-ordinates
c which have been previously calculated using a general density profile
c n(x,y,z) = n(x)*n(y)*n(z),  
c where density in x is given by n(x) = fnx(x,argx1,argx2,argx3,0)
c and integral of the density is given by = fnx(x,argx1,argx2,argx3,1)
c and where density in y is given by n(y) = fny(y,argy1,argy2,argy3,0)
c and integral of the density is given by = fny(y,argy1,argy2,argy3,1)
c and where density in z is given by n(z) = fnz(z,argz1,argz2,argz3,0)
c and integral of the density is given by = fnz(z,argz1,argz2,argz3,1)
c for distributed data with 2D spatial decomposition
c the exact spatial positions are randomized locally.
c particles are not necessarily in the correct processor.
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c nps(m) = starting address of particles in partition m
c fnx/fny/fnz = density and density integral function in x/y/z direction
c argx1,argx2,argx3 = arguments to fnx
c argy1,argy2,argy3 = arguments to fny
c argz1,argz2,argz3 = arguments to fnz
c npx/npy/npz = initial number of particles distributed in x/y/z
c direction
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions
c vranx,vrany,vranz = output arrays for parallel gaussian random numbers
c with zero mean and unit variance
c kstrt = starting data block number
c nvp = number of real or virtual processors
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
c with spatial decomposition
      implicit none
      integer npx, npy, npz, nx, ny, nz, idimp, npmax, mnblok, kstrt
      integer nvp, ndv, nvrp, ipbc, ierr
      integer nps
      real argx1, argx2, argx3, argy1, argy2, argy3, argz1, argz2, argz3
      real part
      double precision vranx, vrany, vranz
      dimension part(idimp,npmax,mnblok)
      dimension nps(mnblok)
      dimension vranx(nvrp,mnblok), vrany(nvrp,mnblok)
      dimension vranz(nvrp,mnblok)
      real fnx, fny, fnz
      external fnx, fny, fnz
c local data
c     integer*8 npxyz, ipp
      double precision dnpxyz
      integer ipp, ks, i, m, n, noff, npd, mdp, ii, jj, nd, id, is
      integer iwork
      real edgelx, edgely, edgelz, anx, any, anz, bnx, bny, bnz
      real xt, yt, zt, xt0, yt0, zt0, x0, y0, z0, xn, yn, zn
      ierr = 0
c particle distribution constants
      ks = kstrt - 2
      dnpxyz = dble(npx)*dble(npy)*dble(npz)
c npd = total number of seeds used
      npd = nvrp*nvp
c ipp = number of particles per seed
      ipp = dnpxyz/dble(npd)
c adjacent mdp processors share the same seed
      mdp = nvp/min0(nvp,ndv)
c check for errors
      if (dnpxyz.ne.(dble(ipp)*dble(npd))) ierr = 1
      call PIMAX(ierr,iwork,1,1)
      if (ierr.gt.0) then
         if (kstrt.eq.1) then
           write (2,*) 'number of seeds not multiple of particle number',
     1, npd, dnpxyz
         endif
         return
      endif
c set boundary values
      edgelx = 0.
      edgely = 0.
      edgelz = 0.
      if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgelz = 1.
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgely = 1.
         edgelz = 0.
      endif
c find normalization for function
      anx = float(nx) - edgelx
      any = float(ny) - edgely
      anz = float(nz) - edgelz
      x0 = fnx(edgelx,argx1,argx2,argx3,1)
      y0 = fny(edgely,argy1,argy2,argy3,1)
      z0 = fnz(edgelz,argz1,argz2,argz3,1)
      bnx = float(npx)/(fnx(anx,argx1,argx2,argx3,1) - x0)
      bny = float(npy)/(fny(any,argy1,argy2,argy3,1) - y0)
      bnz = float(npz)/(fnz(anz,argz1,argz2,argz3,1) - z0)
c randomize co-ordinates
c outer loop over processor blocks which share the same seed
      do 40 ii = 1, mdp
c particles in each block get random numbers from the same seed
      do 30 n = 1, ipp
      call prandom(vranx,kstrt,nvp,nvrp,nd,nvrp,mnblok)
      call prandom(vrany,kstrt,nvp,nvrp,nd,nvrp,mnblok)
      call prandom(vranz,kstrt,nvp,nvrp,nd,nvrp,mnblok)
c inner loop over independent seeds
      do 20 m = 1, mnblok
      noff = nps(m) - 1
      id = (m + ks)/mdp
      is = m + ks - id*mdp + 1
      do 10 jj = 1, nvrp
c i = local particle number belonging to jj seed
      i = n + ipp*(jj - 1)
      xn = vranx(jj,m)
      yn = vrany(jj,m)
      zn = vranz(jj,m)
c nth group keeps nth block of random numbers
      if (ii.eq.is) then
         xt = part(1,i+noff,m)
         xt0 = 1.0/(bnx*fnx(xt,argx1,argx2,argx3,0))
         xn = xt + xt0*(xn - 0.5)
         xn = max(edgelx,min(xn,anx))
         part(1,i+noff,m) = xn
         yt = part(2,i+noff,m)
         yt0 = 1.0/(bny*fny(yt,argy1,argy2,argy3,0))
         yn = yt + yt0*(yn - 0.5)
         yn = max(edgely,min(yn,any))
         part(2,i+noff,m) = yn
         zt = part(3,i+noff,m)
         zt0 = 1.0/(bnz*fnz(zt,argz1,argz2,argz3,0))
         zn = zt + zt0*(zn - 0.5)
         zn = max(edgelz,min(zn,anz))
         part(3,i+noff,m) = zn
      endif
   10 continue
   20 continue
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      function FLDISTR1(x,anlx,anxi,shift,intg)
c this function calculates either a density function or its integral
c for a linear density profile.  Used in initializing particle
c coordinates.
c if intg = 0, n(x) = 1. + anlx*(x*anxi - shift)
c if intg = 1, n(x) = x + .5*anlx*x*(x*anxi - 2.*shift)
      implicit none
      integer intg
      real x, anlx, anxi, shift
c local data
      real FLDISTR1, f
      if (intg.eq.0) then
         f = 1.0 + anlx*(x*anxi - shift)
      else if (intg.eq.1) then
         if (anxi.eq.0.) then
            f = x
         else
            f = x + .5*anlx*x*(x*anxi - 2.*shift)
         endif
      else
         f = -1.0
      endif
      if (f.lt.0.) write (2,*) 'FLDISTR1 Error: f = ', f
      FLDISTR1 = f
      return
      end
c-----------------------------------------------------------------------
      function FSDISTR1(x,ans,dkx,phase,intg)
c this function calculates either a density function or its integral
c for a sinusoidal density profile.  Used in initializing particle
c coordinates.
c if intg = 0, n(x) = 1.0 + ans*sin(dkx*x - phase)
c if intg = 1, n(x) = x - (ans/dkx)*(cos(dkx*x - phase) - cos(phase))
      implicit none
      integer intg
      real x, ans, dkx, phase
c local data
      real FSDISTR1, f
      if (intg.eq.0) then
         f = 1.0 + ans*sin(dkx*x - phase)
      else if (intg.eq.1) then
         if (dkx.eq.0.) then
            f = x - ans*sin(phase)*x
         else
            f = x - (ans/dkx)*(cos(dkx*x - phase) - cos(phase))
         endif
      else
         f = -1.0
      endif
      if (f.lt.0.) write (2,*) 'FSDISTR1 Error: f = ', f
      FSDISTR1 = f
      return
      end
c-----------------------------------------------------------------------
      function FGDISTR1(x,ang,wi,x0,intg)
c this function calculates either a density function or its integral
c for a gaussian density profile.  Used in initializing particle
c coordinates.
c if intg = 0, n(x) = 1.0 + ang*exp(-((x-x0)*wi)**2/2.)
c if intg = 1, n(x) = x + (ang*sqrt(pi/2)/wi)*
c                         (erf((x-x0)*wi/sqrt(2)) + erf(x0*wi/sqrt(2)))
      implicit none
      integer intg
      real x, ang, x0, wi
c local data
      real FGDISTR1, f, sqrt2i, sqtpih, aw, t, erfn
      external erfn
      data sqrt2i, sqtpih /0.7071067811865476,1.253314137397325/
      save sqrt2i, sqtpih
      aw = wi*sqrt2i
      t = (x - x0)*aw
      if (intg.eq.0) then
         if (abs(t).lt.8.) then
            f = 1.0 + ang*exp(-t**2)
         else
            f = 1.0
         endif
      else if (intg.eq.1) then
         if (wi.eq.0.) then
            f = (1.0 + ang)*x
         else
            f = x + (ang*sqtpih/wi)*(erfn(t) + erfn(x0*aw))
         endif
      else
         f = -1.0
      endif
      if (f.lt.0.) write (2,*) 'FGDISTR1 Error: f = ', f
      FGDISTR1 = f
      return
      end
c-----------------------------------------------------------------------
      function FHDISTR1(x,anh,wi,x0,intg)
c this function calculates either a density function or its integral
c for a hyperbolic secant squared density profile.  Used in initializing
c particle coordinates.
c if intg = 0, n(x) = 1.0 + anh*sech((x-x0)*wi)**2
c if intg = 1, n(x) = x + (anh/wi)*(tanh((x-x0)*wi) + tanh(x0*wi))
      implicit none
      integer intg
      real x, anh, x0, wi
c local data
      real FHDISTR1, f, g, t, u
      t = (x - x0)*wi
      if (intg.eq.0) then
         if (abs(t).lt.32.) then
            u = exp(-abs(t))
            f = 1.0 + anh*(2.*u/(1.0 + u*u))**2
         else
            f = 1.0
         endif
      else if (intg.eq.1) then
         if (wi.eq.0.) then
            f = (1.0 + anh)*x
         else
            if (abs(t).lt.32.) then
               u = exp(-abs(t))**2
               f = (1.0 - u)/(1.0 + u)
            else
               f = 1.0
            endif
            if (t.lt.0.) f = -f
            t = x0*wi
            if (abs(t).lt.32.) then
               u = exp(-abs(t))**2
               g = (1.0 - u)/(1.0 + u)
            else
               g = 1.0
            endif
            if (t.lt.0.) g = -g
            f = x + (anh/wi)*(f + g)
         endif
      else
         f = -1.0
      endif
      if (f.lt.0.) write (2,*) 'FHDISTR1 Error: f = ', f
      FHDISTR1 = f
      return
      end
c-----------------------------------------------------------------------
      function FEDISTR1(x,ane,wi,x0,intg)
c this function calculates either a density function or its integral
c for an exponential density profile.  Used in initializing particle
c coordinates.
c if intg = 0, n(x) = 1.0 + ane*exp((x-x0)*wi)
c if intg = 1, n(x) = x + (ane/wi)*(exp((x-x0)*wi) - exp(-x0*wi)
      implicit none
      integer intg
      real x, ane, x0, wi
c local data
      real FEDISTR1, f, t
      t = (x - x0)*wi
      if (intg.eq.0) then
         if (t.gt.(-64.)) then
            f = 1.0 + ane*exp(t)
         else
            f = 1.0
         endif
      else if (intg.eq.1) then
         if (wi.eq.0.) then
            f = (1.0 + ane)*x
         else
            f = x0*wi
            if (f.gt.64) then
               f = 0.
            else
               f = exp(-f)
            endif
            f = x + (ane/wi)*(exp(t) - f)
         endif
      else
         f = -1.0
      endif
      if (f.lt.0.) write (2,*) 'FEDISTR1 Error: f = ', f
      FEDISTR1 = f
      return
      end
c-----------------------------------------------------------------------
      subroutine PVDISTR32(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,
     1npz,idimp,npmax,mblok,nblok,kstrt,nvp,ierr)
c for 3d code, this subroutine calculates initial particle velocities
c with maxwellian velocity with drift for distributed data
c and 2D spatial decomposition
c particles are not necessarily in the correct processor, and the total
c number of particles should be a multiple of the number of processors.
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c npp(m) = number of particles in partition m
c nps(m) = starting address of particles in partition m
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c npx/npy/npz = initial number of particles distributed in x/y/z
c direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mblok/nblok = number of particle partitions in y/z
c kstrt = starting data block number
c nvp = number of real or virtual processors
c ierr = (0,1) = (no,yes) error condition exists
c ranorm = gaussian random number with zero mean and unit variance
      implicit none
      integer npx, npy, npz, idimp, npmax, mblok, nblok, kstrt, nvp
      integer ierr
      real vtx, vty, vtz, vdx, vdy, vdz
      integer npp, nps
      real part
      dimension part(idimp,npmax,mblok*nblok)
      dimension npp(mblok*nblok), nps(mblok*nblok)
c local data
c     integer*8 npxy, npxyz, nppv, joff, koff, i, imin
      double precision dnpxy, dnpxyz, dnpx, dnppv, djoff, dkoff, di
      double precision dimin, dt1
      integer nppv, ks, j, k, l, m, my, mz, moff, mnblok, npt, npxyzp
      integer iwork
      real vxt, vyt, vzt
      double precision ranorm
      double precision sum0, sum1, sum2, sum4, work4
      dimension sum4(4), work4(4)
      ierr = 0
c particle distribution constants
      ks = kstrt - 2
      dnpx = dble(npx)
      dnpxy = dnpx*dble(npy)
      nppv = dnpxy*dble(npz)/dble(nvp)
      nppv = min(nppv,npmax)
      dnppv = dble(nppv)
      dnpxyz = dnppv*dble(nvp)
c check for errors
      if (dnpxyz.ne.(dnpxy*dble(npz))) ierr = 1
      call PIMAX(ierr,iwork,1,1)
      if (ierr.gt.0) then
         ierr = dnpxyz - dnpxy*dble(npz)
         if (kstrt.eq.1) then
            write (2,*) 'particle distribution truncated, np = ', dnpxyz
         endif
         return
      endif
c maxwellian velocity distribution
      do 50 l = 1, npz
      dkoff = dnpxy*dble(l - 1)
      do 40 k = 1, npy
      djoff = dnpx*dble(k - 1) + dkoff
      do 30 j = 1, npx
      di = dble(j) + djoff
c maxwellian velocity distribution
      vxt = vtx*ranorm()
      vyt = vty*ranorm()
      vzt = vtz*ranorm()
      do 20 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 10 my = 1, mblok
      m = my + moff
      dimin = dnppv*dble(m + ks) + 1.0d0
      if ((di.ge.dimin).and.(di.lt.(dimin+dnppv))) then
         npt = npp(m) + 1
         if (npt.le.npmax) then
            part(4,npt,m) = vxt
            part(5,npt,m) = vyt
            part(6,npt,m) = vzt
            npp(m) = npt
         else
            ierr = ierr + 1
         endif
      endif
   10 continue
   20 continue
   30 continue
   40 continue
   50 continue
      npxyzp = 0
c add correct drift
      sum4(1) = 0.0d0
      sum4(2) = 0.0d0
      sum4(3) = 0.0d0
      mnblok = mblok*nblok
      do 70 m = 1, mnblok
      sum0 = 0.0d0
      sum1 = 0.0d0
      sum2 = 0.0d0
      do 60 j = nps(m), npp(m)
      npxyzp = npxyzp + 1
      sum0 = sum0 + part(4,j,m)
      sum1 = sum1 + part(5,j,m)
      sum2 = sum2 + part(6,j,m)
   60 continue
      sum4(1) = sum4(1) + sum0
      sum4(2) = sum4(2) + sum1
      sum4(3) = sum4(3) + sum2
   70 continue
      sum4(4) = npxyzp
      call PDSUM(sum4,work4,4,1)
      dnpxyz = sum4(4)
      call PIMAX(ierr,iwork,1,1)
      dt1 = 1.0d0/dnpxyz
      sum4(1) = dt1*sum4(1) - vdx
      sum4(2) = dt1*sum4(2) - vdy
      sum4(3) = dt1*sum4(3) - vdz
      do 90 m = 1, mnblok
      do 80 j = nps(m), npp(m)
      part(4,j,m) = part(4,j,m) - sum4(1)
      part(5,j,m) = part(5,j,m) - sum4(2)
      part(6,j,m) = part(6,j,m) - sum4(3)
   80 continue
   90 continue
c process errors
      if (ierr.gt.0) then
         if (kstrt.eq.1) then
            write (2,*) 'particle overflow error, ierr = ', ierr
         endif
      else if (dnpxyz.ne.(dnpxy*dble(npz))) then
         if (kstrt.eq.1) then
            write (2,*) 'velocity distribution truncated, np = ', dnpxyz
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PVVISTR32(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,
     1npz,idimp,npmax,mblok,nblok,vranx,vrany,vranz,kstrt,nvp,ndv,nvrp,i
     2err)
c for 3 code, this subroutine calculates initial velocities
c with maxwellian velocity with drift for distributed data
c with 2D spatial decomposition
c on input, the array npp contains the number of particles already
c stored in part, normally zero.
c on output, contains the number of particles stored in part.
c the total number of particles should be a multiple of ndv, if the
c number of processors is <= ndv, or else a multiple of the number of
c processors.
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c npp(m) = number of particles in partition m
c nps(m) = starting address of particles in partition m
c vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
c vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
c npx/npy/npz = initial number of particles distributed in x/y/z
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mblok/nblok = number of particle partitions in y/z
c vranx,vrany,vranz = output arrays for parallel gaussian random numbers
c with zero mean and unit variance
c kstrt = starting data block number
c nvp = number of real or virtual processors
c ndv = total maximum number of random seeds, currently 256
c nvrp = number of parallel seeds per processor
c ierr = (0,1) = (no,yes) error condition exists
c npx*npy*npz should be a multiple of nvrp*nvp
      implicit none
      integer npx, npy, npz, idimp, npmax, mblok, nblok, kstrt, nvp, ndv
      integer nvrp, ierr
      real vtx, vty, vtz, vdx, vdy, vdz
      double precision vranx, vrany, vranz
      integer npp, nps
      real part
      dimension part(idimp,npmax,mblok*nblok)
      dimension npp(mblok*nblok), nps(mblok*nblok)
      dimension vranx(nvrp,mblok*nblok), vrany(nvrp,mblok*nblok)
      dimension vranz(nvrp,mblok*nblok)
c local data
c     integer*8 npxy, npxyz, ipp, nppv
      double precision dnpxy, dnpxyz, dt1
      integer ipp, nppv
      integer ks, npd, mdp, mnblok, moff, ii, n, nd, m, my, mz, id, is
      integer jj, i, j, npt, npxyzp, iwork
      real vxt, vyt, vzt
      double precision sum0, sum1, sum2, sum4, work4
      dimension sum4(4), work4(4)
      ierr = 0
c particle distribution constants
      ks = kstrt - 2
      dnpxy = dble(npx)*dble(npy)
      dnpxyz = dnpxy*dble(npz)
c npd = total number of seeds used
      npd = nvrp*nvp
c ipp = number of particles per seed
      ipp = dnpxyz/dble(npd)
c nppv = number of particles per processor
      nppv = nvrp*ipp
c adjacent mdp processors share the same seed
      mdp = nvp/min0(nvp,ndv)
c check unsupported particle number
      if (dnpxyz.ne.(dble(ipp)*dble(npd))) ierr = 1
      call PIMAX(ierr,iwork,1,1)
      if (ierr.gt.0) then
         if (kstrt.eq.1) then
           write (2,*) 'number of seeds not multiple of particle number' 
     1, npd, dnpxyz
         endif
         return
      endif
      mnblok = mblok*nblok
c outer loop over processor blocks which share the same seed
      do 50 ii = 1, mdp
c particles in each block get random numbers from the same seed
      do 40 n = 1, ipp
      call pranorm(vranx,kstrt,nvp,nvrp,nd,nvrp,mnblok)
      call pranorm(vrany,kstrt,nvp,nvrp,nd,nvrp,mnblok)
      call pranorm(vranz,kstrt,nvp,nvrp,nd,nvrp,mnblok)
c inner loop over independent seeds
      do 30 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 20 my = 1, mblok
      m = my + moff
      id = (m + ks)/mdp
      is = m + ks - id*mdp + 1
      do 10 jj = 1, nvrp
c i = local particle number belonging to jj seed
      i = n + ipp*(jj - 1)
      vxt = vtx*vranx(jj,m)
      vyt = vty*vrany(jj,m)
      vzt = vtz*vranz(jj,m)
c nth group keeps nth block of random numbers
      if (ii.eq.is) then
         npt = nps(m) + i - 1
         if (npt.le.npmax) then
c maxwellian velocity distribution
            part(4,npt,m) = vxt
            part(5,npt,m) = vyt
            part(6,npt,m) = vzt
         else
            ierr = ierr + 1
         endif
      endif
   10 continue
c update particle number
      if (ii.eq.is) npp(m) = npp(m) + nvrp
   20 continue
   30 continue
   40 continue
   50 continue
      npxyzp = 0
c add correct drift
      sum4(1) = 0.0d0
      sum4(2) = 0.0d0
      sum4(3) = 0.0d0
      do 70 m = 1, mnblok
      sum0 = 0.0d0
      sum1 = 0.0d0
      sum2 = 0.0d0
      do 60 j = nps(m), npp(m)
      npxyzp = npxyzp + 1
      sum0 = sum0 + part(4,j,m)
      sum1 = sum1 + part(5,j,m)
      sum2 = sum2 + part(6,j,m)
   60 continue
      sum4(1) = sum4(1) + sum0
      sum4(2) = sum4(2) + sum1
      sum4(3) = sum4(3) + sum2
   70 continue
      sum4(4) = npxyzp
      call PDSUM(sum4,work4,4,1)
      dnpxyz = sum4(4)
      call PIMAX(ierr,iwork,1,1)
      dt1 = 1.0d0/dnpxyz
      sum4(1) = dt1*sum4(1) - vdx
      sum4(2) = dt1*sum4(2) - vdy
      sum4(3) = dt1*sum4(3) - vdz
      do 90 m = 1, mnblok
      do 80 j = nps(m), npp(m)
      part(4,j,m) = part(4,j,m) - sum4(1)
      part(5,j,m) = part(5,j,m) - sum4(2)
      part(6,j,m) = part(6,j,m) - sum4(3)
   80 continue
   90 continue
c process errors
      if (ierr.gt.0) then
         if (kstrt.eq.1) then
            write (2,*) 'particle overflow error, ierr = ', ierr
         endif
      else if (dnpxyz.ne.(dnpxy*dble(npz))) then
         ierr = dnpxyz - dnpxy*dble(npz)
         if (kstrt.eq.1) then
            write (2,*) 'particle distribution truncated, np = ', dnpxyz
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PBDISTR32L(part,bx,by,bz,npp,noff,qbm,nx,ny,nz,idimp,np
     1max,mblok,nblok,nxv,nypmx,nzpmx,idds)
c for 3d code, this subroutine reinterprets current particle positions
c as positions of guiding centers, and calculates the actual particle
c positions, with periodic boundary conditions, for distributed data
c with 2D spatial decomposition
c in converting from guiding center to actual co-ordinates, the
c following equations are used:
c       x(t) = xg(t) - vy(t)/omz
c       y(t) = yg(t) + vx(t)/omz
c       x(t) = xg(t) - (vy(t)*omz - vz(t)*omy)/om**2
c       y(t) = yg(t) - (vz(t)*omx - vx(t)*omz)/om**2
c       z(t) = zg(t) - (vx(t)*omy - vy(t)*omx)/om**2
c where omx = (q/m)*bx(xg(t),yg(t),zg(t)),
c       omy = (q/m)*by(xg(t),yg(t),zg(t)),
c and   omz = (q/m)*bz(xg(t),yg(t),zg(t)),
c bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c bx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*bx(n,m,l)+dx*bx(n+1,m,l))
c                + dy*((1-dx)*bx(n,m+1,l) + dx*bx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*bx(n,m,l+1)+dx*bx(n+1,m,l+1))
c                + dy*((1-dx)*bx(n,m+1,l+1) + dx*bx(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c similarly for by(x,y,z), bz(x,y,z)
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c bx(j,k,l,m) = x component of magnetic field at grid (j,kk,ll)
c by(j,k,l,m) = y component of magnetic field at grid (j,kk,ll)
c bz(j,k,l,m) = z component of magnetic field at grid (j,kk,ll)
c in other words, bx/by/bz is the convolution of the magnetic field
c over the particle shape,
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qbm = particle charge/mass ratio
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mblok/nblok = number of particle partitions in y/z
c nxv = first dimension of field array, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
      dimension part(idimp,npmax,mblok*nblok)
      dimension bx(nxv,nypmx,nzpmx,mblok*nblok)
      dimension by(nxv,nypmx,nzpmx,mblok*nblok)
      dimension bz(nxv,nypmx,nzpmx,mblok*nblok)
      dimension npp(mblok*nblok), noff(idds,mblok*nblok)
      zero = 0.
      anx = float(nx)
      any = float(ny)
      anz = float(nz)
      mnblok = mblok*nblok
c calculate actual position from guiding center
      do 20 m = 1, mnblok
      mnoff = noff(1,m) - 1
      lnoff = noff(2,m) - 1
      do 10 j = 1, npp(m)
c find interpolation weights
      nn = part(1,j,m)
      mm = part(2,j,m)
      ll = part(3,j,m)
      dxp = qbm*(part(1,j,m) - float(nn))
      dyp = part(2,j,m) - float(mm)
      dzp = part(3,j,m) - float(ll)
      nn = nn + 1
      amx = qbm - dxp
      np = nn + 1
      if (np.gt.nx) np = np - nx
      mm = mm - mnoff
      amy = 1. - dyp
      mp = mm + 1
      ll = ll - lnoff
      amz = 1. - dzp
      lp = ll + 1
c find magnetic field
      omx = amz*(amy*(amx*bx(nn,mm,ll,m) + dxp*bx(np,mm,ll,m)) + dyp*(am
     1x*bx(nn,mp,ll,m) + dxp*bx(np,mp,ll,m))) + dzp*(amy*(amx*bx(nn,mm,l
     2p,m) + dxp*bx(np,mm,lp,m)) + dyp*(amx*bx(nn,mp,lp,m) + dxp*bx(np,m
     3p,lp,m)))
      omy = amz*(amy*(amx*by(nn,mm,ll,m) + dxp*by(np,mm,ll,m)) + dyp*(am
     1x*by(nn,mp,ll,m) + dxp*by(np,mp,ll,m))) + dzp*(amy*(amx*by(nn,mm,l
     2p,m) + dxp*by(np,mm,lp,m)) + dyp*(amx*by(nn,mp,lp,m) + dxp*by(np,m
     3p,lp,m)))
      omz = amz*(amy*(amx*bz(nn,mm,ll,m) + dxp*bz(np,mm,ll,m)) + dyp*(am
     1x*bz(nn,mp,ll,m) + dxp*bz(np,mp,ll,m))) + dzp*(amy*(amx*bz(nn,mm,l
     2p,m) + dxp*bz(np,mm,lp,m)) + dyp*(amx*bz(nn,mp,lp,m) + dxp*bz(np,m
     3p,lp,m)))
      at3 = sqrt(omx*omx + omy*omy + omz*omz)
      if (at3.ne.0.) at3 = 1./at3
      at3 = at3*at3
      omxt = omx*at3
      omyt = omy*at3
      omzt = omz*at3
c correct position
      dx = part(1,j,m) - (part(5,j,m)*omzt - part(6,j,m)*omyt)
      dy = part(2,j,m) - (part(6,j,m)*omxt - part(4,j,m)*omzt)
      dz = part(3,j,m) - (part(4,j,m)*omyt - part(5,j,m)*omxt)
c periodic boundary conditions
      nn = abs(dx)/anx
      if (dx.lt.zero) dx = dx + float(nn + 1)*anx
      if (dx.ge.anx) dx = dx - float(nn)*anx
      part(1,j,m) = dx
      mm = abs(dy)/any
      if (dy.lt.zero) dy = dy + float(mm + 1)*any
      if (dy.ge.any) dy = dy - float(mm)*any
      part(2,j,m) = dy
      ll = abs(dz)/anz
      if (dz.lt.zero) dz = dz + float(ll + 1)*anz
      if (dz.ge.anz) dz = dz - float(ll)*anz
      part(3,j,m) = dz
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGBDISTR32L(part,bxyz,npp,noff,qbm,nx,ny,nz,idimp,npmax
     1,mblok,nblok,nxv,nypmx,nzpmx,idds,ipbc)
c for 3d code, this subroutine reinterprets current particle positions
c as positions of guiding centers, and calculates the actual particle
c positions, for distributed data with 2D spatial decomposition
c in converting from guiding center to actual co-ordinates, the
c following equations are used:
c       x(t) = xg(t) - (vy(t)*omz - vz(t)*omy)/om**2
c       y(t) = yg(t) - (vz(t)*omx - vx(t)*omz)/om**2
c       z(t) = zg(t) - (vx(t)*omy - vy(t)*omx)/om**2
c where omx = (q/m)*bx(xg(t),yg(t),zg(t)),
c       omy = (q/m)*by(xg(t),yg(t),zg(t)),
c       omz = (q/m)*bz(xg(t),yg(t),zg(t)),
c and the magnetic field components bxyz(i,x(t),y(t),z(t)) are 
c are approximated by interpolation from the nearest grid points:
c bxyz(i,x,y,z)=(1-dz)*((1-dy)*((1-dx)*bxyz(i,n,m,l)+dx*bxyz(i,n+1,m,l))
c            + dy*((1-dx)*bxyz(i,n,m+1,l) + dx*bxyz(i,n+1,m+1,l)))
c            + dz*((1-dy)*((1-dx)*bxyz(i,n,m,l+1)+dx*bxyz(i,n+1,m,l+1))
c            + dy*((1-dx)*bxyz(i,n,m+1,l+1) + dx*bxyz(i,n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c similarly for by(x,y,z), bz(x,y,z)
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c bxyz(i,j,k,l,m) = i component of magnetic field at grid (j,kk,ll)
c in other words, bxyz is the convolution of the magnetic field
c over the particle shape,
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qbm = particle charge/mass ratio
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mblok/nblok = number of particle partitions in y/z
c nxv = first dimension of field array, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      dimension part(idimp,npmax,mblok*nblok)
      dimension bxyz(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension npp(mblok*nblok), noff(idds,mblok*nblok)
      mnblok = mblok*nblok
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgely = 0.
         edgelz = 0.
         edgerx = float(nx)
         edgery = float(ny)
         edgerz = float(nz)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgelz = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
         edgerz = float(nz-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgely = 1.
         edgelz = 0.
         edgerx = float(nx-1)
         edgery = float(ny-1)
         edgerz = float(nz)
      endif
c calculate actual position from guiding center
      do 30 m = 1, mnblok
      mnoff = noff(1,m) - 1
      lnoff = noff(2,m) - 1
      do 20 j = 1, npp(m)
c find interpolation weights
      nn = part(1,j,m)
      mm = part(2,j,m)
      ll = part(3,j,m)
      dxp = qbm*(part(1,j,m) - float(nn))
      dyp = part(2,j,m) - float(mm)
      dzp = part(3,j,m) - float(ll)
      nn = nn + 1
      amx = qbm - dxp
      amy = 1. - dyp
      np = nn + 1
      mm = mm - mnoff
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + 1
      amx = amx*amy
      amz = 1. - dzp
      ll = ll - lnoff
      amy = dxp*amy
      lp = ll + 1
c find magnetic field
      omx = amz*(amx*bxyz(1,nn,mm,ll,m) + amy*bxyz(1,np,mm,ll,m) + dyp*b
     1xyz(1,nn,mp,ll,m) + dx1*bxyz(1,np,mp,ll,m)) + dzp*(amx*bxyz(1,nn,m
     2m,lp,m) + amy*bxyz(1,np,mm,lp,m) + dyp*bxyz(1,nn,mp,lp,m) + dx1*bx
     3yz(1,np,mp,lp,m))
      omy = amz*(amx*bxyz(2,nn,mm,ll,m) + amy*bxyz(2,np,mm,ll,m) + dyp*b
     1xyz(2,nn,mp,ll,m) + dx1*bxyz(2,np,mp,ll,m)) + dzp*(amx*bxyz(2,nn,m
     2m,lp,m) + amy*bxyz(2,np,mm,lp,m) + dyp*bxyz(2,nn,mp,lp,m) + dx1*bx
     3yz(2,np,mp,lp,m))
      omz = amz*(amx*bxyz(3,nn,mm,ll,m) + amy*bxyz(3,np,mm,ll,m) + dyp*b
     1xyz(3,nn,mp,ll,m) + dx1*bxyz(3,np,mp,ll,m)) + dzp*(amx*bxyz(3,nn,m
     2m,lp,m) + amy*bxyz(3,np,mm,lp,m) + dyp*bxyz(3,nn,mp,lp,m) + dx1*bx
     3yz(3,np,mp,lp,m))
      at3 = sqrt(omx*omx + omy*omy + omz*omz)
      if (at3.ne.0.) at3 = 1./at3
      at3 = at3*at3
      omxt = omx*at3
      omyt = omy*at3
      omzt = omz*at3
c correct position
      dx = part(1,j,m) - (part(5,j,m)*omzt - part(6,j,m)*omyt)
      dy = part(2,j,m) - (part(6,j,m)*omxt - part(4,j,m)*omzt)
      dz = part(3,j,m) - (part(4,j,m)*omyt - part(5,j,m)*omxt)
c periodic boundary conditions
      if (ipbc.eq.1) then
         n = abs(dx)/edgerx
         if (dx.lt.edgelx) dx = dx + float(n + 1)*edgerx
         if (dx.ge.edgerx) dx = dx - float(n)*edgerx
         k = abs(dy)/edgery
         if (dy.lt.edgely) dy = dy + float(k + 1)*edgery
         if (dy.ge.edgery) dy = dy - float(k)*edgery
         l = abs(dz)/edgerz
         if (dz.lt.edgelz) dz = dz + float(l + 1)*edgerz
         if (dz.ge.edgerz) dz = dz - float(l)*edgerz
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx).or.(dy.lt.edgely).or.(dy.g
     1e.edgery).or.(dz.lt.edgelz).or.(dz.ge.edgerz)) then
c try switching vx
            dy = part(2,j,m) - (part(6,j,m)*omxt + part(4,j,m)*omzt)
            dz = part(3,j,m) + (part(4,j,m)*omyt + part(5,j,m)*omxt)
            if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgely).and
     1.(dy.lt.edgery).and.(dz.ge.edgelz).and.(dz.lt.edgerz)) then
               part(4,j,m) = -part(4,j,m)
               go to 10
            endif
c try switching vy
            dx = part(1,j,m) + (part(5,j,m)*omzt + part(6,j,m)*omyt)
            dz = part(3,j,m) - (part(4,j,m)*omyt + part(5,j,m)*omxt)
            if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgely).and
     1.(dy.lt.edgery).and.(dz.ge.edgelz).and.(dz.lt.edgerz)) then
               part(5,j,m) = -part(5,j,m)
               go to 10
            endif
c try switching vz
            dx = part(1,j,m) - (part(5,j,m)*omzt + part(6,j,m)*omyt)
            dy = part(2,j,m) + (part(6,j,m)*omxt + part(4,j,m)*omzt)
            if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgely).and
     1.(dy.lt.edgery).and.(dz.ge.edgelz).and.(dz.lt.edgerz)) then
               part(6,j,m) = -part(6,j,m)
               go to 10
            endif
c try switching both vx and vy
            dx = part(1,j,m) + (part(5,j,m)*omzt + part(6,j,m)*omyt)
            dy = part(2,j,m) - (part(6,j,m)*omxt + part(4,j,m)*omzt)
            dz = part(3,j,m) + (part(4,j,m)*omyt - part(5,j,m)*omxt)
            if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgely).and
     1.(dy.lt.edgery).and.(dz.ge.edgelz).and.(dz.lt.edgerz)) then
               part(4,j,m) = -part(4,j,m)
               part(5,j,m) = -part(5,j,m)
               go to 10
            endif
c try switching both vx and vz
            dx = part(1,j,m) - (part(5,j,m)*omzt + part(6,j,m)*omyt)
            dy = part(2,j,m) + (part(6,j,m)*omxt - part(4,j,m)*omzt)
            dz = part(3,j,m) + (part(4,j,m)*omyt + part(5,j,m)*omxt)
            if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgely).and
     1.(dy.lt.edgery).and.(dz.ge.edgelz).and.(dz.lt.edgerz)) then
               part(4,j,m) = -part(4,j,m)
               part(6,j,m) = -part(6,j,m)
               go to 10
            endif
c try switching both vy and vz
            dx = part(1,j,m) + (part(5,j,m)*omzt - part(6,j,m)*omyt)
            dy = part(2,j,m) + (part(6,j,m)*omxt + part(4,j,m)*omzt)
            dz = part(3,j,m) - (part(4,j,m)*omyt + part(5,j,m)*omxt)
            if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgely).and
     1.(dy.lt.edgery).and.(dz.ge.edgelz).and.(dz.lt.edgerz)) then
               part(5,j,m) = -part(5,j,m)
               part(6,j,m) = -part(6,j,m)
               go to 10
            endif
c try switching vx, vy, vz
            dx = part(1,j,m) + (part(5,j,m)*omzt - part(6,j,m)*omyt)
            dy = part(2,j,m) + (part(6,j,m)*omxt - part(4,j,m)*omzt)
            dz = part(3,j,m) + (part(4,j,m)*omyt - part(5,j,m)*omxt)
            if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgely).and
     1.(dy.lt.edgery).and.(dz.ge.edgelz).and.(dz.lt.edgerz)) then
               part(4,j,m) = -part(4,j,m)
               part(5,j,m) = -part(5,j,m)
               part(6,j,m) = -part(6,j,m)
               go to 10
            endif
c give up if larmor radius is too large
            dx = part(1,j,m)
            dy = part(2,j,m)
            dz = part(3,j,m)
   10       continue
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx).or.(dy.lt.edgely).or.(dy.g
     1e.edgery)) then
            if ((dy.ge.edgely).and.(dy.lt.edgery)) then
c if x co-ordinate only is out of bounds, try switching vy
               dx = part(1,j,m) + (part(5,j,m)*omzt + part(6,j,m)*omyt)
               dz = part(3,j,m) - (part(4,j,m)*omyt + part(5,j,m)*omxt)
               if ((dx.ge.edgelx).and.(dx.lt.edgerx)) then
                  part(5,j,m) = -part(5,j,m)
               else
c then try switching vz
                  dx = part(1,j,m) - (part(5,j,m)*omzt + part(6,j,m)*omy
     1t)
                  dy = part(2,j,m) + (part(6,j,m)*omxt + part(4,j,m)*omz
     1t)
                  if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgel
     1y).and.(dy.lt.edgery)) then
                     part(6,j,m) = -part(6,j,m)
                  else
c otherwise, try switching both vy and vz
                     dx = part(1,j,m) + (part(5,j,m)*omzt - part(6,j,m)*
     1omyt)
                     dy = part(2,j,m) + (part(6,j,m)*omxt + part(4,j,m)*
     1omzt)
                     dz = part(3,j,m) - (part(4,j,m)*omyt + part(5,j,m)*
     1omxt)
                     if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.ed
     1gely).and.(dy.lt.edgery)) then
                        part(5,j,m) = -part(5,j,m)
                        part(6,j,m) = -part(6,j,m)
                     endif
                  endif
               endif
            else if ((dx.ge.edgelx).and.(dx.lt.edgerx)) then
c if y co-ordinate only is out of bounds, try switching vx
               dy = part(2,j,m) - (part(6,j,m)*omxt + part(4,j,m)*omzt)
               dz = part(3,j,m) + (part(4,j,m)*omyt + part(5,j,m)*omxt)
               if ((dy.ge.edgely).and.(dy.lt.edgery)) then
                  part(4,j,m) = -part(4,j,m)
               else
c then try switching vz
                  dx = part(1,j,m) - (part(5,j,m)*omzt + part(6,j,m)*omy
     1t)
                  dy = part(2,j,m) + (part(6,j,m)*omxt + part(4,j,m)*omz
     1t)
                  if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgel
     1y).and.(dy.lt.edgery)) then
                     part(6,j,m) = -part(6,j,m)
                  else
c otherwise, try switching both vx and vz
                     dx = part(1,j,m) - (part(5,j,m)*omzt + part(6,j,m)*
     1omyt)
                     dy = part(2,j,m) + (part(6,j,m)*omxt - part(4,j,m)*
     1omzt)
                     dz = part(3,j,m) + (part(4,j,m)*omyt + part(5,j,m)*
     1omxt)
                     if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.ed
     1gely).and.(dy.lt.edgery)) then
                        part(4,j,m) = -part(4,j,m)
                        part(6,j,m) = -part(6,j,m)
                     endif
                  endif
               endif
            endif
c if both co-ordinates are out of bounds, try switching vx, vy, vz
            if ((dx.lt.edgelx).or.(dx.ge.edgerx).or.(dy.lt.edgely).or.(d
     1y.ge.edgery)) then
               dx = part(1,j,m) + (part(5,j,m)*omzt - part(6,j,m)*omyt)
               dy = part(2,j,m) + (part(6,j,m)*omxt - part(4,j,m)*omzt)
               dz = part(3,j,m) + (part(4,j,m)*omyt - part(5,j,m)*omxt)
               if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgely).
     1and.(dy.lt.edgery)) then
                  part(4,j,m) = -part(4,j,m)
                  part(5,j,m) = -part(5,j,m)
                  part(6,j,m) = -part(6,j,m)
               else
c give up if larmor radius is too large
                  dx = part(1,j,m)
                  dy = part(2,j,m)
                  dz = part(3,j,m)
               endif
            endif
         endif
         l = abs(dz)/edgerz
         if (dz.lt.edgelz) dz = dz + float(l + 1)*edgerz
         if (dz.ge.edgerz) dz = dz - float(l)*edgerz
      endif
c set new position
      part(1,j,m) = dx
      part(2,j,m) = dy
      part(3,j,m) = dz
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGRBDISTR32L(part,bxyz,npp,noff,qbm,ci,nx,ny,nz,idimp,n
     1pmax,mblok,nblok,nxv,nypmx,nzpmx,idds,ipbc)
c for 3d code, this subroutine reinterprets current particle positions
c as positions of guiding centers, and calculates the actual particle
c positions, for relativistic particles, for distributed data
c with 2D spatial decomposition
c in converting from guiding center to actual co-ordinates, the
c following equations are used:
c       x(t) = xg(t) - gami*(py(t)*omz - pz(t)*omy)/om**2
c       y(t) = yg(t) - gami*(pz(t)*omx - px(t)*omz)/om**2
c       z(t) = zg(t) - gami*(px(t)*omy - py(t)*omx)/om**2
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c and omx = (q/m)*bxyz(1,xg(t),yg(t),zg(t)),
c     omy = (q/m)*bxyz(2,xg(t),yg(t),zg(t)),
c     omz = (q/m)*bxyz(3,xg(t),yg(t),zg(t)),
c and the magnetic field components bxyz(i,x(t),y(t),z(t)) are 
c are approximated by interpolation from the nearest grid points:
c bxyz(i,x,y,z)=(1-dz)*((1-dy)*((1-dx)*bxyz(i,n,m,l)+dx*bxyz(i,n+1,m,l))
c            + dy*((1-dx)*bxyz(i,n,m+1,l) + dx*bxyz(i,n+1,m+1,l)))
c            + dz*((1-dy)*((1-dx)*bxyz(i,n,m,l+1)+dx*bxyz(i,n+1,m,l+1))
c            + dy*((1-dx)*bxyz(i,n,m+1,l+1) + dx*bxyz(i,n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c similarly for by(x,y,z), bz(x,y,z)
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = momentum px of particle n in partition m
c part(5,n,m) = momentum py of particle n in partition m
c part(6,n,m) = momentum pz of particle n in partition m
c bxyz(i,j,k,l,m) = i component of magnetic field at grid (j,kk,ll)
c in other words, bxyz is the convolution of the magnetic field
c over the particle shape,
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qbm = particle charge/mass ratio
c ci = reciprical of velocity of light
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mblok/nblok = number of particle partitions in y/z
c nxv = first dimension of field array, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      dimension part(idimp,npmax,mblok*nblok)
      dimension bxyz(3,nxv,nypmx,nzpmx,mblok*nblok)
      dimension npp(mblok*nblok), noff(idds,mblok*nblok)
      ci2 = ci*ci
      mnblok = mblok*nblok
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgely = 0.
         edgelz = 0.
         edgerx = float(nx)
         edgery = float(ny)
         edgerz = float(nz)
      else if (ipbc.eq.2) then
         edgelx = 1.
         edgely = 1.
         edgelz = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
         edgerz = float(nz-1)
      else if (ipbc.eq.3) then
         edgelx = 1.
         edgely = 1.
         edgelz = 0.
         edgerx = float(nx-1)
         edgery = float(ny-1)
         edgerz = float(nz)
      endif
c calculate actual position from guiding center
      do 30 m = 1, mnblok
      mnoff = noff(1,m) - 1
      lnoff = noff(2,m) - 1
      do 20 j = 1, npp(m)
c find interpolation weights
      nn = part(1,j,m)
      mm = part(2,j,m)
      ll = part(3,j,m)
      dxp = qbm*(part(1,j,m) - float(nn))
      dyp = part(2,j,m) - float(mm)
      dzp = part(3,j,m) - float(ll)
      nn = nn + 1
      amx = qbm - dxp
      amy = 1. - dyp
      np = nn + 1
      mm = mm - mnoff
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + 1
      amx = amx*amy
      amz = 1. - dzp
      ll = ll - lnoff
      amy = dxp*amy
      lp = ll + 1
c find inverse gamma
      vx = part(4,j,m)
      vy = part(5,j,m)
      vz = part(6,j,m)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c find magnetic field
      omx = amz*(amx*bxyz(1,nn,mm,ll,m) + amy*bxyz(1,np,mm,ll,m) + dyp*b
     1xyz(1,nn,mp,ll,m) + dx1*bxyz(1,np,mp,ll,m)) + dzp*(amx*bxyz(1,nn,m
     2m,lp,m) + amy*bxyz(1,np,mm,lp,m) + dyp*bxyz(1,nn,mp,lp,m) + dx1*bx
     3yz(1,np,mp,lp,m))
      omy = amz*(amx*bxyz(2,nn,mm,ll,m) + amy*bxyz(2,np,mm,ll,m) + dyp*b
     1xyz(2,nn,mp,ll,m) + dx1*bxyz(2,np,mp,ll,m)) + dzp*(amx*bxyz(2,nn,m
     2m,lp,m) + amy*bxyz(2,np,mm,lp,m) + dyp*bxyz(2,nn,mp,lp,m) + dx1*bx
     3yz(2,np,mp,lp,m))
      omz = amz*(amx*bxyz(3,nn,mm,ll,m) + amy*bxyz(3,np,mm,ll,m) + dyp*b
     1xyz(3,nn,mp,ll,m) + dx1*bxyz(3,np,mp,ll,m)) + dzp*(amx*bxyz(3,nn,m
     2m,lp,m) + amy*bxyz(3,np,mm,lp,m) + dyp*bxyz(3,nn,mp,lp,m) + dx1*bx
     3yz(3,np,mp,lp,m))
      at3 = sqrt(omx*omx + omy*omy + omz*omz)
      if (at3.ne.0.) at3 = 1./at3
      at3 = at3*at3*gami
      omxt = omx*at3
      omyt = omy*at3
      omzt = omz*at3
c correct position
      dx = part(1,j,m) - (part(5,j,m)*omzt - part(6,j,m)*omyt)
      dy = part(2,j,m) - (part(6,j,m)*omxt - part(4,j,m)*omzt)
      dz = part(3,j,m) - (part(4,j,m)*omyt - part(5,j,m)*omxt)
c periodic boundary conditions
      if (ipbc.eq.1) then
         n = abs(dx)/edgerx
         if (dx.lt.edgelx) dx = dx + float(n + 1)*edgerx
         if (dx.ge.edgerx) dx = dx - float(n)*edgerx
         k = abs(dy)/edgery
         if (dy.lt.edgely) dy = dy + float(k + 1)*edgery
         if (dy.ge.edgery) dy = dy - float(k)*edgery
         l = abs(dz)/edgerz
         if (dz.lt.edgelz) dz = dz + float(l + 1)*edgerz
         if (dz.ge.edgerz) dz = dz - float(l)*edgerz
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx).or.(dy.lt.edgely).or.(dy.g
     1e.edgery).or.(dz.lt.edgelz).or.(dz.ge.edgerz)) then
c try switching vx
            dy = part(2,j,m) - (part(6,j,m)*omxt + part(4,j,m)*omzt)
            dz = part(3,j,m) + (part(4,j,m)*omyt + part(5,j,m)*omxt)
            if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgely).and
     1.(dy.lt.edgery).and.(dz.ge.edgelz).and.(dz.lt.edgerz)) then
               part(4,j,m) = -part(4,j,m)
               go to 10
            endif
c try switching vy
            dx = part(1,j,m) + (part(5,j,m)*omzt + part(6,j,m)*omyt)
            dz = part(3,j,m) - (part(4,j,m)*omyt + part(5,j,m)*omxt)
            if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgely).and
     1.(dy.lt.edgery).and.(dz.ge.edgelz).and.(dz.lt.edgerz)) then
               part(5,j,m) = -part(5,j,m)
               go to 10
            endif
c try switching vz
            dx = part(1,j,m) - (part(5,j,m)*omzt + part(6,j,m)*omyt)
            dy = part(2,j,m) + (part(6,j,m)*omxt + part(4,j,m)*omzt)
            if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgely).and
     1.(dy.lt.edgery).and.(dz.ge.edgelz).and.(dz.lt.edgerz)) then
               part(6,j,m) = -part(6,j,m)
               go to 10
            endif
c try switching both vx and vy
            dx = part(1,j,m) + (part(5,j,m)*omzt + part(6,j,m)*omyt)
            dy = part(2,j,m) - (part(6,j,m)*omxt + part(4,j,m)*omzt)
            dz = part(3,j,m) + (part(4,j,m)*omyt - part(5,j,m)*omxt)
            if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgely).and
     1.(dy.lt.edgery).and.(dz.ge.edgelz).and.(dz.lt.edgerz)) then
               part(4,j,m) = -part(4,j,m)
               part(5,j,m) = -part(5,j,m)
               go to 10
            endif
c try switching both vx and vz
            dx = part(1,j,m) - (part(5,j,m)*omzt + part(6,j,m)*omyt)
            dy = part(2,j,m) + (part(6,j,m)*omxt - part(4,j,m)*omzt)
            dz = part(3,j,m) + (part(4,j,m)*omyt + part(5,j,m)*omxt)
            if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgely).and
     1.(dy.lt.edgery).and.(dz.ge.edgelz).and.(dz.lt.edgerz)) then
               part(4,j,m) = -part(4,j,m)
               part(6,j,m) = -part(6,j,m)
               go to 10
            endif
c try switching both vy and vz
            dx = part(1,j,m) + (part(5,j,m)*omzt - part(6,j,m)*omyt)
            dy = part(2,j,m) + (part(6,j,m)*omxt + part(4,j,m)*omzt)
            dz = part(3,j,m) - (part(4,j,m)*omyt + part(5,j,m)*omxt)
            if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgely).and
     1.(dy.lt.edgery).and.(dz.ge.edgelz).and.(dz.lt.edgerz)) then
               part(5,j,m) = -part(5,j,m)
               part(6,j,m) = -part(6,j,m)
               go to 10
            endif
c try switching vx, vy, vz
            dx = part(1,j,m) + (part(5,j,m)*omzt - part(6,j,m)*omyt)
            dy = part(2,j,m) + (part(6,j,m)*omxt - part(4,j,m)*omzt)
            dz = part(3,j,m) + (part(4,j,m)*omyt - part(5,j,m)*omxt)
            if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgely).and
     1.(dy.lt.edgery).and.(dz.ge.edgelz).and.(dz.lt.edgerz)) then
               part(4,j,m) = -part(4,j,m)
               part(5,j,m) = -part(5,j,m)
               part(6,j,m) = -part(6,j,m)
               go to 10
            endif
c give up if larmor radius is too large
            dx = part(1,j,m)
            dy = part(2,j,m)
            dz = part(3,j,m)
   10       continue
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx).or.(dy.lt.edgely).or.(dy.g
     1e.edgery)) then
            if ((dy.ge.edgely).and.(dy.lt.edgery)) then
c if x co-ordinate only is out of bounds, try switching vy
               dx = part(1,j,m) + (part(5,j,m)*omzt + part(6,j,m)*omyt)
               dz = part(3,j,m) - (part(4,j,m)*omyt + part(5,j,m)*omxt)
               if ((dx.ge.edgelx).and.(dx.lt.edgerx)) then
                  part(5,j,m) = -part(5,j,m)
               else
c then try switching vz
                  dx = part(1,j,m) - (part(5,j,m)*omzt + part(6,j,m)*omy
     1t)
                  dy = part(2,j,m) + (part(6,j,m)*omxt + part(4,j,m)*omz
     1t)
                  if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgel
     1y).and.(dy.lt.edgery)) then
                     part(6,j,m) = -part(6,j,m)
                  else
c otherwise, try switching both vy and vz
                     dx = part(1,j,m) + (part(5,j,m)*omzt - part(6,j,m)*
     1omyt)
                     dy = part(2,j,m) + (part(6,j,m)*omxt + part(4,j,m)*
     1omzt)
                     dz = part(3,j,m) - (part(4,j,m)*omyt + part(5,j,m)*
     1omxt)
                     if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.ed
     1gely).and.(dy.lt.edgery)) then
                        part(5,j,m) = -part(5,j,m)
                        part(6,j,m) = -part(6,j,m)
                     endif
                  endif
               endif
            else if ((dx.ge.edgelx).and.(dx.lt.edgerx)) then
c if y co-ordinate only is out of bounds, try switching vx
               dy = part(2,j,m) - (part(6,j,m)*omxt + part(4,j,m)*omzt)
               dz = part(3,j,m) + (part(4,j,m)*omyt + part(5,j,m)*omxt)
               if ((dy.ge.edgely).and.(dy.lt.edgery)) then
                  part(4,j,m) = -part(4,j,m)
               else
c then try switching vz
                  dx = part(1,j,m) - (part(5,j,m)*omzt + part(6,j,m)*omy
     1t)
                  dy = part(2,j,m) + (part(6,j,m)*omxt + part(4,j,m)*omz
     1t)
                  if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgel
     1y).and.(dy.lt.edgery)) then
                     part(6,j,m) = -part(6,j,m)
                  else
c otherwise, try switching both vx and vz
                     dx = part(1,j,m) - (part(5,j,m)*omzt + part(6,j,m)*
     1omyt)
                     dy = part(2,j,m) + (part(6,j,m)*omxt - part(4,j,m)*
     1omzt)
                     dz = part(3,j,m) + (part(4,j,m)*omyt + part(5,j,m)*
     1omxt)
                     if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.ed
     1gely).and.(dy.lt.edgery)) then
                        part(4,j,m) = -part(4,j,m)
                        part(6,j,m) = -part(6,j,m)
                     endif
                  endif
               endif
            endif
c if both co-ordinates are out of bounds, try switching vx, vy, vz
            if ((dx.lt.edgelx).or.(dx.ge.edgerx).or.(dy.lt.edgely).or.(d
     1y.ge.edgery)) then
               dx = part(1,j,m) + (part(5,j,m)*omzt - part(6,j,m)*omyt)
               dy = part(2,j,m) + (part(6,j,m)*omxt - part(4,j,m)*omzt)
               dz = part(3,j,m) + (part(4,j,m)*omyt - part(5,j,m)*omxt)
               if ((dx.ge.edgelx).and.(dx.lt.edgerx).and.(dy.ge.edgely).
     1and.(dy.lt.edgery)) then
                  part(4,j,m) = -part(4,j,m)
                  part(5,j,m) = -part(5,j,m)
                  part(6,j,m) = -part(6,j,m)
               else
c give up if larmor radius is too large
                  dx = part(1,j,m)
                  dy = part(2,j,m)
                  dz = part(3,j,m)
               endif
            endif
         endif
         l = abs(dz)/edgerz
         if (dz.lt.edgelz) dz = dz + float(l + 1)*edgerz
         if (dz.ge.edgerz) dz = dz - float(l)*edgerz
      endif
c set new position
      part(1,j,m) = dx
      part(2,j,m) = dy
      part(3,j,m) = dz
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine FEDGES32(edges,noff,nyzp,fny,argy1,argy2,argy3,fnz,argz
     11,argz2,argz3,ny,nz,nypmin,nypmax,nzpmin,nzpmax,kstrt,nvpy,nvpz,mb
     2lok,nblok,idps,idds,ipbc)
c this subroutines finds new partitions boundaries (edges,noff,nyzp)
c from density integrals given by = fny(z,argy1,argy2,argy3,1) and
c fnz(z,argz1,argz2,argz3,1)
c edges(1,m) = lower boundary in y of particle partition m
c edges(2,m) = upper boundary in y of particle partition m
c edges(3,m) = back boundary in z of particle partition m
c edges(4,m) = front boundary in z of particle partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c fny/fnz = density and density integral function in y/z
c argy1,argy2,argy3 = arguments to fny
c argz1,argz2,argz3 = arguments to fnz
c nypmin/nypmax = minimum/maximum value of nyzp(1,m) in new partition
c nzpmin/nzpmax = minimum/maximum value of nyzp(2,m) in new partition
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c mblok/nblok = number of field partitions in y/z.
c idps = number of partition boundaries
c idds = dimensionality of domain decomposition
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      real edges
      real argy1, argy2, argy3, argz1, argz2, argz3
      integer noff, nyzp
      integer ny, nz, nypmin, nypmax, nzpmin, nzpmax, kstrt, nvpy, nvpz
      integer mblok, nblok, idps, idds, ipbc
      dimension edges(idps,nblok)
      dimension noff(idds,mblok*nblok), nyzp(idds,mblok*nblok)
      real fny, fnz
      external fny, fnz
c local data
      integer js, ks, kl, kr, my, mz, m, moff, mnblok
      real edgely, edgelz, any1, any, anz1, anz, y0, y1, z0, z1
      real anpav, anpl, anpr, sum1, at1, at2
      integer ibflg, iwork
      dimension ibflg(4), iwork(4)
c particle distribution constants
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      mnblok = mblok*nblok
c set boundary values
      edgely = 0.
      edgelz = 0.
      if (ipbc.eq.2) then
         edgely = 1.
         edgelz = 1.
      else if (ipbc.eq.3) then
         edgely = 1.
         edgelz = 0.
      endif
c find normalization for function
      any = real(ny)
      anz = real(nz)
      any1 = any - edgely
      anz1 = anz - edgelz
      y0 = fny(edgely,argy1,argy2,argy3,1)
      z0 = fnz(edgelz,argz1,argz2,argz3,1)
c search for boundaries
      do 60 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 50 my = 1, mblok
      m = my + moff
c anpav = desired number of particles per processor in y
      anpav = (fny(any1,argy1,argy2,argy3,1) - y0)/real(nvpy)
c search for boundaries in y
      kl = my + js
      anpl = real(kl)*anpav
      anpr = real(kl+1)*anpav
      y1 = edgely
      sum1 = 0.
c first find left boundary
   10 at1 = sum1
      sum1 = fny(y1,argy1,argy2,argy3,1) - y0
      y1 = y1 + 1.0
      if ((sum1.lt.anpl).and.(y1.le.any)) go to 10 
      if (sum1.gt.at1) then
         at2 = (y1 - 2.0) + (anpl - at1)/(sum1 - at1)
      else
         at2 = y1 - 1.0
      endif
      edges(1,m) = at2
c set leftmost edge to zero
      if (kl.eq.0) edges(1,m) = 0.
c then find right boundary
   20 at1 = sum1
      sum1 = fny(y1,argy1,argy2,argy3,1) - y0
      y1 = y1 + 1.0
      if ((sum1.lt.anpr).and.(y1.le.any)) go to 20
      at2 = (y1 - 2.0) + (anpr - at1)/(sum1 - at1)
      edges(2,m) = at2
c set rightmost edge to ny
      if ((kl+1).eq.nvpy) edges(2,m) = any
c anpav = desired number of particles per processor in z
      anpav = (fnz(anz1,argz1,argz2,argz3,1) - z0)/real(nvpz)
c search for boundaries in z
      kl = mz + ks
      anpl = real(kl)*anpav
      anpr = real(kl+1)*anpav
      z1 = edgelz
      sum1 = 0.
c first find left boundary
   30 at1 = sum1
      sum1 = fnz(z1,argz1,argz2,argz3,1) - z0
      z1 = z1 + 1.0
      if ((sum1.lt.anpl).and.(z1.le.anz)) go to 30 
      if (sum1.gt.at1) then
         at2 = (z1 - 2.0) + (anpl - at1)/(sum1 - at1)
      else
         at2 = z1 - 1.0
      endif
      edges(3,m) = at2
c set leftmost edge to zero
      if (kl.eq.0) edges(3,m) = 0.
c then find right boundary
   40 at1 = sum1
      sum1 = fnz(z1,argz1,argz2,argz3,1) - z0
      z1 = z1 + 1.0
      if ((sum1.lt.anpr).and.(z1.le.anz)) go to 40
      at2 = (z1 - 2.0) + (anpr - at1)/(sum1 - at1)
      edges(4,m) = at2
c set rightmost edge to nz
      if ((kl+1).eq.nvpz) edges(4,m) = anz
   50 continue
   60 continue
c calculate number of grids and offsets in new partitions
      do 70 m = 1, mnblok
      kl = edges(1,m) + .5
      noff(1,m) = kl
      kr = edges(2,m) + .5
      nyzp(1,m) = kr - kl
      edges(1,m) = real(kl)
      edges(2,m) = real(kr)
      kl = edges(3,m) + .5
      noff(2,m) = kl
      kr = edges(4,m) + .5
      nyzp(2,m) = kr - kl
      edges(3,m) = real(kl)
      edges(4,m) = real(kr)
   70 continue
c find minimum and maximum partition size
      nypmin = nyzp(1,1)
      nypmax = nyzp(1,1)
      nzpmin = nyzp(2,1)
      nzpmax = nyzp(2,1)
      do 80 m = 1, mnblok
      nypmin = min0(nypmin,nyzp(1,m))
      nypmax = max0(nypmax,nyzp(1,m))
      nzpmin = min0(nzpmin,nyzp(2,m))
      nzpmax = max0(nzpmax,nyzp(2,m))
   80 continue
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
      function randum()
c this is a version of the random number generator dprandom due to
c c. bingham and the yale computer center, producing numbers
c in the interval (0,1).  written for the sun by viktor k. decyk, ucla
      integer r1,r2
      double precision randum,h1l,h1u,r0,r3,asc,bsc
      save r1,r2,h1l,h1u
      data r1,r2 /1271199957,1013501921/
      data h1l,h1u /65533.0d0,32767.0d0/
      isc = 65536
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
      randum = (dble(r1) + dble(r2)*asc)*asc
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
c-----------------------------------------------------------------------
      subroutine prandom(vran,kstrt,nvp,nvrp,ndp,nvrd,nblok)
c this program calculates nvrp random numbers for nvp processors from a 
c uniform distribution, producing numbers in the interval (0,1)
c written for the ibm by viktor k. decyk, ucla
c parallel version
c each random number is generated from a separate seed, which is just
c the normal seed calculated every 100 million numbers apart.
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
      integer r1,r2
      integer r1a,r1b,r1c,r1d,r2a,r2b,r2c,r2d
      double precision vran,h1l,h1u,r0,r3,asc,bsc
      dimension vran(nvrd,nblok)
      dimension r1(ndv), r2(ndv)
      dimension r1a(ndvb), r1b(ndvb), r1c(ndvb), r1d(ndvb)
      dimension r2a(ndvb), r2b(ndvb), r2c(ndvb), r2d(ndvb)
      equivalence (r1a(1),r1(1)), (r1b(1),r1(ndvb+1))
      equivalence (r1c(1),r1(2*ndvb+1)), (r1d(1),r1(3*ndvb+1))
      equivalence (r2a(1),r2(1)), (r2b(1),r2(ndvb+1))
      equivalence (r2c(1),r2(2*ndvb+1)), (r2d(1),r2(3*ndvb+1))
      save h1l,h1u,r1,r2
      data h1l,h1u /65533.0d0,32767.0d0/
      data r1a /1614663166,1129329378,526800021,776978983,2108255405,139
     17073720,831792859,2081030824,998880178,46863628,1925407687,1594703
     2607,307275694,850067710,1175916538,1456798645,936529216,1762783151
     3,471085076,1480506753,1061464328,1905157566,312592563,1837762811,9
     460479143,2098259915,352902777,1448633027,1404933564,1105984119,207
     59618565,665669498,2048536551,466887265,762418167,1810781119,196421
     67033,736093225,1022311664,1122169906,180994817,2022163567,18864270
     73,28620148,515818799,229001172,1720861301,1391402020,1373019636,47
     87623031,713009217,1451207906,576147951,463161977,1610460051,110244
     91551,199812544,1729155833,1566920971,2135121674,664100360,16144837
     a35,1718721866,1329301939/
      data r1b /1271447397,232687987,1279235311,1729360107,743509626,100
     18844974,223231130,295967057,2053001188,479017831,1971230956,132785
     29589,75663429,1712029246,1041477252,1966325799,1454488636,19614304
     384,766673410,1322010584,1514270826,1792306888,945952775,1011786552
     4,1978345837,1323181498,529233202,1400657637,768044007,762885963,30
     51571619,2138630017,1129495672,947431450,1272551290,550099371,14330
     683839,1729395911,918185432,760656388,5282396,367311347,160239069,1
     7349055786,2107977967,652469309,1755616295,179054527,1147987993,560
     8889925,95812952,1378733667,688768121,1380127404,1071140367,8449837
     933,2102381231,1803690001,1132759789,1812761173,811401052,615394580
     a,714455184,359573731/
      data r1c /1858168607,1521204822,474536859,974736896,1114310296,801
     1200246,1832712363,2063859275,849994344,731349012,1064506466,780680
     2804,1680867118,1530444240,575550303,6394604,957440394,1896974923,8
     348517698,296583298,1891896861,1693480484,979996076,418089926,12100
     459141,335938044,1364841020,33372505,425678949,1625600626,208556353
     59,237374938,2082340698,1400966182,1301717327,1946914281,1700397830
     6,959550392,81736978,2142584359,1198745608,1152306633,792325244,848
     7077619,841757441,773013240,1192873003,1667056045,1923639440,139011
     81718,1414781121,1582603830,862850389,1441126705,203364062,84032877
     9,476390896,1064504187,1897171744,1851237362,1926739714,488998756,1
     a932474921,111695461/
      data r1d /963757354,570544522,265037977,799660391,1346458889,11816
     142383,1906989838,1618194839,347499325,1748230931,283825706,1165976
     2470,27462760,1786557364,1393598092,1621652483,1329282346,144004840
     34,721467605,690775934,320143906,2016181203,956443044,732611239,161
     43258736,80903311,1790834073,707070495,1724865332,880405482,4571221
     532,1006551685,348518543,1698123395,854765942,1845325945,891960480,
     6981543163,1202170526,1648924828,276573191,1026555889,1016009069,18
     785978509,211667813,2071877637,1648093826,1162635757,1288904536,688
     8436696,379796093,54402138,1371679875,1053662729,1696335358,1643010
     9934,426964867,455972152,643781030,1315875809,1062175994,569057285,
     a545792532,187863609/
      data r2a /1071956961,1978709985,512535521,2042142689,1198822369,12
     103800033,983333857,1611165665,2013553633,1116755937,2141998049,186
     28054497,1368667105,1717577697,1841044449,665325537,1411646433,8587
     381665,80473057,150462433,2142491617,687851489,1155251169,323465185
     4,1413719009,1204787169,770411489,1184333793,1372812257,262105057,1
     5073437665,585584609,2019771361,7288801,2064329697,674701281,120711
     62673,440338401,1595603937,1451683809,1082319841,1561253857,1814744
     7033,769048545,1645392865,1222551521,574266337,774279137,748848097,
     81571715041,21654497,1467375585,540169185,461260769,156908513,70085
     94241,1019356129,38672353,980028385,622198753,38925281,303949793,34
     a3530465,1231409121/
      data r2b /1893843937,1257093089,394898401,381001697,141661153,7506
     118593,1134132193,218460129,1224827873,932009953,413748193,74378441
     27,848376801,1801267169,381230049,1956974561,1159791585,1210906593,
     31036577761,1710546913,11588577,1308411873,232307681,4501473,169873
     45073,2093783009,115903457,1133805537,1926263777,1419536353,6873650
     589,803491809,694174689,1433155553,1946692577,1161043937,149951457,
     62134640609,1746402273,58978273,293594081,1376507873,86494177,17922
     762113,1125102561,1306240993,1261935585,2065928161,496993249,192383
     89969,977759201,879976417,556749793,1081821153,1381448673,381890529
     9,1304372193,927668193,325520353,571670497,592376801,1461381089,210
     a4941537,1449316321/
      data r2c /568247265,535476193,277261281,867344353,1231983585,29743
     17153,1284930529,973238241,436102113,747263969,832981985,1766997985
     2,328086497,1884956641,1068899297,1101139937,907936737,1563031521,1
     3992682465,1123147745,28169185,1928972257,1456847841,1833021409,198
     43751137,835295201,1608879073,1083277281,332231649,429484001,301292
     5513,1021399009,1516061665,711538657,1829055457,1647386593,12402738
     689,1681459169,1897200609,813756385,1652351969,1191761889,505727969
     7,667992033,604812257,1389930465,1949604833,1210093537,245138401,12
     88481249,1933863905,292577249,573330401,1702381537,458505185,629268
     917,1589388257,1816664033,1818495969,521142241,1145828321,471328737
     a,1718868961,1667223521/
      data r2d /1390134241,1961342945,159624161,1353687009,174822369,199
     11739361,1435728865,1728016353,1794860001,562517985,1252215777,6427
     227905,1955279841,1968646113,1756568545,245305313,656081889,1915156
     3449,801303521,535748577,44749793,402048993,533904353,1514057697,12
     41283553,1724291041,954371041,1032749025,885683169,1586915297,20627
     503585,1239306209,190464993,2137405409,1711418337,2133729249,183112
     6673,1228277729,2047998945,1568534497,863626209,1007015905,92496176
     71,1691205601,84521953,1473619937,489790433,354258913,2140767201,48
     80606177,742484961,1852661729,589911009,175458273,1683045345,189144
     96753,1874404321,558176225,1163987937,470613985,1699279841,16287600
     a33,1332796385,1885130721/
      ks = kstrt - 2
      ndvp = (ndv - 1)/nvp + 1
      ndp = min0(ndvp,nvrp)
      mdp = nvp/min0(nvp,ndv)
      do 20 k = 1, nblok
      id = (k + ks)/mdp
      do 10 j = 1, ndp
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
      vran(j,k) = (dble(r1(l)) + dble(r2(l))*asc)*asc
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      function erfn(x)
c this function calculates the real error function, according to the
c formulae given in Abramowitz and Stegun, Handbook of Mathematical
c Functions, p. 299.  Error is < 1.5 x 10-7.
      implicit none
      real x
c local data
      real erfn, p, a1, a2, a3, a4, a5, t, f
      data p, a1, a2 /0.3275911,0.254829592,-0.284496736/
      data a3, a4, a5 /1.421413741,-1.453152027,1.061405429/
      save p, a1, a2, a3, a4, a5
      f = abs(x)
      t = 1.0/(1.0 + p*f)
      if (f.le.8.) then
         erfn = 1.0 - t*(a1 + t*(a2 + t*(a3 + t*(a4 + t*a5))))*exp(-x*x)
      else
         erfn = 1.0
      endif
      if (x.lt.0.) erfn = -erfn
      return
      end
c-----------------------------------------------------------------------
      function e1ln(x)
c this function calculates the sum of the exponential integral and the
c natural logarithm, according to the formulae given in Abramowitz and
c Stegun, Handbook of Mathematical Functions, p. 231.
c Error is < 2.0 x 10-7.
      implicit none
      real x
c local data
      real e1ln, a0, a1, a2, a3, a4, a5, b1, b2, b3, b4, c1, c2, c3, c4
      data a0, a1, a2 /-0.57721566,0.99999193,-0.24991055/
      data a3, a4, a5 /0.05519968,-0.00976004,0.00107857/
      data b1, b2, b3 /8.5733287401,18.0590169730,8.6347608925/
      data c1, c2, c3 /9.5733223454,25.6329561486,21.0996530827/
      data b4, c4 /0.2677737343,3.9584969228/
      save 
      if (x.le.1.0) then
         e1ln = a0 + x*(a1 + x*(a2 + x*(a3 + x*(a4 + x*a5))))
      else if (x.lt.50.0) then
         e1ln = alog(x) + (exp(-x)/x)*((b4 + x*(b3 + x*(b2 + x*(b1 + x))
     1))/(c4 + x*(c3 + x*(c2 + x*(c1 + x)))))
      else
         e1ln = alog(x)
      endif
      return
      end
