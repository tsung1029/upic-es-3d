c modified by Jay Fahlen from PVVISTR32 to do random postions of
c particles distributed randomly in space, but all on the correct
c processor
      subroutine UNIF_RNDM_INIT(part,npp,nps,xs,ys,zs,lx,ly,lz,npx,npy,
     1npz,idimp,npmax,mblok,nblok,posx,posy,posz,kstrt,nvp,ndv,nvrp,
     2ierr)
c for 3 code, this subroutine calculates initial positions
c  for distributed data
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
c xs,ys,zs is the start position of this processor
c lx,ly,lz is the range over which the positions are distributed
c npx/npy/npz = initial number of particles distributed in x/y/z
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mblok/nblok = number of particle partitions in y/z
c posx,posy,posz = output arrays for parallel random numbers
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
      real lx, ly, lz, xs, ys, zs
      double precision posx, posy, posz
      integer npp, nps
      real part
      dimension part(idimp,npmax,mblok*nblok)
      dimension npp(mblok*nblok), nps(mblok*nblok)
      dimension posx(nvrp,mblok*nblok), posy(nvrp,mblok*nblok)
      dimension posz(nvrp,mblok*nblok)
c local data
c     integer*8 npxy, npxyz, ipp, nppv
      double precision dnpxy, dnpxyz, dt1
      integer ipp, nppv
      integer ks, npd, mdp, mnblok, moff, ii, n, nd, m, my, mz, id, is
      integer jj, i, j, npt, npxyzp, iwork
      real px, py, pz
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
      call prandom(posx,kstrt,nvp,nvrp,nd,nvrp,mnblok)
      call prandom(posy,kstrt,nvp,nvrp,nd,nvrp,mnblok)
      call prandom(posz,kstrt,nvp,nvrp,nd,nvrp,mnblok)
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
      px = lx*posx(jj,m) + xs
      py = ly*posy(jj,m) + ys
      pz = lz*posz(jj,m) + zs
c nth group keeps nth block of random numbers
      if (ii.eq.is) then
         npt = nps(m) + i - 1
         if (npt.le.npmax) then
            part(1,npt,m) = px
            part(2,npt,m) = py
            part(3,npt,m) = pz
         else
            ierr = ierr + 1
         endif
      endif
   10 continue
c update particle number
c      if (ii.eq.is) npp(m) = npp(m) + nvrp
   20 continue
   30 continue
   40 continue
   50 continue
c process errors
      if (ierr.gt.0) then
         if (kstrt.eq.1) then
            write (2,*) 'particle overflow error, ierr = ', ierr
         endif
      endif
      return
      end
