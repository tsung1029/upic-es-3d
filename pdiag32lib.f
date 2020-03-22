c 3d parallel PIC library for diagnostics
c with 2D domain decomposition
c written by viktor k. decyk, ucla
c copyright 1994, regents of the university of california
c update: march 28, 2008
c-----------------------------------------------------------------------
      subroutine PVDIST32(part,fv,fvm,npp,idimp,npmax,mnblok,nmv,nmvf)
c for 3d code, this subroutine calculates 3d velocity distribution
c and velocity moments, with 2D domain decomposition
c input: all, output: fv
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c fv = distribution function, number of particles in each velocity range
c maximum velocity (used for scaling) is contained in first element fv.
c vdrift for i-th dimension is contained in fvm(1,i)
c vth for i-th dimension is contained in fvm(2,i)
c entropy for i-th dimension is contained in fvm(3,i), defined to be:
c s/k = -sum(f(v)/np)*log(f(v)/(np*delta_v)).  Assumes that distribution is
c uniform in space and distributions in each dimension are independent.
c npp(m) = number of particles in partition m
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions
c the number of velocity bins used is 2*nmv + 1, nmvf >= 2*nmv+2
      implicit none
      integer idimp, npmax, mnblok, nmv, nmvf
      integer npp
      real part, fv, fvm
      dimension part(idimp,npmax,mnblok)
      dimension fv(nmvf,3,mnblok), fvm(3,3,mnblok)
      dimension npp(mnblok)
c local data
      integer j, m, nvx, nvy, nvz
      real anmv, svx, svy, svz
      double precision sumvx, sumvy, sumvz, sumvx2, sumvy2, sumvz2, anp
      double precision sum7, work7
      dimension sum7(7), work7(7)
      anmv = real(nmv)
      do 40 m = 1, mnblok
      svx = anmv/fv(1,1,m)
      svy = anmv/fv(1,2,m)
      svz = anmv/fv(1,3,m)
c zero out distribution
      do 10 j = 2, nmvf
      fv(j,1,m) = 0.
      fv(j,2,m) = 0.
      fv(j,3,m) = 0.
   10 continue
c count particles in each velocity region
      anp = 0.0d0
      anmv = anmv + 2.5
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvz = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      sumvz2 = 0.0d0
      do 20 j = 1, npp(m)
      anp = anp + 1.
      nvx = part(4,j,m)*svx + anmv
      sumvx = sumvx + part(4,j,m)
      sumvx2 = sumvx2 + part(4,j,m)**2
      nvy = part(5,j,m)*svy + anmv
      sumvy = sumvy + part(5,j,m)
      sumvy2 = sumvy2 + part(5,j,m)**2
      nvz = part(6,j,m)*svz + anmv
      sumvz = sumvz + part(6,j,m)
      sumvz2 = sumvz2 + part(6,j,m)**2
      if ((nvx.ge.2).and.(nvx.le.nmvf)) fv(nvx,1,m) = fv(nvx,1,m) + 1.
      if ((nvy.ge.2).and.(nvy.le.nmvf)) fv(nvy,2,m) = fv(nvy,2,m) + 1.
      if ((nvz.ge.2).and.(nvz.le.nmvf)) fv(nvz,3,m) = fv(nvz,3,m) + 1.
   20 continue
      sum7(1) = sumvx
      sum7(2) = sumvy
      sum7(3) = sumvz
      sum7(4) = sumvx2
      sum7(5) = sumvy2
      sum7(6) = sumvz2
      sum7(7) = anp
      call PDSUM(sum7,work7,7,1)
      sumvx = sum7(1)
      sumvy = sum7(2)
      sumvz = sum7(3)
      sumvx2 = sum7(4)
      sumvy2 = sum7(5)
      sumvz2 = sum7(6)
      anp = sum7(7)
c calculate velocity moments
      if (anp.ne.0.0d0) anp = 1.0d0/anp
      sumvx = sumvx*anp
      fvm(1,1,m) = sumvx
      fvm(2,1,m) = dsqrt(sumvx2*anp - sumvx**2)
      sumvy = sumvy*anp
      fvm(1,2,m) = sumvy
      fvm(2,2,m) = dsqrt(sumvy2*anp - sumvy**2)
      sumvz = sumvz*anp
      fvm(1,3,m) = sumvz
      fvm(2,3,m) = dsqrt(sumvz2*anp - sumvz**2)
c calculate entropy
      sumvx = 0.0d0
      sumvy = 0.0d0
      sumvz = 0.0d0
      sumvx2 = 0.0d0
      sumvy2 = 0.0d0
      sumvz2 = 0.0d0
      do 30 j = 2, nmvf
      if (fv(j,1,m).gt.0.) then
         sumvx = sumvx + fv(j,1,m)
         sumvx2 = sumvx2 + fv(j,1,m)*dlog(dble(fv(j,1,m)*svx))
      endif
      if (fv(j,2,m).gt.0.) then
         sumvy = sumvy + fv(j,2,m)
         sumvy2 = sumvy2 + fv(j,2,m)*dlog(dble(fv(j,2,m)*svy))
      endif
      if (fv(j,3,m).gt.0.) then
         sumvz = sumvz + fv(j,3,m)
         sumvz2 = sumvz2 + fv(j,3,m)*dlog(dble(fv(j,3,m)*svz))
      endif
   30 continue
      if (sumvx.gt.0.0d0) sumvx = -sumvx2/sumvx + dlog(sumvx)
      if (sumvy.gt.0.0d0) sumvy = -sumvy2/sumvy + dlog(sumvy)
      if (sumvz.gt.0.0d0) sumvz = -sumvz2/sumvz + dlog(sumvz)
      fvm(3,1,m) = sumvx
      fvm(3,2,m) = sumvy
      fvm(3,3,m) = sumvz
   40 continue
      return
      end
