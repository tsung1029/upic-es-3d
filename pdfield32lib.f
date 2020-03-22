c 3d parallel PIC library for solving field equations with dirichlet
c boundary conditions and 2D domain decomposition
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: september 27, 2005
c-----------------------------------------------------------------------
      subroutine PLCGUARD32X(fxyz,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
c this subroutine replicates vector field so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation in the x direction
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c nx = system length in x direction
c nxe = first dimension of input array fxyz, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
c mnblok = number of particle partitions
c for distributed data with 2D decomposition
      implicit none
      real fxyz
      integer nyzp, nx, nxe, nypmx, nzpmx, idds, mnblok
      dimension fxyz(3,nxe,nypmx,nzpmx,mnblok), nyzp(idds,mnblok)
c local data
      integer i, k, l, m, nyp3, nzp3
      do 40 m = 1, mnblok
      nyp3 = nyzp(1,m) + 3
      nzp3 = nyzp(2,m) + 3
      do 30 l = 1, nzp3
      do 20 k = 1, nyp3
      do 10 i = 1, 3
      fxyz(i,1,k,l,m) = 2.*fxyz(i,2,k,l,m) - fxyz(i,3,k,l,m)
      fxyz(i,nx+3,k,l,m) = 2.*fxyz(i,nx+2,k,l,m) - fxyz(i,nx+1,k,l,m)
   10 continue
   20 continue
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLDGUARD32X(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
c this subroutine replicates scalar field so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation in the x direction
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c nx = system length in x direction
c nxe = first dimension of input array q, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
c mnblok = number of particle partitions
c for distributed data with 2D decomposition
      implicit none
      real q
      integer nyzp, nx, nxe, nypmx, nzpmx, idds, mnblok
      dimension q(nxe,nypmx,nzpmx,mnblok), nyzp(idds,mnblok)
c local data
      integer k, l, m, nyp3, nzp3
      do 30 m = 1, mnblok
      nyp3 = nyzp(1,m) + 3
      nzp3 = nyzp(2,m) + 3
      do 20 l = 1, nzp3
      do 10 k = 1, nyp3
      q(1,k,l,m) = 2.*q(2,k,l,m) - q(3,k,l,m)
      q(nx+3,k,l,m) = 2.*q(nx+2,k,l,m) - q(nx+1,k,l,m)
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLSCGUARD32X(cu,kstrt,nvpy,nvpz,noff,nyzp,xj0,yj0,zj0,n
     1x,ny,nz,ngx,ngy,ngz,nxe,nypmx,nzpmx,idds,mblok,nblok)
c initialize extended non-periodic field
c cu(i,j+1,k,l,m) = ith component of current density at grid point 
c (j,kk,ll), where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c xj0/yj0/zj0 = initialization constants in x/y/z direction
c nx/ny/nz = system length in x/y/z direction
c ngx/ngy/ngz = (0,1) = number of grid cells away from edge
c nxe = first dimension of current array, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
c mblok/nblok = number of particle partitions in y/z
c quadratic interpolation, for distributed data with 2D decomposition
      implicit none
      real cu, xj0, yj0, zj0
      integer noff, nyzp
      integer kstrt, nvpy, nvpz, nx, ny, nz, ngx, ngy, ngz, nxe
      integer nypmx, nzpmx, idds, mblok, nblok
      dimension cu(3,nxe,nypmx,nzpmx,mblok*nblok)
      dimension noff(idds,mblok*nblok), nyzp(idds,mblok*nblok)
      integer i, j, k, l, m, my, mz, js, ks, moff, kk, ll, nyp3, nzp3
      integer nxg, nx1, nx3
      real chx, chy, chz
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1).or.(ngz.l
     1t.0).or.(ngz.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx1 = nx + 1
      nx3 = nx + 3
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      chx = .5*xj0
      chy = .5*yj0
      chz = .5*zj0
      do 350 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 340 my = 1, mblok
      m = my + moff
      nyp3 = nyzp(1,m) + 3
      nzp3 = nyzp(2,m) + 3
c handle first grid point in y
      if ((my+js).eq.0) then
         do 30 l = 1, nzp3
         do 20 j = 1, nx3
         do 10 i = 1, 3
         cu(i,j,2,l,m) = 0.
   10    continue
   20    continue
   30    continue
      endif
c handle first grid point in z
      if ((mz+ks).eq.0) then
         do 60 k = 1, nyp3
         do 50 j = 1, nx3
         do 40 i = 1, 3
         cu(i,j,k,2,m) = 0.
   40    continue
   50    continue
   60    continue
      endif
c handle interior grid points in z
      do 230 l = 1, nyzp(2,m)
      ll = l + noff(2,m)
c guard cells in y
      do 80 j = 1, nx3
      do 70 i = 1, 3
      cu(i,j,1,l+1,m) = 0.
      cu(i,j,nyp3-1,l+1,m) = 0.
      cu(i,j,nyp3,l+1,m) = 0.
   70 continue
   80 continue
      if ((ll.ge.(ngz+2)).and.(ll.le.(nz-ngz))) then
c handle interior grid points in y
         do 130 k = 1, nyzp(1,m)
         kk = k + noff(1,m)
         if ((kk.ge.(ngy+2)).and.(kk.le.(ny-ngy))) then
            do 90 j = 2, nxg
            cu(1,j+ngx+1,k+1,l+1,m) = xj0
            cu(2,j+ngx+1,k+1,l+1,m) = yj0
            cu(3,j+ngx+1,k+1,l+1,m) = zj0
   90       continue
            do 100 i = 1, 3
            cu(i,1,k+1,l+1,m) = 0.
            cu(i,2,k+1,l+1,m) = 0.
            cu(i,nx+2,k+1,l+1,m) = 0.
            cu(i,nx+3,k+1,l+1,m) = 0.
  100       continue
            cu(1,ngx+2,k+1,l+1,m) = chx
            cu(2,ngx+2,k+1,l+1,m) = chy
            cu(3,ngx+2,k+1,l+1,m) = chz
            cu(1,nx-ngx+2,k+1,l+1,m) = chx
            cu(2,nx-ngx+2,k+1,l+1,m) = chy
            cu(3,nx-ngx+2,k+1,l+1,m) = chz
         else if ((kk.eq.(ngy+1)).or.(kk.eq.(ny-ngy+1))) then
            do 110 j = 2, nxg
            cu(1,j+ngx+1,k+1,l+1,m) = chx
            cu(2,j+ngx+1,k+1,l+1,m) = chy
            cu(3,j+ngx+1,k+1,l+1,m) = chz
  110       continue
            do 120 i = 1, 3
            cu(i,1,k+1,l+1,m) = 0.
            cu(i,2,k+1,l+1,m) = 0.
            cu(i,nx+2,k+1,l+1,m) = 0.
            cu(i,nx+3,k+1,l+1,m) = 0.
  120       continue
            cu(1,ngx+2,k+1,l+1,m) = .5*chx
            cu(2,ngx+2,k+1,l+1,m) = .5*chy
            cu(3,ngx+2,k+1,l+1,m) = .5*chz
            cu(1,nx-ngx+2,k+1,l+1,m) = .5*chx
            cu(2,nx-ngx+2,k+1,l+1,m) = .5*chy
            cu(3,nx-ngx+2,k+1,l+1,m) = .5*chz
         endif
  130    continue
c handle last grid point in y
         if (((nyzp(1,m)+noff(1,m)).lt.(ny-ngy+1)).and.((my+js).eq.(nvpy
     1-1))) then
            do 140 j = 2, nxg
            cu(1,j+ngx+1,nyp3-ngy-1,l+1,m) = chx
            cu(2,j+ngx+1,nyp3-ngy-1,l+1,m) = chy
            cu(3,j+ngx+1,nyp3-ngy-1,l+1,m) = chz
  140       continue
            do 150 i = 1, 3
            cu(i,1,nyp3-ngy-1,l+1,m) = 0.
            cu(i,2,nyp3-ngy-1,l+1,m) = 0.
            cu(i,nx+2,nyp3-ngy-1,l+1,m) = 0.
            cu(i,nx+3,nyp3-ngy-1,l+1,m) = 0.
  150       continue
            cu(1,ngx+2,nyp3-ngy-1,l+1,m) = .5*chx
            cu(2,ngx+2,nyp3-ngy-1,l+1,m) = .5*chy
            cu(3,ngx+2,nyp3-ngy-1,l+1,m) = .5*chz
            cu(1,nx-ngx+2,nyp3-ngy-1,l+1,m) = .5*chx
            cu(2,nx-ngx+2,nyp3-ngy-1,l+1,m) = .5*chy
            cu(3,nx-ngx+2,nyp3-ngy-1,l+1,m) = .5*chz
         endif
      else if ((ll.eq.(ngz+1)).or.(ll.eq.(nz-ngz+1))) then
         do 200 k = 1, nyzp(1,m)
         kk = k + noff(1,m)
         if ((kk.ge.(ngy+2)).and.(kk.le.(ny-ngy))) then
            do 160 j = 2, nxg
            cu(1,j+ngx+1,k+1,l+1,m) = chx
            cu(2,j+ngx+1,k+1,l+1,m) = chy
            cu(3,j+ngx+1,k+1,l+1,m) = chz
  160       continue
            do 170 i = 1, 3
            cu(i,1,k+1,l+1,m) = 0.
            cu(i,2,k+1,l+1,m) = 0.
            cu(i,nx+2,k+1,l+1,m) = 0.
            cu(i,nx+3,k+1,l+1,m) = 0.
  170       continue
            cu(1,ngx+2,k+1,l+1,m) = .5*chx
            cu(2,ngx+2,k+1,l+1,m) = .5*chy
            cu(3,ngx+2,k+1,l+1,m) = .5*chz
            cu(1,nx-ngx+2,k+1,l+1,m) = .5*chx
            cu(2,nx-ngx+2,k+1,l+1,m) = .5*chy
            cu(3,nx-ngx+2,k+1,l+1,m) = .5*chz
         else if ((kk.eq.(ngy+1)).or.(kk.eq.(ny-ngy+1))) then
            do 180 j = 2, nxg
            cu(1,j+ngx+1,k+1,l+1,m) = .5*chx
            cu(2,j+ngx+1,k+1,l+1,m) = .5*chy
            cu(3,j+ngx+1,k+1,l+1,m) = .5*chz
  180       continue
            do 190 i = 1, 3
            cu(i,1,k+1,l+1,m) = 0.
            cu(i,2,k+1,l+1,m) = 0.
            cu(i,nx+2,k+1,l+1,m) = 0.
            cu(i,nx+3,k+1,l+1,m) = 0.
  190       continue
            cu(1,ngx+2,k+1,l+1,m) = .25*chx
            cu(2,ngx+2,k+1,l+1,m) = .25*chy
            cu(3,ngx+2,k+1,l+1,m) = .25*chz
            cu(1,nx-ngx+2,k+1,l+1,m) = .25*chx
            cu(2,nx-ngx+2,k+1,l+1,m) = .25*chy
            cu(3,nx-ngx+2,k+1,l+1,m) = .25*chz
         endif
  200    continue
c handle last grid point in y
         if (((nyzp(1,m)+noff(1,m)).lt.(ny-ngy+1)).and.((my+js).eq.(nvpy
     1-1))) then
            do 210 j = 2, nxg
            cu(1,j+ngx+1,nyp3-ngy-1,l+1,m) = .5*chx
            cu(2,j+ngx+1,nyp3-ngy-1,l+1,m) = .5*chy
            cu(3,j+ngx+1,nyp3-ngy-1,l+1,m) = .5*chz
  210       continue
            do 220 i = 1, 3
            cu(i,1,nyp3-ngy-1,l+1,m) = 0.
            cu(i,2,nyp3-ngy-1,l+1,m) = 0.
            cu(i,nx+2,nyp3-ngy-1,l+1,m) = 0.
            cu(i,nx+3,nyp3-ngy-1,l+1,m) = 0.
  220       continue
            cu(1,ngx+2,nyp3-ngy-1,l+1,m) = .25*chx
            cu(2,ngx+2,nyp3-ngy-1,l+1,m) = .25*chy
            cu(3,ngx+2,nyp3-ngy-1,l+1,m) = .25*chz
            cu(1,nx-ngx+2,nyp3-ngy-1,l+1,m) = .25*chx
            cu(2,nx-ngx+2,nyp3-ngy-1,l+1,m) = .25*chy
            cu(3,nx-ngx+2,nyp3-ngy-1,l+1,m) = .25*chz
         endif
      endif
  230 continue
c guard cells in z
      do 260 k = 1, nyp3
      do 250 j = 1, nx3
      do 240 i = 1, 3
      cu(i,j,k,1,m) = 0.
      cu(i,j,k,nzp3-1,m) = 0.
      cu(i,j,k,nzp3,m) = 0.
  240 continue 
  250 continue
  260 continue
c handle last grid point in z
      if (((nyzp(2,m)+noff(2,m)).lt.(nz-ngz+1)).and.((mz+ks).eq.(nvpz-1)
     1)) then
         do 310 k = 1, nyzp(1,m)
         kk = k + noff(1,m)
         if ((kk.ge.(ngy+2)).and.(kk.le.(ny-ngy))) then
            do 270 j = 2, nxg
            cu(1,j+ngx+1,k+1,nzp3-ngz-1,m) = chx
            cu(2,j+ngx+1,k+1,nzp3-ngz-1,m) = chy
            cu(3,j+ngx+1,k+1,nzp3-ngz-1,m) = chz
  270       continue
            do 280 i = 1, 3
            cu(i,1,k+1,nzp3-ngz-1,m) = 0.
            cu(i,2,k+1,nzp3-ngz-1,m) = 0.
            cu(i,nx+2,k+1,nzp3-ngz-1,m) = 0.
            cu(i,nx+3,k+1,nzp3-ngz-1,m) = 0.
  280       continue
            cu(1,ngx+2,k+1,nzp3-ngz-1,m) = .5*chx
            cu(2,ngx+2,k+1,nzp3-ngz-1,m) = .5*chy
            cu(3,ngx+2,k+1,nzp3-ngz-1,m) = .5*chz
            cu(1,nx-ngx+2,k+1,nzp3-ngz-1,m) = .5*chx
            cu(2,nx-ngx+2,k+1,nzp3-ngz-1,m) = .5*chy
            cu(3,nx-ngx+2,k+1,nzp3-ngz-1,m) = .5*chz
         else if ((kk.eq.(ngy+1)).or.(kk.eq.(ny-ngy+1))) then
            do 290 j = 2, nxg
            cu(1,j+ngx+1,k+1,nzp3-ngz-1,m) = .5*chx
            cu(2,j+ngx+1,k+1,nzp3-ngz-1,m) = .5*chy
            cu(3,j+ngx+1,k+1,nzp3-ngz-1,m) = .5*chz
  290       continue
            do 300 i = 1, 3
            cu(i,1,k+1,nzp3-ngz-1,m) = 0.
            cu(i,2,k+1,nzp3-ngz-1,m) = 0.
            cu(i,nx+2,k+1,nzp3-ngz-1,m) = 0.
            cu(i,nx+3,k+1,nzp3-ngz-1,m) = 0.
  300       continue
            cu(1,ngx+2,k+1,nzp3-ngz-1,m) = .25*chx
            cu(2,ngx+2,k+1,nzp3-ngz-1,m) = .25*chy
            cu(3,ngx+2,k+1,nzp3-ngz-1,m) = .25*chz
            cu(1,nx-ngx+2,k+1,nzp3-ngz-1,m) = .25*chx
            cu(2,nx-ngx+2,k+1,nzp3-ngz-1,m) = .25*chy
            cu(3,nx-ngx+2,k+1,nzp3-ngz-1,m) = .25*chz
         endif
  310    continue
c handle last grid point in y
         if (((nyzp(1,m)+noff(1,m)).lt.(ny-ngy+1)).and.((my+js).eq.(nvpy
     1-1))) then
            do 320 j = 2, nxg
            cu(1,j+ngx+1,nyp3-ngy-1,nzp3-ngz-1,m) = .5*chx
            cu(2,j+ngx+1,nyp3-ngy-1,nzp3-ngz-1,m) = .5*chy
            cu(3,j+ngx+1,nyp3-ngy-1,nzp3-ngz-1,m) = .5*chz
  320       continue
            do 330 i = 1, 3
            cu(i,1,nyp3-ngy-1,nzp3-ngz-1,m) = 0.
            cu(i,2,nyp3-ngy-1,nzp3-ngz-1,m) = 0.
            cu(i,nx+2,nyp3-ngy-1,nzp3-ngz-1,m) = 0.
            cu(i,nx+3,nyp3-ngy-1,nzp3-ngz-1,m) = 0.
  330       continue
            cu(1,ngx+2,nyp3-ngy-1,nzp3-ngz-1,m) = .25*chx
            cu(2,ngx+2,nyp3-ngy-1,nzp3-ngz-1,m) = .25*chy
            cu(3,ngx+2,nyp3-ngy-1,nzp3-ngz-1,m) = .25*chz
            cu(1,nx-ngx+2,nyp3-ngy-1,nzp3-ngz-1,m) = .25*chx
            cu(2,nx-ngx+2,nyp3-ngy-1,nzp3-ngz-1,m) = .25*chy
            cu(3,nx-ngx+2,nyp3-ngy-1,nzp3-ngz-1,m) = .25*chz
         endif
      endif
  340 continue
  350 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLSGUARD32X(q,kstrt,nvpy,nvpz,noff,nyzp,qi0,nx,ny,nz,ng
     1x,ngy,ngz,nxe,nypmx,nzpmx,idds,mblok,nblok)
c initialize extended non-periodic scalar field
c q(j+1,k,l,m) = charge density at grid point (j,kk,ll), 
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c qi0 = initialization constant
c nx/ny/nz = system length in x/y/z direction
c ngx/ngy/ngz = (0,1) = number of grid cells away from edge
c nxe = first dimension of charge array, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
c mblok/nblok = number of particle partitions in y/z
c quadratic interpolation, for distributed data with 2D decomposition
      implicit none
      real q, qi0
      integer noff, nyzp
      integer kstrt, nvpy, nvpz, nx, ny, nz, ngx, ngy, ngz, nxe
      integer nypmx, nzpmx, idds, mblok, nblok
      dimension q(nxe,nypmx,nzpmx,mblok*nblok)
      dimension noff(idds,mblok*nblok), nyzp(idds,mblok*nblok)
      integer j, k, l, m, my, mz, js, ks, moff, kk, ll, nyp3, nzp3, nxg
      integer nx3
      real qh
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1).or.(ngz.l
     1t.0).or.(ngz.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx3 = nx + 3
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      qh = .5*qi0
      do 220 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 210 my = 1, mblok
      m = my + moff
      nyp3 = nyzp(1,m) + 3
      nzp3 = nyzp(2,m) + 3
c handle first grid point in y
      if ((my+js).eq.0) then
         do 20 l = 1, nzp3
         do 10 j = 1, nx3
         q(j,2,l,m) = 0.
   10    continue
   20    continue
      endif
c handle first grid point in z
      if ((mz+ks).eq.0) then
         do 40 k = 1, nyp3
         do 30 j = 1, nx3
         q(j,k,2,m) = 0.
   30    continue
   40    continue
      endif
c handle interior grid points in z
      do 140 l = 1, nyzp(2,m)
      ll = l + noff(2,m)
c guard cells in y
      do 50 j = 1, nx3
      q(j,1,l+1,m) = 0.
      q(j,nyp3-1,l+1,m) = 0.
      q(j,nyp3,l+1,m) = 0.
   50 continue
      if ((ll.ge.(ngz+2)).and.(ll.le.(nz-ngz))) then
c handle interior grid points in y
         do 80 k = 1, nyzp(1,m)
         kk = k + noff(1,m)
         if ((kk.ge.(ngy+2)).and.(kk.le.(ny-ngy))) then
            do 60 j = 2, nxg
            q(j+ngx+1,k+1,l+1,m) = qi0
   60       continue
            q(1,k+1,l+1,m) = 0.
            q(2,k+1,l+1,m) = 0.
            q(nx+2,k+1,l+1,m) = 0.
            q(nx+3,k+1,l+1,m) = 0.
            q(ngx+2,k+1,l+1,m) = qh
            q(nx-ngx+2,k+1,l+1,m) = qh
         else if ((kk.eq.(ngy+1)).or.(kk.eq.(ny-ngy+1))) then
            do 70 j = 2, nxg
            q(j+ngx+1,k+1,l+1,m) = qh
   70       continue
            q(1,k+1,l+1,m) = 0.
            q(2,k+1,l+1,m) = 0.
            q(nx+2,k+1,l+1,m) = 0.
            q(nx+3,k+1,l+1,m) = 0.
            q(ngx+2,k+1,l+1,m) = .5*qh
            q(nx-ngx+2,k+1,l+1,m) = .5*qh
         endif
   80    continue
c handle last grid point in y
         if (((nyzp(1,m)+noff(1,m)).lt.(ny-ngy+1)).and.((my+js).eq.(nvpy
     1-1))) then
            do 90 j = 2, nxg
            q(j+ngx+1,nyp3-ngy-1,l+1,m) = qh
   90       continue
            q(1,nyp3-ngy-1,l+1,m) = 0.
            q(2,nyp3-ngy-1,l+1,m) = 0.
            q(nx+2,nyp3-ngy-1,l+1,m) = 0.
            q(nx+3,nyp3-ngy-1,l+1,m) = 0.
            q(ngx+2,nyp3-ngy-1,l+1,m) = .5*qh
            q(nx-ngx+2,nyp3-ngy-1,l+1,m) = .5*qh
         endif
      else if ((ll.eq.(ngz+1)).or.(ll.eq.(nz-ngz+1))) then
         do 120 k = 1, nyzp(1,m)
         kk = k + noff(1,m)
         if ((kk.ge.(ngy+2)).and.(kk.le.(ny-ngy))) then
            do 100 j = 2, nxg
            q(j+ngx+1,k+1,l+1,m) = qh
  100       continue
            q(1,k+1,l+1,m) = 0.
            q(2,k+1,l+1,m) = 0.
            q(nx+2,k+1,l+1,m) = 0.
            q(nx+3,k+1,l+1,m) = 0.
            q(ngx+2,k+1,l+1,m) = .5*qh
            q(nx-ngx+2,k+1,l+1,m) = .5*qh
         else if ((kk.eq.(ngy+1)).or.(kk.eq.(ny-ngy+1))) then
            do 110 j = 2, nxg
            q(j+ngx+1,k+1,l+1,m) = .5*qh
  110       continue
            q(1,k+1,l+1,m) = 0.
            q(2,k+1,l+1,m) = 0.
            q(nx+2,k+1,l+1,m) = 0.
            q(nx+3,k+1,l+1,m) = 0.
            q(ngx+2,k+1,l+1,m) = .25*qh
            q(nx-ngx+2,k+1,l+1,m) = .25*qh
         endif
  120    continue
c handle last grid point in y
         if (((nyzp(1,m)+noff(1,m)).lt.(ny-ngy+1)).and.((my+js).eq.(nvpy
     1-1))) then
            do 130 j = 2, nxg
            q(j+ngx+1,nyp3-ngy-1,l+1,m) = .5*qh
  130       continue
            q(1,nyp3-ngy-1,l+1,m) = 0.
            q(2,nyp3-ngy-1,l+1,m) = 0.
            q(nx+2,nyp3-ngy-1,l+1,m) = 0.
            q(nx+3,nyp3-ngy-1,l+1,m) = 0.
            q(ngx+2,nyp3-ngy-1,l+1,m) = .25*qh
            q(nx-ngx+2,nyp3-ngy-1,l+1,m) = .25*qh
         endif
      endif
  140 continue
c guard cells in z
      do 160 k = 1, nyp3
      do 150 j = 1, nx3
      q(j,k,1,m) = 0.
      q(j,k,nzp3-1,m) = 0.
      q(j,k,nzp3,m) = 0.
  150 continue
  160 continue
c handle last grid point in z
      if (((nyzp(2,m)+noff(2,m)).lt.(nz-ngz+1)).and.((mz+ks).eq.(nvpz-1)
     1)) then
         do 190 k = 1, nyzp(1,m)
         kk = k + noff(1,m)
         if ((kk.ge.(ngy+2)).and.(kk.le.(ny-ngy))) then
            do 170 j = 2, nxg
            q(j+ngx+1,k+1,nzp3-ngz-1,m) = qh
  170       continue
            q(1,k+1,nzp3-ngz-1,m) = 0.
            q(2,k+1,nzp3-ngz-1,m) = 0.
            q(nx+2,k+1,nzp3-ngz-1,m) = 0.
            q(nx+3,k+1,nzp3-ngz-1,m) = 0.
            q(ngx+2,k+1,nzp3-ngz-1,m) = .5*qh
            q(nx-ngx+2,k+1,nzp3-ngz-1,m) = .5*qh
         else if ((kk.eq.(ngy+1)).or.(kk.eq.(ny-ngy+1))) then
            do 180 j = 2, nxg
            q(j+ngx+1,k+1,nzp3-ngz-1,m) = .5*qh
  180       continue
            q(1,k+1,nzp3-ngz-1,m) = 0.
            q(2,k+1,nzp3-ngz-1,m) = 0.
            q(nx+2,k+1,nzp3-ngz-1,m) = 0.
            q(nx+3,k+1,nzp3-ngz-1,m) = 0.
            q(ngx+2,k+1,nzp3-ngz-1,m) = .25*qh
            q(nx-ngx+2,k+1,nzp3-ngz-1,m) = .25*qh
         endif
  190    continue
c handle last grid point in y
         if (((nyzp(1,m)+noff(1,m)).lt.(ny-ngy+1)).and.((my+js).eq.(nvpy
     1-1))) then
            do 200 j = 2, nxg
            q(j+ngx+1,nyp3-ngy-1,nzp3-ngz-1,m) = .5*qh
  200       continue
            q(1,nyp3-ngy-1,nzp3-ngz-1,m) = 0.
            q(2,nyp3-ngy-1,nzp3-ngz-1,m) = 0.
            q(nx+2,nyp3-ngy-1,nzp3-ngz-1,m) = 0.
            q(nx+3,nyp3-ngy-1,nzp3-ngz-1,m) = 0.
            q(ngx+2,nyp3-ngy-1,nzp3-ngz-1,m) = .25*qh
            q(nx-ngx+2,nyp3-ngy-1,nzp3-ngz-1,m) = .25*qh
         endif
      endif
  210 continue
  220 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLACGUARD32X(cu,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
c this subroutine adds adds up guard cells so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation in the x direction
c for vector field
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c nx = system length in x direction
c nxe = first dimension of input array q, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition, including guard cells
c idds = dimensionality of domain decomposition
c mnblok = number of particle partitions
c for distributed data with 2D decomposition
      implicit none
      real cu
      integer nyzp, nx, nxe, nypmx, nzpmx, idds, mnblok
      dimension cu(3,nxe,nypmx,nzpmx,mnblok), nyzp(idds,mnblok)
c local data
      integer i, k, l, m, nyp3, nzp3
c add up guard cells
      do 40 m = 1, mnblok
      nyp3 = nyzp(1,m) + 3
      nzp3 = nyzp(2,m) + 3
      do 30 l = 1, nzp3
      do 20 k = 1, nyp3
      do 10 i = 1, 3
      cu(i,2,k,l,m) = cu(i,2,k,l,m) + 2.*cu(i,1,k,l,m)
      cu(i,3,k,l,m) = cu(i,3,k,l,m) - cu(i,1,k,l,m)
      cu(i,nx+1,k,l,m) = cu(i,nx+1,k,l,m) - cu(i,nx+3,k,l,m)
      cu(i,nx+2,k,l,m) = cu(i,nx+2,k,l,m) + 2.*cu(i,nx+3,k,l,m)
      cu(i,1,k,l,m) = 0.
      cu(i,nx+3,k,l,m) = 0.
   10 continue
   20 continue
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLAGUARD32X(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
c this subroutine adds adds up guard cells so as to disable
c quadratic interpolation within a half a cell of the edges,
c and reduce it to linear interpolation in the x direction
c for scalar field
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c nx = system length in x direction
c nxe = first dimension of input array q, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition, including guard cells
c idds = dimensionality of domain decomposition
c mnblok = number of particle partitions
c for distributed data with 2D decomposition
      implicit none
      real q
      integer nyzp, nx, nxe, nypmx, nzpmx, idds, mnblok
      dimension q(nxe,nypmx,nzpmx,mnblok), nyzp(idds,mnblok)
c local data
      integer k, l, m, nyp3, nzp3
c add up guard cells
      do 30 m = 1, mnblok
      nyp3 = nyzp(1,m) + 3
      nzp3 = nyzp(2,m) + 3
      do 20 l = 1, nzp3
      do 10 k = 1, nyp3
      q(2,k,l,m) = q(2,k,l,m) + 2.*q(1,k,l,m)
      q(3,k,l,m) = q(3,k,l,m) - q(1,k,l,m)
      q(nx+1,k,l,m) = q(nx+1,k,l,m) - q(nx+3,k,l,m)
      q(nx+2,k,l,m) = q(nx+2,k,l,m) + 2.*q(nx+3,k,l,m)
      q(1,k,l,m) = 0.
      q(nx+3,k,l,m) = 0.
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLSCGUARD32XL(cu,kstrt,nvpy,nvpz,noff,nyzp,xj0,yj0,zj0,
     1nx,ny,nz,ngx,ngy,ngz,nxe,nypmx,nzpmx,idds,mblok,nblok)
c initialize extended non-periodic field
c cu(i,j,k,l,m) = ith component of current density at grid point
c (j,kk,ll), where kk = k + noff(1,m), and ll = l + noff(2,m)
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c xj0/yj0/zj0 = initialization constants in x/y/z direction
c nx/ny/nz = system length in x/y/z direction
c ngx/ngy/ngz = (0,1) = number of grid cells away from edge
c nxe = first dimension of current array, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
c mblok/nblok = number of particle partitions in y/z
c linear interpolation, for distributed data with 2D decomposition
      implicit none
      real cu, xj0, yj0, zj0
      integer noff, nyzp
      integer kstrt, nvpy, nvpz, nx, ny, nz, ngx, ngy, ngz, nxe
      integer nypmx, nzpmx, idds, mblok, nblok
      dimension cu(3,nxe,nypmx,nzpmx,mblok*nblok)
      dimension noff(idds,mblok*nblok), nyzp(idds,mblok*nblok)
      integer i, j, k, l, m, my, mz, js, ks, moff, kk, ll, nyp1, nzp1
      integer nxg, nx1
      real chx, chy, chz
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1).or.(ngz.l
     1t.0).or.(ngz.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx1 = nx + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      chx = .5*xj0
      chy = .5*yj0
      chz = .5*zj0
      do 350 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 340 my = 1, mblok
      m = my + moff
      nyp1 = nyzp(1,m) + 1
      nzp1 = nyzp(2,m) + 1
c handle first grid point in y
      if ((my+js).eq.0) then
         do 30 l = 1, nzp1
         do 20 j = 1, nx1
         do 10 i = 1, 3
         cu(i,j,1,l,m) = 0.
   10    continue
   20    continue
   30    continue
      endif
c handle first grid point in z
      if ((mz+ks).eq.0) then
         do 60 k = 1, nyp1
         do 50 j = 1, nx1
         do 40 i = 1, 3
         cu(i,j,k,1,m) = 0.
   40    continue
   50    continue
   60    continue
      endif
c handle interior grid points in z
      do 230 l = 1, nyzp(2,m)
      ll = l + noff(2,m)
c guard cells in y
      do 80 j = 1, nx1
      do 70 i = 1, 3
      cu(i,j,nyp1,l,m) = 0.
   70 continue
   80 continue
      if ((ll.ge.(ngz+2)).and.(ll.le.(nz-ngz))) then
c handle interior grid points in y
         do 130 k = 1, nyzp(1,m)
         kk = k + noff(1,m)
         if ((kk.ge.(ngy+2)).and.(kk.le.(ny-ngy))) then
            do 90 j = 2, nxg
            cu(1,j+ngx,k,l,m) = xj0
            cu(2,j+ngx,k,l,m) = yj0
            cu(3,j+ngx,k,l,m) = zj0
   90       continue
            do 100 i = 1, 3
            cu(i,1,k,l,m) = 0.
            cu(i,nx+1,k,l,m) = 0.
  100       continue
            cu(1,ngx+1,k,l,m) = chx
            cu(2,ngx+1,k,l,m) = chy
            cu(3,ngx+1,k,l,m) = chz
            cu(1,nx-ngx+1,k,l,m) = chx
            cu(2,nx-ngx+1,k,l,m) = chy
            cu(3,nx-ngx+1,k,l,m) = chz
         else if ((kk.eq.(ngy+1)).or.(kk.eq.(ny-ngy+1))) then
            do 110 j = 2, nxg
            cu(1,j+ngx,k,l,m) = chx
            cu(2,j+ngx,k,l,m) = chy
            cu(3,j+ngx,k,l,m) = chz
  110       continue
            do 120 i = 1, 3
            cu(i,1,k,l,m) = 0.
            cu(i,nx+1,k,l,m) = 0.
  120       continue
            cu(1,ngx+1,k,l,m) = .5*chx
            cu(2,ngx+1,k,l,m) = .5*chy
            cu(3,ngx+1,k,l,m) = .5*chz
            cu(1,nx-ngx+1,k,l,m) = .5*chx
            cu(2,nx-ngx+1,k,l,m) = .5*chy
            cu(3,nx-ngx+1,k,l,m) = .5*chz
         endif
  130    continue
c handle last grid point in y
         if (((nyzp(1,m)+noff(1,m)).lt.(ny-ngy+1)).and.((my+js).eq.(nvpy
     1-1))) then
            do 140 j = 2, nxg
            cu(1,j+ngx,nyp1-ngy,l,m) = chx
            cu(2,j+ngx,nyp1-ngy,l,m) = chy
            cu(3,j+ngx,nyp1-ngy,l,m) = chz
  140       continue
            do 150 i = 1, 3
            cu(i,1,nyp1-ngy,l,m) = 0.
            cu(i,nx+1,nyp1-ngy,l,m) = 0.
  150       continue
            cu(1,ngx+1,nyp1-ngy,l,m) = .5*chx
            cu(2,ngx+1,nyp1-ngy,l,m) = .5*chy
            cu(3,ngx+1,nyp1-ngy,l,m) = .5*chz
            cu(1,nx-ngx+1,nyp1-ngy,l,m) = .5*chx
            cu(2,nx-ngx+1,nyp1-ngy,l,m) = .5*chy
            cu(3,nx-ngx+1,nyp1-ngy,l,m) = .5*chz
         endif
      else if ((ll.eq.(ngz+1)).or.(ll.eq.(nz-ngz+1))) then
         do 200 k = 1, nyzp(1,m)
         kk = k + noff(1,m)
         if ((kk.ge.(ngy+2)).and.(kk.le.(ny-ngy))) then
            do 160 j = 2, nxg
            cu(1,j+ngx,k,l,m) = chx
            cu(2,j+ngx,k,l,m) = chy
            cu(3,j+ngx,k,l,m) = chz
  160       continue
            do 170 i = 1, 3
            cu(i,1,k,l,m) = 0.
            cu(i,nx+1,k,l,m) = 0.
  170       continue
            cu(1,ngx+1,k,l,m) = .5*chx
            cu(2,ngx+1,k,l,m) = .5*chy
            cu(3,ngx+1,k,l,m) = .5*chz
            cu(1,nx-ngx+1,k,l,m) = .5*chx
            cu(2,nx-ngx+1,k,l,m) = .5*chy
            cu(3,nx-ngx+1,k,l,m) = .5*chz
         else if ((kk.eq.(ngy+1)).or.(kk.eq.(ny-ngy+1))) then
            do 180 j = 2, nxg
            cu(1,j+ngx,k,l,m) = .5*chx
            cu(2,j+ngx,k,l,m) = .5*chy
            cu(3,j+ngx,k,l,m) = .5*chz
  180       continue
            do 190 i = 1, 3
            cu(i,1,k,l,m) = 0.
            cu(i,nx+1,k,l,m) = 0.
  190       continue
            cu(1,ngx+1,k,l,m) = .25*chx
            cu(2,ngx+1,k,l,m) = .25*chy
            cu(3,ngx+1,k,l,m) = .25*chz
            cu(1,nx-ngx+1,k,l,m) = .25*chx
            cu(2,nx-ngx+1,k,l,m) = .25*chy
            cu(3,nx-ngx+1,k,l,m) = .25*chz
         endif
  200    continue
c handle last grid point in y
         if (((nyzp(1,m)+noff(1,m)).lt.(ny-ngy+1)).and.((my+js).eq.(nvpy
     1-1))) then
            do 210 j = 2, nxg
            cu(1,j+ngx,nyp1-ngy,l,m) = .5*chx
            cu(2,j+ngx,nyp1-ngy,l,m) = .5*chy
            cu(3,j+ngx,nyp1-ngy,l,m) = .5*chz
  210       continue
            do 220 i = 1, 3
            cu(i,1,nyp1-ngy,l,m) = 0.
            cu(i,nx+1,nyp1-ngy,l,m) = 0.
  220       continue
            cu(1,ngx+1,nyp1-ngy,l,m) = .25*chx
            cu(2,ngx+1,nyp1-ngy,l,m) = .25*chy
            cu(3,ngx+1,nyp1-ngy,l,m) = .25*chz
            cu(1,nx-ngx+1,nyp1-ngy,l,m) = .25*chx
            cu(2,nx-ngx+1,nyp1-ngy,l,m) = .25*chy
            cu(3,nx-ngx+1,nyp1-ngy,l,m) = .25*chz
         endif
      endif
  230 continue
c guard cells in z
      do 260 k = 1, nyp1
      do 250 j = 1, nx1
      do 240 i = 1, 3
      cu(i,j,k,nzp1,m) = 0.
  240 continue 
  250 continue
  260 continue
c handle last grid point in z
      if (((nyzp(2,m)+noff(2,m)).lt.(nz-ngz+1)).and.((mz+ks).eq.(nvpz-1)
     1)) then
         do 310 k = 1, nyzp(1,m)
         kk = k + noff(1,m)
         if ((kk.ge.(ngy+2)).and.(kk.le.(ny-ngy))) then
            do 270 j = 2, nxg
            cu(1,j+ngx,k,nzp1-ngz,m) = chx
            cu(2,j+ngx,k,nzp1-ngz,m) = chy
            cu(3,j+ngx,k,nzp1-ngz,m) = chz
  270       continue
            do 280 i = 1, 3
            cu(i,1,k,nzp1-ngz,m) = 0.
            cu(i,nx+1,k,nzp1-ngz,m) = 0.
  280       continue
            cu(1,ngx+1,k,nzp1-ngz,m) = .5*chx
            cu(2,ngx+1,k,nzp1-ngz,m) = .5*chy
            cu(3,ngx+1,k,nzp1-ngz,m) = .5*chz
            cu(1,nx-ngx+1,k,nzp1-ngz,m) = .5*chx
            cu(2,nx-ngx+1,k,nzp1-ngz,m) = .5*chy
            cu(3,nx-ngx+1,k,nzp1-ngz,m) = .5*chz
         else if ((kk.eq.(ngy+1)).or.(kk.eq.(ny-ngy+1))) then
            do 290 j = 2, nxg
            cu(1,j+ngx,k,nzp1-ngz,m) = .5*chx
            cu(2,j+ngx,k,nzp1-ngz,m) = .5*chy
            cu(3,j+ngx,k,nzp1-ngz,m) = .5*chz
  290       continue
            do 300 i = 1, 3
            cu(i,1,k,nzp1-ngz,m) = 0.
            cu(i,nx+1,k,nzp1-ngz,m) = 0.
  300       continue
            cu(1,ngx+1,k,nzp1-ngz,m) = .25*chx
            cu(2,ngx+1,k,nzp1-ngz,m) = .25*chy
            cu(3,ngx+1,k,nzp1-ngz,m) = .25*chz
            cu(1,nx-ngx+1,k,nzp1-ngz,m) = .25*chx
            cu(2,nx-ngx+1,k,nzp1-ngz,m) = .25*chy
            cu(3,nx-ngx+1,k,nzp1-ngz,m) = .25*chz
         endif
  310    continue
c handle last grid point in y
         if (((nyzp(1,m)+noff(1,m)).lt.(ny-ngy+1)).and.((my+js).eq.(nvpy
     1-1))) then
            do 320 j = 2, nxg
            cu(1,j+ngx,nyp1-ngy,nzp1-ngz,m) = .5*chx
            cu(2,j+ngx,nyp1-ngy,nzp1-ngz,m) = .5*chy
            cu(3,j+ngx,nyp1-ngy,nzp1-ngz,m) = .5*chz
  320       continue
            do 330 i = 1, 3
            cu(i,1,nyp1-ngy,nzp1-ngz,m) = 0.
            cu(i,nx+1,nyp1-ngy,nzp1-ngz,m) = 0.
  330       continue
            cu(1,ngx+1,nyp1-ngy,nzp1-ngz,m) = .25*chx
            cu(2,ngx+1,nyp1-ngy,nzp1-ngz,m) = .25*chy
            cu(3,ngx+1,nyp1-ngy,nzp1-ngz,m) = .25*chz
            cu(1,nx-ngx+1,nyp1-ngy,nzp1-ngz,m) = .25*chx
            cu(2,nx-ngx+1,nyp1-ngy,nzp1-ngz,m) = .25*chy
            cu(3,nx-ngx+1,nyp1-ngy,nzp1-ngz,m) = .25*chz
         endif
      endif
  340 continue
  350 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PLSGUARD32XL(q,kstrt,nvpy,nvpz,noff,nyzp,qi0,nx,ny,nz,n
     1gx,ngy,ngz,nxe,nypmx,nzpmx,idds,mblok,nblok)
c initialize extended non-periodic field
c q(j,k,l,m) = charge density at grid point (j,kk,ll),
c where kk = k + noff(1,m), and ll = l + noff(2,m)
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c qi0 = initialization constant
c nx/ny/nz = system length in x/y/z direction
c ngx/ngy/ngz = (0,1) = number of grid cells away from edge
c nxe = first dimension of current array, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
c mblok/nblok = number of particle partitions in y/z
c linear interpolation, for distributed data with 2D decomposition
      implicit none
      real q, qi0
      integer noff, nyzp
      integer kstrt, nvpy, nvpz, nx, ny, nz, ngx, ngy, ngz, nxe
      integer nypmx, nzpmx, idds, mblok, nblok
      dimension q(nxe,nypmx,nzpmx,mblok*nblok)
      dimension noff(idds,mblok*nblok), nyzp(idds,mblok*nblok)
      integer j, k, l, m, my, mz, js, ks, moff, kk, ll, nyp1, nzp1, nxg
      integer nx1
      real qh
      if ((ngx.lt.0).or.(ngx.gt.1).or.(ngy.lt.0).or.(ngy.gt.1).or.(ngz.l
     1t.0).or.(ngz.gt.1)) return
c initialize extended field, with zero in the edges
      nxg = nx - 2*ngx
      nx1 = nx + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      qh = .5*qi0
      do 220 mz = 1, nblok
      moff = mblok*(mz - 1)
      do 210 my = 1, mblok
      m = my + moff
      nyp1 = nyzp(1,m) + 1
      nzp1 = nyzp(2,m) + 1
c handle first grid point in y
      if ((my+js).eq.0) then
         do 20 l = 1, nzp1
         do 10 j = 1, nx1
         q(j,1,l,m) = 0.
   10    continue
   20    continue
      endif
c handle first grid point in z
      if ((mz+ks).eq.0) then
         do 40 k = 1, nyp1
         do 30 j = 1, nx1
         q(j,k,1,m) = 0.
   30    continue
   40    continue
      endif
c handle interior grid points in z
      do 140 l = 1, nyzp(2,m)
      ll = l + noff(2,m)
c guard cells in y
      do 50 j = 1, nx1
      q(j,nyp1,l,m) = 0.
   50 continue
      if ((ll.ge.(ngz+2)).and.(ll.le.(nz-ngz))) then
c handle interior grid points in y
         do 80 k = 1, nyzp(1,m)
         kk = k + noff(1,m)
         if ((kk.ge.(ngy+2)).and.(kk.le.(ny-ngy))) then
            do 60 j = 2, nxg
            q(j+ngx,k,l,m) = qi0
   60       continue
            q(1,k,l,m) = 0.
            q(nx+1,k,l,m) = 0.
            q(ngx+1,k,l,m) = qh
            q(nx-ngx+1,k,l,m) = qh
         else if ((kk.eq.(ngy+1)).or.(kk.eq.(ny-ngy+1))) then
            do 70 j = 2, nxg
            q(j+ngx,k,l,m) = qh
   70       continue
            q(1,k,l,m) = 0.
            q(nx+1,k,l,m) = 0.
            q(ngx+1,k,l,m) = .5*qh
            q(nx-ngx+1,k,l,m) = .5*qh
         endif
   80    continue
c handle last grid point in y
         if (((nyzp(1,m)+noff(1,m)).lt.(ny-ngy+1)).and.((my+js).eq.(nvpy
     1-1))) then
            do 90 j = 2, nxg
            q(j+ngx,nyp1-ngy,l,m) = qh
   90       continue
            q(1,nyp1-ngy,l,m) = 0.
            q(nx+1,nyp1-ngy,l,m) = 0.
            q(ngx+1,nyp1-ngy,l,m) = .5*qh
            q(nx-ngx+1,nyp1-ngy,l,m) = .5*qh
         endif
      else if ((ll.eq.(ngz+1)).or.(ll.eq.(nz-ngz+1))) then
         do 120 k = 1, nyzp(1,m)
         kk = k + noff(1,m)
         if ((kk.ge.(ngy+2)).and.(kk.le.(ny-ngy))) then
            do 100 j = 2, nxg
            q(j+ngx,k,l,m) = qh
  100       continue
            q(1,k,l,m) = 0.
            q(nx+1,k,l,m) = 0.
            q(ngx+1,k,l,m) = .5*qh
            q(nx-ngx+1,k,l,m) = .5*qh
         else if ((kk.eq.(ngy+1)).or.(kk.eq.(ny-ngy+1))) then
            do 110 j = 2, nxg
            q(j+ngx,k,l,m) = .5*qh
  110       continue
            q(1,k,l,m) = 0.
            q(nx+1,k,l,m) = 0.
            q(ngx+1,k,l,m) = .25*qh
            q(nx-ngx+1,k,l,m) = .25*qh
         endif
  120    continue
c handle last grid point in y
         if (((nyzp(1,m)+noff(1,m)).lt.(ny-ngy+1)).and.((my+js).eq.(nvpy
     1-1))) then
            do 130 j = 2, nxg
            q(j+ngx,nyp1-ngy,l,m) = .5*qh
  130       continue
            q(1,nyp1-ngy,l,m) = 0.
            q(nx+1,nyp1-ngy,l,m) = 0.
            q(ngx+1,nyp1-ngy,l,m) = .25*qh
            q(nx-ngx+1,nyp1-ngy,l,m) = .25*qh
         endif
      endif
  140 continue
c guard cells in z
      do 160 k = 1, nyp1
      do 150 j = 1, nx1
      q(j,k,nzp1,m) = 0.
  150 continue
  160 continue
c handle last grid point in z
      if (((nyzp(2,m)+noff(2,m)).lt.(nz-ngz+1)).and.((mz+ks).eq.(nvpz-1)
     1)) then
         do 190 k = 1, nyzp(1,m)
         kk = k + noff(1,m)
         if ((kk.ge.(ngy+2)).and.(kk.le.(ny-ngy))) then
            do 170 j = 2, nxg
            q(j+ngx,k,nzp1-ngz,m) = qh
  170       continue
            q(1,k,nzp1-ngz,m) = 0.
            q(nx+1,k,nzp1-ngz,m) = 0.
            q(ngx+1,k,nzp1-ngz,m) = .5*qh
            q(nx-ngx+1,k,nzp1-ngz,m) = .5*qh
         else if ((kk.eq.(ngy+1)).or.(kk.eq.(ny-ngy+1))) then
            do 180 j = 2, nxg
            q(j+ngx,k,nzp1-ngz,m) = .5*qh
  180       continue
            q(1,k,nzp1-ngz,m) = 0.
            q(nx+1,k,nzp1-ngz,m) = 0.
            q(ngx+1,k,nzp1-ngz,m) = .25*qh
            q(nx-ngx+1,k,nzp1-ngz,m) = .25*qh
         endif
  190    continue
c handle last grid point in y
         if (((nyzp(1,m)+noff(1,m)).lt.(ny-ngy+1)).and.((my+js).eq.(nvpy
     1-1))) then
            do 200 j = 2, nxg
            q(j+ngx,nyp1-ngy,nzp1-ngz,m) = .5*qh
  200       continue
            q(1,nyp1-ngy,nzp1-ngz,m) = 0.
            q(nx+1,nyp1-ngy,nzp1-ngz,m) = 0.
            q(ngx+1,nyp1-ngy,nzp1-ngz,m) = .25*qh
            q(nx-ngx+1,nyp1-ngy,nzp1-ngz,m) = .25*qh
         endif
      endif
  210 continue
  220 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISDX32(q,fx,fy,fz,isign,ffd,ax,ay,az,affp,we,nx,ny,n
     1z,kstrt,nz2d,kxyp2,kyzp2,j2blok,m2blok,nzd)
c this subroutine solves 3d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides smoothing a smoothing function,
c with dirichlet boundary conditions (zero potential)
c for distributed data with 2D spatial decomposition
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0, input: isign,ax,ay,az,affp,nx,ny,nz,kstrt,kxyp2,kyzp2,
c                       j2blok,m2blok,nzd
c output: ffd
c for isign = -1, input: q,ffd,isign,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2,
c                        j2blok,m2blok,nzd
c output: fx,fy,fz,we
c approximate flop count is: 12*nxc*nyc*nzc
c for isign = 1, input: q,ffd,isign,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2,
c                       j2blok,m2blok,nzd
c output: fx,we
c approximate flop count is: 6*nxc*nyc*nzc
c for isign = 2, input: q,ffd,isign,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2,
c                       j2blok,m2blok,nzd
c output: fy
c approximate flop count is: 2*nxc*nyc*nzc
c where nxc = (nx-1)/nvpy, nyc = (ny-1)/nvpz, nzc = nz - 1, and
c nvpy/nvpz = number of procs in y/z
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky,kz) = -sqrt(-1)*kx*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c fy(kx,ky,kz) = -sqrt(-1)*ky*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c fz(kx,ky,kz) = -sqrt(-1)*kz*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2)
c if isign = 1, potential is calculated using the equation:
c fx(kx,ky,kz) = g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz)
c if isign = 2, smoothing is calculated using the equation:
c fy(kx,ky,kz) = q(kx,ky,kz)*s(kx,ky,kz)
c q(l,j,k,m) = complex charge density for fourier mode jj-1,kk-1,l-1
c fx(l,j,k,m) = x component of force/charge
c fy(l,j,k,m) = y component of force/charge
c fz(l,j,k,m) = z component of force/charge
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp2*(mx - 1) and
c kk = k + kyzp2*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c j2blok/m2blok = number of field partitions in x/y
c kxyp2/kyzp2 = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c if isign = 0, form factor array is prepared
c aimag(ffd(l,j,k,m)) = finite-size particle shape factor s
c real(ffd(l,j,k,m)) = potential green's function g
c for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp2*(mx - 1) and
c kk = k + kyzp2*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c ax/ay/az = half-width of particle in x/y/z direction
c affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*
c    |q(kx,ky,kz)*s(kx,ky,kz)|**2)
c nx/ny/nz = system length in x/y/z direction
c nz2d = first dimension of field arrays, must be >= 2*nz
c nzd = first dimension of form factor array, must be >= nz
      double precision wp
      complex q, fx, fy, fz, ffd, zero
      dimension q(nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension fx(nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension fy(nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension fz(nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension ffd(nzd,kxyp2,kyzp2,j2blok*m2blok)
      nyy = ny + ny
      nzz = nz + nz
      nz2 = nzz + 2
      kxb = nx/kxyp2
      kyzb = nyy/kyzp2
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(nyy)
      dnz = 6.28318530717959/float(nzz)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 60
      if (kstrt.gt.(kxb*kyzb)) return
c prepare form factor array
      do 50 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 40 mx = 1, j2blok
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 30 k = 1, kyzp2
      k1 = k + koff
      if (k1.gt.ny) k1 = k1 - nyy
      dky = dny*float(k1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 20 j = 1, kxyp2
      dkx = dnx*float(j + joff)
      at3 = dkx*dkx + at1
      at4 = (dkx*ax)**2 + at2
      do 10 l = 1, nz
      dkz = dnz*float(l - 1)
      at5 = dkz*dkz + at3
      at6 = exp(-.5*((dkz*az)**2 + at4))
      if (at5.eq.0.) then
         ffd(l,j,k,m) = cmplx(affp,1.)
      else
         ffd(l,j,k,m) = cmplx(affp*at6/at5,at6)
      endif
   10 continue
   20 continue
   30 continue
   40 continue
   50 continue
      return
   60 if (isign.gt.0) go to 180
c calculate force/charge and sum field energy
      wp = 0.0d0
      if (kstrt.gt.(kxb*kyzb)) go to 170
      do 160 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 150 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 100 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         if (k1.gt.ny) k1 = k1 - nyy
         dky = dny*float(k1)
         do 80 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 70 l = 2, nz
            l1 = nz2 - l
            at1 = real(ffd(l,j,k,m))*aimag(ffd(l,j,k,m))
            at4 = aimag(q(l,j,k,m))*at1
            at2 = dkx*at4
            at3 = dky*at4
            at4 = dnz*float(l - 1)*at4
            fx(l,j,k,m) = cmplx(at2,0.)
            fx(l1,j,k,m) = cmplx(-at2,0.)
            fy(l,j,k,m) = cmplx(at3,0.)
            fy(l1,j,k,m) = cmplx(-at3,0.)
            fz(l,j,k,m) = cmplx(at4,0.)
            fz(l1,j,k,m) = cmplx(at4,0.)
            wp = wp + at1*aimag(q(l,j,k,m))**2
   70       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            fx(1,j,k,m) = zero
            fx(l1,j,k,m) = zero
            fy(1,j,k,m) = zero
            fy(l1,j,k,m) = zero
            fz(1,j,k,m) = zero
            fz(l1,j,k,m) = zero
         endif
   80    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 90 l = 1, nzz
            fx(l,1,k,m) = zero
            fy(l,1,k,m) = zero
            fz(l,1,k,m) = zero
   90       continue
         endif
      endif
  100 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
      if (n2.eq.0) then
         do 120 j = 1, kxyp2
         do 110 l = 1, nzz
         fx(l,j,1,m) = zero
         fy(l,j,1,m) = zero
         fz(l,j,1,m) = zero
  110    continue
  120    continue
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 140 j = 1, kxyp2
         do 130 l = 1, nzz
         fx(l,j,k1,m) = zero
         fy(l,j,k1,m) = zero
         fz(l,j,k1,m) = zero
  130    continue
  140    continue
      endif
  150 continue
  160 continue
  170 continue
      we = 2.0*float(nx*ny*nz)*wp
      return
c calculate potential and sum field energy
  180 if (isign.gt.1) go to 300
c calculate force/charge and sum field energy
      wp = 0.0d0
      if (kstrt.gt.(kxb*kyzb)) go to 290
      do 280 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 270 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 220 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         do 200 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 190 l = 2, nz
            l1 = nz2 - l
            at2 = real(ffd(l,j,k,m))
            at1 = at2*aimag(ffd(l,j,k,m))
            at2 = at2*aimag(q(l,j,k,m))
            fx(l,j,k,m) = cmplx(0.,at2)
            fx(l1,j,k,m) = cmplx(0.,-at2)
            wp = wp + at1*aimag(q(l,j,k,m))**2
  190       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            fx(1,j,k,m) = zero
            fx(l1,j,k,m) = zero
         endif
  200    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 210 l = 1, nzz
            fx(l,1,k,m) = zero
  210       continue
         endif
      endif
  220 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
      if (n2.eq.0) then
         do 240 j = 1, kxyp2
         do 230 l = 1, nzz
         fx(l,j,1,m) = zero
  230    continue
  240    continue
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 260 j = 1, kxyp2
         do 250 l = 1, nzz
         fx(l,j,k1,m) = zero
  250    continue
  260    continue
      endif
  270 continue
  280 continue
  290 continue
      we = 2.0*float(nx*ny*nz)*wp
      return
c calculate smoothing
  300 if (kstrt.gt.(kxb*kyzb)) go to 410
      do 400 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 390 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 340 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         do 320 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 310 l = 2, nz
            l1 = nz2 - l
            at1 = aimag(ffd(l,j,k,m))*aimag(q(l,j,k,m))
            fy(l,j,k,m) = cmplx(0.,at1)
            fy(l1,j,k,m) = cmplx(0.,-at1)
  310       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            fy(1,j,k,m) = zero
            fy(l1,j,k,m) = zero
         endif
  320    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 330 l = 1, nzz
            fy(l,1,k,m) = zero
  330       continue
         endif
      endif
  340 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
      if (n2.eq.0) then
         do 360 j = 1, kxyp2
         do 350 l = 1, nzz
         fy(l,j,1,m) = zero
  350    continue
  360    continue
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 380 j = 1, kxyp2
         do 370 l = 1, nzz
         fy(l,j,k1,m) = zero
  370    continue
  380    continue
      endif
  390 continue
  400 continue
  410 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISD32(q,fx,fy,fz,isign,ffd,ax,ay,az,affp,we,nx,ny,nz
     1,kstrt,nzv,kxyp2,kyzp2,j2blok,m2blok,nzd)
c this subroutine solves 3d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides smoothing a smoothing function,
c with dirichlet boundary conditions (zero potential)
c for distributed data with 2D spatial decomposition
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate sin-cos, or cos-sin transform
c for isign = 0, input: isign,ax,ay,az,affp,nx,ny,nz,kstrt,kxyp2,kyzp2,
c                       j2blok,m2blok,nzd
c output: ffd
c for isign = -1, input: q,ffd,isign,nx,ny,nz,kstrt,nzv,kxyp2,kyzp2,
c                        j2blok,m2blok,nzd
c output: fx,fy,fz,we
c approximate flop count is: 10*nxc*nyc*nzc
c for isign = 1, input: q,ffd,isign,nx,ny,nz,kstrt,nzv,kxyp2,kyzp2,
c                       j2blok,m2blok,nzd
c output: fx,we
c approximate flop count is: 5*nxc*nyc*nzc
c for isign = 2, input: q,ffd,isign,nx,ny,nz,kstrt,nzv,kxyp2,kyzp2,
c                       j2blok,m2blok,nzd
c output: fy
c approximate flop count is: 1*nxc*nyc*nzc
c where nxc = (nx-1)/nvpy, nyc = (ny-1)/nvpz, nzc = nz - 1, and
c nvpy/nvpz = number of procs in y/z
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky,kz) = -sqrt(-1)*kx*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c fy(kx,ky,kz) = -sqrt(-1)*ky*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c fz(kx,ky,kz) = -sqrt(-1)*kz*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2)
c if isign = 1, potential is calculated using the equation:
c fx(kx,ky,kz) = g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz)
c if isign = 2, smoothing is calculated using the equation:
c fy(kx,ky,kz) = q(kx,ky,kz)*s(kx,ky,kz)
c q(l,j,k,m) = complex charge density for fourier mode jj-1,kk-1,l-1
c fx(l,j,k,m) = x component of force/charge
c fy(l,j,k,m) = y component of force/charge
c fz(l,j,k,m) = z component of force/charge
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp2*(mx - 1) and
c kk = k + kyzp2*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c j2blok/m2blok = number of field partitions in x/y
c kxyp2/kyzp2 = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c if isign = 0, form factor array is prepared
c aimag(ffd(l,j,k,m)) = finite-size particle shape factor s
c real(ffd(l,j,k,m)) = potential green's function g
c for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp2*(mx - 1) and
c kk = k + kyzp2*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c ax/ay/az = half-width of particle in x/y/z direction
c affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*
c    |q(kx,ky,kz)*s(kx,ky,kz)|**2)
c nx/ny/nz = system length in x/y/z direction
c nzv = first dimension of field arrays, must be >= nz
c nzd = first dimension of form factor array, must be >= nz
      double precision wp
      complex ffd
      dimension q(nzv,kxyp2,kyzp2,j2blok*m2blok)
      dimension fx(nzv,kxyp2,kyzp2,j2blok*m2blok)
      dimension fy(nzv,kxyp2,kyzp2,j2blok*m2blok)
      dimension fz(nzv,kxyp2,kyzp2,j2blok*m2blok)
      dimension ffd(nzd,kxyp2,kyzp2,j2blok*m2blok)
      nyy = ny + ny
      kxb = nx/kxyp2
      kyzb = nyy/kyzp2
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(nyy)
      dnz = 6.28318530717959/float(nz + nz)
      if (isign.ne.0) go to 60
      if (kstrt.gt.(kxb*kyzb)) return
c prepare form factor array
      do 50 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 40 mx = 1, j2blok
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 30 k = 1, kyzp2
      k1 = k + koff
      if (k1.gt.ny) k1 = k1 - nyy
      dky = dny*float(k1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 20 j = 1, kxyp2
      dkx = dnx*float(j + joff)
      at3 = dkx*dkx + at1
      at4 = (dkx*ax)**2 + at2
      do 10 l = 1, nz
      dkz = dnz*float(l - 1)
      at5 = dkz*dkz + at3
      at6 = exp(-.5*((dkz*az)**2 + at4))
      if (at5.eq.0.) then
         ffd(l,j,k,m) = cmplx(affp,1.)
      else
         ffd(l,j,k,m) = cmplx(affp*at6/at5,at6)
      endif
   10 continue
   20 continue
   30 continue
   40 continue
   50 continue
      return
   60 if (isign.gt.0) go to 180
c calculate force/charge and sum field energy
      wp = 0.0d0
      if (kstrt.gt.(kxb*kyzb)) go to 170
      do 160 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 150 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 100 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         if (k1.gt.ny) k1 = k1 - nyy
         dky = dny*float(k1)
         do 80 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 70 l = 2, nz
            at1 = real(ffd(l,j,k,m))*aimag(ffd(l,j,k,m))
            at4 = q(l,j,k,m)*at1
            at2 = dkx*at4
            at3 = dky*at4
            at4 = dnz*float(l - 1)*at4
            fx(l,j,k,m) = at2
            fy(l,j,k,m) = at3
            fz(l,j,k,m) = at4
            wp = wp + at1*q(l,j,k,m)**2
   70       continue
c mode numbers kz = 0, nz
            fx(1,j,k,m) = 0.
            fy(1,j,k,m) = 0.
            fz(1,j,k,m) = 0.
         endif
   80    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 90 l = 1, nz
            fx(l,1,k,m) = 0.
            fy(l,1,k,m) = 0.
            fz(l,1,k,m) = 0.
   90       continue
         endif
      endif
  100 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
      if (n2.eq.0) then
         do 120 j = 1, kxyp2
         do 110 l = 1, nz
         fx(l,j,1,m) = 0.
         fy(l,j,1,m) = 0.
         fz(l,j,1,m) = 0.
  110    continue
  120    continue
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 140 j = 1, kxyp2
         do 130 l = 1, nz
         fx(l,j,k1,m) = 0.
         fy(l,j,k1,m) = 0.
         fz(l,j,k1,m) = 0.
  130    continue
  140    continue
      endif
  150 continue
  160 continue
  170 continue
      we = 2.0*float(nx*ny*nz)*wp
      return
c calculate potential and sum field energy
  180 if (isign.gt.1) go to 300
c calculate force/charge and sum field energy
      wp = 0.0d0
      if (kstrt.gt.(kxb*kyzb)) go to 290
      do 280 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 270 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 220 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         do 200 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 190 l = 2, nz
            at2 = real(ffd(l,j,k,m))
            at1 = at2*aimag(ffd(l,j,k,m))
            fx(l,j,k,m) = at2*q(l,j,k,m)
            wp = wp + at1*q(l,j,k,m)**2
  190       continue
c mode numbers kz = 0, nz
            fx(1,j,k,m) = 0.
         endif
  200    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 210 l = 1, nz
            fx(l,1,k,m) = 0.
  210       continue
         endif
      endif
  220 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
      if (n2.eq.0) then
         do 240 j = 1, kxyp2
         do 230 l = 1, nz
         fx(l,j,1,m) = 0.
  230    continue
  240    continue
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 260 j = 1, kxyp2
         do 250 l = 1, nz
         fx(l,j,k1,m) = 0.
  250    continue
  260    continue
      endif
  270 continue
  280 continue
  290 continue
      we = 2.0*float(nx*ny*nz)*wp
      return
c calculate smoothing
  300 if (kstrt.gt.(kxb*kyzb)) go to 410
      do 400 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 390 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 340 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         do 320 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 310 l = 2, nz
            fy(l,j,k,m) = aimag(ffd(l,j,k,m))*q(l,j,k,m)
  310       continue
c mode numbers kz = 0, nz
            fy(1,j,k,m) = 0.
         endif
  320    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 330 l = 1, nz
            fy(l,1,k,m) = 0.
  330       continue
         endif
      endif
  340 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
      if (n2.eq.0) then
         do 360 j = 1, kxyp2
         do 350 l = 1, nz
         fy(l,j,1,m) = 0.
  350    continue
  360    continue
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 380 j = 1, kxyp2
         do 370 l = 1, nz
         fy(l,j,k1,m) = 0.
  370    continue
  380    continue
      endif
  390 continue
  400 continue
  410 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISDX332(q,fxyz,isign,ffd,ax,ay,az,affp,we,nx,ny,nz,k
     1strt,nz2d,kxyp2,kyzp2,j2blok,m2blok,nzd)
c this subroutine solves 3d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with dirichlet boundary conditions (zero potential)
c for distributed data with 2D spatial decomposition
c for isign = 0, input: isign,ax,ay,az,affp,nx,ny,nz,kstrt,kxyp2,kyzp2,
c                       j2blok,m2blok,nzd
c output: ffd
c for isign =/ 0, input: q,ffd,isign,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2,
c                        j2blok,m2blok,nzd
c output: fxyz,we
c approximate flop count is: 12*nxc*nyc*nzc
c where nxc = (nx-1)/nvpy, nyc = (ny-1)/nvpz, nzc = nz - 1, and
c nvpy/nvpz = number of procs in y/z
c the equation used is:
c fx(kx,ky,kz) = -sqrt(-1)*kx*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c fy(kx,ky,kz) = -sqrt(-1)*ky*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c fz(kx,ky,kz) = -sqrt(-1)*kz*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2)
c q(l,j,k,m) = complex charge density for fourier mode jj-1,kk-1,l-1
c fxyz(1,l,j,k,m) = x component of force/charge
c fxyz(2,l,j,k,m) = y component of force/charge
c fxyz(3,l,j,k,m) = z component of force/charge
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp2*(mx - 1) and
c kk = k + kyzp2*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c j2blok/m2blok = number of field partitions in x/y
c kxyp2/kyzp2 = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c if isign = 0, form factor array is prepared
c aimag(ffd(l,j,k,m)) = finite-size particle shape factor s
c real(ffd(l,j,k,m)) = potential green's function g
c for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp2*(mx - 1) and
c kk = k + kyzp2*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c ax/ay/az = half-width of particle in x/y/z direction
c affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*
c    |q(kx,ky,kz)*s(kx,ky,kz)|**2)
c nx/ny/nz = system length in x/y/z direction
c nz2d = first dimension of field arrays, must be >= 2*nz
c nzd = first dimension of form factor array, must be >= nz
      double precision wp
      complex q, fxyz, ffd, zero
      dimension q(nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension fxyz(3,nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension ffd(nzd,kxyp2,kyzp2,j2blok*m2blok)
      nyy = ny + ny
      nzz = nz + nz
      nz2 = nzz + 2
      kxb = nx/kxyp2
      kyzb = nyy/kyzp2
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(nyy)
      dnz = 6.28318530717959/float(nzz)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 60
      if (kstrt.gt.(kxb*kyzb)) return
c prepare form factor array
      do 50 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 40 mx = 1, j2blok
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 30 k = 1, kyzp2
      k1 = k + koff
      if (k1.gt.ny) k1 = k1 - nyy
      dky = dny*float(k1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 20 j = 1, kxyp2
      dkx = dnx*float(j + joff)
      at3 = dkx*dkx + at1
      at4 = (dkx*ax)**2 + at2
      do 10 l = 1, nz
      dkz = dnz*float(l - 1)
      at5 = dkz*dkz + at3
      at6 = exp(-.5*((dkz*az)**2 + at4))
      if (at5.eq.0.) then
         ffd(l,j,k,m) = cmplx(0.,1.)
      else
         ffd(l,j,k,m) = cmplx(affp*at6/at5,at6)
      endif
   10 continue
   20 continue
   30 continue
   40 continue
   50 continue
      return
c calculate force/charge and sum field energy
   60 wp = 0.0d0
      if (kstrt.gt.(kxb*kyzb)) go to 200
      do 190 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 180 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 100 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         if (k1.gt.ny) k1 = k1 - nyy
         dky = dny*float(k1)
         do 80 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 70 l = 2, nz
            l1 = nz2 - l
            at1 = real(ffd(l,j,k,m))*aimag(ffd(l,j,k,m))
            at4 = aimag(q(l,j,k,m))*at1
            at2 = dkx*at4
            at3 = dky*at4
            at4 = dnz*float(l - 1)*at4
            fxyz(1,l,j,k,m) = cmplx(at2,0.)
            fxyz(2,l,j,k,m) = cmplx(at3,0.)
            fxyz(3,l,j,k,m) = cmplx(at4,0.)
            fxyz(1,l1,j,k,m) = cmplx(-at2,0.)
            fxyz(2,l1,j,k,m) = cmplx(-at3,0.)
            fxyz(3,l1,j,k,m) = cmplx(at4,0.)
            wp = wp + at1*aimag(q(l,j,k,m))**2
   70       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            fxyz(1,1,j,k,m) = zero
            fxyz(2,1,j,k,m) = zero
            fxyz(3,1,j,k,m) = zero
            fxyz(1,l1,j,k,m) = zero
            fxyz(2,l1,j,k,m) = zero
            fxyz(3,l1,j,k,m) = zero
         endif
   80    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 90 l = 1, nzz
            fxyz(1,l,1,k,m) = zero
            fxyz(2,l,1,k,m) = zero
            fxyz(3,l,1,k,m) = zero
   90       continue
         endif
      endif
  100 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
      if (n2.eq.0) then
         do 120 j = 1, kxyp2
         do 110 l = 1, nzz
         fxyz(1,l,j,1,m) = zero
         fxyz(2,l,j,1,m) = zero
         fxyz(3,l,j,1,m) = zero
  110    continue
  120    continue
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 160 j = 1, kxyp2
         do 150 l = 1, nzz
         fxyz(1,l,j,k1,m) = zero
         fxyz(2,l,j,k1,m) = zero
         fxyz(3,l,j,k1,m) = zero
  150    continue
  160    continue
      endif
  180 continue
  190 continue
  200 continue
      we = 2.0*float(nx*ny*nz)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISD332(q,fxyz,isign,ffd,ax,ay,az,affp,we,nx,ny,nz,ks
     1trt,nzv,kxyp2,kyzp2,j2blok,m2blok,nzd)
c this subroutine solves 3d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with dirichlet boundary conditions (zero potential)
c for distributed data with 2D spatial decomposition
c for isign = 0, input: isign,ax,ay,az,affp,nx,ny,nz,kstrt,kxyp2,kyzp2,
c                       j2blok,m2blok,nzd
c output: ffd
c for isign =/ 0, input: q,ffd,isign,nx,ny,nz,kstrt,nzv,kxyp2,kyzp2,
c                        j2blok,m2blok,nzd
c output: fxyz,we
c approximate flop count is: 10*nxc*nyc*nzc
c where nxc = (nx-1)/nvpy, nyc = (ny-1)/nvpz, nzc = nz - 1, and
c nvpy/nvpz = number of procs in y/z
c the equation used is:
c fx(kx,ky,kz) = -sqrt(-1)*kx*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c fy(kx,ky,kz) = -sqrt(-1)*ky*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c fz(kx,ky,kz) = -sqrt(-1)*kz*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2)
c q(l,j,k,m) = complex charge density for fourier mode jj-1,kk-1,l-1
c fxyz(1,l,j,k,m) = x component of force/charge
c fxyz(2,l,j,k,m) = y component of force/charge
c fxyz(3,l,j,k,m) = z component of force/charge
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp2*(mx - 1) and
c kk = k + kyzp2*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c j2blok/m2blok = number of field partitions in x/y
c kxyp2/kyzp2 = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c if isign = 0, form factor array is prepared
c aimag(ffd(l,j,k,m)) = finite-size particle shape factor s
c real(ffd(l,j,k,m)) = potential green's function g
c for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp2*(mx - 1) and
c kk = k + kyzp2*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c ax/ay/az = half-width of particle in x/y/z direction
c affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*
c    |q(kx,ky,kz)*s(kx,ky,kz)|**2)
c nx/ny/nz = system length in x/y/z direction
c nzv = first dimension of field arrays, must be >= nz
c nzd = first dimension of form factor array, must be >= nz
      double precision wp
      complex ffd
      dimension q(nzv,kxyp2,kyzp2,j2blok*m2blok)
      dimension fxyz(3,nzv,kxyp2,kyzp2,j2blok*m2blok)
      dimension ffd(nzd,kxyp2,kyzp2,j2blok*m2blok)
      nyy = ny + ny
      kxb = nx/kxyp2
      kyzb = nyy/kyzp2
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(nyy)
      dnz = 6.28318530717959/float(nz + nz)
      if (isign.ne.0) go to 60
      if (kstrt.gt.(kxb*kyzb)) return
c prepare form factor array
      do 50 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 40 mx = 1, j2blok
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 30 k = 1, kyzp2
      k1 = k + koff
      if (k1.gt.ny) k1 = k1 - nyy
      dky = dny*float(k1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 20 j = 1, kxyp2
      dkx = dnx*float(j + joff)
      at3 = dkx*dkx + at1
      at4 = (dkx*ax)**2 + at2
      do 10 l = 1, nz
      dkz = dnz*float(l - 1)
      at5 = dkz*dkz + at3
      at6 = exp(-.5*((dkz*az)**2 + at4))
      if (at5.eq.0.) then
         ffd(l,j,k,m) = cmplx(0.,1.)
      else
         ffd(l,j,k,m) = cmplx(affp*at6/at5,at6)
      endif
   10 continue
   20 continue
   30 continue
   40 continue
   50 continue
      return
c calculate force/charge and sum field energy
   60 wp = 0.0d0
      if (kstrt.gt.(kxb*kyzb)) go to 200
      do 190 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 180 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 100 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         if (k1.gt.ny) k1 = k1 - nyy
         dky = dny*float(k1)
         do 80 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 70 l = 2, nz
            at1 = real(ffd(l,j,k,m))*aimag(ffd(l,j,k,m))
            at4 = q(l,j,k,m)*at1
            at2 = dkx*at4
            at3 = dky*at4
            at4 = dnz*float(l - 1)*at4
            fxyz(1,l,j,k,m) = at2
            fxyz(2,l,j,k,m) = at3
            fxyz(3,l,j,k,m) = at4
            wp = wp + at1*q(l,j,k,m)**2
   70       continue
c mode numbers kz = 0, nz
            fxyz(1,1,j,k,m) = 0.
            fxyz(2,1,j,k,m) = 0.
            fxyz(3,1,j,k,m) = 0.
         endif
   80    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 90 l = 1, nz
            fxyz(1,l,1,k,m) = 0.
            fxyz(2,l,1,k,m) = 0.
            fxyz(3,l,1,k,m) = 0.
   90       continue
         endif
      endif
  100 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
      if (n2.eq.0) then
         do 120 j = 1, kxyp2
         do 110 l = 1, nz
         fxyz(1,l,j,1,m) = 0.
         fxyz(2,l,j,1,m) = 0.
         fxyz(3,l,j,1,m) = 0.
  110    continue
  120    continue
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 160 j = 1, kxyp2
         do 150 l = 1, nz
         fxyz(1,l,j,k1,m) = 0.
         fxyz(2,l,j,k1,m) = 0.
         fxyz(3,l,j,k1,m) = 0.
  150    continue
  160    continue
      endif
  180 continue
  190 continue
  200 continue
      we = 2.0*float(nx*ny*nz)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PCUPERPDX32(cu,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2,j2blok,m
     12blok)
c this subroutine calculates the transverse current in fourier space
c with dirichlet boundary conditions (zero potential).
c for distributed data with 2D spatial decomposition
c input: all output: cu
c approximate flop count is:
c 19*nxc*nyc*nzc and (nx/2)*nyc*nzc divides
c where nxc = (nx-1)/nvpy, nyc = (ny-1)/nvpz, nzc = nz - 1, and
c nvpy/nvpz = number of procs in y/z
c the transverse current is calculated using the equation:
c cux(kx,ky,kz) = cux(kx,ky,kz) - kx*(kx*cux(kx,ky,kz)+ky*cuy(kx,ky,kz)+
c                                 kz*cuz(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c cuy(kx,ky,kz) = cuy(kx,ky,kz) - ky*(kx*cux(kx,ky,kz)+ky*cuy(kx,ky,kz)+
c                                 kz*cuz(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c cuz(kx,ky,kz) = cuz(kx,ky,kz) - kz*(kx*cux(kx,ky,kz)+ky*cuy(kx,ky,kz)+
c                                 kz*cuz(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers
c nx/ny/nz = system length in x/y/z direction
c j2blok/m2blok = number of field partitions in x/y
c kxyp2/kyzp2 = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c nz2d = first dimension of field arrays, must be >= 2*nz
      complex cu, zero
      dimension cu(3,nz2d,kxyp2,kyzp2,j2blok*m2blok)
      nyy = ny + ny
      nzz = nz + nz
      nz2 = nzz + 2
      kxb = nx/kxyp2
      kyzb = nyy/kyzp2
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(nyy)
      dnz = 6.28318530717959/float(nzz)
      zero = cmplx(0.,0.)
c calculate transverse part of current
      if (kstrt.gt.(kxb*kyzb)) return
      do 120 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 110 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 50 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         if (k1.gt.ny) k1 = k1 - nyy
         dky = dny*float(k1)
         dky2 = dky*dky
         do 20 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         dkxy2 = dkx*dkx + dky2
         if ((j+joff).gt.0) then
            do 10 l = 2, nz
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            at1 = 1./(dkz*dkz + dkxy2)
            at3 = real(cu(1,l,j,k,m))
            at4 = real(cu(2,l,j,k,m))
            at5 = real(cu(3,l,j,k,m))
            at2 = at1*(dkx*at3 + dky*at4 + dkz*at5)
            at3 = at3 - dkx*at2
            at4 = at4 - dky*at2
            at5 = at5 - dkz*at2
            cu(1,l,j,k,m) = cmplx(at3,0.)
            cu(2,l,j,k,m) = cmplx(at4,0.)
            cu(3,l,j,k,m) = cmplx(at5,0.)
            cu(1,l1,j,k,m) = cmplx(-at3,0.)
            cu(2,l1,j,k,m) = cmplx(-at4,0.)
            cu(3,l1,j,k,m) = cmplx(at5,0.)
   10       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            cu(1,1,j,k,m) = zero
            cu(2,1,j,k,m) = zero
            cu(1,l1,j,k,m) = zero
            cu(2,l1,j,k,m) = zero
            cu(3,l1,j,k,m) = zero
         endif
   20    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nz
               l1 = nz2 - l
               cu(2,l,1,k,m) = zero
               cu(3,l,1,k,m) = zero
               cu(2,l1,1,k,m) = zero
               cu(3,l1,1,k,m) = zero
   30          continue
c mode numbers kz = 0, nz
               l1 = nz + 1
               cu(1,1,1,k,m) = zero
               cu(2,1,1,k,m) = zero
               cu(1,l1,1,k,m) = zero
               cu(2,l1,1,k,m) = zero
               cu(3,l1,1,k,m) = zero
c throw away kx = nx
            else
               do 40 l = 1, nzz
               cu(1,l,1,k,m) = zero
               cu(2,l,1,k,m) = zero
               cu(3,l,1,k,m) = zero
   40          continue
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 70 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 60 l = 2, nz
            l1 = nz2 - l
            cu(1,l,j,1,m) = zero
            cu(3,l,j,1,m) = zero
            cu(1,l1,j,1,m) = zero
            cu(3,l1,j,1,m) = zero
   60       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            cu(1,1,j,1,m) = zero
            cu(2,1,j,1,m) = zero
            cu(3,1,j,1,m) = zero
            cu(1,l1,j,1,m) = zero
            cu(2,l1,j,1,m) = zero
            cu(3,l1,j,1,m) = zero
         endif
   70    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 80 l = 2, nz
            l1 = nz2 - l
            cu(1,l,1,1,m) = zero
            cu(3,l,1,1,m) = zero
            cu(1,l1,1,1,m) = zero
            cu(2,l1,1,1,m) = zero
            cu(3,l1,1,1,m) = zero
   80       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            cu(1,1,1,1,m) = zero
            cu(2,1,1,1,m) = zero
            cu(3,1,1,1,m) = zero
            cu(1,l1,1,1,m) = zero
            cu(2,l1,1,1,m) = zero
            cu(3,l1,1,1,m) = zero
         endif
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 100 j = 1, kxyp2
         do 90 l = 1, nzz
         cu(1,l,j,k1,m) = zero
         cu(2,l,j,k1,m) = zero
         cu(3,l,j,k1,m) = zero
   90    continue
  100    continue
      endif
  110 continue
  120 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCUPERPD32(cu,nx,ny,nz,kstrt,nzv,kxyp2,kyzp2,j2blok,m2b
     1lok)
c this subroutine calculates the transverse current in fourier space
c with dirichlet boundary conditions (zero potential).
c for distributed data with 2D spatial decomposition
c input: all output: cu
c approximate flop count is:
c 16*nxc*nyc*nzc and (nx/2)*nyc*nzc divides
c where nxc = (nx-1)/nvpy, nyc = (ny-1)/nvpz, nzc = nz - 1, and
c nvpy/nvpz = number of procs in y/z
c the transverse current is calculated using the equation:
c cux(kx,ky,kz) = cux(kx,ky,kz) - kx*(kx*cux(kx,ky,kz)+ky*cuy(kx,ky,kz)+
c                                 kz*cuz(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c cuy(kx,ky,kz) = cuy(kx,ky,kz) - ky*(kx*cux(kx,ky,kz)+ky*cuy(kx,ky,kz)+
c                                 kz*cuz(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c cuz(kx,ky,kz) = cuz(kx,ky,kz) - kz*(kx*cux(kx,ky,kz)+ky*cuy(kx,ky,kz)+
c                                 kz*cuz(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers
c nx/ny/nz = system length in x/y/z direction
c j2blok/m2blok = number of field partitions in x/y
c kxyp2/kyzp2 = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c nzv = first dimension of field arrays, must be >= 2*nz
      dimension cu(3,nzv,kxyp2,kyzp2,j2blok*m2blok)
      nyy = ny + ny
      kxb = nx/kxyp2
      kyzb = nyy/kyzp2
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(nyy)
      dnz = 6.28318530717959/float(nz + nz)
c calculate transverse part of current
      if (kstrt.gt.(kxb*kyzb)) return
      do 110 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 100 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 50 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         if (k1.gt.ny) k1 = k1 - nyy
         dky = dny*float(k1)
         dky2 = dky*dky
         do 20 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         dkxy2 = dkx*dkx + dky2
         if ((j+joff).gt.0) then
            do 10 l = 2, nz
            dkz = dnz*float(l - 1)
            at1 = 1./(dkz*dkz + dkxy2)
            at3 = cu(1,l,j,k,m)
            at4 = cu(2,l,j,k,m)
            at5 = cu(3,l,j,k,m)
            at2 = at1*(dkx*at3 + dky*at4 + dkz*at5)
            cu(1,l,j,k,m) = at3 - dkx*at2
            cu(2,l,j,k,m) = at4 - dky*at2
            cu(3,l,j,k,m) = at5 - dkz*at2
   10       continue
c mode numbers kz = 0, nz
            cu(1,1,j,k,m) = 0.
            cu(2,1,j,k,m) = 0.
         endif
   20    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nz
               cu(2,l,1,k,m) = 0.
               cu(3,l,1,k,m) = 0.
   30          continue
c mode numbers kz = 0, nz
               cu(1,1,1,k,m) = 0.
               cu(2,1,1,k,m) = 0.
c throw away kx = nx
            else
               do 40 l = 1, nz
               cu(1,l,1,k,m) = 0.
               cu(2,l,1,k,m) = 0.
               cu(3,l,1,k,m) = 0.
   40          continue
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 70 j = 1, kxyp2
         do 60 l = 2, nz
         cu(1,l,j,1,m) = 0.
         cu(3,l,j,1,m) = 0.
   60    continue
c mode numbers kz = 0, nz
         cu(1,1,j,1,m) = 0.
         cu(2,1,j,1,m) = 0.
         cu(3,1,j,1,m) = 0.
   70    continue
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 90 j = 1, kxyp2
         do 80 l = 1, nz
         cu(1,l,j,k1,m) = 0.
         cu(2,l,j,k1,m) = 0.
         cu(3,l,j,k1,m) = 0.
   80    continue
   90    continue
      endif
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPOISDX332(cu,bxyz,isign,ffd,ax,ay,az,affp,ci,wm,nx,ny
     1,nz,kstrt,nz2d,kxyp2,kyzp2,j2blok,m2blok,nzd)
c this subroutine solves 3d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with dirichlet boundary conditions (zero potential)
c for distributed data with 2D spatial decomposition
c for isign = 0, output: ffd
c input: isign,ax,ay,az,affp,nx,ny,nz,kstrt,kxyp2,kyzp2,2blok,m2blok,nzd
c for isign = -1, output: bxyz,wm
c input: cu,ffd,isign,ci,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2,j2blok,m2blok,
c        nzd
c approximate flop count is: 
c 24*nxc*nyc*nzc + 11*(nxc*nyc + nxc*nzc + nyc*nzc)
c input: cu,ffd,isign,ci,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2,2blok,m2blok,
c        nzd
c approximate flop count is:
c 14*nxc*nyc*nzc + 7*(nxc*nyc + nxc*nzc + nyc*nzc)
c for isign = 2, output: bxyz
c input: cu,ffd,isign,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2,2blok,m2blok,nzd
c approximate flop count is:
c 5*nxc*nyc*nzc + 1*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx-1)/nvpy, nyc = (ny-1)/nvpz, nzc = nz - 1, and
c nvpy/nvpz = number of procs in y/z
c if isign = 0, form factor array is prepared
c if (isign < 0) the magnetic field is calculated using the equations:
c bx(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
c                (ky*cuz(kx,ky,kz)-kz*cuy(kx,ky,kz))*s(kx,ky,kz),
c by(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
c                (kz*cux(kx,ky,kz)-kx*cuz(kx,ky,kz))*s(kx,ky,kz),
c bz(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
c                (kx*cuy(kx,ky,kz)-ky*cux(kx,ky,kz))*s(kx,ky,kz),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2)
c if isign = 1, vector potential is calculated using the equation:
c bx(kx,ky,kz) = ci*ci*g(kx,ky,kz)*cux(kx,ky,kz)
c by(kx,ky,kz) = ci*ci*g(kx,ky,kz)*cuy(kx,ky,kz)
c bz(kx,ky,kz) = ci*ci*g(kx,ky,kz)*cuz(kx,ky,kz)
c if isign = 2, smoothing is calculated using the equation:
c bx(kx,ky,kz) = cux(kx,ky,kz)*s(kx,ky,kz)
c by(kx,ky,kz) = cuy(kx,ky,kz)*s(kx,ky,kz)
c bz(kx,ky,kz) = cuz(kx,ky,kz)*s(kx,ky,kz)
c cu(l,j,k,m) = complex current density for fourier mode jj-1,kk-1,l-1
c bxyz(1,l,j,k,m) = x component of complex magnetic field
c bxyz(2,l,j,k,m) = y component of complex magnetic field
c bxyz(3,l,j,k,m) = z component of complex magnetic field
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp2*(mx - 1) and
c kk = k + kyzp2*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c aimag(ffd(l,j,k,m)) = finite-size particle shape factor s
c real(ffd(l,j,k,m)) = potential green's function g
c ax/ay/az = half-width of particle in x/y/z direction
c affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny/nz = system length in x/y/z direction
c j2blok/m2blok = number of field partitions in x/y
c kxyp2/kyzp2 = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c nz2d = first dimension of field arrays, must be >= 2*nz
c nzd = first dimension of form factor array, must be >= nz
      double precision wp
      complex cu, bxyz, ffd, zero
      dimension cu(3,nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension bxyz(3,nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension ffd(nzd,kxyp2,kyzp2,j2blok*m2blok)
      nyy = ny + ny
      nzz = nz + nz
      nz2 = nzz + 2
      kxb = nx/kxyp2
      kyzb = nyy/kyzp2
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(nyy)
      dnz = 6.28318530717959/float(nzz)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 60
      if (kstrt.gt.(kxb*kyzb)) return
c prepare form factor array
      do 50 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 40 mx = 1, j2blok
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 30 k = 1, kyzp2
      k1 = k + koff
      if (k1.gt.ny) k1 = k1 - nyy
      dky = dny*float(k1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 20 j = 1, kxyp2
      dkx = dnx*float(j + joff)
      at3 = dkx*dkx + at1
      at4 = (dkx*ax)**2 + at2
      do 10 l = 1, nz
      dkz = dnz*float(l - 1)
      at5 = dkz*dkz + at3
      at6 = exp(-.5*((dkz*az)**2 + at4))
      if (at5.eq.0.) then
         ffd(l,j,k,m) = cmplx(affp,1.)
      else
         ffd(l,j,k,m) = cmplx(affp*at6/at5,at6)
      endif
   10 continue
   20 continue
   30 continue
   40 continue
   50 continue
      return
   60 if (isign.gt.0) go to 200
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.(kxb*kyzb)) go to 190
      do 180 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 170 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 110 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         if (k1.gt.ny) k1 = k1 - nyy
         dky = dny*float(k1)
         do 80 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 70 l = 2, nz
            l1 = nz2 - l
            at1 = ci2*real(ffd(l,j,k,m))*aimag(ffd(l,j,k,m))
            at2 = dkx*at1
            at3 = dky*at1
            at4 = dnz*float(l - 1)*at1
            at5 = real(cu(3,l,j,k,m))
            at6 = real(cu(2,l,j,k,m))
            at7 = real(cu(1,l,j,k,m))
            at8 = at3*at5 - at4*at6
            at5 = at4*at7 - at2*at5
            at2 = at2*at6 - at3*at7
            bxyz(1,l,j,k,m) = cmplx(0.,at8)
            bxyz(2,l,j,k,m) = cmplx(0.,at5)
            bxyz(3,l,j,k,m) = cmplx(0.,at2)
            bxyz(1,l1,j,k,m) = cmplx(0.,at8)
            bxyz(2,l1,j,k,m) = cmplx(0.,at5)
            bxyz(3,l1,j,k,m) = cmplx(0.,-at2)
            wp = wp + at1*(real(cu(1,l,j,k,m))**2 + real(cu(2,l,j,k,m))*
     1*2 + real(cu(3,l,j,k,m))**2)
   70       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            at1 = ci2*real(ffd(1,j,k,m))*aimag(ffd(1,j,k,m))
            at2 = dkx*at1
            at3 = dky*at1
            at1 = .5*at1
            at5 = real(cu(3,1,j,k,m))
            at8 = at3*at5
            at5 = -at2*at5
            bxyz(1,1,j,k,m) = cmplx(0.,at8)
            bxyz(2,1,j,k,m) = cmplx(0.,at5)
            bxyz(3,1,j,k,m) = zero
            bxyz(1,l1,j,k,m) = zero
            bxyz(2,l1,j,k,m) = zero
            bxyz(3,l1,j,k,m) = zero
            wp = wp + at1*real(cu(3,1,j,k,m))**2
         endif
   80    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 90 l = 2, nz
               l1 = nz2 - l
               at1 = ci2*real(ffd(l,1,k,m))*aimag(ffd(l,1,k,m))
               at3 = dky*at1
               at4 = dnz*float(l - 1)*at1
               at7 = real(cu(1,l,1,k,m))
               at5 = at4*at7
               at2 = -at3*at7
               bxyz(1,l,1,k,m) = zero
               bxyz(2,l,1,k,m) = cmplx(0.,at5)
               bxyz(3,l,1,k,m) = cmplx(0.,at2)
               bxyz(1,l1,1,k,m) = zero
               bxyz(2,l1,1,k,m) = cmplx(0.,at5)
               bxyz(3,l1,1,k,m) = cmplx(0.,-at2)
               wp = wp + at1*real(cu(1,l,1,k,m))**2
   90          continue
c mode numbers kz = 0, nz
               l1 = nz + 1
               bxyz(1,1,1,k,m) = zero
               bxyz(2,1,1,k,m) = zero
               bxyz(3,1,1,k,m) = zero
               bxyz(1,l1,1,k,m) = zero
               bxyz(2,l1,1,k,m) = zero
               bxyz(3,l1,1,k,m) = zero
c throw away kx = nx
            else
               do 100 l = 1, nzz
               bxyz(1,l,1,k,m) = zero
               bxyz(2,l,1,k,m) = zero
               bxyz(3,l,1,k,m) = zero
  100          continue
            endif
         endif
      endif
  110 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 130 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 120 l = 2, nz
            l1 = nz2 - l
            at1 = ci2*real(ffd(l,j,1,m))*aimag(ffd(l,j,1,m))
            at2 = dkx*at1
            at4 = dnz*float(l - 1)*at1
            at6 = real(cu(2,l,j,1,m))
            at8 = -at4*at6
            at2 = at2*at6
            bxyz(1,l,j,1,m) = cmplx(0.,at8)
            bxyz(2,l,j,1,m) = zero
            bxyz(3,l,j,1,m) = cmplx(0.,at2)
            bxyz(1,l1,j,1,m) = cmplx(0.,at8)
            bxyz(2,l1,j,1,m) = zero
            bxyz(3,l1,j,1,m) = cmplx(0.,-at2)
            wp = wp + at1*real(cu(2,l,j,1,m))**2
  120       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            bxyz(1,1,j,1,m) = zero
            bxyz(2,1,j,1,m) = zero
            bxyz(3,1,j,1,m) = zero
            bxyz(1,l1,j,1,m) = zero
            bxyz(2,l1,j,1,m) = zero
            bxyz(3,l1,j,1,m) = zero
         endif
  130    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 140 l = 1, nzz
            bxyz(1,l,1,1,m) = zero
            bxyz(2,l,1,1,m) = zero
            bxyz(3,l,1,1,m) = zero
  140       continue
         endif
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 160 j = 1, kxyp2
         do 150 l = 1, nzz
         bxyz(1,l,j,k1,m) = zero
         bxyz(2,l,j,k1,m) = zero
         bxyz(3,l,j,k1,m) = zero
  150    continue
  160    continue
      endif
  170 continue
  180 continue
  190 continue
      wm = 2.0*float(nx*ny*nz)*wp
      return
c calculate vector potential and sum field energy
  200 if (isign.gt.1) go to 340
      wp = 0.0d0
      if (kstrt.gt.(kxb*kyzb)) go to 330
      do 320 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 310 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 250 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         if (k1.gt.ny) k1 = k1 - nyy
         do 220 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 210 l = 2, nz
            l1 = nz2 - l
            at2 = ci2*real(ffd(l,j,k,m))
            at1 = at2*aimag(ffd(l,j,k,m))
            at5 = at2*real(cu(1,l,j,k,m))
            at6 = at2*real(cu(2,l,j,k,m))
            at7 = at2*real(cu(3,l,j,k,m))
            bxyz(1,l,j,k,m) = cmplx(at5,0.)
            bxyz(2,l,j,k,m) = cmplx(at6,0.)
            bxyz(3,l,j,k,m) = cmplx(at7,0.)
            bxyz(1,l1,j,k,m) = cmplx(-at5,0.)
            bxyz(2,l1,j,k,m) = cmplx(-at6,0.)
            bxyz(3,l1,j,k,m) = cmplx(at7,0.)
            wp = wp + at1*(real(cu(1,l,j,k,m))**2 + real(cu(2,l,j,k,m))*
     1*2 + real(cu(3,l,j,k,m))**2)
  210       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            at2 = ci2*real(ffd(1,j,k,m))
            at1 = .5*at2*aimag(ffd(1,j,k,m))
            at7 = at2*real(cu(3,1,j,k,m))
            bxyz(1,1,j,k,m) = zero
            bxyz(2,1,j,k,m) = zero
            bxyz(3,1,j,k,m) = cmplx(at7,0.)
            bxyz(1,l1,j,k,m) = zero
            bxyz(2,l1,j,k,m) = zero
            bxyz(3,l1,j,k,m) = zero
            wp = wp + at1*real(cu(3,1,j,k,m))**2
         endif
  220    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 230 l = 2, nz
               l1 = nz2 - l
               at2 = ci2*real(ffd(l,1,k,m))
               at1 = at2*aimag(ffd(l,1,k,m))
               at5 = at2*real(cu(1,l,1,k,m))
               bxyz(1,l,1,k,m) = cmplx(at5,0.)
               bxyz(2,l,1,k,m) = zero
               bxyz(3,l,1,k,m) = zero
               bxyz(1,l1,1,k,m) = cmplx(-at5,0.)
               bxyz(2,l1,1,k,m) = zero
               bxyz(3,l1,1,k,m) = zero
               wp = wp + at1*real(cu(1,l,1,k,m))**2
  230          continue
c mode numbers kz = 0, nz
               l1 = nz + 1
               bxyz(1,1,1,k,m) = zero
               bxyz(2,1,1,k,m) = zero
               bxyz(3,1,1,k,m) = zero
               bxyz(1,l1,1,k,m) = zero
               bxyz(2,l1,1,k,m) = zero
               bxyz(3,l1,1,k,m) = zero
c throw away kx = nx/2
            else
               do 240 l = 1, nzz
               bxyz(1,l,1,k,m) = zero
               bxyz(2,l,1,k,m) = zero
               bxyz(3,l,1,k,m) = zero
  240          continue
            endif
         endif
      endif
  250 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 270 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 260 l = 2, nz
            l1 = nz2 - l
            at2 = ci2*real(ffd(l,j,1,m))
            at1 = at2*aimag(ffd(l,j,1,m))
            at6 = at2*real(cu(2,l,j,1,m))
            bxyz(1,l,j,1,m) = zero
            bxyz(2,l,j,1,m) = cmplx(at6,0.)
            bxyz(3,l,j,1,m) = zero
            bxyz(1,l1,j,1,m) = zero
            bxyz(2,l1,j,1,m) = cmplx(-at6,0.)
            bxyz(3,l1,j,1,m) = zero
            wp = wp + at1*real(cu(2,l,j,1,m))**2
  260       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            bxyz(1,1,j,1,m) = zero
            bxyz(2,1,j,1,m) = zero
            bxyz(3,1,j,1,m) = zero
            bxyz(1,l1,j,1,m) = zero
            bxyz(2,l1,j,1,m) = zero
            bxyz(3,l1,j,1,m) = zero
         endif
  270    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 280 l = 1, nzz
            bxyz(1,l,1,1,m) = zero
            bxyz(2,l,1,1,m) = zero
            bxyz(3,l,1,1,m) = zero
  280       continue
         endif
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 300 j = 1, kxyp2
         do 290 l = 1, nzz
         bxyz(1,l,j,k1,m) = zero
         bxyz(2,l,j,k1,m) = zero
         bxyz(3,l,j,k1,m) = zero
  290    continue
  300    continue
      endif
  310 continue
  320 continue
  330 continue
      wm = 2.0*float(nx*ny*nz)*wp
      return
c calculate smoothing
  340 if (kstrt.gt.(kxb*kyzb)) go to 470
      do 460 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 450 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 390 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         if (k1.gt.ny) k1 = k1 - nyy
         do 360 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 350 l = 2, nz
            l1 = nz2 - l
            at1 = aimag(ffd(l,j,k,m))
            at5 = at1*real(cu(1,l,j,k,m))
            at6 = at1*real(cu(2,l,j,k,m))
            at7 = at1*real(cu(3,l,j,k,m))
            bxyz(1,l,j,k,m) = cmplx(at5,0.)
            bxyz(2,l,j,k,m) = cmplx(at6,0.)
            bxyz(3,l,j,k,m) = cmplx(at7,0.)
            bxyz(1,l1,j,k,m) = cmplx(-at5,0.)
            bxyz(2,l1,j,k,m) = cmplx(-at6,0.)
            bxyz(3,l1,j,k,m) = cmplx(at7,0.)
  350       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            at1 = aimag(ffd(1,j,k,m))*real(cu(3,1,j,k,m))
            bxyz(1,1,j,k,m) = zero
            bxyz(2,1,j,k,m) = zero
            bxyz(3,1,j,k,m) = cmplx(at1,0.)
            bxyz(1,l1,j,k,m) = zero
            bxyz(2,l1,j,k,m) = zero
            bxyz(3,l1,j,k,m) = zero
         endif
  360    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 370 l = 2, nz
               l1 = nz2 - l
               at1 = aimag(ffd(l,1,k,m))*real(cu(1,l,1,k,m))
               bxyz(1,l,1,k,m) = cmplx(at1,0.)
               bxyz(2,l,1,k,m) = zero
               bxyz(3,l,1,k,m) = zero
               bxyz(1,l1,1,k,m) = cmplx(-at1,0.)
               bxyz(2,l1,1,k,m) = zero
               bxyz(3,l1,1,k,m) = zero
  370          continue
c mode numbers kz = 0, nz
               l1 = nz + 1
               bxyz(1,1,1,k,m) = zero
               bxyz(2,1,1,k,m) = zero
               bxyz(3,1,1,k,m) = zero
               bxyz(1,l1,1,k,m) = zero
               bxyz(2,l1,1,k,m) = zero
               bxyz(3,l1,1,k,m) = zero
c throw away kx = nx
            else
               do 380 l = 1, nzz
               bxyz(1,l,1,k,m) = zero
               bxyz(2,l,1,k,m) = zero
               bxyz(3,l,1,k,m) = zero
  380          continue
            endif
         endif
      endif
  390 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 410 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 400 l = 2, nz
            l1 = nz2 - l
            at1 = aimag(ffd(l,j,1,m))*real(cu(2,l,j,1,m))
            bxyz(1,l,j,1,m) = zero
            bxyz(2,l,j,1,m) = cmplx(at1,0.)
            bxyz(3,l,j,1,m) = zero
            bxyz(1,l1,j,1,m) = zero
            bxyz(2,l1,j,1,m) = cmplx(-at1,0.)
            bxyz(3,l1,j,1,m) = zero
  400       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            bxyz(1,1,j,1,m) = zero
            bxyz(2,1,j,1,m) = zero
            bxyz(3,1,j,1,m) = zero
            bxyz(1,l1,j,1,m) = zero
            bxyz(2,l1,j,1,m) = zero
            bxyz(3,l1,j,1,m) = zero
         endif
  410    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 420 l = 1, nzz
            bxyz(1,l,1,1,m) = zero
            bxyz(2,l,1,1,m) = zero
            bxyz(3,l,1,1,m) = zero
  420       continue
         endif
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 440 j = 1, kxyp2
         do 430 l = 1, nzz
         bxyz(1,l,j,k1,m) = zero
         bxyz(2,l,j,k1,m) = zero
         bxyz(3,l,j,k1,m) = zero
  430    continue
  440    continue
      endif
  450 continue
  460 continue
  470 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPOISD332(cu,bxyz,isign,ffd,ax,ay,az,affp,ci,wm,nx,ny,
     1nz,kstrt,nzv,kxyp2,kyzp2,j2blok,m2blok,nzd)
c this subroutine solves 3d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with dirichlet boundary conditions (zero potential)
c for distributed data with 2D spatial decomposition
c for isign = 0, output: ffd
c input: isign,ax,ay,az,affp,nx,ny,nz,kstrt,kxyp2,kyzp2,2blok,m2blok,nzd
c for isign = -1, output: bxyz,wm
c input: cu,ffd,isign,ci,nx,ny,nz,kstrt,nzv,kxyp2,kyzp2,j2blok,m2blok,
c        nzd
c approximate flop count is: 
c 23*nxc*nyc*nzc + 11*(nxc*nyc + nxc*nzc + nyc*nzc)
c input: cu,ffd,isign,ci,nx,ny,nz,kstrt,nzv,kxyp2,kyzp2,2blok,m2blok,
c        nzd
c approximate flop count is:
c 12*nxc*nyc*nzc + 7*(nxc*nyc + nxc*nzc + nyc*nzc)
c for isign = 2, output: bxyz
c input: cu,ffd,isign,nx,ny,nz,kstrt,nzv,kxyp2,kyzp2,2blok,m2blok,nzd
c approximate flop count is:
c 3*nxc*nyc*nzc + 1*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx-1)/nvpy, nyc = (ny-1)/nvpz, nzc = nz - 1, and
c nvpy/nvpz = number of procs in y/z
c if isign = 0, form factor array is prepared
c if (isign < 0) the magnetic field is calculated using the equations:
c bx(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
c                (ky*cuz(kx,ky,kz)-kz*cuy(kx,ky,kz))*s(kx,ky,kz),
c by(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
c                (kz*cux(kx,ky,kz)-kx*cuz(kx,ky,kz))*s(kx,ky,kz),
c bz(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
c                (kx*cuy(kx,ky,kz)-ky*cux(kx,ky,kz))*s(kx,ky,kz),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2)
c if isign = 1, vector potential is calculated using the equation:
c bx(kx,ky,kz) = ci*ci*g(kx,ky,kz)*cux(kx,ky,kz)
c by(kx,ky,kz) = ci*ci*g(kx,ky,kz)*cuy(kx,ky,kz)
c bz(kx,ky,kz) = ci*ci*g(kx,ky,kz)*cuz(kx,ky,kz)
c if isign = 2, smoothing is calculated using the equation:
c bx(kx,ky,kz) = cux(kx,ky,kz)*s(kx,ky,kz)
c by(kx,ky,kz) = cuy(kx,ky,kz)*s(kx,ky,kz)
c bz(kx,ky,kz) = cuz(kx,ky,kz)*s(kx,ky,kz)
c cu(l,j,k,m) = complex current density for fourier mode jj-1,kk-1,l-1
c bxyz(1,l,j,k,m) = x component of complex magnetic field
c bxyz(2,l,j,k,m) = y component of complex magnetic field
c bxyz(3,l,j,k,m) = z component of complex magnetic field
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp2*(mx - 1) and
c kk = k + kyzp2*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c aimag(ffd(l,j,k,m)) = finite-size particle shape factor s
c real(ffd(l,j,k,m)) = potential green's function g
c ax/ay/az = half-width of particle in x/y/z direction
c affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny/nz = system length in x/y/z direction
c j2blok/m2blok = number of field partitions in x/y
c kxyp2/kyzp2 = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c nzv = first dimension of field arrays, must be >= nz
c nzd = first dimension of form factor array, must be >= nz
      double precision wp
      complex ffd
      dimension cu(3,nzv,kxyp2,kyzp2,j2blok*m2blok)
      dimension bxyz(3,nzv,kxyp2,kyzp2,j2blok*m2blok)
      dimension ffd(nzd,kxyp2,kyzp2,j2blok*m2blok)
      nyy = ny + ny
      kxb = nx/kxyp2
      kyzb = nyy/kyzp2
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(nyy)
      dnz = 6.28318530717959/float(nz + nz)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 60
      if (kstrt.gt.(kxb*kyzb)) return
c prepare form factor array
      do 50 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 40 mx = 1, j2blok
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 30 k = 1, kyzp2
      k1 = k + koff
      if (k1.gt.ny) k1 = k1 - nyy
      dky = dny*float(k1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 20 j = 1, kxyp2
      dkx = dnx*float(j + joff)
      at3 = dkx*dkx + at1
      at4 = (dkx*ax)**2 + at2
      do 10 l = 1, nz
      dkz = dnz*float(l - 1)
      at5 = dkz*dkz + at3
      at6 = exp(-.5*((dkz*az)**2 + at4))
      if (at5.eq.0.) then
         ffd(l,j,k,m) = cmplx(affp,1.)
      else
         ffd(l,j,k,m) = cmplx(affp*at6/at5,at6)
      endif
   10 continue
   20 continue
   30 continue
   40 continue
   50 continue
      return
   60 if (isign.gt.0) go to 200
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.(kxb*kyzb)) go to 190
      do 180 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 170 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 110 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         if (k1.gt.ny) k1 = k1 - nyy
         dky = dny*float(k1)
         do 80 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 70 l = 2, nz
            at1 = ci2*real(ffd(l,j,k,m))*aimag(ffd(l,j,k,m))
            at2 = dkx*at1
            at3 = dky*at1
            at4 = dnz*float(l - 1)*at1
            at5 = cu(3,l,j,k,m)
            at6 = cu(2,l,j,k,m)
            at7 = cu(1,l,j,k,m)
            bxyz(1,l,j,k,m) = at3*at5 - at4*at6
            bxyz(2,l,j,k,m) = at4*at7 - at2*at5
            bxyz(3,l,j,k,m) = at2*at6 - at3*at7
            wp = wp + at1*(cu(1,l,j,k,m)**2 + cu(2,l,j,k,m)**2 + cu(3,l,
     1j,k,m)**2)
   70       continue
c mode numbers kz = 0, nz
            at1 = ci2*real(ffd(1,j,k,m))*aimag(ffd(1,j,k,m))
            at2 = dkx*at1
            at3 = dky*at1
            at1 = .5*at1
            at5 = cu(3,1,j,k,m)
            bxyz(1,1,j,k,m) = at3*at5
            bxyz(2,1,j,k,m) = -at2*at5
            bxyz(3,1,j,k,m) = 0.
            wp = wp + at1*cu(3,1,j,k,m)**2
         endif
   80    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 90 l = 2, nz
               at1 = ci2*real(ffd(l,1,k,m))*aimag(ffd(l,1,k,m))
               at3 = dky*at1
               at4 = dnz*float(l - 1)*at1
               at7 = cu(1,l,1,k,m)
               bxyz(1,l,1,k,m) = 0.
               bxyz(2,l,1,k,m) = at4*at7
               bxyz(3,l,1,k,m) = -at3*at7
               wp = wp + at1*cu(1,l,1,k,m)**2
   90          continue
c mode numbers kz = 0, nz
               bxyz(1,1,1,k,m) = 0.
               bxyz(2,1,1,k,m) = 0.
               bxyz(3,1,1,k,m) = 0.
c throw away kx = nx
            else
               do 100 l = 1, nz
               bxyz(1,l,1,k,m) = 0.
               bxyz(2,l,1,k,m) = 0.
               bxyz(3,l,1,k,m) = 0.
  100          continue
            endif
         endif
      endif
  110 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 130 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 120 l = 2, nz
            at1 = ci2*real(ffd(l,j,1,m))*aimag(ffd(l,j,1,m))
            at2 = dkx*at1
            at4 = dnz*float(l - 1)*at1
            at6 = cu(2,l,j,1,m)
            bxyz(1,l,j,1,m) = -at4*at6
            bxyz(2,l,j,1,m) = 0.
            bxyz(3,l,j,1,m) = at2*at6
            wp = wp + at1*cu(2,l,j,1,m)**2
  120       continue
c mode numbers kz = 0, nz
            bxyz(1,1,j,1,m) = 0.
            bxyz(2,1,j,1,m) = 0.
            bxyz(3,1,j,1,m) = 0.
         endif
  130    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 140 l = 1, nz
            bxyz(1,l,1,1,m) = 0.
            bxyz(2,l,1,1,m) = 0.
            bxyz(3,l,1,1,m) = 0.
  140       continue
         endif
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 160 j = 1, kxyp2
         do 150 l = 1, nz
         bxyz(1,l,j,k1,m) = 0.
         bxyz(2,l,j,k1,m) = 0.
         bxyz(3,l,j,k1,m) = 0.
  150    continue
  160    continue
      endif
  170 continue
  180 continue
  190 continue
      wm = 2.0*float(nx*ny*nz)*wp
      return
c calculate vector potential and sum field energy
  200 if (isign.gt.1) go to 340
      wp = 0.0d0
      if (kstrt.gt.(kxb*kyzb)) go to 330
      do 320 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 310 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 250 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         if (k1.gt.ny) k1 = k1 - nyy
         do 220 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 210 l = 2, nz
            at2 = ci2*real(ffd(l,j,k,m))
            at1 = at2*aimag(ffd(l,j,k,m))
            bxyz(1,l,j,k,m) = at2*cu(1,l,j,k,m)
            bxyz(2,l,j,k,m) = at2*cu(2,l,j,k,m)
            bxyz(3,l,j,k,m) = at2*cu(3,l,j,k,m)
            wp = wp + at1*(cu(1,l,j,k,m)**2 + cu(2,l,j,k,m)**2 + cu(3,l,
     1j,k,m)**2)
  210       continue
c mode numbers kz = 0, nz
            at2 = ci2*real(ffd(1,j,k,m))
            at1 = .5*at2*aimag(ffd(1,j,k,m))
            bxyz(1,1,j,k,m) = 0.
            bxyz(2,1,j,k,m) = 0.
            bxyz(3,1,j,k,m) = at2*cu(3,1,j,k,m)
            wp = wp + at1*cu(3,1,j,k,m)**2
         endif
  220    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 230 l = 2, nz
               at2 = ci2*real(ffd(l,1,k,m))
               at1 = at2*aimag(ffd(l,1,k,m))
               bxyz(1,l,1,k,m) = at2*cu(1,l,1,k,m)
               bxyz(2,l,1,k,m) = 0.
               bxyz(3,l,1,k,m) = 0.
               wp = wp + at1*cu(1,l,1,k,m)**2
  230          continue
c mode numbers kz = 0, nz
               bxyz(1,1,1,k,m) = 0.
               bxyz(2,1,1,k,m) = 0.
               bxyz(3,1,1,k,m) = 0.
c throw away kx = nx/2
            else
               do 240 l = 1, nz
               bxyz(1,l,1,k,m) = 0.
               bxyz(2,l,1,k,m) = 0.
               bxyz(3,l,1,k,m) = 0.
  240          continue
            endif
         endif
      endif
  250 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 270 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 260 l = 2, nz
            at2 = ci2*real(ffd(l,j,1,m))
            at1 = at2*aimag(ffd(l,j,1,m))
            bxyz(1,l,j,1,m) = 0.
            bxyz(2,l,j,1,m) = at2*cu(2,l,j,1,m)
            bxyz(3,l,j,1,m) = 0.
            wp = wp + at1*cu(2,l,j,1,m)**2
  260       continue
c mode numbers kz = 0, nz
            bxyz(1,1,j,1,m) = 0.
            bxyz(2,1,j,1,m) = 0.
            bxyz(3,1,j,1,m) = 0.
         endif
  270    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 280 l = 1, nz
            bxyz(1,l,1,1,m) = 0.
            bxyz(2,l,1,1,m) = 0.
            bxyz(3,l,1,1,m) = 0.
  280       continue
         endif
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 300 j = 1, kxyp2
         do 290 l = 1, nz
         bxyz(1,l,j,k1,m) = 0.
         bxyz(2,l,j,k1,m) = 0.
         bxyz(3,l,j,k1,m) = 0.
  290    continue
  300    continue
      endif
  310 continue
  320 continue
  330 continue
      wm = 2.0*float(nx*ny*nz)*wp
      return
c calculate smoothing
  340 if (kstrt.gt.(kxb*kyzb)) go to 470
      do 460 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 450 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 390 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         if (k1.gt.ny) k1 = k1 - nyy
         do 360 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 350 l = 2, nz
            at1 = aimag(ffd(l,j,k,m))
            bxyz(1,l,j,k,m) = at1*cu(1,l,j,k,m)
            bxyz(2,l,j,k,m) = at1*cu(2,l,j,k,m)
            bxyz(3,l,j,k,m) = at1*cu(3,l,j,k,m)
  350       continue
c mode numbers kz = 0, nz
            bxyz(1,1,j,k,m) = 0.
            bxyz(2,1,j,k,m) = 0.
            bxyz(3,1,j,k,m) = aimag(ffd(1,j,k,m))*cu(3,1,j,k,m)
         endif
  360    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 370 l = 2, nz
               bxyz(1,l,1,k,m) = aimag(ffd(l,1,k,m))*cu(1,l,1,k,m)
               bxyz(2,l,1,k,m) = 0.
               bxyz(3,l,1,k,m) = 0.
  370          continue
c mode numbers kz = 0, nz
               bxyz(1,1,1,k,m) = 0.
               bxyz(2,1,1,k,m) = 0.
               bxyz(3,1,1,k,m) = 0.
c throw away kx = nx
            else
               do 380 l = 1, nz
               bxyz(1,l,1,k,m) = 0.
               bxyz(2,l,1,k,m) = 0.
               bxyz(3,l,1,k,m) = 0.
  380          continue
            endif
         endif
      endif
  390 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 410 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 400 l = 2, nz
            bxyz(1,l,j,1,m) = 0.
            bxyz(2,l,j,1,m) = aimag(ffd(l,j,1,m))*cu(2,l,j,1,m)
            bxyz(3,l,j,1,m) = 0.
  400       continue
c mode numbers kz = 0, nz
            bxyz(1,1,j,1,m) = 0.
            bxyz(2,1,j,1,m) = 0.
            bxyz(3,1,j,1,m) = 0.
         endif
  410    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 420 l = 1, nz
            bxyz(1,l,1,1,m) = 0.
            bxyz(2,l,1,1,m) = 0.
            bxyz(3,l,1,1,m) = 0.
  420       continue
         endif
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 440 j = 1, kxyp2
         do 430 l = 1, nz
         bxyz(1,l,j,k1,m) = 0.
         bxyz(2,l,j,k1,m) = 0.
         bxyz(3,l,j,k1,m) = 0.
  430    continue
  440    continue
      endif
  450 continue
  460 continue
  470 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine IPBPOISDX332(cu,bxyz,ffd,ci,wm,nx,ny,nz,kstrt,nz2d,kxyp
     12,kyzp2,j2blok,m2blok,nzd)
c this subroutine solves 3d poisson's equation in fourier space for
c magnetic field with dirichlet boundary conditions (zero potential)
c for distributed data with 2D spatial decomposition
c input: cu,ffd,ci,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2,j2blok,m2blok,nzd
c output: bxyz, wm
c approximate flop count is:
c 24*nxc*nyc*nzc + 11*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx-1)/nvpy, nyc = (ny-1)/nvpz, nzc = nz - 1, and
c nvpy/nvpz = number of procs in y/z
c the magnetic field is calculated using the equations:
c bx(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
c                (ky*cuz(kx,ky,kz)-kz*cuy(kx,ky,kz)),
c by(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
c                (kz*cux(kx,ky,kz)-kx*cuz(kx,ky,kz)),
c bz(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
c                (kx*cuy(kx,ky,kz)-ky*cux(kx,ky,kz)),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2)
c cu(l,j,k,m) = complex current density for fourier mode jj-1,kk-1,l-1
c bxyz(i,l,j,k,m) = i component of complex magnetic field
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp2*(mx - 1) and
c kk = k + kyzp2*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c aimag(ffd(l,j,k,m)) = finite-size particle shape factor s
c real(ffd(l,j,k,m)) = potential green's function g
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny/nz = system length in x/y/z direction
c j2blok/m2blok = number of field partitions in x/y
c kxyp2/kyzp2 = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c nz2d = first dimension of field arrays, must be >= 2*nz
c nzd = first dimension of form factor array, must be >= nz
      double precision wp
      complex cu, bxyz, ffd, zero
      dimension cu(3,nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension bxyz(3,nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension ffd(nzd,kxyp2,kyzp2,j2blok*m2blok)
      nyy = ny + ny
      nzz = nz + nz
      nz2 = nzz + 2
      kxb = nx/kxyp2
      kyzb = nyy/kyzp2
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(nyy)
      dnz = 6.28318530717959/float(nzz)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.(kxb*kyzb)) go to 130
      do 120 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 110 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 50 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         if (k1.gt.ny) k1 = k1 - nyy
         dky = dny*float(k1)
         do 20 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 10 l = 2, nz
            l1 = nz2 - l
            at1 = ci2*real(ffd(l,j,k,m))
            at2 = dkx*at1
            at3 = dky*at1
            at4 = dnz*float(l - 1)*at1
            at1 = at1*aimag(ffd(l,j,k,m))
            at5 = real(cu(3,l,j,k,m))
            at6 = real(cu(2,l,j,k,m))
            at7 = real(cu(1,l,j,k,m))
            at8 = at3*at5 - at4*at6
            at5 = at4*at7 - at2*at5
            at2 = at2*at6 - at3*at7
            bxyz(1,l,j,k,m) = cmplx(0.,at8)
            bxyz(2,l,j,k,m) = cmplx(0.,at5)
            bxyz(3,l,j,k,m) = cmplx(0.,at2)
            bxyz(1,l1,j,k,m) = cmplx(0.,at8)
            bxyz(2,l1,j,k,m) = cmplx(0.,at5)
            bxyz(3,l1,j,k,m) = cmplx(0.,-at2)
            wp = wp + at1*(real(cu(1,l,j,k,m))**2 + real(cu(2,l,j,k,m))*
     1*2 + real(cu(3,l,j,k,m))**2)
   10       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            at1 = ci2*real(ffd(1,j,k,m))
            at2 = dkx*at1
            at3 = dky*at1
            at1 = .5*at1*aimag(ffd(1,j,k,m))
            at5 = real(cu(3,1,j,k,m))
            at8 = at3*at5
            at5 = -at2*at5
            bxyz(1,1,j,k,m) = cmplx(0.,at8)
            bxyz(2,1,j,k,m) = cmplx(0.,at5)
            bxyz(3,1,j,k,m) = zero
            bxyz(1,l1,j,k,m) = zero
            bxyz(2,l1,j,k,m) = zero
            bxyz(3,l1,j,k,m) = zero
            wp = wp + at1*real(cu(3,1,j,k,m))**2
         endif
   20    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nz
               l1 = nz2 - l
               at1 = ci2*real(ffd(l,1,k,m))
               at3 = dky*at1
               at4 = dnz*float(l - 1)*at1
               at1 = at1*aimag(ffd(l,1,k,m))
               at7 = real(cu(1,l,1,k,m))
               at5 = at4*at7
               at2 = -at3*at7
               bxyz(1,l,1,k,m) = zero
               bxyz(2,l,1,k,m) = cmplx(0.,at5)
               bxyz(3,l,1,k,m) = cmplx(0.,at2)
               bxyz(1,l1,1,k,m) = zero
               bxyz(2,l1,1,k,m) = cmplx(0.,at5)
               bxyz(3,l1,1,k,m) = cmplx(0.,-at2)
               wp = wp + at1*real(cu(1,l,1,k,m))**2
   30          continue
c mode numbers kz = 0, nz
               l1 = nz + 1
               bxyz(1,1,1,k,m) = zero
               bxyz(2,1,1,k,m) = zero
               bxyz(3,1,1,k,m) = zero
               bxyz(1,l1,1,k,m) = zero
               bxyz(2,l1,1,k,m) = zero
               bxyz(3,l1,1,k,m) = zero
c throw away kx = nx
            else
               do 40 l = 1, nzz
               bxyz(1,l,1,k,m) = zero
               bxyz(2,l,1,k,m) = zero
               bxyz(3,l,1,k,m) = zero
   40          continue
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 70 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 60 l = 2, nz
            l1 = nz2 - l
            at1 = ci2*real(ffd(l,j,1,m))
            at2 = dkx*at1
            at4 = dnz*float(l - 1)*at1
            at1 = at1*aimag(ffd(l,j,1,m))
            at6 = real(cu(2,l,j,1,m))
            at8 = -at4*at6
            at2 = at2*at6
            bxyz(1,l,j,1,m) = cmplx(0.,at8)
            bxyz(2,l,j,1,m) = zero
            bxyz(3,l,j,1,m) = cmplx(0.,at2)
            bxyz(1,l1,j,1,m) = cmplx(0.,at8)
            bxyz(2,l1,j,1,m) = zero
            bxyz(3,l1,j,1,m) = cmplx(0.,-at2)
            wp = wp + at1*real(cu(2,l,j,1,m))**2
   60       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            bxyz(1,1,j,1,m) = zero
            bxyz(2,1,j,1,m) = zero
            bxyz(3,1,j,1,m) = zero
            bxyz(1,l1,j,1,m) = zero
            bxyz(2,l1,j,1,m) = zero
            bxyz(3,l1,j,1,m) = zero
         endif
   70    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 80 l = 1, nzz
            bxyz(1,l,1,1,m) = zero
            bxyz(2,l,1,1,m) = zero
            bxyz(3,l,1,1,m) = zero
   80       continue
         endif
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 100 j = 1, kxyp2
         do 90 l = 1, nzz
         bxyz(1,l,j,k1,m) = zero
         bxyz(2,l,j,k1,m) = zero
         bxyz(3,l,j,k1,m) = zero
   90    continue
  100    continue
      endif
  110 continue
  120 continue
  130 continue
      wm = 2.0*float(nx*ny*nz)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine IPBPOISD332(cu,bxyz,ffd,ci,wm,nx,ny,nz,kstrt,nzv,kxyp2,
     1kyzp2,j2blok,m2blok,nzd)
c this subroutine solves 3d poisson's equation in fourier space for
c magnetic field with dirichlet boundary conditions (zero potential)
c for distributed data with 2D spatial decomposition
c input: cu,ffd,ci,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2,j2blok,m2blok,nzd
c output: bxyz, wm
c approximate flop count is:
c 23*nxc*nyc*nzc + 11*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx-1)/nvpy, nyc = (ny-1)/nvpz, nzc = nz - 1, and
c nvpy/nvpz = number of procs in y/z
c the magnetic field is calculated using the equations:
c bx(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
c                (ky*cuz(kx,ky,kz)-kz*cuy(kx,ky,kz)),
c by(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
c                (kz*cux(kx,ky,kz)-kx*cuz(kx,ky,kz)),
c bz(kx,ky,kz) = ci*ci*sqrt(-1)*g(kx,ky,kz)*
c                (kx*cuy(kx,ky,kz)-ky*cux(kx,ky,kz)),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2)
c cu(l,j,k,m) = complex current density for fourier mode jj-1,kk-1,l-1
c bxyz(i,l,j,k,m) = i component of complex magnetic field
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp2*(mx - 1) and
c kk = k + kyzp2*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c aimag(ffd(l,j,k,m)) = finite-size particle shape factor s
c real(ffd(l,j,k,m)) = potential green's function g
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny/nz = system length in x/y/z direction
c j2blok/m2blok = number of field partitions in x/y
c kxyp2/kyzp2 = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c nzv = first dimension of field arrays, must be >= nz
c nzd = first dimension of form factor array, must be >= nz
      double precision wp
      complex ffd
      dimension cu(3,nzv,kxyp2,kyzp2,j2blok*m2blok)
      dimension bxyz(3,nzv,kxyp2,kyzp2,j2blok*m2blok)
      dimension ffd(nzd,kxyp2,kyzp2,j2blok*m2blok)
      nyy = ny + ny
      kxb = nx/kxyp2
      kyzb = nyy/kyzp2
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(nyy)
      dnz = 6.28318530717959/float(nz + nz)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.(kxb*kyzb)) go to 130
      do 120 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 110 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 50 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         if (k1.gt.ny) k1 = k1 - nyy
         dky = dny*float(k1)
         do 20 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 10 l = 2, nz
            at1 = ci2*real(ffd(l,j,k,m))
            at2 = dkx*at1
            at3 = dky*at1
            at4 = dnz*float(l - 1)*at1
            at1 = at1*aimag(ffd(l,j,k,m))
            at5 = cu(3,l,j,k,m)
            at6 = cu(2,l,j,k,m)
            at7 = cu(1,l,j,k,m)
            bxyz(1,l,j,k,m) = at3*at5 - at4*at6
            bxyz(2,l,j,k,m) = at4*at7 - at2*at5
            bxyz(3,l,j,k,m) = at2*at6 - at3*at7
            wp = wp + at1*(cu(1,l,j,k,m)**2 + cu(2,l,j,k,m)**2 + cu(3,l,
     1j,k,m)**2)
   10       continue
c mode numbers kz = 0, nz
            at1 = ci2*real(ffd(1,j,k,m))
            at2 = dkx*at1
            at3 = dky*at1
            at1 = .5*at1*aimag(ffd(1,j,k,m))
            at5 = cu(3,1,j,k,m)
            bxyz(1,1,j,k,m) = at3*at5
            bxyz(2,1,j,k,m) = -at2*at5
            bxyz(3,1,j,k,m) = 0.
            wp = wp + at1*cu(3,1,j,k,m)**2
         endif
   20    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nz
               at1 = ci2*real(ffd(l,1,k,m))
               at3 = dky*at1
               at4 = dnz*float(l - 1)*at1
               at1 = at1*aimag(ffd(l,1,k,m))
               at7 = cu(1,l,1,k,m)
               bxyz(1,l,1,k,m) = 0.
               bxyz(2,l,1,k,m) = at4*at7
               bxyz(3,l,1,k,m) = -at3*at7
               wp = wp + at1*cu(1,l,1,k,m)**2
   30          continue
c mode numbers kz = 0, nz
               bxyz(1,1,1,k,m) = 0.
               bxyz(2,1,1,k,m) = 0.
               bxyz(3,1,1,k,m) = 0.
c throw away kx = nx
            else
               do 40 l = 1, nz
               bxyz(1,l,1,k,m) = 0.
               bxyz(2,l,1,k,m) = 0.
               bxyz(3,l,1,k,m) = 0.
   40          continue
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 70 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 60 l = 2, nz
            at1 = ci2*real(ffd(l,j,1,m))
            at2 = dkx*at1
            at4 = dnz*float(l - 1)*at1
            at1 = at1*aimag(ffd(l,j,1,m))
            at6 = cu(2,l,j,1,m)
            bxyz(1,l,j,1,m) = -at4*at6
            bxyz(2,l,j,1,m) = 0.
            bxyz(3,l,j,1,m) = at2*at6
            wp = wp + at1*cu(2,l,j,1,m)**2
   60       continue
c mode numbers kz = 0, nz
            bxyz(1,1,j,1,m) = 0.
            bxyz(2,1,j,1,m) = 0.
            bxyz(3,1,j,1,m) = 0.
         endif
   70    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 80 l = 1, nz
            bxyz(1,l,1,1,m) = 0.
            bxyz(2,l,1,1,m) = 0.
            bxyz(3,l,1,1,m) = 0.
   80       continue
         endif
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 100 j = 1, kxyp2
         do 90 l = 1, nz
         bxyz(1,l,j,k1,m) = 0.
         bxyz(2,l,j,k1,m) = 0.
         bxyz(3,l,j,k1,m) = 0.
   90    continue
  100    continue
      endif
  110 continue
  120 continue
  130 continue
      wm = 2.0*float(nx*ny*nz)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PMAXWELDX32(exyz,bxyz,cu,ffd,affp,ci,dt,wf,wm,nx,ny,nz,
     1kstrt,nz2d,kxyp2,kyzp2,j2blok,m2blok,nzd)
c this subroutine solves 3d maxwell's equation in fourier space for
c transverse electric and magnetic fields with dirichlet boundary
c conditions (zero potential) for distributed data
c with 2D spatial decomposition
c input: all, output: wf, wm, exyz, bxyz
c approximate flop count is:
c 70*nxc*nyc*nzc + 30*(nxc*nyc + nxc*nzc + nyc*nzc)
c plus nxc*nyc*nzc divides
c where nxc = (nx-1)/nvp, nyc = ny - 1, nzc = nz - 1, and
c nvp = number of procs
c the magnetic field is first updated half a step using the equations:
c bx(kx,ky,kz) = bx(kx,ky,kz) - .5*dt*sqrt(-1)*
c                (ky*ez(kx,ky,kz)-kz*ey(kx,ky,kz))
c by(kx,ky,kz) = by(kx,ky,kz) - .5*dt*sqrt(-1)*
c               (kz*ex(kx,ky,kz)-kx*ez(kx,ky,kz))
c bz(kx,ky,kz) = bz(kx,ky,kz) - .5*dt*sqrt(-1)*
c               (kx*ey(kx,ky,kz)-ky*ex(kx,ky,kz))
c the electric field is then updated a whole step using the equations:
c ex(kx,ky,kz) = ex(kx,ky,kz) + c2*dt*sqrt(-1)*
c  (ky*bz(kx,ky,kz)-kz*by(kx,ky,kz)) - affp*dt*cux(kx,ky,kz)*s(kx,ky,kz)
c ey(kx,ky,kz) = ey(kx,ky,kz) + c2*dt*sqrt(-1)*
c  (kz*bx(kx,ky,kz)-kx*bz(kx,ky,kz)) - affp*dt*cuy(kx,ky,kz)*s(kx,ky,kz)
c ez(kx,ky,kz) = ez(kx,ky,kz) + c2*dt*sqrt(-1)*
c  (kx*by(kx,ky,kz)-ky*bx(kx,ky,kz)) - affp*dt*cuz(kx,ky,kz)*s(kx,ky,kz)
c the magnetic field is finally updated the remaining half step with
c the new electric field and the previous magnetic field equations.
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, c2 = 1./(ci*ci)
c and s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)
c j,k,l = fourier mode numbers
c and similarly for bx, by, bz.
c exyz(i,l,j,k,m) = i component of complex transverse electric field
c bxyz(i,l,j,k,m) = i component of complex magnetic field
c cu(i,l,j,k,m) = i component of complex current density
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp2*(mx - 1) and
c kk = k + kyzp*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c aimag(ffd(l,j,k,m)) = finite-size particle shape factor s,
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2)
c affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny*nz**sum((1/affp)*|exyz(kx,ky,kz)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny*nz**sum((c2/affp)*|bxyz(kx,ky,kz)|**2)
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c kxyp2/kyzp2 = number of complex grids in each field partition in
c x/y direction
c nz2d = first dimension of field arrays, must be >= 2*nz
c j2blok/m2blok = number of field partitions in x/y
c nzd = first dimension of form factor array, must be >= nz
      double precision wp, ws
      complex exyz, bxyz, cu, ffd, zero
      dimension exyz(3,nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension bxyz(3,nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension cu(3,nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension ffd(nzd,kxyp2,kyzp2,j2blok*m2blok)
      if (ci.le.0.) return
      nyy = ny + ny
      nzz = nz + nz
      nz2 = nzz + 2
      kxb = nx/kxyp2
      kyzb = nyy/kyzp2
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(nyy)
      dnz = 6.28318530717959/float(nzz)
      dth = .5*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      adt = affp*dt
      zero = cmplx(0.,0.)
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
c calculate the electromagnetic fields
      if (kstrt.gt.(kxb*kyzb)) go to 130
      do 120 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 110 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 50 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         if (k1.gt.ny) k1 = k1 - nyy
         dky = dny*float(k1)
         do 20 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 10 l = 2, nz
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            afdt = adt*aimag(ffd(l,j,k,m))
c update magnetic field half time step, ky > 0, kz > 0
            at7 = real(exyz(1,l,j,k,m))
            at8 = real(exyz(2,l,j,k,m))
            at9 = real(exyz(3,l,j,k,m))
            at4 = aimag(bxyz(1,l,j,k,m)) + dth*(dkz*at8 - dky*at9)
            at5 = aimag(bxyz(2,l,j,k,m)) + dth*(dkx*at9 - dkz*at7)
            at6 = aimag(bxyz(3,l,j,k,m)) + dth*(dky*at7 - dkx*at8)
c update electric field whole time step
            at7 = at7 + cdt*(dkz*at5 - dky*at6) - afdt*real(cu(1,l,j,k,m
     1))
            at8 = at8 + cdt*(dkx*at6 - dkz*at4) - afdt*real(cu(2,l,j,k,m
     1))
            at9 = at9 + cdt*(dky*at4 - dkx*at5) - afdt*real(cu(3,l,j,k,m
     1))
c update magnetic field half time step and store electric field
            at4 = at4 + dth*(dkz*at8 - dky*at9)
            at5 = at5 + dth*(dkx*at9 - dkz*at7)
            at6 = at6 + dth*(dky*at7 - dkx*at8)
            ws = ws + anorm*(at7*at7 + at8*at8 + at9*at9)
            wp = wp + anorm*(at4*at4 + at5*at5 + at6*at6)
            exyz(1,l,j,k,m) = cmplx(at7,0.)
            exyz(2,l,j,k,m) = cmplx(at8,0.)
            exyz(3,l,j,k,m) = cmplx(at9,0.)
            bxyz(1,l,j,k,m) = cmplx(0.,at4)
            bxyz(2,l,j,k,m) = cmplx(0.,at5)
            bxyz(3,l,j,k,m) = cmplx(0.,at6)
c update magnetic field half time step, ky > 0, kz < 0
            exyz(1,l1,j,k,m) = cmplx(-at7,0.)
            exyz(2,l1,j,k,m) = cmplx(-at8,0.)
            exyz(3,l1,j,k,m) = cmplx(at9,0.)
            bxyz(1,l1,j,k,m) = cmplx(0.,at4)
            bxyz(2,l1,j,k,m) = cmplx(0.,at5)
            bxyz(3,l1,j,k,m) = cmplx(0.,-at6)
   10       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            afdt = adt*aimag(ffd(1,j,k,m))
c update magnetic field half time step, ky > 0
            at9 = real(exyz(3,1,j,k,m))
            at4 = aimag(bxyz(1,1,j,k,m)) - dth*(dky*at9)
            at5 = aimag(bxyz(2,1,j,k,m)) + dth*(dkx*at9)
c update electric field whole time step
            at9 = at9 + cdt*(dky*at4 - dkx*at5) - afdt*real(cu(3,1,j,k,m
     1))
c update magnetic field half time step and store electric field
            at4 = at4 - dth*(dky*at9)
            at5 = at5 + dth*(dkx*at9)
            ws = ws + .5*anorm*(at9*at9)
            wp = wp + .5*anorm*(at4*at4 + at5*at5)
            exyz(1,1,j,k,m) = zero
            exyz(2,1,j,k,m) = zero
            exyz(3,1,j,k,m) = cmplx(at9,0.)
            bxyz(1,1,j,k,m) = cmplx(0.,at4)
            bxyz(2,1,j,k,m) = cmplx(0.,at5)
            bxyz(3,1,j,k,m) = zero
            bxyz(1,l1,j,k,m) = zero
            bxyz(2,l1,j,k,m) = zero
            bxyz(3,l1,j,k,m) = zero
            exyz(1,l1,j,k,m) = zero
            exyz(2,l1,j,k,m) = zero
            exyz(3,l1,j,k,m) = zero
         endif
   20    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nz
               l1 = nz2 - l
               dkz = dnz*float(l - 1)
               afdt = adt*aimag(ffd(l,1,k,m))
c update magnetic field half time step, kz > 0
               at7 = real(exyz(1,l,1,k,m))
               at5 = aimag(bxyz(2,l,1,k,m)) - dth*(dkz*at7)
               at6 = aimag(bxyz(3,l,1,k,m)) + dth*(dky*at7)
c update electric field whole time step
               at7 = at7 + cdt*(dkz*at5 - dky*at6) - afdt*real(cu(1,l,1,
     1k,m))
c update magnetic field half time step and store electric field
               at5 = at5 - dth*(dkz*at7)
               at6 = at6 + dth*(dky*at7) 
               ws = ws + anorm*(at7*at7)
               wp = wp + anorm*(at5*at5 + at6*at6)
               exyz(1,l,1,k,m) = cmplx(at7,0.)
               exyz(2,l,1,k,m) = zero
               exyz(3,l,1,k,m) = zero
               bxyz(1,l,1,k,m) = zero
               bxyz(2,l,1,k,m) = cmplx(0.,at5)
               bxyz(3,l,1,k,m) = cmplx(0.,at6)
               exyz(1,l1,1,k,m) = cmplx(-at7,0.)
               exyz(2,l1,1,k,m) = zero
               exyz(3,l1,1,k,m) = zero
               bxyz(1,l1,1,k,m) = zero
               bxyz(2,l1,1,k,m) = cmplx(0.,at5)
               bxyz(3,l1,1,k,m) = cmplx(0.,-at6)
   30          continue
c mode numbers kz = 0, nz
               l1 = nz + 1
               exyz(1,1,1,k,m) = zero
               exyz(2,1,1,k,m) = zero
               exyz(3,1,1,k,m) = zero
               bxyz(1,1,1,k,m) = zero
               bxyz(2,1,1,k,m) = zero
               bxyz(3,1,1,k,m) = zero
               bxyz(1,l1,1,k,m) = zero
               bxyz(2,l1,1,k,m) = zero
               bxyz(3,l1,1,k,m) = zero
               exyz(1,l1,1,k,m) = zero
               exyz(2,l1,1,k,m) = zero
               exyz(3,l1,1,k,m) = zero
c throw away kx = nx
            else
               do 40 l = 1, nzz
               bxyz(1,l,1,k,m) = zero
               bxyz(2,l,1,k,m) = zero
               bxyz(3,l,1,k,m) = zero
               exyz(1,l,1,k,m) = zero
               exyz(2,l,1,k,m) = zero
               exyz(3,l,1,k,m) = zero
   40          continue
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 70 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 60 l = 2, nz
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            afdt = adt*aimag(ffd(l,j,1,m))
            at8 = real(exyz(2,l,j,1,m))
            at4 = aimag(bxyz(1,l,j,1,m)) + dth*(dkz*at8)
            at6 = aimag(bxyz(3,l,j,1,m)) - dth*(dkx*at8)
c update electric field whole time step
            at8 = at8 + cdt*(dkx*at6 - dkz*at4) - afdt*real(cu(2,l,j,1,m
     1))
c update magnetic field half time step and store electric field
            at4 = at4 + dth*(dkz*at8)
            at6 = at6 - dth*(dkx*at8)
            ws = ws + anorm*(at8*at8)
            wp = wp + anorm*(at4*at4 + at6*at6)
            exyz(1,l,j,1,m) = zero
            exyz(2,l,j,1,m) = cmplx(at8,0.)
            exyz(3,l,j,1,m) = zero
            bxyz(1,l,j,1,m) = cmplx(0.,at4)
            bxyz(2,l,j,1,m) = zero
            bxyz(3,l,j,1,m) = cmplx(0.,at6)
c update magnetic field half time step, kz < 0
            exyz(1,l1,j,1,m) = zero
            exyz(2,l1,j,1,m) = cmplx(-at8,0.)
            exyz(3,l1,j,1,m) = zero
            bxyz(1,l1,j,1,m) = cmplx(0.,at4)
            bxyz(2,l1,j,1,m) = zero
            bxyz(3,l1,j,1,m) = cmplx(0.,-at6)
   60       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            exyz(1,1,j,1,m) = zero
            exyz(2,1,j,1,m) = zero
            exyz(3,1,j,1,m) = zero
            bxyz(1,1,j,1,m) = zero
            bxyz(2,1,j,1,m) = zero
            bxyz(3,1,j,1,m) = zero
            bxyz(1,l1,j,1,m) = zero
            bxyz(2,l1,j,1,m) = zero
            bxyz(3,l1,j,1,m) = zero
            exyz(1,l1,j,1,m) = zero
            exyz(2,l1,j,1,m) = zero
            exyz(3,l1,j,1,m) = zero
         endif
   70    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 80 l = 1, nzz
            exyz(1,l,1,1,m) = zero
            exyz(2,l,1,1,m) = zero
            exyz(3,l,1,1,m) = zero
            bxyz(1,l,1,1,m) = zero
            bxyz(2,l,1,1,m) = zero
            bxyz(3,l,1,1,m) = zero
   80       continue
         endif
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 100 j = 1, kxyp2
         do 90 l = 1, nzz
         exyz(1,l,j,k1,m) = zero
         exyz(2,l,j,k1,m) = zero
         exyz(3,l,j,k1,m) = zero
         bxyz(1,l,j,k1,m) = zero
         bxyz(2,l,j,k1,m) = zero
         bxyz(3,l,j,k1,m) = zero
   90    continue
  100    continue
      endif
  110 continue
  120 continue
  130 continue      
      wf = 2.0*float(nx*ny*nz)*ws
      wm = 2.0*float(nx*ny*nz)*c2*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PMAXWELD32(exyz,bxyz,cu,ffd,affp,ci,dt,wf,wm,nx,ny,nz,k
     1strt,nzv,kxyp2,kyzp2,j2blok,m2blok,nzd)
c this subroutine solves 3d maxwell's equation in fourier space for
c transverse electric and magnetic fields with dirichlet boundary
c conditions (zero potential) for distributed data
c with 2D spatial decomposition
c input: all, output: wf, wm, exyz, bxyz
c approximate flop count is:
c 62*nxc*nyc*nzc + 29*(nxc*nyc + nxc*nzc + nyc*nzc)
c plus nxc*nyc*nzc divides
c where nxc = (nx-1)/nvp, nyc = ny - 1, nzc = nz - 1, and
c nvp = number of procs
c the magnetic field is first updated half a step using the equations:
c bx(kx,ky,kz) = bx(kx,ky,kz) - .5*dt*sqrt(-1)*
c                (ky*ez(kx,ky,kz)-kz*ey(kx,ky,kz))
c by(kx,ky,kz) = by(kx,ky,kz) - .5*dt*sqrt(-1)*
c               (kz*ex(kx,ky,kz)-kx*ez(kx,ky,kz))
c bz(kx,ky,kz) = bz(kx,ky,kz) - .5*dt*sqrt(-1)*
c               (kx*ey(kx,ky,kz)-ky*ex(kx,ky,kz))
c the electric field is then updated a whole step using the equations:
c ex(kx,ky,kz) = ex(kx,ky,kz) + c2*dt*sqrt(-1)*
c  (ky*bz(kx,ky,kz)-kz*by(kx,ky,kz)) - affp*dt*cux(kx,ky,kz)*s(kx,ky,kz)
c ey(kx,ky,kz) = ey(kx,ky,kz) + c2*dt*sqrt(-1)*
c  (kz*bx(kx,ky,kz)-kx*bz(kx,ky,kz)) - affp*dt*cuy(kx,ky,kz)*s(kx,ky,kz)
c ez(kx,ky,kz) = ez(kx,ky,kz) + c2*dt*sqrt(-1)*
c  (kx*by(kx,ky,kz)-ky*bx(kx,ky,kz)) - affp*dt*cuz(kx,ky,kz)*s(kx,ky,kz)
c the magnetic field is finally updated the remaining half step with
c the new electric field and the previous magnetic field equations.
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, c2 = 1./(ci*ci)
c and s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)
c j,k,l = fourier mode numbers
c and similarly for bx, by, bz.
c exyz(i,l,j,k,m) = i component of complex transverse electric field
c bxyz(i,l,j,k,m) = i component of complex magnetic field
c cu(i,l,j,k,m) = i component of complex current density
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp2*(mx - 1) and
c kk = k + kyzp*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c aimag(ffd(l,j,k,m)) = finite-size particle shape factor s,
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2)
c affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c ci = reciprical of velocity of light
c dt = time interval between successive calculations
c transverse electric field energy is also calculated, using
c wf = nx*ny*nz**sum((1/affp)*|exyz(kx,ky,kz)|**2)
c magnetic field energy is also calculated, using
c wm = nx*ny*nz**sum((c2/affp)*|bxyz(kx,ky,kz)|**2)
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c kxyp2/kyzp2 = number of complex grids in each field partition in
c x/y direction
c nzv = first dimension of field arrays, must be >= nz
c j2blok/m2blok = number of field partitions in x/y
c nzd = first dimension of form factor array, must be >= nz
      double precision wp, ws
      complex ffd
      dimension exyz(3,nzv,kxyp2,kyzp2,j2blok*m2blok)
      dimension bxyz(3,nzv,kxyp2,kyzp2,j2blok*m2blok)
      dimension cu(3,nzv,kxyp2,kyzp2,j2blok*m2blok)
      dimension ffd(nzd,kxyp2,kyzp2,j2blok*m2blok)
      if (ci.le.0.) return
      nyy = ny + ny
      kxb = nx/kxyp2
      kyzb = nyy/kyzp2
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(nyy)
      dnz = 6.28318530717959/float(nz + nz)
      dth = .5*dt
      c2 = 1./(ci*ci)
      cdt = c2*dt
      adt = affp*dt
      anorm = 1.0/affp
c update electromagnetic field and sum field energies
      ws = 0.0d0
      wp = 0.0d0
c calculate the electromagnetic fields
      if (kstrt.gt.(kxb*kyzb)) go to 130
      do 120 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 110 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 50 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         if (k1.gt.ny) k1 = k1 - nyy
         dky = dny*float(k1)
         do 20 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 10 l = 2, nz
            dkz = dnz*float(l - 1)
            afdt = adt*aimag(ffd(l,j,k,m))
c update magnetic field half time step, ky > 0, kz > 0
            at7 = exyz(1,l,j,k,m)
            at8 = exyz(2,l,j,k,m)
            at9 = exyz(3,l,j,k,m)
            at4 = bxyz(1,l,j,k,m) + dth*(dkz*at8 - dky*at9)
            at5 = bxyz(2,l,j,k,m) + dth*(dkx*at9 - dkz*at7)
            at6 = bxyz(3,l,j,k,m) + dth*(dky*at7 - dkx*at8)
c update electric field whole time step
            at7 = at7 + cdt*(dkz*at5 - dky*at6) - afdt*cu(1,l,j,k,m)
            at8 = at8 + cdt*(dkx*at6 - dkz*at4) - afdt*cu(2,l,j,k,m)
            at9 = at9 + cdt*(dky*at4 - dkx*at5) - afdt*cu(3,l,j,k,m)
c update magnetic field half time step and store electric field
            at4 = at4 + dth*(dkz*at8 - dky*at9)
            at5 = at5 + dth*(dkx*at9 - dkz*at7)
            at6 = at6 + dth*(dky*at7 - dkx*at8)
            ws = ws + anorm*(at7*at7 + at8*at8 + at9*at9)
            wp = wp + anorm*(at4*at4 + at5*at5 + at6*at6)
            exyz(1,l,j,k,m) = at7
            exyz(2,l,j,k,m) = at8
            exyz(3,l,j,k,m) = at9
            bxyz(1,l,j,k,m) = at4
            bxyz(2,l,j,k,m) = at5
            bxyz(3,l,j,k,m) = at6
   10       continue
c mode numbers kz = 0, nz
            afdt = adt*aimag(ffd(1,j,k,m))
c update magnetic field half time step, ky > 0
            at9 = exyz(3,1,j,k,m)
            at4 = bxyz(1,1,j,k,m) - dth*(dky*at9)
            at5 = bxyz(2,1,j,k,m) + dth*(dkx*at9)
c update electric field whole time step
            at9 = at9 + cdt*(dky*at4 - dkx*at5) - afdt*cu(3,1,j,k,m)
c update magnetic field half time step and store electric field
            at4 = at4 - dth*(dky*at9)
            at5 = at5 + dth*(dkx*at9)
            ws = ws + .5*anorm*(at9*at9)
            wp = wp + .5*anorm*(at4*at4 + at5*at5)
            exyz(1,1,j,k,m) = 0.
            exyz(2,1,j,k,m) = 0.
            exyz(3,1,j,k,m) = at9
            bxyz(1,1,j,k,m) = at4
            bxyz(2,1,j,k,m) = at5
            bxyz(3,1,j,k,m) = 0.
         endif
   20    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nz
               dkz = dnz*float(l - 1)
               afdt = adt*aimag(ffd(l,1,k,m))
c update magnetic field half time step, kz > 0
               at7 = exyz(1,l,1,k,m)
               at5 = bxyz(2,l,1,k,m) - dth*(dkz*at7)
               at6 = bxyz(3,l,1,k,m) + dth*(dky*at7)
c update electric field whole time step
               at7 = at7 + cdt*(dkz*at5 - dky*at6) - afdt*cu(1,l,1,k,m)
c update magnetic field half time step and store electric field
               at5 = at5 - dth*(dkz*at7)
               at6 = at6 + dth*(dky*at7) 
               ws = ws + anorm*(at7*at7)
               wp = wp + anorm*(at5*at5 + at6*at6)
               exyz(1,l,1,k,m) = at7
               exyz(2,l,1,k,m) = 0.
               exyz(3,l,1,k,m) = 0.
               bxyz(1,l,1,k,m) = 0.
               bxyz(2,l,1,k,m) = at5
               bxyz(3,l,1,k,m) = at6
   30          continue
c mode numbers kz = 0, nz
               exyz(1,1,1,k,m) = 0.
               exyz(2,1,1,k,m) = 0.
               exyz(3,1,1,k,m) = 0.
               bxyz(1,1,1,k,m) = 0.
               bxyz(2,1,1,k,m) = 0.
               bxyz(3,1,1,k,m) = 0.
c throw away kx = nx
            else
               do 40 l = 1, nz
               bxyz(1,l,1,k,m) = 0.
               bxyz(2,l,1,k,m) = 0.
               bxyz(3,l,1,k,m) = 0.
               exyz(1,l,1,k,m) = 0.
               exyz(2,l,1,k,m) = 0.
               exyz(3,l,1,k,m) = 0.
   40          continue
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 70 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 60 l = 2, nz
            dkz = dnz*float(l - 1)
            afdt = adt*aimag(ffd(l,j,1,m))
            at8 = exyz(2,l,j,1,m)
            at4 = bxyz(1,l,j,1,m) + dth*(dkz*at8)
            at6 = bxyz(3,l,j,1,m) - dth*(dkx*at8)
c update electric field whole time step
            at8 = at8 + cdt*(dkx*at6 - dkz*at4) - afdt*cu(2,l,j,1,m)
c update magnetic field half time step and store electric field
            at4 = at4 + dth*(dkz*at8)
            at6 = at6 - dth*(dkx*at8)
            ws = ws + anorm*(at8*at8)
            wp = wp + anorm*(at4*at4 + at6*at6)
            exyz(1,l,j,1,m) = 0.
            exyz(2,l,j,1,m) = at8
            exyz(3,l,j,1,m) = 0.
            bxyz(1,l,j,1,m) = at4
            bxyz(2,l,j,1,m) = 0.
            bxyz(3,l,j,1,m) = at6
   60       continue
c mode numbers kz = 0, nz
            exyz(1,1,j,1,m) = 0.
            exyz(2,1,j,1,m) = 0.
            exyz(3,1,j,1,m) = 0.
            bxyz(1,1,j,1,m) = 0.
            bxyz(2,1,j,1,m) = 0.
            bxyz(3,1,j,1,m) = 0.
         endif
   70    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 80 l = 1, nz
            exyz(1,l,1,1,m) = 0.
            exyz(2,l,1,1,m) = 0.
            exyz(3,l,1,1,m) = 0.
            bxyz(1,l,1,1,m) = 0.
            bxyz(2,l,1,1,m) = 0.
            bxyz(3,l,1,1,m) = 0.
   80       continue
         endif
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 100 j = 1, kxyp2
         do 90 l = 1, nz
         exyz(1,l,j,k1,m) = 0.
         exyz(2,l,j,k1,m) = 0.
         exyz(3,l,j,k1,m) = 0.
         bxyz(1,l,j,k1,m) = 0.
         bxyz(2,l,j,k1,m) = 0.
         bxyz(3,l,j,k1,m) = 0.
   90    continue
  100    continue
      endif
  110 continue
  120 continue
  130 continue      
      wf = 2.0*float(nx*ny*nz)*ws
      wm = 2.0*float(nx*ny*nz)*c2*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PDMFIELDD32(q3,q,nx,ny,nz,kstrt,nz2d,nzv,kxyp2,kyzp2,j2
     1blok,m2blok)
c this subroutine copies the charge density into a smaller array
      implicit none
      integer nx, ny, nz, kstrt, nz2d, nzv, kxyp2, kyzp2, j2blok, m2blok
      complex q3
      real q
      dimension q3(nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension q(nzv,kxyp2,kyzp2,j2blok*m2blok)
      integer j, k, l, m, kxb, kyzb
      kxb = nx/kxyp2
      kyzb = (ny + ny)/kyzp2
      if (kstrt.gt.(kxb*kyzb)) return
      do 40 m = 1, j2blok*m2blok
      do 30 k = 1, kyzp2
      do 20 j = 1, kxyp2
      do 10 l = 1, nz
      q(l,j,k,m) = aimag(q3(l,j,k,m))
   10 continue
   20 continue
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCMFIELDD32(cu3,cu,nx,ny,nz,kstrt,nz2d,nzv,kxyp2,kyzp2,
     1j2blok,m2blok)
c this subroutine copies the current density into a smaller array
      implicit none
      integer nx, ny, nz, kstrt, nz2d, nzv, kxyp2, kyzp2, j2blok, m2blok
      complex cu3
      real cu
      dimension cu3(3,nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension cu(3,nzv,kxyp2,kyzp2,j2blok*m2blok)
      integer j, k, l, m, kxb, kyzb
      kxb = nx/kxyp2
      kyzb = (ny + ny)/kyzp2
      if (kstrt.gt.(kxb*kyzb)) return
      do 40 m = 1, j2blok*m2blok
      do 30 k = 1, kyzp2
      do 20 j = 1, kxyp2
      do 10 l = 1, nz
      cu(1,l,j,k,m) = real(cu3(1,l,j,k,m))
      cu(2,l,j,k,m) = real(cu3(2,l,j,k,m))
      cu(3,l,j,k,m) = real(cu3(3,l,j,k,m))
   10 continue
   20 continue
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PEMFIELDD32(fxyz,exyz,ffd,isign,nx,ny,nz,kstrt,nz2d,nzv
     1,kxyp2,kyzp2,j2blok,m2blok,nzd)
c this subroutine either adds complex vector fields if isign > 0
c or copies complex vector fields if isign <= 0
c adds image charges appropriate for electric field if isign >= 0
c or appropriate for magnetic field if isign < 0
c includes additional smoothing for isign /= 0
      implicit none
      integer isign, nx, ny, nz, kstrt, nz2d, nzv, kxyp2, kyzp2
      integer j2blok, m2blok, nzd
      complex fxyz, ffd
      real exyz
      dimension fxyz(3,nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension exyz(3,nzv,kxyp2,kyzp2,j2blok*m2blok)
      dimension ffd(nzd,kxyp2,kyzp2,j2blok*m2blok)
      complex zero
      integer j, k, l, m, mx, my, nyy, nzz, nz2, kxb, kyzb, js, ks
      integer joff, koff, moff, k1, l1, n1, n2
      real at1, at2, at3, at4
      nyy = ny + ny
      nzz = nz + nz
      nz2 = nzz + 2
      kxb = nx/kxyp2
      kyzb = nyy/kyzp2
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      zero = cmplx(0.,0.)
      if (kstrt.gt.(kxb*kyzb)) return
c add the fields
      if (isign.gt.0) then
         do 90 my = 1, m2blok
         koff = kyzp2*(my + ks) - 1
         moff = j2blok*(my - 1)
         do 80 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
         joff = kxyp2*(mx + js) - 1
         m = mx + moff
         do 40 k = 1, kyzp2
         k1 = k + koff
         if ((k1.gt.0).and.(k1.ne.ny)) then
            do 20 j = 1, kxyp2
            if ((j+joff).gt.0) then
               do 10 l = 2, nz
               l1 = nz2 - l
               at1 = aimag(ffd(l,j,k,m))
               at2 = exyz(1,l,j,k,m)*at1
               at3 = exyz(2,l,j,k,m)*at1
               at4 = exyz(3,l,j,k,m)*at1
               fxyz(1,l,j,k,m) = fxyz(1,l,j,k,m) + cmplx(at2,0.)
               fxyz(2,l,j,k,m) = fxyz(2,l,j,k,m) + cmplx(at3,0.)
               fxyz(3,l,j,k,m) = fxyz(3,l,j,k,m) + cmplx(at4,0.)
               fxyz(1,l1,j,k,m) = fxyz(1,l1,j,k,m) + cmplx(-at2,0.)
               fxyz(2,l1,j,k,m) = fxyz(2,l1,j,k,m) + cmplx(-at3,0.)
               fxyz(3,l1,j,k,m) = fxyz(3,l1,j,k,m) + cmplx(at4,0.)
   10          continue
c mode numbers kz = 0, nz
               l1 = nz + 1
               at1 = aimag(ffd(1,j,k,m))
               at2 = exyz(1,1,j,k,m)*at1
               at3 = exyz(2,1,j,k,m)*at1
               at4 = exyz(3,1,j,k,m)*at1
               fxyz(1,1,j,k,m) = fxyz(1,1,j,k,m) + cmplx(at2,0.)
               fxyz(2,1,j,k,m) = fxyz(2,1,j,k,m) + cmplx(at3,0.)
               fxyz(3,1,j,k,m) = fxyz(3,1,j,k,m) + cmplx(at4,0.)
            endif
   20       continue
c mode numbers kx = 0, nx
            n1 = joff + 1
            if (n1.eq.0) then
c keep kx = 0
               if (k1.gt.0) then
                  do 30 l = 2, nz
                  l1 = nz2 - l
                  at1 = aimag(ffd(l,1,k,m))
                  at2 = exyz(1,l,1,k,m)*at1
                  at3 = exyz(2,l,1,k,m)*at1
                  at4 = exyz(3,l,1,k,m)*at1
                  fxyz(1,l,1,k,m) = fxyz(1,l,1,k,m) + cmplx(at2,0.)
                  fxyz(2,l,1,k,m) = fxyz(2,l,1,k,m) + cmplx(at3,0.)
                  fxyz(3,l,1,k,m) = fxyz(3,l,1,k,m) + cmplx(at4,0.)
                  fxyz(1,l1,1,k,m) = fxyz(1,l1,1,k,m) + cmplx(-at2,0.)
                  fxyz(2,l1,1,k,m) = fxyz(2,l1,1,k,m) + cmplx(-at3,0.)
                  fxyz(3,l1,1,k,m) = fxyz(3,l1,1,k,m) + cmplx(at4,0.)
   30             continue
c mode numbers kz = 0, nz
                  l1 = nz + 1
                  at1 = aimag(ffd(1,1,k,m))
                  at2 = exyz(1,1,1,k,m)*at1
                  at3 = exyz(2,1,1,k,m)*at1
                  at4 = exyz(3,1,1,k,m)*at1
                  fxyz(1,1,1,k,m) = fxyz(1,1,1,k,m) + cmplx(at2,0.)
                  fxyz(2,1,1,k,m) = fxyz(2,1,1,k,m) + cmplx(at3,0.)
                  fxyz(3,1,1,k,m) = fxyz(3,1,1,k,m) + cmplx(at4,0.)
               endif
            endif
         endif
   40    continue
c mode numbers ky = 0, ny
         n2 = koff + 1
c keep ky = 0
         if (n2.eq.0) then
            do 60 j = 1, kxyp2
            if ((j+joff).gt.0) then
               do 50 l = 2, nz
               l1 = nz2 - l
               at1 = aimag(ffd(l,j,1,m))
               at2 = exyz(1,l,j,1,m)*at1
               at3 = exyz(2,l,j,1,m)*at1
               at4 = exyz(3,l,j,1,m)*at1
               fxyz(1,l,j,1,m) = fxyz(1,l,j,1,m) + cmplx(at2,0.)
               fxyz(2,l,j,1,m) = fxyz(2,l,j,1,m) + cmplx(at3,0.)
               fxyz(3,l,j,1,m) = fxyz(3,l,j,1,m) + cmplx(at4,0.)
               fxyz(1,l1,j,1,m) = fxyz(1,l1,j,1,m) + cmplx(-at2,0.)
               fxyz(2,l1,j,1,m) = fxyz(2,l1,j,1,m) + cmplx(-at3,0.)
               fxyz(3,l1,j,1,m) = fxyz(3,l1,j,1,m) + cmplx(at4,0.)
   50          continue
c mode numbers kz = 0, nz
               l1 = nz + 1
               at1 = aimag(ffd(1,j,1,m))
               at2 = exyz(1,1,j,1,m)*at1
               at3 = exyz(2,1,j,1,m)*at1
               at4 = exyz(3,1,j,1,m)*at1
               fxyz(1,1,j,1,m) = fxyz(1,1,j,1,m) + cmplx(at2,0.)
               fxyz(2,1,j,1,m) = fxyz(2,1,j,1,m) + cmplx(at3,0.)
               fxyz(3,1,j,1,m) = fxyz(3,1,j,1,m) + cmplx(at4,0.)
            endif
   60       continue
c mode numbers kx = 0, nx
            n1 = joff + 1
            if (n1.eq.0) then
               do 70 l = 2, nz
               l1 = nz2 - l
               at1 = aimag(ffd(l,1,1,m))
               at2 = exyz(1,l,1,1,m)*at1
               at3 = exyz(2,l,1,1,m)*at1
               at4 = exyz(3,l,1,1,m)*at1
               fxyz(1,l,1,1,m) = fxyz(1,l,1,1,m) + cmplx(at2,0.)
               fxyz(2,l,1,1,m) = fxyz(2,l,1,1,m) + cmplx(at3,0.)
               fxyz(3,l,1,1,m) = fxyz(3,l,1,1,m) + cmplx(at4,0.)
   70          continue
               l1 = nz + 1
               at1 = aimag(ffd(1,1,1,m))
               at2 = exyz(1,1,1,1,m)*at1
               at3 = exyz(2,1,1,1,m)*at1
               at4 = exyz(3,1,1,1,m)*at1
               fxyz(1,1,1,1,m) = fxyz(1,1,1,1,m) + cmplx(at2,0.)
               fxyz(2,1,1,1,m) = fxyz(2,1,1,1,m) + cmplx(at3,0.)
               fxyz(3,1,1,1,m) = fxyz(3,1,1,1,m) + cmplx(at4,0.)
            endif
         endif
   80    continue
   90    continue
c copy the fields
      else if (isign.lt.0) then
         do 210 my = 1, m2blok
         koff = kyzp2*(my + ks) - 1
         moff = j2blok*(my - 1)
         do 200 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
         joff = kxyp2*(mx + js) - 1
         m = mx + moff
         do 140 k = 1, kyzp2
         k1 = k + koff
         if ((k1.gt.0).and.(k1.ne.ny)) then
            do 110 j = 1, kxyp2
            if ((j+joff).gt.0) then
               do 100 l = 2, nz
               l1 = nz2 - l
               at1 = aimag(ffd(l,j,k,m))
               at2 = exyz(1,l,j,k,m)*at1
               at3 = exyz(2,l,j,k,m)*at1
               at4 = exyz(3,l,j,k,m)*at1
               fxyz(1,l,j,k,m) = cmplx(0.,at2)
               fxyz(2,l,j,k,m) = cmplx(0.,at3)
               fxyz(3,l,j,k,m) = cmplx(0.,at4)
               fxyz(1,l1,j,k,m) = cmplx(0.,at2)
               fxyz(2,l1,j,k,m) = cmplx(0.,at3)
               fxyz(3,l1,j,k,m) = cmplx(0.,-at4)
  100          continue
c mode numbers kz = 0, nz
               l1 = nz + 1
               at1 = aimag(ffd(1,j,k,m))
               at2 = exyz(1,1,j,k,m)*at1
               at3 = exyz(2,1,j,k,m)*at1
               at4 = exyz(3,1,j,k,m)*at1
               fxyz(1,1,j,k,m) = cmplx(0.,at2)
               fxyz(2,1,j,k,m) = cmplx(0.,at3)
               fxyz(3,1,j,k,m) = cmplx(0.,at4)
               fxyz(1,l1,j,k,m) = zero
               fxyz(2,l1,j,k,m) = zero
               fxyz(3,l1,j,k,m) = zero
            endif
  110       continue
c mode numbers kx = 0, nx
            n1 = joff + 1
            if (n1.eq.0) then
c keep kx = 0
               if (k1.gt.0) then
                  do 120 l = 2, nz
                  l1 = nz2 - l
                  at1 = aimag(ffd(l,1,k,m))
                  at2 = exyz(1,l,1,k,m)*at1
                  at3 = exyz(2,l,1,k,m)*at1
                  at4 = exyz(3,l,1,k,m)*at1
                  fxyz(1,l,1,k,m) = cmplx(0.,at2)
                  fxyz(2,l,1,k,m) = cmplx(0.,at3)
                  fxyz(3,l,1,k,m) = cmplx(0.,at4)
                  fxyz(1,l1,1,k,m) = cmplx(0.,at2)
                  fxyz(2,l1,1,k,m) = cmplx(0.,at3)
                  fxyz(3,l1,1,k,m) = cmplx(0.,-at4)
  120             continue
c mode numbers kz = 0, nz
                  l1 = nz + 1
                  at1 = aimag(ffd(1,1,k,m))
                  at2 = exyz(1,1,1,k,m)*at1
                  at3 = exyz(2,1,1,k,m)*at1
                  at4 = exyz(3,1,1,k,m)*at1
                  fxyz(1,1,1,k,m) = cmplx(0.,at2)
                  fxyz(2,1,1,k,m) = cmplx(0.,at3)
                  fxyz(3,1,1,k,m) = cmplx(0.,at4)
                  fxyz(1,l1,1,k,m) = zero
                  fxyz(2,l1,1,k,m) = zero
                  fxyz(3,l1,1,k,m) = zero
c throw away kx = nx
               else
                  do 130 l = 1, nzz
                  fxyz(1,l,1,k,m) = zero
                  fxyz(2,l,1,k,m) = zero
                  fxyz(3,l,1,k,m) = zero
  130             continue
               endif
            endif
         endif
  140    continue
c mode numbers ky = 0, ny
         n2 = koff + 1
c keep ky = 0
         if (n2.eq.0) then
            do 160 j = 1, kxyp2
            if ((j+joff).gt.0) then
               do 150 l = 2, nz
               l1 = nz2 - l
               at1 = aimag(ffd(l,j,1,m))
               at2 = exyz(1,l,j,1,m)*at1
               at3 = exyz(2,l,j,1,m)*at1
               at4 = exyz(3,l,j,1,m)*at1
               fxyz(1,l,j,1,m) = cmplx(0.,at2)
               fxyz(2,l,j,1,m) = cmplx(0.,at3)
               fxyz(3,l,j,1,m) = cmplx(0.,at4)
               fxyz(1,l1,j,1,m) = cmplx(0.,at2)
               fxyz(2,l1,j,1,m) = cmplx(0.,at3)
               fxyz(3,l1,j,1,m) = cmplx(0.,-at4)
  150          continue
c mode numbers kz = 0, nz
               l1 = nz + 1
               at1 = aimag(ffd(1,j,1,m))
               at2 = exyz(1,1,j,1,m)*at1
               at3 = exyz(2,1,j,1,m)*at1
               at4 = exyz(3,1,j,1,m)*at1
               fxyz(1,1,j,1,m) = cmplx(0.,at2)
               fxyz(2,1,j,1,m) = cmplx(0.,at3)
               fxyz(3,1,j,1,m) = cmplx(0.,at4)
               fxyz(1,l1,j,1,m) = zero
               fxyz(2,l1,j,1,m) = zero
               fxyz(3,l1,j,1,m) = zero
            endif
  160       continue
c mode numbers kx = 0, nx
            n1 = joff + 1
            if (n1.eq.0) then
               do 170 l = 2, nz
               l1 = nz2 - l
               at1 = aimag(ffd(l,1,1,m))
               at2 = exyz(1,l,1,1,m)*at1
               at3 = exyz(2,l,1,1,m)*at1
               at4 = exyz(3,l,1,1,m)*at1
               fxyz(1,l,1,1,m) = cmplx(0.,at2)
               fxyz(2,l,1,1,m) = cmplx(0.,at3)
               fxyz(3,l,1,1,m) = cmplx(0.,at4)
               fxyz(1,l1,1,1,m) = zero
               fxyz(2,l1,1,1,m) = zero
               fxyz(3,l1,1,1,m) = zero
  170          continue
               l1 = nz + 1
               at1 = aimag(ffd(1,1,1,m))
               at2 = exyz(1,1,1,1,m)*at1
               at3 = exyz(2,1,1,1,m)*at1
               at4 = exyz(3,1,1,1,m)*at1
               fxyz(1,1,1,1,m) = cmplx(0.,at2)
               fxyz(2,1,1,1,m) = cmplx(0.,at3)
               fxyz(3,1,1,1,m) = cmplx(0.,at4)
               fxyz(1,l1,1,1,m) = zero
               fxyz(2,l1,1,1,m) = zero
               fxyz(3,l1,1,1,m) = zero
            endif
         endif
c throw away ky = ny
         k1 = (ny/kyzp2)*kyzp2
         if (n2.eq.k1) then
            k1 = ny - k1 + 1
            do 190 j = 1, kxyp2
            do 180 l = 1, nzz
            fxyz(1,l,j,k1,m) = zero
            fxyz(2,l,j,k1,m) = zero
            fxyz(3,l,j,k1,m) = zero
  180       continue
  190       continue
         endif
  200    continue
  210    continue
c copy the electric fields
      else
         do 330 my = 1, m2blok
         koff = kyzp2*(my + ks) - 1
         moff = j2blok*(my - 1)
         do 320 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
         joff = kxyp2*(mx + js) - 1
         m = mx + moff
         do 260 k = 1, kyzp2
         k1 = k + koff
         if ((k1.gt.0).and.(k1.ne.ny)) then
            do 230 j = 1, kxyp2
            if ((j+joff).gt.0) then
               do 220 l = 2, nz
               l1 = nz2 - l
               at2 = exyz(1,l,j,k,m)
               at3 = exyz(2,l,j,k,m)
               at4 = exyz(3,l,j,k,m)
               fxyz(1,l,j,k,m) = cmplx(at2,0.)
               fxyz(2,l,j,k,m) = cmplx(at3,0.)
               fxyz(3,l,j,k,m) = cmplx(at4,0.)
               fxyz(1,l1,j,k,m) = cmplx(-at2,0.)
               fxyz(2,l1,j,k,m) = cmplx(-at3,0.)
               fxyz(3,l1,j,k,m) = cmplx(at4,0.)
  220          continue
c mode numbers kz = 0, nz
               l1 = nz + 1
               at2 = exyz(1,1,j,k,m)
               at3 = exyz(2,1,j,k,m)
               at4 = exyz(3,1,j,k,m)
               fxyz(1,1,j,k,m) = cmplx(at2,0.)
               fxyz(2,1,j,k,m) = cmplx(at3,0.)
               fxyz(3,1,j,k,m) = cmplx(at4,0.)
               fxyz(1,l1,j,k,m) = zero
               fxyz(2,l1,j,k,m) = zero
               fxyz(3,l1,j,k,m) = zero
            endif
  230       continue
c mode numbers kx = 0, nx
            n1 = joff + 1
            if (n1.eq.0) then
c keep kx = 0
               if (k1.gt.0) then
                  do 240 l = 2, nz
                  l1 = nz2 - l
                  at2 = exyz(1,l,1,k,m)
                  at3 = exyz(2,l,1,k,m)
                  at4 = exyz(3,l,1,k,m)
                  fxyz(1,l,1,k,m) = cmplx(at2,0.)
                  fxyz(2,l,1,k,m) = cmplx(at3,0.)
                  fxyz(3,l,1,k,m) = cmplx(at4,0.)
                  fxyz(1,l1,1,k,m) = cmplx(-at2,0.)
                  fxyz(2,l1,1,k,m) = cmplx(-at3,0.)
                  fxyz(3,l1,1,k,m) = cmplx(at4,0.)
  240             continue
c mode numbers kz = 0, nz
                  l1 = nz + 1
                  at2 = exyz(1,1,1,k,m)
                  at3 = exyz(2,1,1,k,m)
                  at4 = exyz(3,1,1,k,m)
                  fxyz(1,1,1,k,m) = cmplx(at2,0.)
                  fxyz(2,1,1,k,m) = cmplx(at3,0.)
                  fxyz(3,1,1,k,m) = cmplx(at4,0.)
                  fxyz(1,l1,1,k,m) = zero
                  fxyz(2,l1,1,k,m) = zero
                  fxyz(3,l1,1,k,m) = zero
c throw away kx = nx
               else
                  do 250 l = 1, nzz
                  fxyz(1,l,1,k,m) = zero
                  fxyz(2,l,1,k,m) = zero
                  fxyz(3,l,1,k,m) = zero
  250             continue
               endif
            endif
         endif
  260    continue
c mode numbers ky = 0, ny
         n2 = koff + 1
c keep ky = 0
         if (n2.eq.0) then
            do 280 j = 1, kxyp2
            if ((j+joff).gt.0) then
               do 270 l = 2, nz
               l1 = nz2 - l
               at2 = exyz(1,l,j,1,m)
               at3 = exyz(2,l,j,1,m)
               at4 = exyz(3,l,j,1,m)
               fxyz(1,l,j,1,m) = cmplx(at2,0.)
               fxyz(2,l,j,1,m) = cmplx(at3,0.)
               fxyz(3,l,j,1,m) = cmplx(at4,0.)
               fxyz(1,l1,j,1,m) = cmplx(-at2,0.)
               fxyz(2,l1,j,1,m) = cmplx(-at3,0.)
               fxyz(3,l1,j,1,m) = cmplx(at4,0.)
  270          continue
c mode numbers kz = 0, nz
               l1 = nz + 1
               at2 = exyz(1,1,j,1,m)
               at3 = exyz(2,1,j,1,m)
               at4 = exyz(3,1,j,1,m)
               fxyz(1,1,j,1,m) = cmplx(at2,0.)
               fxyz(2,1,j,1,m) = cmplx(at3,0.)
               fxyz(3,1,j,1,m) = cmplx(at4,0.)
               fxyz(1,l1,j,1,m) = zero
               fxyz(2,l1,j,1,m) = zero
               fxyz(3,l1,j,1,m) = zero
            endif
  280       continue
c mode numbers kx = 0, nx
            n1 = joff + 1
            if (n1.eq.0) then
               do 290 l = 2, nz
               l1 = nz2 - l
               at2 = exyz(1,l,1,1,m)
               at3 = exyz(2,l,1,1,m)
               at4 = exyz(3,l,1,1,m)
               fxyz(1,l,1,1,m) = cmplx(at2,0.)
               fxyz(2,l,1,1,m) = cmplx(at3,0.)
               fxyz(3,l,1,1,m) = cmplx(at4,0.)
               fxyz(1,l1,1,1,m) = zero
               fxyz(2,l1,1,1,m) = zero
               fxyz(3,l1,1,1,m) = zero
  290          continue
               l1 = nz + 1
               at2 = exyz(1,1,1,1,m)
               at3 = exyz(2,1,1,1,m)
               at4 = exyz(3,1,1,1,m)
               fxyz(1,1,1,1,m) = cmplx(at2,0.)
               fxyz(2,1,1,1,m) = cmplx(at3,0.)
               fxyz(3,1,1,1,m) = cmplx(at4,0.)
               fxyz(1,l1,1,1,m) = zero
               fxyz(2,l1,1,1,m) = zero
               fxyz(3,l1,1,1,m) = zero
            endif
         endif
c throw away ky = ny
         k1 = (ny/kyzp2)*kyzp2
         if (n2.eq.k1) then
            k1 = ny - k1 + 1
            do 310 j = 1, kxyp2
            do 300 l = 1, nzz
            fxyz(1,l,j,k1,m) = zero
            fxyz(2,l,j,k1,m) = zero
            fxyz(3,l,j,k1,m) = zero
  300       continue
  310       continue
         endif
  320    continue
  330    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPMFIELDD32(pot3,pot,nx,ny,nz,kstrt,nz2d,nzv,kxyp2,kyzp
     12,j2blok,m2blok)
c copies image charges appropriate for potential
      implicit none
      integer nx, ny, nz, kstrt, nz2d, nzv, kxyp2, kyzp2, j2blok, m2blok
      complex pot3
      real pot
      dimension pot3(nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension pot(nzv,kxyp2,kyzp2,j2blok*m2blok)
      complex zero
      integer j, k, l, m, mx, my, nyy, nzz, nz2, kxb, kyzb, js, ks
      integer joff, koff, moff, k1, l1, n1, n2
      real at1
      nyy = ny + ny
      nzz = nz + nz
      nz2 = nzz + 2
      kxb = nx/kxyp2
      kyzb = nyy/kyzp2
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      zero = cmplx(0.,0.)
      if (kstrt.gt.(kxb*kyzb)) return
      do 120 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 110 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 50 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         do 20 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 10 l = 2, nz
            l1 = nz2 - l
            at1 = pot(l,j,k,m)
            pot3(l,j,k,m) = cmplx(0.,at1)
            pot3(l1,j,k,m) = cmplx(0.,-at1)
   10       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            at1 = pot(1,j,k,m)
            pot3(1,j,k,m) = cmplx(0.,at1)
            pot3(l1,j,k,m) = zero
         endif
   20    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nz
               l1 = nz2 - l
               at1 = pot(l,1,k,m)
               pot3(l,1,k,m) = cmplx(0.,at1)
               pot3(l1,1,k,m) = cmplx(0.,-at1)
   30          continue
c mode numbers kz = 0, nz
               l1 = nz + 1
               at1 = pot(1,1,k,m)
               pot3(1,1,k,m) = cmplx(0.,at1)
               pot3(l1,1,k,m) = zero
c throw away kx = nx
            else
               do 40 l = 1, nzz
               pot3(l,1,k,m) = zero
   40          continue
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 70 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 60 l = 2, nz
            l1 = nz2 - l
            at1 = pot(l,j,1,m)
            pot3(l,j,1,m) = cmplx(0.,at1)
            pot3(l1,j,1,m) = cmplx(0.,-at1)
   60       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            at1 = pot(1,j,1,m)
            pot3(1,j,1,m) = cmplx(0.,at1)
            pot3(l1,j,1,m) = zero
         endif
   70    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 80 l = 2, nz
            l1 = nz2 - l
            at1 = pot(l,1,1,m)
            pot3(l,1,1,m) = cmplx(0.,at1)
            pot3(l1,1,1,m) = zero
   80       continue
            l1 = nz + 1
            at1 = pot(1,1,1,m)
            pot3(1,1,1,m) = cmplx(0.,at1)
            pot3(l1,1,1,m) = zero
         endif
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 100 j = 1, kxyp2
         do 90 l = 1, nzz
         pot3(l,j,k1,m) = zero
   90    continue
  100    continue
      endif
  110 continue
  120 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCPFIELDD32(fxyz,exyz,nx,ny,nz,kstrt,nz2d,nzv,kxyp2,kyz
     1p2,j2blok,m2blok)
c this subroutine copies complex vector fields and adds image charges
c appropriate for electric field
      implicit none
      integer nx, ny, nz, kstrt, nz2d, nzv, kxyp2, kyzp2, j2blok, m2blok
      complex fxyz
      real exyz
      dimension fxyz(3,nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension exyz(3,nzv,kxyp2,kyzp2,j2blok*m2blok)
c local data
      integer isign, nzd
      complex ffd
      dimension ffd(1,1,1,1)
      isign = 0
      nzd = 1
      call PEMFIELDD32(fxyz,exyz,ffd,isign,nx,ny,nz,kstrt,nz2d,nzv,kxyp2
     1,kyzp2,j2blok,m2blok,nzd)
      return
      end
c-----------------------------------------------------------------------
      subroutine PAVPOTDX332(bxyz,axyz,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2,j
     12blok,m2blok)
c this subroutine calculates 3d vector potential from magnetic field
c in fourier space with dirichlet boundary conditions (zero potential)
c for distributed data with 2D spatial decomposition
c input: bxyz,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2,j2blok,m2blok
c output: axyz
c approximate flop count is:
c 18*nxc*nyc*nzc + 5*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx-1)/nvpy, nyc = (ny-1)/nvpz, nzc = nz - 1, and
c nvpy/nvpz = number of procs in y/z
c the vector potential is calculated using the equations:
c ax(kx,ky,kz) = sqrt(-1)*
c                (ky*bz(kx,ky,kz)-kz*by(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c ay(kx,ky,kz) = sqrt(-1)*
c                (kz*bx(kx,ky,kz)-kx*bz(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c az(kx,ky,kz) = sqrt(-1)*
c                (kx*by(kx,ky,kz)-ky*bx(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c bxyz(i,l,j,k,m) = i component of complex magnetic field
c axyz(i,l,j,k,m) = i component of complex vector potential
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp2*(mx - 1) and
c kk = k + kyzp2*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c nx/ny/nz = system length in x/y/z direction
c j2blok/m2blok = number of field partitions in x/y
c kxyp2/kyzp2 = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c nz2d = first dimension of field arrays, must be >= 2*nz
      complex bxyz, axyz, zero
      dimension bxyz(3,nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension axyz(3,nz2d,kxyp2,kyzp2,j2blok*m2blok)
      nyy = ny + ny
      nzz = nz + nz
      nz2 = nzz + 2
      kxb = nx/kxyp2
      kyzb = nyy/kyzp2
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(nyy)
      dnz = 6.28318530717959/float(nzz)
      zero = cmplx(0.,0.)
c calculate vector potential
      if (kstrt.gt.(kxb*kyzb)) go to 130
      do 120 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 110 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 50 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         if (k1.gt.ny) k1 = k1 - nyy
         dky = dny*float(k1)
         dky2 = dky*dky
         do 20 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         dkxy2 = dkx*dkx + dky2
         if ((j+joff).gt.0) then
            do 10 l = 2, nz
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            at1 = 1./(dkz*dkz + dkxy2)
            at2 = dkx*at1
            at3 = dky*at1
            at4 = dkz*at1
            at5 = aimag(bxyz(3,l,j,k,m))
            at6 = aimag(bxyz(2,l,j,k,m))
            at7 = aimag(bxyz(1,l,j,k,m))
            at8 = at4*at6 - at3*at5
            at5 = at2*at5 - at4*at7
            at2 = at3*at7 - at2*at6
            axyz(1,l,j,k,m) = cmplx(at8,0.)
            axyz(2,l,j,k,m) = cmplx(at5,0.)
            axyz(3,l,j,k,m) = cmplx(at2,0.)
            axyz(1,l1,j,k,m) = cmplx(-at8,0.)
            axyz(2,l1,j,k,m) = cmplx(-at5,0.)
            axyz(3,l1,j,k,m) = cmplx(at2,0.)
   10       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            at1 = 1./dkxy2
            at2 = dkx*at1
            at3 = dky*at1
            at2 = at3*aimag(bxyz(1,1,j,k,m)) - at2*aimag(bxyz(2,1,j,k,m)
     1)
            axyz(1,1,j,k,m) = zero
            axyz(2,1,j,k,m) = zero
            axyz(3,1,j,k,m) = cmplx(at2,0.)
            axyz(1,l1,j,k,m) = zero
            axyz(2,l1,j,k,m) = zero
            axyz(3,l1,j,k,m) = zero
         endif
   20    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nz
               l1 = nz2 - l
               dkz = dnz*float(l - 1)
               at1 = 1./(dkz*dkz + dky2)
               at3 = dky*at1
               at4 = dkz*at1
               at8 = at4*aimag(bxyz(2,l,1,k,m)) - at3*aimag(bxyz(3,l,1,k
     1,m))
               axyz(1,l,1,k,m) = cmplx(at8,0.)
               axyz(2,l,1,k,m) = zero
               axyz(3,l,1,k,m) = zero
               axyz(1,l1,1,k,m) = cmplx(-at8,0.)
               axyz(2,l1,1,k,m) = zero
               axyz(3,l1,1,k,m) = zero
   30          continue
c mode numbers kz = 0, nz
               l1 = nz + 1
               axyz(1,1,1,k,m) = zero
               axyz(2,1,1,k,m) = zero
               axyz(3,1,1,k,m) = zero
               axyz(1,l1,1,k,m) = zero
               axyz(2,l1,1,k,m) = zero
               axyz(3,l1,1,k,m) = zero
c throw away kx = nx/2
            else
               do 40 l = 1, nzz
               axyz(1,l,1,k,m) = zero
               axyz(2,l,1,k,m) = zero
               axyz(3,l,1,k,m) = zero
   40          continue
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 70 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         dkx2 = dkx*dkx
         if ((j+joff).gt.0) then
            do 60 l = 2, nz
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            at1 = 1./(dkz*dkz + dkx2)
            at2 = dkx*at1
            at4 = dkz*at1
            at5 = at2*aimag(bxyz(3,l,j,1,m)) - at4*aimag(bxyz(1,l,j,1,m)
     1)
            axyz(1,l,j,1,m) = zero
            axyz(2,l,j,1,m) = cmplx(at5,0.)
            axyz(3,l,j,1,m) = zero
            axyz(1,l1,j,1,m) = zero
            axyz(2,l1,j,1,m) = cmplx(-at5,0.)
            axyz(3,l1,j,1,m) = zero
   60       continue
c mode numbers kz = 0, nz
            l1 = nz + 1
            axyz(1,1,j,1,m) = zero
            axyz(2,1,j,1,m) = zero
            axyz(3,1,j,1,m) = zero
            axyz(1,l1,j,1,m) = zero
            axyz(2,l1,j,1,m) = zero
            axyz(3,l1,j,1,m) = zero
         endif
   70    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 80 l = 1, nzz
            axyz(1,l,1,1,m) = zero
            axyz(2,l,1,1,m) = zero
            axyz(3,l,1,1,m) = zero
   80       continue
         endif
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 100 j = 1, kxyp2
         do 90 l = 1, nzz
         axyz(1,l,j,k1,m) = zero
         axyz(2,l,j,k1,m) = zero
         axyz(3,l,j,k1,m) = zero
   90    continue
  100    continue
      endif
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PAVPOTD332(bxyz,axyz,nx,ny,nz,kstrt,nzv,kxyp2,kyzp2,j2b
     1lok,m2blok)
c this subroutine calculates 3d vector potential from magnetic field
c in fourier space with dirichlet boundary conditions (zero potential)
c for distributed data with 2D spatial decomposition
c input: bxyz,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2,j2blok,m2blok
c output: axyz
c approximate flop count is:
c 16*nxc*nyc*nzc + 5*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx-1)/nvpy, nyc = (ny-1)/nvpz, nzc = nz - 1, and
c nvpy/nvpz = number of procs in y/z
c the vector potential is calculated using the equations:
c ax(kx,ky,kz) = sqrt(-1)*
c                (ky*bz(kx,ky,kz)-kz*by(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c ay(kx,ky,kz) = sqrt(-1)*
c                (kz*bx(kx,ky,kz)-kx*bz(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c az(kx,ky,kz) = sqrt(-1)*
c                (kx*by(kx,ky,kz)-ky*bx(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c bxyz(i,l,j,k,m) = i component of complex magnetic field
c axyz(i,l,j,k,m) = i component of complex vector potential
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp2*(mx - 1) and
c kk = k + kyzp2*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c nx/ny/nz = system length in x/y/z direction
c j2blok/m2blok = number of field partitions in x/y
c kxyp2/kyzp2 = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c nzv = first dimension of field arrays, must be >= nz
      dimension bxyz(3,nzv,kxyp2,kyzp2,j2blok*m2blok)
      dimension axyz(3,nzv,kxyp2,kyzp2,j2blok*m2blok)
      nyy = ny + ny
      kxb = nx/kxyp2
      kyzb = nyy/kyzp2
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx + nx)
      dny = 6.28318530717959/float(nyy)
      dnz = 6.28318530717959/float(nz + nz)
c calculate vector potential
      if (kstrt.gt.(kxb*kyzb)) go to 130
      do 120 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 110 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      m = mx + moff
      do 50 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         if (k1.gt.ny) k1 = k1 - nyy
         dky = dny*float(k1)
         dky2 = dky*dky
         do 20 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         dkxy2 = dkx*dkx + dky2
         if ((j+joff).gt.0) then
            do 10 l = 2, nz
            dkz = dnz*float(l - 1)
            at1 = 1./(dkz*dkz + dkxy2)
            at2 = dkx*at1
            at3 = dky*at1
            at4 = dkz*at1
            at5 = bxyz(3,l,j,k,m)
            at6 = bxyz(2,l,j,k,m)
            at7 = bxyz(1,l,j,k,m)
            axyz(1,l,j,k,m) = at4*at6 - at3*at5
            axyz(2,l,j,k,m) = at2*at5 - at4*at7
            axyz(3,l,j,k,m) = at3*at7 - at2*at6
   10       continue
c mode numbers kz = 0, nz
            at1 = 1./dkxy2
            at2 = dkx*at1
            at3 = dky*at1
            axyz(1,1,j,k,m) = 0.
            axyz(2,1,j,k,m) = 0.
            axyz(3,1,j,k,m) = at3*bxyz(1,1,j,k,m) - at2*bxyz(2,1,j,k,m)
         endif
   20    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nz
               dkz = dnz*float(l - 1)
               at1 = 1./(dkz*dkz + dky2)
               at3 = dky*at1
               at4 = dkz*at1
               axyz(1,l,1,k,m) = at4*bxyz(2,l,1,k,m) - at3*bxyz(3,l,1,k,
     1m)
               axyz(2,l,1,k,m) = 0.
               axyz(3,l,1,k,m) = 0.
   30          continue
c mode numbers kz = 0, nz
               axyz(1,1,1,k,m) = 0.
               axyz(2,1,1,k,m) = 0.
               axyz(3,1,1,k,m) = 0.
c throw away kx = nx/2
            else
               do 40 l = 1, nz
               axyz(1,l,1,k,m) = 0.
               axyz(2,l,1,k,m) = 0.
               axyz(3,l,1,k,m) = 0.
   40          continue
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 70 j = 1, kxyp2
         dkx = dnx*float(j + joff)
         dkx2 = dkx*dkx
         if ((j+joff).gt.0) then
            do 60 l = 2, nz
            dkz = dnz*float(l - 1)
            at1 = 1./(dkz*dkz + dkx2)
            at2 = dkx*at1
            at4 = dkz*at1
            axyz(1,l,j,1,m) = 0.
            axyz(2,l,j,1,m) = at2*bxyz(3,l,j,1,m) - at4*bxyz(1,l,j,1,m)
            axyz(3,l,j,1,m) = 0.
   60       continue
c mode numbers kz = 0, nz
            axyz(1,1,j,1,m) = 0.
            axyz(2,1,j,1,m) = 0.
            axyz(3,1,j,1,m) = 0.
         endif
   70    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 80 l = 1, nz
            axyz(1,l,1,1,m) = 0.
            axyz(2,l,1,1,m) = 0.
            axyz(3,l,1,1,m) = 0.
   80       continue
         endif
      endif
c throw away ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 100 j = 1, kxyp2
         do 90 l = 1, nz
         axyz(1,l,j,k1,m) = 0.
         axyz(2,l,j,k1,m) = 0.
         axyz(3,l,j,k1,m) = 0.
   90    continue
  100    continue
      endif
  110 continue
  120 continue
  130 continue
      return
      end
