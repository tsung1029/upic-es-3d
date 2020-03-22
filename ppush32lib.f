c 3d parallel PIC library for pushing particles and depositing charge
c with 2D domain decomposition
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: august 27, 2005
c-----------------------------------------------------------------------
      subroutine PDOST32(part,q,npp,noff,qm,nx,idimp,npmax,mnblok,nxv,ny
     1pmx,nzpmx,idds)
c for 3d code, this subroutine calculates particle charge density
c using second-order spline interpolation, periodic boundaries
c and distributed data with 2D spatial decomposition
c baseline scalar distributed version
c 118 flops/particle, 30 loads, 27 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m,l)=qm*(.75-dx**2)*(.75-dy**2)*(.75-dz**2)
c q(n+1,m,l)=.5*qm*((.5+dx)**2)*(.75-dy**2)*(.75-dz**2)
c q(n-1,m,l)=.5*qm*((.5-dx)**2)*(.75-dy**2)*(.75-dz**2)
c q(n,m+1,l)=.5*qm*(.75-dx**2)*(.5+dy)**2*(.75-dz**2)
c q(n+1,m+1,l)=.25*qm*((.5+dx)**2)*((.5+dy)**2)*(.75-dz**2)
c q(n-1,m+1,l)=.25*qm*((.5-dx)**2)*((.5+dy)**2)*(.75-dz**2)
c q(n,m-1,l)=.5*qm*(.75-dx**2)*((.5-dy)**2)*(.75-dz**2)
c q(n+1,m-1,l)=.25*qm*((.5+dx)**2)*((.5-dy)**2)*(.75-dz**2)
c q(n-1,m-1,l)=.25*qm*((.5-dx)**2)*((.5-dy)**2)*(.75-dz**2)
c q(n,m,l+1)=.5*qm*(.75-dx**2)*(.75-dy**2)*((.5+dz)**2)
c q(n+1,m,l+1)=.25*qm*((.5+dx)**2)*(.75-dy**2)*((.5+dz)**2)
c q(n-1,m,l+1)=.25*qm*((.5-dx)**2)*(.75-dy**2)*((.5+dz)**2)
c q(n,m+1,l+1)=.25*qm*(.75-dx**2)*(.5+dy)**2*((.5+dz)**2)
c q(n+1,m+1,l+1)=.125*qm*((.5+dx)**2)*((.5+dy)**2)*((.5+dz)**2)
c q(n-1,m+1,l+1)=.125*qm*((.5-dx)**2)*((.5+dy)**2)*((.5+dz)**2)
c q(n,m-1,l+1)=.25*qm*(.75-dx**2)*((.5-dy)**2)*((.5+dz)**2)
c q(n+1,m-1,l+1)=.125*qm*((.5+dx)**2)*((.5-dy)**2)*((.5+dz)**2)
c q(n-1,m-1,l+1)=.125*qm*((.5-dx)**2)*((.5-dy)**2)*((.5+dz)**2)
c q(n,m,l-1)=.5*qm*(.75-dx**2)*(.75-dy**2)*((.5-dz)**2)
c q(n+1,m,l-1)=.25*qm*((.5+dx)**2)*(.75-dy**2)*((.5-dz)**2)
c q(n-1,m,l-1)=.25*qm*((.5-dx)**2)*(.75-dy**2)*((.5-dz)**2)
c q(n,m+1,l-1)=.25*qm*(.75-dx**2)*(.5+dy)**2*((.5-dz)**2)
c q(n+1,m+1,l-1)=.125*qm*((.5+dx)**2)*((.5+dy)**2)*((.5-dz)**2)
c q(n-1,m+1,l-1)=.125*qm*((.5-dx)**2)*((.5+dy)**2)*((.5-dz)**2)
c q(n,m-1,l-1)=.25*qm*(.75-dx**2)*((.5-dy)**2)*((.5-dz)**2)
c q(n+1,m-1,l-1)=.125*qm*((.5+dx)**2)*((.5-dy)**2)*((.5-dz)**2)
c q(n-1,m-1,l-1)=.125*qm*((.5-dx)**2)*((.5-dy)**2)*((.5-dz)**2)
c where n,m,l = nearest grid points and dx = x-n, dy = y-m, dz = z-l
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c q(j,k,l,m) = charge density at grid point (j,kk,ll),
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qm = charge on particle, in units of e
c nx = system length in x direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first dimension of charge array, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
      dimension part(idimp,npmax,mnblok), q(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      qmh = .5*qm
      do 20 m = 1, mnblok
      mnoff = noff(1,m) - 1
      lnoff = noff(2,m) - 1
c find interpolation weights
      do 10 j = 1, npp(m)
      nn = part(1,j,m) + .5
      mm = part(2,j,m) + .5
      ll = part(3,j,m) + .5
      dx = part(1,j,m) - float(nn)
      dy = part(2,j,m) - float(mm)
      dz = part(3,j,m) - float(ll)
      nl = nn
      if (nl.lt.1) nl = nl + nx
      amx = qm*(.75 - dx*dx)
      ml = mm - mnoff
      dxl = qmh*(.5 - dx)**2
      lm = ll - lnoff
      dxp = qmh*(.5 + dx)**2
      nn = nn + 1
      if (nn.gt.nx) nn = nn - nx
      dyl = .5*(.5 - dy)**2
      np = nn + 1
      if (np.gt.nx) np = np - nx
      amy = .75 - dy*dy
      mm = ml + 1
      dyp = .5*(.5 + dy)**2
      mp = mm + 1
      dzl = .5*(.5 - dz)**2
      ll = lm + 1
      amz = .75 - dz*dz
      lp = ll + 1
      dzp = .5*(.5 + dz)**2
c deposit charge
      q(nl,ml,lm,m) = q(nl,ml,lm,m) + dxl*dyl*dzl
      q(nn,ml,lm,m) = q(nn,ml,lm,m) + amx*dyl*dzl
      q(np,ml,lm,m) = q(np,ml,lm,m) + dxp*dyl*dzl
      q(nl,mm,lm,m) = q(nl,mm,lm,m) + dxl*amy*dzl
      q(nn,mm,lm,m) = q(nn,mm,lm,m) + amx*amy*dzl
      q(np,mm,lm,m) = q(np,mm,lm,m) + dxp*amy*dzl
      q(nl,mp,lm,m) = q(nl,mp,lm,m) + dxl*dyp*dzl
      q(nn,mp,lm,m) = q(nn,mp,lm,m) + amx*dyp*dzl
      q(np,mp,lm,m) = q(np,mp,lm,m) + dxp*dyp*dzl
      q(nl,ml,ll,m) = q(nl,ml,ll,m) + dxl*dyl*amz
      q(nn,ml,ll,m) = q(nn,ml,ll,m) + amx*dyl*amz
      q(np,ml,ll,m) = q(np,ml,ll,m) + dxp*dyl*amz
      q(nl,mm,ll,m) = q(nl,mm,ll,m) + dxl*amy*amz
      q(nn,mm,ll,m) = q(nn,mm,ll,m) + amx*amy*amz
      q(np,mm,ll,m) = q(np,mm,ll,m) + dxp*amy*amz
      q(nl,mp,ll,m) = q(nl,mp,ll,m) + dxl*dyp*amz
      q(nn,mp,ll,m) = q(nn,mp,ll,m) + amx*dyp*amz
      q(np,mp,ll,m) = q(np,mp,ll,m) + dxp*dyp*amz
      q(nl,ml,lp,m) = q(nl,ml,lp,m) + dxl*dyl*dzp
      q(nn,ml,lp,m) = q(nn,ml,lp,m) + amx*dyl*dzp
      q(np,ml,lp,m) = q(np,ml,lp,m) + dxp*dyl*dzp
      q(nl,mm,lp,m) = q(nl,mm,lp,m) + dxl*amy*dzp
      q(nn,mm,lp,m) = q(nn,mm,lp,m) + amx*amy*dzp
      q(np,mm,lp,m) = q(np,mm,lp,m) + dxp*amy*dzp
      q(nl,mp,lp,m) = q(nl,mp,lp,m) + dxl*dyp*dzp
      q(nn,mp,lp,m) = q(nn,mp,lp,m) + amx*dyp*dzp
      q(np,mp,lp,m) = q(np,mp,lp,m) + dxp*dyp*dzp
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGPOST32(part,q,npp,noff,qm,idimp,npmax,mnblok,nxv,nypm
     1x,nzpmx,idds)
c for 3d code, this subroutine calculates particle charge density
c using second-order spline interpolation, and distributed data
c with 2D spatial decomposition
c scalar version using guard cells, for distributed data
c 100 flops/particle, 30 loads, 27 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m,l)=qm*(.75-dx**2)*(.75-dy**2)*(.75-dz**2)
c q(n+1,m,l)=.5*qm*((.5+dx)**2)*(.75-dy**2)*(.75-dz**2)
c q(n-1,m,l)=.5*qm*((.5-dx)**2)*(.75-dy**2)*(.75-dz**2)
c q(n,m+1,l)=.5*qm*(.75-dx**2)*(.5+dy)**2*(.75-dz**2)
c q(n+1,m+1,l)=.25*qm*((.5+dx)**2)*((.5+dy)**2)*(.75-dz**2)
c q(n-1,m+1,l)=.25*qm*((.5-dx)**2)*((.5+dy)**2)*(.75-dz**2)
c q(n,m-1,l)=.5*qm*(.75-dx**2)*((.5-dy)**2)*(.75-dz**2)
c q(n+1,m-1,l)=.25*qm*((.5+dx)**2)*((.5-dy)**2)*(.75-dz**2)
c q(n-1,m-1,l)=.25*qm*((.5-dx)**2)*((.5-dy)**2)*(.75-dz**2)
c q(n,m,l+1)=.5*qm*(.75-dx**2)*(.75-dy**2)*((.5+dz)**2)
c q(n+1,m,l+1)=.25*qm*((.5+dx)**2)*(.75-dy**2)*((.5+dz)**2)
c q(n-1,m,l+1)=.25*qm*((.5-dx)**2)*(.75-dy**2)*((.5+dz)**2)
c q(n,m+1,l+1)=.25*qm*(.75-dx**2)*(.5+dy)**2*((.5+dz)**2)
c q(n+1,m+1,l+1)=.125*qm*((.5+dx)**2)*((.5+dy)**2)*((.5+dz)**2)
c q(n-1,m+1,l+1)=.125*qm*((.5-dx)**2)*((.5+dy)**2)*((.5+dz)**2)
c q(n,m-1,l+1)=.25*qm*(.75-dx**2)*((.5-dy)**2)*((.5+dz)**2)
c q(n+1,m-1,l+1)=.125*qm*((.5+dx)**2)*((.5-dy)**2)*((.5+dz)**2)
c q(n-1,m-1,l+1)=.125*qm*((.5-dx)**2)*((.5-dy)**2)*((.5+dz)**2)
c q(n,m,l-1)=.5*qm*(.75-dx**2)*(.75-dy**2)*((.5-dz)**2)
c q(n+1,m,l-1)=.25*qm*((.5+dx)**2)*(.75-dy**2)*((.5-dz)**2)
c q(n-1,m,l-1)=.25*qm*((.5-dx)**2)*(.75-dy**2)*((.5-dz)**2)
c q(n,m+1,l-1)=.25*qm*(.75-dx**2)*(.5+dy)**2*((.5-dz)**2)
c q(n+1,m+1,l-1)=.125*qm*((.5+dx)**2)*((.5+dy)**2)*((.5-dz)**2)
c q(n-1,m+1,l-1)=.125*qm*((.5-dx)**2)*((.5+dy)**2)*((.5-dz)**2)
c q(n,m-1,l-1)=.25*qm*(.75-dx**2)*((.5-dy)**2)*((.5-dz)**2)
c q(n+1,m-1,l-1)=.125*qm*((.5+dx)**2)*((.5-dy)**2)*((.5-dz)**2)
c q(n-1,m-1,l-1)=.125*qm*((.5-dx)**2)*((.5-dy)**2)*((.5-dz)**2)
c where n,m,l = nearest grid points and dx = x-n, dy = y-m, dz = z-l
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c q(j+1,k,l,m) = charge density at grid point (j,kk,ll),
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qm = charge on particle, in units of e
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first dimension of charge array, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
      dimension part(idimp,npmax,mnblok), q(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      qmh = .5*qm
      do 20 m = 1, mnblok
      mnoff = noff(1,m) - 1
      lnoff = noff(2,m) - 1
c find interpolation weights
      do 10 j = 1, npp(m)
      nn = part(1,j,m) + .5
      mm = part(2,j,m) + .5
      ll = part(3,j,m) + .5
      dx = part(1,j,m) - float(nn)
      dy = part(2,j,m) - float(mm)
      dz = part(3,j,m) - float(ll)
      nl = nn + 1
      amx = qm*(.75 - dx*dx)
      ml = mm - mnoff
      dxl = qmh*(.5 - dx)**2
      lm = ll - lnoff
      dxp = qmh*(.5 + dx)**2
      nn = nl + 1
      dyl = .5*(.5 - dy)**2
      np = nl + 2
      amy = .75 - dy*dy
      mm = ml + 1
      dyp = .5*(.5 + dy)**2
      mp = ml + 2
      dx1 = dxl*dyl
      dx2 = amx*dyl
      dyl = dxp*dyl
      dzl = .5*(.5 - dz)**2
      ll = lm + 1
      dy1 = dxl*amy
      dy2 = amx*amy
      amy = dxp*amy
      amz = .75 - dz*dz
      lp = lm + 2
      dxl = dxl*dyp
      amx = amx*dyp
      dyp = dxp*dyp
      dzp = .5*(.5 + dz)**2
c deposit charge
      q(nl,ml,lm,m) = q(nl,ml,lm,m) + dx1*dzl
      q(nn,ml,lm,m) = q(nn,ml,lm,m) + dx2*dzl
      q(np,ml,lm,m) = q(np,ml,lm,m) + dyl*dzl
      q(nl,mm,lm,m) = q(nl,mm,lm,m) + dy1*dzl
      q(nn,mm,lm,m) = q(nn,mm,lm,m) + dy2*dzl
      q(np,mm,lm,m) = q(np,mm,lm,m) + amy*dzl
      q(nl,mp,lm,m) = q(nl,mp,lm,m) + dxl*dzl
      q(nn,mp,lm,m) = q(nn,mp,lm,m) + amx*dzl
      q(np,mp,lm,m) = q(np,mp,lm,m) + dyp*dzl
      q(nl,ml,ll,m) = q(nl,ml,ll,m) + dx1*amz
      q(nn,ml,ll,m) = q(nn,ml,ll,m) + dx2*amz
      q(np,ml,ll,m) = q(np,ml,ll,m) + dyl*amz
      q(nl,mm,ll,m) = q(nl,mm,ll,m) + dy1*amz
      q(nn,mm,ll,m) = q(nn,mm,ll,m) + dy2*amz
      q(np,mm,ll,m) = q(np,mm,ll,m) + amy*amz
      q(nl,mp,ll,m) = q(nl,mp,ll,m) + dxl*amz
      q(nn,mp,ll,m) = q(nn,mp,ll,m) + amx*amz
      q(np,mp,ll,m) = q(np,mp,ll,m) + dyp*amz
      q(nl,ml,lp,m) = q(nl,ml,lp,m) + dx1*dzp
      q(nn,ml,lp,m) = q(nn,ml,lp,m) + dx2*dzp
      q(np,ml,lp,m) = q(np,ml,lp,m) + dyl*dzp
      q(nl,mm,lp,m) = q(nl,mm,lp,m) + dy1*dzp
      q(nn,mm,lp,m) = q(nn,mm,lp,m) + dy2*dzp
      q(np,mm,lp,m) = q(np,mm,lp,m) + amy*dzp
      q(nl,mp,lp,m) = q(nl,mp,lp,m) + dxl*dzp
      q(nn,mp,lp,m) = q(nn,mp,lp,m) + amx*dzp
      q(np,mp,lp,m) = q(np,mp,lp,m) + dyp*dzp
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSPOST32(part,q,npp,noff,qm,idimp,npmax,mnblok,nxv,nyp
     1mx,nxyzp,idds)
c for 3d code, this subroutine calculates particle charge density
c using second-order spline interpolation, and distributed data
c with 2D spatial decomposition
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 100 flops/particle, 30 loads, 27 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m,l)=qm*(.75-dx**2)*(.75-dy**2)*(.75-dz**2)
c q(n+1,m,l)=.5*qm*((.5+dx)**2)*(.75-dy**2)*(.75-dz**2)
c q(n-1,m,l)=.5*qm*((.5-dx)**2)*(.75-dy**2)*(.75-dz**2)
c q(n,m+1,l)=.5*qm*(.75-dx**2)*(.5+dy)**2*(.75-dz**2)
c q(n+1,m+1,l)=.25*qm*((.5+dx)**2)*((.5+dy)**2)*(.75-dz**2)
c q(n-1,m+1,l)=.25*qm*((.5-dx)**2)*((.5+dy)**2)*(.75-dz**2)
c q(n,m-1,l)=.5*qm*(.75-dx**2)*((.5-dy)**2)*(.75-dz**2)
c q(n+1,m-1,l)=.25*qm*((.5+dx)**2)*((.5-dy)**2)*(.75-dz**2)
c q(n-1,m-1,l)=.25*qm*((.5-dx)**2)*((.5-dy)**2)*(.75-dz**2)
c q(n,m,l+1)=.5*qm*(.75-dx**2)*(.75-dy**2)*((.5+dz)**2)
c q(n+1,m,l+1)=.25*qm*((.5+dx)**2)*(.75-dy**2)*((.5+dz)**2)
c q(n-1,m,l+1)=.25*qm*((.5-dx)**2)*(.75-dy**2)*((.5+dz)**2)
c q(n,m+1,l+1)=.25*qm*(.75-dx**2)*(.5+dy)**2*((.5+dz)**2)
c q(n+1,m+1,l+1)=.125*qm*((.5+dx)**2)*((.5+dy)**2)*((.5+dz)**2)
c q(n-1,m+1,l+1)=.125*qm*((.5-dx)**2)*((.5+dy)**2)*((.5+dz)**2)
c q(n,m-1,l+1)=.25*qm*(.75-dx**2)*((.5-dy)**2)*((.5+dz)**2)
c q(n+1,m-1,l+1)=.125*qm*((.5+dx)**2)*((.5-dy)**2)*((.5+dz)**2)
c q(n-1,m-1,l+1)=.125*qm*((.5-dx)**2)*((.5-dy)**2)*((.5+dz)**2)
c q(n,m,l-1)=.5*qm*(.75-dx**2)*(.75-dy**2)*((.5-dz)**2)
c q(n+1,m,l-1)=.25*qm*((.5+dx)**2)*(.75-dy**2)*((.5-dz)**2)
c q(n-1,m,l-1)=.25*qm*((.5-dx)**2)*(.75-dy**2)*((.5-dz)**2)
c q(n,m+1,l-1)=.25*qm*(.75-dx**2)*(.5+dy)**2*((.5-dz)**2)
c q(n+1,m+1,l-1)=.125*qm*((.5+dx)**2)*((.5+dy)**2)*((.5-dz)**2)
c q(n-1,m+1,l-1)=.125*qm*((.5-dx)**2)*((.5+dy)**2)*((.5-dz)**2)
c q(n,m-1,l-1)=.25*qm*(.75-dx**2)*((.5-dy)**2)*((.5-dz)**2)
c q(n+1,m-1,l-1)=.125*qm*((.5+dx)**2)*((.5-dy)**2)*((.5-dz)**2)
c q(n-1,m-1,l-1)=.125*qm*((.5-dx)**2)*((.5-dy)**2)*((.5-dz)**2)
c where n,m,l = nearest grid points and dx = x-n, dy = y-m, dz = z-l
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c q(n,m) = charge density at grid point (j,kk,ll),
c where n = j + nxv*kk + nxv*nypmx*ll
c and kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qm = charge on particle, in units of e
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first virtual dimension of charge array, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nxyzp = dimension of charge array, must be >= nxv*nypmx*nzpmx
c idds = dimensionality of domain decomposition
      dimension part(idimp,npmax,mnblok), q(nxyzp,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      qmh = .5*qm
      nxvy = nxv*nypmx
      do 20 m = 1, mnblok
      if (npp(m).lt.1) go to 20
      mnoff = noff(1,m)
      lnoff = noff(2,m)
c begin first particle
      nnn = part(1,1,m) + .5
      mmm = part(2,1,m) + .5
      lll = part(3,1,m) + .5
      dxn = part(1,1,m) - float(nnn)
      dyn = part(2,1,m) - float(mmm)
      dzn = part(3,1,m) - float(lll)
      mmm = mmm - mnoff
      lll = lll - lnoff
c find interpolation weights
      do 10 j = 2, npp(m)
      nn = nnn + 1
      mm = nxv*mmm + nxvy*lll
      nnn = part(1,j,m) + .5
      mmm = part(2,j,m) + .5
      lll = part(3,j,m) + .5
      dx = dxn
      dy = dyn
      dz = dzn
      dxn = part(1,j,m) - float(nnn)
      dyn = part(2,j,m) - float(mmm)
      dzn = part(3,j,m) - float(lll)
      amx = qm*(.75 - dx*dx)
      dxl = qmh*(.5 - dx)**2
      dxp = qmh*(.5 + dx)**2
      ml = mm + nn
      dyl = .5*(.5 - dy)**2
      amy = .75 - dy*dy
      dyp = .5*(.5 + dy)**2
      mn = ml + nxv
      dx1 = dxl*dyl
      dx2 = amx*dyl
      dyl = dxp*dyl
      dzl = .5*(.5 - dz)**2
      dy1 = dxl*amy
      dy2 = amx*amy
      amy = dxp*amy
      amz = .75 - dz*dz
      dxl = dxl*dyp
      amx = amx*dyp
      dyp = dxp*dyp
      dzp = .5*(.5 + dz)**2
      mp = mn + nxv
      mmm = mmm - mnoff
      lll = lll - lnoff
c deposit charge
      dx = q(ml,m) + dx1*dzl
      dy = q(ml+1,m) + dx2*dzl
      dz = q(ml+2,m) + dyl*dzl
      dx3 = q(mn,m) + dy1*dzl
      dy3 = q(mn+1,m) + dy2*dzl
      dz3 = q(mn+2,m) + amy*dzl
      dxp = q(mp,m) + dxl*dzl
      dy4 = q(mp+1,m) + amx*dzl
      dzl = q(mp+2,m) + dyp*dzl
      q(ml,m) = dx
      q(ml+1,m) = dy
      q(ml+2,m) = dz
      q(mn,m) = dx3
      q(mn+1,m) = dy3
      q(mn+2,m) = dz3
      q(mp,m) = dxp
      q(mp+1,m) = dy4
      q(mp+2,m) = dzl
      ml = ml + nxvy
      mn = mn + nxvy
      mp = mp + nxvy
      dx = q(ml,m) + dx1*amz
      dy = q(ml+1,m) + dx2*amz
      dz = q(ml+2,m) + dyl*amz
      dx3 = q(mn,m) + dy1*amz
      dy3 = q(mn+1,m) + dy2*amz
      dz3 = q(mn+2,m) + amy*amz
      dxp = q(mp,m) + dxl*amz
      dzl = q(mp+1,m) + amx*amz
      amz = q(mp+2,m) + dyp*amz
      q(ml,m) = dx
      q(ml+1,m) = dy
      q(ml+2,m) = dz
      q(mn,m) = dx3
      q(mn+1,m) = dy3
      q(mn+2,m) = dz3
      q(mp,m) = dxp
      q(mp+1,m) = dzl
      q(mp+2,m) = amz
      ml = ml + nxvy
      mn = mn + nxvy
      mp = mp + nxvy
      dx = q(ml,m) + dx1*dzp
      dy = q(ml+1,m) + dx2*dzp
      dz = q(ml+2,m) + dyl*dzp
      dx3 = q(mn,m) + dy1*dzp
      dy3 = q(mn+1,m) + dy2*dzp
      dz3 = q(mn+2,m) + amy*dzp
      dxp = q(mp,m) + dxl*dzp
      dzl = q(mp+1,m) + amx*dzp
      amz = q(mp+2,m) + dyp*dzp
      q(ml,m) = dx
      q(ml+1,m) = dy
      q(ml+2,m) = dz
      q(mn,m) = dx3
      q(mn+1,m) = dy3
      q(mn+2,m) = dz3
      q(mp,m) = dxp
      q(mp+1,m) = dzl
      q(mp+2,m) = amz
   10 continue
c deposit charge for last particle
      nn = nnn + 1
      mm = nxv*mmm + nxvy*lll
      amx = qm*(.75 - dxn*dxn)
      dxl = qmh*(.5 - dxn)**2
      dxp = qmh*(.5 + dxn)**2
      ml = mm + nn
      dyl = .5*(.5 - dyn)**2
      amy = .75 - dyn*dyn
      dyp = .5*(.5 + dyn)**2
      mn = ml + nxv
      dx1 = dxl*dyl
      dx2 = amx*dyl
      dyl = dxp*dyl
      dzl = .5*(.5 - dzn)**2
      dy1 = dxl*amy
      dy2 = amx*amy
      amy = dxp*amy
      amz = .75 - dzn*dzn
      dxl = dxl*dyp
      amx = amx*dyp
      dyp = dxp*dyp
      dzp = .5*(.5 + dzn)**2
      mp = mn + nxv
c deposit charge
      q(ml,m) = q(ml,m) + dx1*dzl
      q(ml+1,m) = q(ml+1,m) + dx2*dzl
      q(ml+2,m) = q(ml+2,m) + dyl*dzl
      q(mn,m) = q(mn,m) + dy1*dzl
      q(mn+1,m) = q(mn+1,m) + dy2*dzl
      q(mn+2,m) = q(mn+2,m) + amy*dzl
      q(mp,m) = q(mp,m) + dxl*dzl
      q(mp+1,m) = q(mp+1,m) + amx*dzl
      q(mp+2,m) = q(mp+2,m) + dyp*dzl
      ml = ml + nxvy
      mn = mn + nxvy
      mp = mp + nxvy
      q(ml,m) = q(ml,m) + dx1*amz
      q(ml+1,m) = q(ml+1,m) + dx2*amz
      q(ml+2,m) = q(ml+2,m) + dyl*amz
      q(mn,m) = q(mn,m) + dy1*amz
      q(mn+1,m) = q(mn+1,m) + dy2*amz
      q(mn+2,m) = q(mn+2,m) + amy*amz
      q(mp,m) = q(mp,m) + dxl*amz
      q(mp+1,m) = q(mp+1,m) + amx*amz
      q(mp+2,m) = q(mp+2,m) + dyp*amz
      ml = ml + nxvy
      mn = mn + nxvy
      mp = mp + nxvy
      q(ml,m) = q(ml,m) + dx1*dzp
      q(ml+1,m) = q(ml+1,m) + dx2*dzp
      q(ml+2,m) = q(ml+2,m) + dyl*dzp
      q(mn,m) = q(mn,m) + dy1*dzp
      q(mn+1,m) = q(mn+1,m) + dy2*dzp
      q(mn+2,m) = q(mn+2,m) + amy*dzp
      q(mp,m) = q(mp,m) + dxl*dzp
      q(mp+1,m) = q(mp+1,m) + amx*dzp
      q(mp+2,m) = q(mp+2,m) + dyp*dzp
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSOST32X(part,q,npp,noff,nn,amxyz,qm,nx,idimp,npmax,mnb
     1lok,nxv,nypmx,nxvyzp,idds,npd,n27)
c for 3d code, this subroutine calculates particle charge density
c using second-order spline interpolation, periodic boundaries,
c with short vectors over independent weights, and distributed data,
c with 2D spatial decomposition
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c vectorized distributed version with 1d addressing
c 118 flops/particle, 30 loads, 27 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m,l)=qm*(.75-dx**2)*(.75-dy**2)*(.75-dz**2)
c q(n+1,m,l)=.5*qm*((.5+dx)**2)*(.75-dy**2)*(.75-dz**2)
c q(n-1,m,l)=.5*qm*((.5-dx)**2)*(.75-dy**2)*(.75-dz**2)
c q(n,m+1,l)=.5*qm*(.75-dx**2)*(.5+dy)**2*(.75-dz**2)
c q(n+1,m+1,l)=.25*qm*((.5+dx)**2)*((.5+dy)**2)*(.75-dz**2)
c q(n-1,m+1,l)=.25*qm*((.5-dx)**2)*((.5+dy)**2)*(.75-dz**2)
c q(n,m-1,l)=.5*qm*(.75-dx**2)*((.5-dy)**2)*(.75-dz**2)
c q(n+1,m-1,l)=.25*qm*((.5+dx)**2)*((.5-dy)**2)*(.75-dz**2)
c q(n-1,m-1,l)=.25*qm*((.5-dx)**2)*((.5-dy)**2)*(.75-dz**2)
c q(n,m,l+1)=.5*qm*(.75-dx**2)*(.75-dy**2)*((.5+dz)**2)
c q(n+1,m,l+1)=.25*qm*((.5+dx)**2)*(.75-dy**2)*((.5+dz)**2)
c q(n-1,m,l+1)=.25*qm*((.5-dx)**2)*(.75-dy**2)*((.5+dz)**2)
c q(n,m+1,l+1)=.25*qm*(.75-dx**2)*(.5+dy)**2*((.5+dz)**2)
c q(n+1,m+1,l+1)=.125*qm*((.5+dx)**2)*((.5+dy)**2)*((.5+dz)**2)
c q(n-1,m+1,l+1)=.125*qm*((.5-dx)**2)*((.5+dy)**2)*((.5+dz)**2)
c q(n,m-1,l+1)=.25*qm*(.75-dx**2)*((.5-dy)**2)*((.5+dz)**2)
c q(n+1,m-1,l+1)=.125*qm*((.5+dx)**2)*((.5-dy)**2)*((.5+dz)**2)
c q(n-1,m-1,l+1)=.125*qm*((.5-dx)**2)*((.5-dy)**2)*((.5+dz)**2)
c q(n,m,l-1)=.5*qm*(.75-dx**2)*(.75-dy**2)*((.5-dz)**2)
c q(n+1,m,l-1)=.25*qm*((.5+dx)**2)*(.75-dy**2)*((.5-dz)**2)
c q(n-1,m,l-1)=.25*qm*((.5-dx)**2)*(.75-dy**2)*((.5-dz)**2)
c q(n,m+1,l-1)=.25*qm*(.75-dx**2)*(.5+dy)**2*((.5-dz)**2)
c q(n+1,m+1,l-1)=.125*qm*((.5+dx)**2)*((.5+dy)**2)*((.5-dz)**2)
c q(n-1,m+1,l-1)=.125*qm*((.5-dx)**2)*((.5+dy)**2)*((.5-dz)**2)
c q(n,m-1,l-1)=.25*qm*(.75-dx**2)*((.5-dy)**2)*((.5-dz)**2)
c q(n+1,m-1,l-1)=.125*qm*((.5+dx)**2)*((.5-dy)**2)*((.5-dz)**2)
c q(n-1,m-1,l-1)=.125*qm*((.5-dx)**2)*((.5-dy)**2)*((.5-dz)**2)
c where n,m,l = nearest grid points and dx = x-n, dy = y-m, dz = z-l
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c q(j,k,l,m) = charge density at grid point (j,kk,ll),
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nn = scratch address array for vectorized charge deposition
c amxyz = scratch weight array for vectorized charge deposition
c qm = charge on particle, in units of e
c nx = system length in x direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first virtual dimension of charge array, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nxvyzp = nxv*nypmx*nzpmx, first dimension of charge array
c idds = dimensionality of domain decomposition
c npd = size of scratch buffers for vectorized push/charge deposition
c n27 = number of independent weights
c version with spatial decomposition and 1d addressing scheme
      dimension part(idimp,npmax,mnblok), q(nxvyzp,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      dimension nn(n27,npd,mnblok), amxyz(n27,npd,mnblok)
      nxvy = nxv*nypmx
c parallel loop
      do 50 m = 1, mnblok
      mnoff = noff(1,m) - 1
      lnoff = noff(2,m) - 1
      npb = npd
      if (npp(m).gt.npd) then
         ipp = float(npp(m) - 1)/float(npd) + 1.
      else
         ipp = 1
      endif
c outer loop over blocks of particles
      do 40 j = 1, ipp
      jb = (j - 1)*npd
      if (j.ge.ipp) npb = npp(m) - (ipp - 1)*npd
c find interpolation weights
      do 10 i = 1, npb
      n = part(1,i+jb,m) + .5
      dxl = part(1,i+jb,m) - float(n)
      n = n + 1
      if (n.gt.nx) n = n - nx
      amx = qm*(.75 - dxl*dxl)
      np = n + 1
      if (np.gt.nx) np = np - nx
      dxp = .5*qm*(.5 + dxl)**2
      nl = n - 1
      if (nl.lt.1) nl = nl + nx
      dxl = .5*qm*(.5 - dxl)**2
      mm = part(2,i+jb,m) + .5
      dyl = part(2,i+jb,m) - float(mm)
      mm = nxv*(mm - mnoff)
      amy = .75 - dyl*dyl
      mp = mm + nxv
      dyp = .5*(.5 + dyl)**2
      ml = mm - nxv
      dyl = .5*(.5 - dyl)**2
      l = part(3,i+jb,m) + .5
      dzl = part(3,i+jb,m) - float(l)
      l = nxvy*(l - lnoff)
      amz = .75 - dzl*dzl
      lp = l + nxvy
      dzp = .5*(.5 + dzl)**2
      lm = l - nxvy
      dzl = .5*(.5 - dzl)**2
      nn(1,i,m) = n + mm + l
      nn(2,i,m) = np + mm + l
      nn(3,i,m) = nl + mm + l
      nn(4,i,m) = n + mp + l
      nn(5,i,m) = np + mp + l
      nn(6,i,m) = nl + mp + l
      nn(7,i,m) = n + ml + l
      nn(8,i,m) = np + ml + l
      nn(9,i,m) = nl + ml + l
      nn(10,i,m) = n + mm + lp
      nn(11,i,m) = np + mm + lp
      nn(12,i,m) = nl + mm + lp
      nn(13,i,m) = n + mp + lp
      nn(14,i,m) = np + mp + lp
      nn(15,i,m) = nl + mp + lp
      nn(16,i,m) = n + ml + lp
      nn(17,i,m) = np + ml + lp
      nn(18,i,m) = nl + ml + lp
      nn(19,i,m) = n + mm + lm
      nn(20,i,m) = np + mm + lm
      nn(21,i,m) = nl + mm + lm
      nn(22,i,m) = n + mp + lm
      nn(23,i,m) = np + mp + lm
      nn(24,i,m) = nl + mp + lm
      nn(25,i,m) = n + ml + lm
      nn(26,i,m) = np + ml + lm
      nn(27,i,m) = nl + ml + lm
      amxyz(1,i,m) = amx*amy*amz
      amxyz(2,i,m) = dxp*amy*amz
      amxyz(3,i,m) = dxl*amy*amz
      amxyz(4,i,m) = amx*dyp*amz
      amxyz(5,i,m) = dxp*dyp*amz
      amxyz(6,i,m) = dxl*dyp*amz
      amxyz(7,i,m) = amx*dyl*amz
      amxyz(8,i,m) = dxp*dyl*amz
      amxyz(9,i,m) = dxl*dyl*amz
      amxyz(10,i,m) = amx*amy*dzp
      amxyz(11,i,m) = dxp*amy*dzp
      amxyz(12,i,m) = dxl*amy*dzp
      amxyz(13,i,m) = amx*dyp*dzp
      amxyz(14,i,m) = dxp*dyp*dzp
      amxyz(15,i,m) = dxl*dyp*dzp
      amxyz(16,i,m) = amx*dyl*dzp
      amxyz(17,i,m) = dxp*dyl*dzp
      amxyz(18,i,m) = dxl*dyl*dzp
      amxyz(19,i,m) = amx*amy*dzl
      amxyz(20,i,m) = dxp*amy*dzl
      amxyz(21,i,m) = dxl*amy*dzl
      amxyz(22,i,m) = amx*dyp*dzl
      amxyz(23,i,m) = dxp*dyp*dzl
      amxyz(24,i,m) = dxl*dyp*dzl
      amxyz(25,i,m) = amx*dyl*dzl
      amxyz(26,i,m) = dxp*dyl*dzl
      amxyz(27,i,m) = dxl*dyl*dzl
   10 continue
c deposit charge
      do 30 i = 1, npb
      do 20 k = 1, 27
      q(nn(k,i,m),m) = q(nn(k,i,m),m) + amxyz(k,i,m)
   20 continue
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSOST32X(part,q,npp,noff,nn,amxyz,qm,idimp,npmax,mnblo
     1k,nxv,nypmx,nxvyzp,idds,npd,n27)
c for 3d code, this subroutine calculates particle charge density
c using second-order spline interpolation,
c with short vectors over independent weights, and distributed data,
c with 2D spatial decomposition
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c vectorized version with guard cells and 1d addressing
c 100 flops/particle, 30 loads, 27 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m,l)=qm*(.75-dx**2)*(.75-dy**2)*(.75-dz**2)
c q(n+1,m,l)=.5*qm*((.5+dx)**2)*(.75-dy**2)*(.75-dz**2)
c q(n-1,m,l)=.5*qm*((.5-dx)**2)*(.75-dy**2)*(.75-dz**2)
c q(n,m+1,l)=.5*qm*(.75-dx**2)*(.5+dy)**2*(.75-dz**2)
c q(n+1,m+1,l)=.25*qm*((.5+dx)**2)*((.5+dy)**2)*(.75-dz**2)
c q(n-1,m+1,l)=.25*qm*((.5-dx)**2)*((.5+dy)**2)*(.75-dz**2)
c q(n,m-1,l)=.5*qm*(.75-dx**2)*((.5-dy)**2)*(.75-dz**2)
c q(n+1,m-1,l)=.25*qm*((.5+dx)**2)*((.5-dy)**2)*(.75-dz**2)
c q(n-1,m-1,l)=.25*qm*((.5-dx)**2)*((.5-dy)**2)*(.75-dz**2)
c q(n,m,l+1)=.5*qm*(.75-dx**2)*(.75-dy**2)*((.5+dz)**2)
c q(n+1,m,l+1)=.25*qm*((.5+dx)**2)*(.75-dy**2)*((.5+dz)**2)
c q(n-1,m,l+1)=.25*qm*((.5-dx)**2)*(.75-dy**2)*((.5+dz)**2)
c q(n,m+1,l+1)=.25*qm*(.75-dx**2)*(.5+dy)**2*((.5+dz)**2)
c q(n+1,m+1,l+1)=.125*qm*((.5+dx)**2)*((.5+dy)**2)*((.5+dz)**2)
c q(n-1,m+1,l+1)=.125*qm*((.5-dx)**2)*((.5+dy)**2)*((.5+dz)**2)
c q(n,m-1,l+1)=.25*qm*(.75-dx**2)*((.5-dy)**2)*((.5+dz)**2)
c q(n+1,m-1,l+1)=.125*qm*((.5+dx)**2)*((.5-dy)**2)*((.5+dz)**2)
c q(n-1,m-1,l+1)=.125*qm*((.5-dx)**2)*((.5-dy)**2)*((.5+dz)**2)
c q(n,m,l-1)=.5*qm*(.75-dx**2)*(.75-dy**2)*((.5-dz)**2)
c q(n+1,m,l-1)=.25*qm*((.5+dx)**2)*(.75-dy**2)*((.5-dz)**2)
c q(n-1,m,l-1)=.25*qm*((.5-dx)**2)*(.75-dy**2)*((.5-dz)**2)
c q(n,m+1,l-1)=.25*qm*(.75-dx**2)*(.5+dy)**2*((.5-dz)**2)
c q(n+1,m+1,l-1)=.125*qm*((.5+dx)**2)*((.5+dy)**2)*((.5-dz)**2)
c q(n-1,m+1,l-1)=.125*qm*((.5-dx)**2)*((.5+dy)**2)*((.5-dz)**2)
c q(n,m-1,l-1)=.25*qm*(.75-dx**2)*((.5-dy)**2)*((.5-dz)**2)
c q(n+1,m-1,l-1)=.125*qm*((.5+dx)**2)*((.5-dy)**2)*((.5-dz)**2)
c q(n-1,m-1,l-1)=.125*qm*((.5-dx)**2)*((.5-dy)**2)*((.5-dz)**2)
c where n,m,l = nearest grid points and dx = x-n, dy = y-m, dz = z-l
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c q(j+1,k,l,m) = charge density at grid point (j,kk,ll),
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nn = scratch address array for vectorized charge deposition
c amxyz = scratch weight array for vectorized charge deposition
c qm = charge on particle, in units of e
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first virtual dimension of charge array, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nxvyzp = nxv*nypmx*nzpmx, first dimension of charge array
c idds = dimensionality of domain decomposition
c npd = size of scratch buffers for vectorized push/charge deposition
c n27 = number of independent weights
      dimension part(idimp,npmax,mnblok), q(nxvyzp,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      dimension nn(n27,npd,mnblok), amxyz(n27,npd,mnblok)
      nxvy = nxv*nypmx
      qmh = .5*qm
c parallel loop
      do 50 m = 1, mnblok
      mnoff = noff(1,m)
      lnoff = noff(2,m)
      npb = npd
      if (npp(m).gt.npd) then
         ipp = float(npp(m) - 1)/float(npd) + 1.
      else
         ipp = 1
      endif
c outer loop over blocks of particles
      do 40 j = 1, ipp
      jb = (j - 1)*npd
      if (j.ge.ipp) npb = npp(m) - (ipp - 1)*npd
c find interpolation weights
      do 10 i = 1, npb
      n = part(1,i+jb,m) + .5
      mm = part(2,i+jb,m) + .5
      ll = part(3,i+jb,m) + .5
      dx = part(1,i+jb,m) - float(n)
      dy = part(2,i+jb,m) - float(mm)
      dz = part(3,i+jb,m) - float(ll)
      n = n + 1
      mm = nxv*(mm - mnoff) + nxvy*(ll - lnoff)
      amx = qm*(.75 - dx*dx)
      dxl = qmh*(.5 - dx)**2
      dxp = qmh*(.5 + dx)**2
      ml = mm + n
      dyl = .5*(.5 - dy)**2
      amy = .75 - dy*dy
      dyp = .5*(.5 + dy)**2
      mn = ml + nxv
      dx1 = dxl*dyl
      dx2 = amx*dyl
      dyl = dxp*dyl
      dzl = .5*(.5 - dz)**2
      dy1 = dxl*amy
      dy2 = amx*amy
      amy = dxp*amy
      amz = .75 - dz*dz
      dxl = dxl*dyp
      amx = amx*dyp
      dyp = dxp*dyp
      dzp = .5*(.5 + dz)**2
      mp = mn + nxv
      nn(1,i,m) = ml
      nn(2,i,m) = ml + 1
      nn(3,i,m) = ml + 2
      nn(4,i,m) = mn
      nn(5,i,m) = mn + 1
      nn(6,i,m) = mn + 2
      nn(7,i,m) = mp
      nn(8,i,m) = mp + 1
      nn(9,i,m) = mp + 2
      ml = ml + nxvy
      mn = mn + nxvy
      mp = mp + nxvy
      nn(10,i,m) = ml
      nn(11,i,m) = ml + 1
      nn(12,i,m) = ml + 2
      nn(13,i,m) = mn
      nn(14,i,m) = mn + 1
      nn(15,i,m) = mn + 2
      nn(16,i,m) = mp
      nn(17,i,m) = mp + 1
      nn(18,i,m) = mp + 2
      ml = ml + nxvy
      mn = mn + nxvy
      mp = mp + nxvy
      nn(19,i,m) = ml
      nn(20,i,m) = ml + 1
      nn(21,i,m) = ml + 2
      nn(22,i,m) = mn
      nn(23,i,m) = mn + 1
      nn(24,i,m) = mn + 2
      nn(25,i,m) = mp
      nn(26,i,m) = mp + 1
      nn(27,i,m) = mp + 2
c deposit charge
      amxyz(1,i,m) = dx1*dzl
      amxyz(2,i,m) = dx2*dzl
      amxyz(3,i,m) = dyl*dzl
      amxyz(4,i,m) = dy1*dzl
      amxyz(5,i,m) = dy2*dzl
      amxyz(6,i,m) = amy*dzl
      amxyz(7,i,m) = dxl*dzl
      amxyz(8,i,m) = amx*dzl
      amxyz(9,i,m) = dyp*dzl
      amxyz(10,i,m) = dx1*amz
      amxyz(11,i,m) = dx2*amz
      amxyz(12,i,m) = dyl*amz
      amxyz(13,i,m) = dy1*amz
      amxyz(14,i,m) = dy2*amz
      amxyz(15,i,m) = amy*amz
      amxyz(16,i,m) = dxl*amz
      amxyz(17,i,m) = amx*amz
      amxyz(18,i,m) = dyp*amz
      amxyz(19,i,m) = dx1*dzp
      amxyz(20,i,m) = dx2*dzp
      amxyz(21,i,m) = dyl*dzp
      amxyz(22,i,m) = dy1*dzp
      amxyz(23,i,m) = dy2*dzp
      amxyz(24,i,m) = amy*dzp
      amxyz(25,i,m) = dxl*dzp
      amxyz(26,i,m) = amx*dzp
      amxyz(27,i,m) = dyp*dzp
   10 continue
c deposit charge
      do 30 i = 1, npb
cdir$ ivdep
      do 20 k = 1, 27
      q(nn(k,i,m),m) = q(nn(k,i,m),m) + amxyz(k,i,m)
   20 continue
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDOST32L(part,q,npp,noff,qm,nx,idimp,npmax,mnblok,nxv,n
     1ypmx,nzpmx,idds)
c for 3d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries
c and distributed data with 2D spatial decomposition
c baseline scalar distributed version
c 37 flops/particle, 11 loads, 8 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m,l)=qm*(1.-dx)*(1.-dy)*(1.-dz)
c q(n+1,m,l)=qm*dx*(1.-dy)*(1.-dz)
c q(n,m+1,l)=qm*(1.-dx)*dy*(1.-dz)
c q(n+1,m+1,l)=qm*dx*dy*(1.-dz)
c q(n,m,l+1)=qm*(1.-dx)*(1.-dy)*dz
c q(n+1,m,l+1)=qm*dx*(1.-dy)*dz
c q(n,m+1,l+1)=qm*(1.-dx)*dy*dz
c q(n+1,m+1,l+1)=qm*dx*dy*dz
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c q(j,k,l,m) = charge density at grid point (j,kk,ll),
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qm = charge on particle, in units of e
c nx = system length in x direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first dimension of charge array, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
      dimension part(idimp,npmax,mnblok), q(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      do 20 m = 1, mnblok
      mnoff = noff(1,m) - 1
      lnoff = noff(2,m) - 1
c find interpolation weights
      do 10 j = 1, npp(m)
      nn = part(1,j,m)
      mm = part(2,j,m)
      ll = part(3,j,m)
      dxp = qm*(part(1,j,m) - float(nn))
      dyp = part(2,j,m) - float(mm)
      dzp = part(3,j,m) - float(ll)
      nn = nn + 1
      amx = qm - dxp
      np = nn + 1
      if (np.gt.nx) np = np - nx
      mm = mm - mnoff
      amy = 1. - dyp
      mp = mm + 1
      ll = ll - lnoff
      amz = 1. - dzp
      lp = ll + 1
c deposit charge
      q(nn,mm,ll,m) = q(nn,mm,ll,m) + amx*amy*amz
      q(np,mm,ll,m) = q(np,mm,ll,m) + dxp*amy*amz
      q(nn,mp,ll,m) = q(nn,mp,ll,m) + amx*dyp*amz
      q(np,mp,ll,m) = q(np,mp,ll,m) + dxp*dyp*amz
      q(nn,mm,lp,m) = q(nn,mm,lp,m) + amx*amy*dzp
      q(np,mm,lp,m) = q(np,mm,lp,m) + dxp*amy*dzp
      q(nn,mp,lp,m) = q(nn,mp,lp,m) + amx*dyp*dzp
      q(np,mp,lp,m) = q(np,mp,lp,m) + dxp*dyp*dzp
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGPOST32L(part,q,npp,noff,qm,idimp,npmax,mnblok,nxv,nyp
     1mx,nzpmx,idds)
c for 3d code, this subroutine calculates particle charge density
c using first-order linear interpolation, and distributed data
c with 2D spatial decomposition
c scalar version using guard cells, for distributed data
c 33 flops/particle, 11 loads, 8 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m,l)=qm*(1.-dx)*(1.-dy)*(1.-dz)
c q(n+1,m,l)=qm*dx*(1.-dy)*(1.-dz)
c q(n,m+1,l)=qm*(1.-dx)*dy*(1.-dz)
c q(n+1,m+1,l)=qm*dx*dy*(1.-dz)
c q(n,m,l+1)=qm*(1.-dx)*(1.-dy)*dz
c q(n+1,m,l+1)=qm*dx*(1.-dy)*dz
c q(n,m+1,l+1)=qm*(1.-dx)*dy*dz
c q(n+1,m+1,l+1)=qm*dx*dy*dz
c where n,m,l = nearest grid points and dx = x-n, dy = y-m, dz = z-l
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c q(j,k,l,m) = charge density at grid point (j,kk,ll),
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qm = charge on particle, in units of e
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first dimension of charge array, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
      dimension part(idimp,npmax,mnblok), q(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      do 20 m = 1, mnblok
      mnoff = noff(1,m) - 1
      lnoff = noff(2,m) - 1
c find interpolation weights
      do 10 j = 1, npp(m)
      nn = part(1,j,m)
      mm = part(2,j,m)
      ll = part(3,j,m)
      dxp = qm*(part(1,j,m) - float(nn))
      dyp = part(2,j,m) - float(mm)
      dzp = part(3,j,m) - float(ll)
      nn = nn + 1
      amx = qm - dxp
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
c deposit charge
      q(nn,mm,ll,m) = q(nn,mm,ll,m) + amx*amz
      q(np,mm,ll,m) = q(np,mm,ll,m) + amy*amz
      q(nn,mp,ll,m) = q(nn,mp,ll,m) + dyp*amz
      q(np,mp,ll,m) = q(np,mp,ll,m) + dx1*amz
      q(nn,mm,lp,m) = q(nn,mm,lp,m) + amx*dzp
      q(np,mm,lp,m) = q(np,mm,lp,m) + amy*dzp
      q(nn,mp,lp,m) = q(nn,mp,lp,m) + dyp*dzp
      q(np,mp,lp,m) = q(np,mp,lp,m) + dx1*dzp
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSPOST32L(part,q,npp,noff,qm,idimp,npmax,mnblok,nxv,ny
     1pmx,nxyzp,idds)
c for 3d code, this subroutine calculates particle charge density
c using first-order linear interpolation, and distributed data
c with 2D spatial decomposition
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 33 flops/particle, 11 loads, 8 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m,l)=qm*(1.-dx)*(1.-dy)*(1.-dz)
c q(n+1,m,l)=qm*dx*(1.-dy)*(1.-dz)
c q(n,m+1,l)=qm*(1.-dx)*dy*(1.-dz)
c q(n+1,m+1,l)=qm*dx*dy*(1.-dz)
c q(n,m,l+1)=qm*(1.-dx)*(1.-dy)*dz
c q(n+1,m,l+1)=qm*dx*(1.-dy)*dz
c q(n,m+1,l+1)=qm*(1.-dx)*dy*dz
c q(n+1,m+1,l+1)=qm*dx*dy*dz
c where n,m,l = nearest grid points and dx = x-n, dy = y-m, dz = z-l
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c q(j,k,l,m) = charge density at grid point (j,kk,ll),
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qm = charge on particle, in units of e
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first virtual dimension of charge array, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nxyzp = dimension of charge array, must be >= nxv*nypmx*nzpmx
c idds = dimensionality of domain decomposition
      dimension part(idimp,npmax,mnblok), q(nxyzp,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      nxvy = nxv*nypmx
      do 20 m = 1, mnblok
      if (npp(m).lt.1) go to 20
      mnoff = noff(1,m)
      lnoff = noff(2,m)
c begin first particle
      nnn = part(1,1,m)
      mmm = part(2,1,m)
      lll = part(3,1,m)
      dxn = part(1,1,m) - float(nnn)
      dyn = part(2,1,m) - float(mmm)
      dzn = part(3,1,m) - float(lll)
      mmm = mmm - mnoff
      lll = lll - lnoff
c find interpolation weights
      do 10 j = 2, npp(m)
      nn = nnn + 1
      mm = nxv*mmm + nxvy*lll
      nnn = part(1,j,m)
      mmm = part(2,j,m)
      lll = part(3,j,m)
      dxp = qm*dxn
      dyp = dyn
      dzp = dzn
      dxn = part(1,j,m) - float(nnn)
      dyn = part(2,j,m) - float(mmm)
      dzn = part(3,j,m) - float(lll)
      amx = qm - dxp
      amy = 1. - dyp
      mm = mm + nn
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + nxv
      amx = amx*amy
      amz = 1. - dzp
      ll = mm + nxvy
      amy = dxp*amy
      lp = mp + nxvy
      mmm = mmm - mnoff
      lll = lll - lnoff
c deposit charge
      dxp = q(mm,m) + amx*amz
      dx2 = q(mm+1,m) + amy*amz
      dx3 = q(mp,m) + dyp*amz
      amz = q(mp+1,m) + dx1*amz
      amx = q(ll,m) + amx*dzp
      amy = q(ll+1,m) + amy*dzp
      dyp = q(lp,m) + dyp*dzp
      dzp = q(lp+1,m) + dx1*dzp
      q(mm,m) = dxp
      q(mm+1,m) = dx2
      q(mp,m) = dx3
      q(mp+1,m) = amz
      q(ll,m) = amx
      q(ll+1,m) = amy
      q(lp,m) = dyp
      q(lp+1,m) = dzp
   10 continue
c deposit charge for last particle
      nn = nnn + 1
      mm = nxv*mmm + nxvy*lll
      dxp = qm*dxn
      amx = qm - dxp
      amy = 1. - dyn
      mm = mm + nn
      dx1 = dxp*dyn
      dyp = amx*dyn
      mp = mm + nxv
      amx = amx*amy
      amz = 1. - dzn
      ll = mm + nxvy
      amy = dxp*amy
      lp = mp + nxvy
c deposit charge
      q(mm,m) = q(mm,m) + amx*amz
      q(mm+1,m) = q(mm+1,m) + amy*amz
      q(mp,m) = q(mp,m) + dyp*amz
      q(mp+1,m) = q(mp+1,m) + dx1*amz
      q(ll,m) = q(ll,m) + amx*dzn
      q(ll+1,m) = q(ll+1,m) + amy*dzn
      q(lp,m) = q(lp,m) + dyp*dzn
      q(lp+1,m) = q(lp+1,m) + dx1*dzn
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSOST32XL(part,q,npp,noff,nn,amxyz,qm,nx,idimp,npmax,mn
     1blok,nxv,nypmx,nxvyzp,idds,npd,ieight)
c for 3d code, this subroutine calculates particle charge density
c using first-order linear interpolation, periodic boundaries,
c with 2D spatial decomposition
c with short vectors over independent weights, and distributed data,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c vectorized distributed version with 1d addressing
c 37 flops/particle, 11 loads, 8 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m,l)=qm*(1.-dx)*(1.-dy)*(1.-dz)
c q(n+1,m,l)=qm*dx*(1.-dy)*(1.-dz)
c q(n,m+1,l)=qm*(1.-dx)*dy*(1.-dz)
c q(n+1,m+1,l)=qm*dx*dy*(1.-dz)
c q(n,m,l+1)=qm*(1.-dx)*(1.-dy)*dz
c q(n+1,m,l+1)=qm*dx*(1.-dy)*dz
c q(n,m+1,l+1)=qm*(1.-dx)*dy*dz
c q(n+1,m+1,l+1)=qm*dx*dy*dz
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c q(j,k,l,m) = charge density at grid point (j,kk,ll),
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nn = scratch address array for vectorized charge deposition
c amxyz = scratch weight array for vectorized charge deposition
c qm = charge on particle, in units of e
c nx = system length in x direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first virtual dimension of charge array, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nxvyzp = nxv*nypmx*nzpmx, first dimension of charge array
c idds = dimensionality of domain decomposition
c npd = size of scratch buffers for vectorized push/charge deposition
c ieight = number of independent weights
      dimension part(idimp,npmax,mnblok), q(nxvyzp,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      dimension nn(ieight,npd,mnblok), amxyz(ieight,npd,mnblok)
      nxvy = nxv*nypmx
c parallel loop
      do 50 m = 1, mnblok
      mnoff = noff(1,m)
      lnoff = noff(2,m)
      npb = npd
      if (npp(m).gt.npd) then
         ipp = float(npp(m) - 1)/float(npd) + 1.
      else
         ipp = 1
      endif
c outer loop over blocks of particles
      do 40 j = 1, ipp
      jb = (j - 1)*npd
      if (j.ge.ipp) npb = npp(m) - (ipp - 1)*npd
c find interpolation weights
      do 10 i = 1, npb
      n = part(1,i+jb,m)
      dxp = qm*(part(1,i+jb,m) - float(n))
      n = n + 1
      amx = qm - dxp
      np = n + 1
      if (np.gt.nx) np = np - nx
      mm = part(2,i+jb,m)
      dyp = part(2,i+jb,m) - float(mm)
      mm = nxv*(mm - mnoff)
      amy = 1. - dyp
      mp = mm + nxv
      l = part(3,i+jb,m)
      dzp = part(3,i+jb,m) - float(l)
      l = nxvy*(l - lnoff)
      amz = 1. - dzp
      lp = l + nxvy
      nn(1,i,m) = n + mm + l
      nn(2,i,m) = np + mm + l
      nn(3,i,m) = n + mp + l
      nn(4,i,m) = np + mp + l
      nn(5,i,m) = n + mm + lp
      nn(6,i,m) = np + mm + lp
      nn(7,i,m) = n + mp + lp
      nn(8,i,m) = np + mp + lp
      amxyz(1,i,m) = amx*amy*amz
      amxyz(2,i,m) = dxp*amy*amz
      amxyz(3,i,m) = amx*dyp*amz
      amxyz(4,i,m) = dxp*dyp*amz
      amxyz(5,i,m) = amx*amy*dzp
      amxyz(6,i,m) = dxp*amy*dzp
      amxyz(7,i,m) = amx*dyp*dzp
      amxyz(8,i,m) = dxp*dyp*dzp
   10 continue
c deposit charge
      do 30 i = 1, npb
cdir$ ivdep
      do 20 k = 1, 8
      q(nn(k,i,m),m) = q(nn(k,i,m),m) + amxyz(k,i,m)
   20 continue
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSOST32XL(part,q,npp,noff,nn,amxyz,qm,idimp,npmax,mnbl
     1ok,nxv,nypmx,nxvyzp,idds,npd,ieight)
c for 3d code, this subroutine calculates particle charge density
c using first-order linear interpolation,
c with short vectors over independent weights, and distributed data,
c with 2D spatial decomposition
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c vectorized distributed version with guard cells and 1d addressing
c 33 flops/particle, 11 loads, 8 stores
c input: all, output: q
c charge density is approximated by values at the nearest grid points
c q(n,m,l)=qm*(1.-dx)*(1.-dy)*(1.-dz)
c q(n+1,m,l)=qm*dx*(1.-dy)*(1.-dz)
c q(n,m+1,l)=qm*(1.-dx)*dy*(1.-dz)
c q(n+1,m+1,l)=qm*dx*dy*(1.-dz)
c q(n,m,l+1)=qm*(1.-dx)*(1.-dy)*dz
c q(n+1,m,l+1)=qm*dx*(1.-dy)*dz
c q(n,m+1,l+1)=qm*(1.-dx)*dy*dz
c q(n+1,m+1,l+1)=qm*dx*dy*dz
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c q(n,m) = charge density at grid point (j,kk,ll),
c where n = j + nxv*(kk - 1) + nxv*nypmx*(ll - 1), 
c kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nn = scratch address array for vectorized charge deposition
c amxyz = scratch weight array for vectorized charge deposition
c qm = charge on particle, in units of e
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first virtual dimension of charge array, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nxvyzp = nxv*nypmx*nzpmx, first dimension of charge array
c idds = dimensionality of domain decomposition
c npd = size of scratch buffers for vectorized push/charge deposition
c ieight = number of independent weights
      dimension part(idimp,npmax,mnblok), q(nxvyzp,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      dimension nn(ieight,npd,mnblok), amxyz(ieight,npd,mnblok)
      nxvy = nxv*nypmx
c parallel loop
      do 50 m = 1, mnblok
      mnoff = noff(1,m)
      lnoff = noff(2,m)
      npb = npd
      if (npp(m).gt.npd) then
         ipp = float(npp(m) - 1)/float(npd) + 1.
      else
         ipp = 1
      endif
c outer loop over blocks of particles
      do 40 j = 1, ipp
      jb = (j - 1)*npd
      if (j.ge.ipp) npb = npp(m) - (ipp - 1)*npd
c find interpolation weights
      do 10 i = 1, npb
      n = part(1,i+jb,m)
      mm = part(2,i+jb,m)
      ll = part(3,i+jb,m)
      dxp = qm*(part(1,i+jb,m) - float(n))
      dyp = part(2,i+jb,m) - float(mm)
      dzp = part(3,i+jb,m) - float(ll)
      n = n + 1
      mm = nxv*(mm - mnoff) + nxvy*(ll - lnoff)
      amx = qm - dxp
      amy = 1. - dyp
      mm = mm + n
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + nxv
      amx = amx*amy
      amz = 1. - dzp
      ll = mm + nxvy
      amy = dxp*amy
      lp = mp + nxvy
      nn(1,i,m) = mm
      nn(2,i,m) = mm + 1
      nn(3,i,m) = mp
      nn(4,i,m) = mp + 1
      nn(5,i,m) = ll
      nn(6,i,m) = ll + 1
      nn(7,i,m) = lp
      nn(8,i,m) = lp + 1
      amxyz(1,i,m) = amx*amz
      amxyz(2,i,m) = amy*amz
      amxyz(3,i,m) = dyp*amz
      amxyz(4,i,m) = dx1*amz
      amxyz(5,i,m) = amx*dzp
      amxyz(6,i,m) = amy*dzp
      amxyz(7,i,m) = dyp*dzp
      amxyz(8,i,m) = dx1*dzp
   10 continue
c deposit charge
      do 30 i = 1, npb
cdir$ ivdep
      do 20 k = 1, 8
      q(nn(k,i,m),m) = q(nn(k,i,m),m) + amxyz(k,i,m)
   20 continue
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPUSH32(part,fx,fy,fz,npp,noff,qbm,dt,ek,nx,idimp,npmax
     1,mnblok,nxv,nypmx,nzpmx,idds)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with periodic boundary conditions,
c for distributed data with 2D spatial decomposition
c baseline scalar distributed version
c 254 flops/particle, 87 loads, 6 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
c z(t+dt) = z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (.75-dz**2)*((.75-dy**2)*((.75-dx**2)*fx(n,m,l)+
c (.5*(.5+dx)**2)*fx(n+1,m,l)+(.5*(.5-dx)**2)*fx(n-1,m,l)) +
c (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1,l)+
c (.5*(.5+dx)**2)*fx(n+1,m+1,l)+(.5*(.5-dx)**2)*fx(n-1,m+1,l)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1,l)+
c (.5*(.5+dx)**2)*fx(n+1,m-1,l)+(.5*(.5-dx)**2)*fx(n-1,m-1,l))) +
c (.5*(.5+dz)**2)*((.75-dy**2)*((.75-dx**2)*fx(n,m,l+1)+
c (.5*(.5+dx)**2)*fx(n+1,m,l+1)+(.5*(.5-dx)**2)*fx(n-1,m,l+1)) +
c (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1,l+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1,l+1)+(.5*(.5-dx)**2)*fx(n-1,m+1,l+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1,l+1)+
c (.5*(.5+dx)**2)*fx(n+1,m-1,l+1)+(.5*(.5-dx)**2)*fx(n-1,m-1,l+1)))
c (.5*(.5-dz)**2)*((.75-dy**2)*((.75-dx**2)*fx(n,m,l-1)+
c (.5*(.5+dx)**2)*fx(n+1,m,l-1)+(.5*(.5-dx)**2)*fx(n-1,m,l-1)) +
c (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1,l-1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1,l-1)+(.5*(.5-dx)**2)*fx(n-1,m+1,l-1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1,l-1)+
c (.5*(.5+dx)**2)*fx(n+1,m-1,l-1)+(.5*(.5-dx)**2)*fx(n-1,m-1,l-1)))
c where n,m,l = nearest grid points and dx = x-n, dy = y-m, dz = z-l
c and similarly for fy(x,y,z) and fz(x,y,z)
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c fx(j,k,l,m) = x component of force/charge at grid (j,kk,ll)
c fy(j,k,l,m) = y component of force/charge at grid (j,kk,ll)
c fz(j,k,l,m) = z component of force/charge at grid (j,kk,ll)
c in other words, fx/fy/fz are the convolutions of the electric field
c over the particle shape,
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
c (vz(t+dt/2)+vz(t-dt/2))**2)
c nx = system length in x direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
      double precision sum1
      dimension part(idimp,npmax,mnblok)
      dimension fx(nxv,nypmx,nzpmx,mnblok), fy(nxv,nypmx,nzpmx,mnblok)
      dimension fz(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      zero = 0.
      anx = float(nx)
      qtm = qbm*dt
      sum1 = 0.0d0
      do 20 m = 1, mnblok
      mnoff = noff(1,m) - 1
      lnoff = noff(2,m) - 1
c find interpolation weights
      do 10 j = 1, npp(m)
      nn = part(1,j,m) + .5
      mm = part(2,j,m) + .5
      ll = part(3,j,m) + .5
      dx = part(1,j,m) - float(nn)
      dy = part(2,j,m) - float(mm)
      dz = part(3,j,m) - float(ll)
      nl = nn
      if (nl.lt.1) nl = nl + nx
      amx = .75 - dx*dx
      ml = mm - mnoff
      dxl = .5*(.5 - dx)**2
      lm = ll - lnoff
      dxp = .5*(.5 + dx)**2
      nn = nn + 1
      if (nn.gt.nx) nn = nn - nx
      dyl = .5*(.5 - dy)**2
      np = nn + 1
      if (np.gt.nx) np = np - nx
      amy = .75 - dy*dy
      mm = ml + 1
      dyp = .5*(.5 + dy)**2
      mp = mm + 1
      dzl = .5*(.5 - dz)**2
      ll = lm + 1
      amz = .75 - dz*dz
      lp = ll + 1
      dzp = .5*(.5 + dz)**2
c find acceleration
      dx = dzl*(dyl*(dxl*fx(nl,ml,lm,m) + amx*fx(nn,ml,lm,m) + dxp*fx(np
     1,ml,lm,m)) + amy*(dxl*fx(nl,mm,lm,m) + amx*fx(nn,mm,lm,m) + dxp*fx
     2(np,mm,lm,m)) + dyp*(dxl*fx(nl,mp,lm,m) + amx*fx(nn,mp,lm,m) + dxp
     3*fx(np,mp,lm,m))) + amz*(dyl*(dxl*fx(nl,ml,ll,m) + amx*fx(nn,ml,ll
     4,m) + dxp*fx(np,ml,ll,m)) + amy*(dxl*fx(nl,mm,ll,m) + amx*fx(nn,mm
     5,ll,m) + dxp*fx(np,mm,ll,m)) + dyp*(dxl*fx(nl,mp,ll,m) + amx*fx(nn
     6,mp,ll,m) + dxp*fx(np,mp,ll,m))) + dzp*(dyl*(dxl*fx(nl,ml,lp,m) + 
     7amx*fx(nn,ml,lp,m) + dxp*fx(np,ml,lp,m)) + amy*(dxl*fx(nl,mm,lp,m) 
     8 + amx*fx(nn,mm,lp,m) + dxp*fx(np,mm,lp,m)) + dyp*(dxl*fx(nl,mp,lp
     9,m) + amx*fx(nn,mp,lp,m) + dxp*fx(np,mp,lp,m)))
      dy = dzl*(dyl*(dxl*fy(nl,ml,lm,m) + amx*fy(nn,ml,lm,m) + dxp*fy(np
     1,ml,lm,m)) + amy*(dxl*fy(nl,mm,lm,m) + amx*fy(nn,mm,lm,m) + dxp*fy
     2(np,mm,lm,m)) + dyp*(dxl*fy(nl,mp,lm,m) + amx*fy(nn,mp,lm,m) + dxp
     3*fy(np,mp,lm,m))) + amz*(dyl*(dxl*fy(nl,ml,ll,m) + amx*fy(nn,ml,ll
     4,m) + dxp*fy(np,ml,ll,m)) + amy*(dxl*fy(nl,mm,ll,m) + amx*fy(nn,mm
     5,ll,m) + dxp*fy(np,mm,ll,m)) + dyp*(dxl*fy(nl,mp,ll,m) + amx*fy(nn
     6,mp,ll,m) + dxp*fy(np,mp,ll,m))) + dzp*(dyl*(dxl*fy(nl,ml,lp,m) + 
     7amx*fy(nn,ml,lp,m) + dxp*fy(np,ml,lp,m)) + amy*(dxl*fy(nl,mm,lp,m) 
     8 + amx*fy(nn,mm,lp,m) + dxp*fy(np,mm,lp,m)) + dyp*(dxl*fy(nl,mp,lp
     9,m) + amx*fy(nn,mp,lp,m) + dxp*fy(np,mp,lp,m)))
      dz = dzl*(dyl*(dxl*fz(nl,ml,lm,m) + amx*fz(nn,ml,lm,m) + dxp*fz(np
     1,ml,lm,m)) + amy*(dxl*fz(nl,mm,lm,m) + amx*fz(nn,mm,lm,m) + dxp*fz
     2(np,mm,lm,m)) + dyp*(dxl*fz(nl,mp,lm,m) + amx*fz(nn,mp,lm,m) + dxp
     3*fz(np,mp,lm,m))) + amz*(dyl*(dxl*fz(nl,ml,ll,m) + amx*fz(nn,ml,ll
     4,m) + dxp*fz(np,ml,ll,m)) + amy*(dxl*fz(nl,mm,ll,m) + amx*fz(nn,mm
     5,ll,m) + dxp*fz(np,mm,ll,m)) + dyp*(dxl*fz(nl,mp,ll,m) + amx*fz(nn
     6,mp,ll,m) + dxp*fz(np,mp,ll,m))) + dzp*(dyl*(dxl*fz(nl,ml,lp,m) + 
     7amx*fz(nn,ml,lp,m) + dxp*fz(np,ml,lp,m)) + amy*(dxl*fz(nl,mm,lp,m) 
     8 + amx*fz(nn,mm,lp,m) + dxp*fz(np,mm,lp,m)) + dyp*(dxl*fz(nl,mp,lp
     9,m) + amx*fz(nn,mp,lp,m) + dxp*fz(np,mp,lp,m)))
c new velocity
      dx = part(4,j,m) + qtm*dx
      dy = part(5,j,m) + qtm*dy
      dz = part(6,j,m) + qtm*dz
c average kinetic energy
      sum1 = sum1 + (dx + part(4,j,m))**2 + (dy + part(5,j,m))**2 + (dz 
     1+ part(6,j,m))**2
      part(4,j,m) = dx
      part(5,j,m) = dy
      part(6,j,m) = dz
c new position
      dx = part(1,j,m) + dx*dt
      dy = part(2,j,m) + dy*dt
      dz = part(3,j,m) + dz*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j,m) = dx
      part(2,j,m) = dy
      part(3,j,m) = dz
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGPUSH32(part,fxyz,npp,noff,qbm,dt,ek,nx,ny,nz,idimp,np
     1max,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, for distributed data
c with 2D spatial decomposition
c scalar version using guard cells,
c 238 flops/particle, 87 loads, 6 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
c z(t+dt) = z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (.75-dz**2)*((.75-dy**2)*((.75-dx**2)*fx(n,m,l)+
c (.5*(.5+dx)**2)*fx(n+1,m,l)+(.5*(.5-dx)**2)*fx(n-1,m,l)) +
c (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1,l)+
c (.5*(.5+dx)**2)*fx(n+1,m+1,l)+(.5*(.5-dx)**2)*fx(n-1,m+1,l)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1,l)+
c (.5*(.5+dx)**2)*fx(n+1,m-1,l)+(.5*(.5-dx)**2)*fx(n-1,m-1,l))) +
c (.5*(.5+dz)**2)*((.75-dy**2)*((.75-dx**2)*fx(n,m,l+1)+
c (.5*(.5+dx)**2)*fx(n+1,m,l+1)+(.5*(.5-dx)**2)*fx(n-1,m,l+1)) +
c (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1,l+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1,l+1)+(.5*(.5-dx)**2)*fx(n-1,m+1,l+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1,l+1)+
c (.5*(.5+dx)**2)*fx(n+1,m-1,l+1)+(.5*(.5-dx)**2)*fx(n-1,m-1,l+1)))
c (.5*(.5-dz)**2)*((.75-dy**2)*((.75-dx**2)*fx(n,m,l-1)+
c (.5*(.5+dx)**2)*fx(n+1,m,l-1)+(.5*(.5-dx)**2)*fx(n-1,m,l-1)) +
c (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1,l-1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1,l-1)+(.5*(.5-dx)**2)*fx(n-1,m+1,l-1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1,l-1)+
c (.5*(.5+dx)**2)*fx(n+1,m-1,l-1)+(.5*(.5-dx)**2)*fx(n-1,m-1,l-1)))
c where n,m,l = nearest grid points and dx = x-n, dy = y-m, dz = z-l
c and similarly for fy(x,y,z) and fz(x,y,z)
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c fxyz(1,j+1,k,l,m) = x component of force/charge at grid (j,kk,ll)
c fxyz(2,j+1,k,l,m) = y component of force/charge at grid (j,kk,ll)
c fxyz(3,j+1,k,l,m) = z component of force/charge at grid (j,kk,ll)
c in other words, fxyz are the convolutions of the electric field
c over the particle shape,
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
c (vz(t+dt/2)+vz(t-dt/2))**2)
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first dimension of field array, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      double precision sum1
      dimension part(idimp,npmax,mnblok)
      dimension fxyz(3,nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      qtm = qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
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
         edgerx = float(nx-1)
         edgery = float(ny-1)
      endif
      do 20 m = 1, mnblok
      mnoff = noff(1,m) - 1
      lnoff = noff(2,m) - 1
c find interpolation weights
      do 10 j = 1, npp(m)
      nn = part(1,j,m) + .5
      mm = part(2,j,m) + .5
      ll = part(3,j,m) + .5
      dx = part(1,j,m) - float(nn)
      dy = part(2,j,m) - float(mm)
      dz = part(3,j,m) - float(ll)
      nl = nn + 1
      amx = .75 - dx*dx
      ml = mm - mnoff
      dxl = .5*(.5 - dx)**2
      lm = ll - lnoff
      dxp = .5*(.5 + dx)**2
      nn = nl + 1
      dyl = .5*(.5 - dy)**2
      np = nl + 2
      amy = .75 - dy*dy
      mm = ml + 1
      dyp = .5*(.5 + dy)**2
      mp = ml + 2
      dx1 = dxl*dyl
      dx2 = amx*dyl
      dyl = dxp*dyl
      dzl = .5*(.5 - dz)**2
      ll = lm + 1
      dy1 = dxl*amy
      dy2 = amx*amy
      amy = dxp*amy
      amz = .75 - dz*dz
      lp = lm + 2
      dxl = dxl*dyp
      amx = amx*dyp
      dyp = dxp*dyp
      dzp = .5*(.5 + dz)**2
c find acceleration
      dx = dzl*(dx1*fxyz(1,nl,ml,lm,m) + dx2*fxyz(1,nn,ml,lm,m) + dyl*fx
     1yz(1,np,ml,lm,m) + dy1*fxyz(1,nl,mm,lm,m) + dy2*fxyz(1,nn,mm,lm,m)
     2 + amy*fxyz(1,np,mm,lm,m) + dxl*fxyz(1,nl,mp,lm,m) + amx*fxyz(1,nn
     3,mp,lm,m) + dyp*fxyz(1,np,mp,lm,m))
      dy = dzl*(dx1*fxyz(2,nl,ml,lm,m) + dx2*fxyz(2,nn,ml,lm,m) + dyl*fx
     1yz(2,np,ml,lm,m) + dy1*fxyz(2,nl,mm,lm,m) + dy2*fxyz(2,nn,mm,lm,m)
     2 + amy*fxyz(2,np,mm,lm,m) + dxl*fxyz(2,nl,mp,lm,m) + amx*fxyz(2,nn
     3,mp,lm,m) + dyp*fxyz(2,np,mp,lm,m))
      dz = dzl*(dx1*fxyz(3,nl,ml,lm,m) + dx2*fxyz(3,nn,ml,lm,m) + dyl*fx
     1yz(3,np,ml,lm,m) + dy1*fxyz(3,nl,mm,lm,m) + dy2*fxyz(3,nn,mm,lm,m)
     2 + amy*fxyz(3,np,mm,lm,m) + dxl*fxyz(3,nl,mp,lm,m) + amx*fxyz(3,nn
     3,mp,lm,m) + dyp*fxyz(3,np,mp,lm,m))
      dx = dx + amz*(dx1*fxyz(1,nl,ml,ll,m) + dx2*fxyz(1,nn,ml,ll,m) + d
     1yl*fxyz(1,np,ml,ll,m) + dy1*fxyz(1,nl,mm,ll,m) + dy2*fxyz(1,nn,mm,
     2ll,m) + amy*fxyz(1,np,mm,ll,m) + dxl*fxyz(1,nl,mp,ll,m) + amx*fxyz
     3(1,nn,mp,ll,m) + dyp*fxyz(1,np,mp,ll,m))
      dy = dy + amz*(dx1*fxyz(2,nl,ml,ll,m) + dx2*fxyz(2,nn,ml,ll,m) + d
     1yl*fxyz(2,np,ml,ll,m) + dy1*fxyz(2,nl,mm,ll,m) + dy2*fxyz(2,nn,mm,
     2ll,m) + amy*fxyz(2,np,mm,ll,m) + dxl*fxyz(2,nl,mp,ll,m) + amx*fxyz
     3(2,nn,mp,ll,m) + dyp*fxyz(2,np,mp,ll,m))
      dz = dz + amz*(dx1*fxyz(3,nl,ml,ll,m) + dx2*fxyz(3,nn,ml,ll,m) + d
     1yl*fxyz(3,np,ml,ll,m) + dy1*fxyz(3,nl,mm,ll,m) + dy2*fxyz(3,nn,mm,
     2ll,m) + amy*fxyz(3,np,mm,ll,m) + dxl*fxyz(3,nl,mp,ll,m) + amx*fxyz
     3(3,nn,mp,ll,m) + dyp*fxyz(3,np,mp,ll,m))
      dx = dx + dzp*(dx1*fxyz(1,nl,ml,lp,m) + dx2*fxyz(1,nn,ml,lp,m) + d
     1yl*fxyz(1,np,ml,lp,m) + dy1*fxyz(1,nl,mm,lp,m) + dy2*fxyz(1,nn,mm,
     2lp,m) + amy*fxyz(1,np,mm,lp,m) + dxl*fxyz(1,nl,mp,lp,m) + amx*fxyz
     3(1,nn,mp,lp,m) + dyp*fxyz(1,np,mp,lp,m))
      dy = dy + dzp*(dx1*fxyz(2,nl,ml,lp,m) + dx2*fxyz(2,nn,ml,lp,m) + d
     1yl*fxyz(2,np,ml,lp,m) + dy1*fxyz(2,nl,mm,lp,m) + dy2*fxyz(2,nn,mm,
     2lp,m) + amy*fxyz(2,np,mm,lp,m) + dxl*fxyz(2,nl,mp,lp,m) + amx*fxyz
     3(2,nn,mp,lp,m) + dyp*fxyz(2,np,mp,lp,m))
      dz = dz + dzp*(dx1*fxyz(3,nl,ml,lp,m) + dx2*fxyz(3,nn,ml,lp,m) + d
     1yl*fxyz(3,np,ml,lp,m) + dy1*fxyz(3,nl,mm,lp,m) + dy2*fxyz(3,nn,mm,
     2lp,m) + amy*fxyz(3,np,mm,lp,m) + dxl*fxyz(3,nl,mp,lp,m) + amx*fxyz
     3(3,nn,mp,lp,m) + dyp*fxyz(3,np,mp,lp,m))
c new velocity
      dx = part(4,j,m) + qtm*dx
      dy = part(5,j,m) + qtm*dy
      dz = part(6,j,m) + qtm*dz
c average kinetic energy
      sum1 = sum1 + (dx + part(4,j,m))**2 + (dy + part(5,j,m))**2 + (dz 
     1+ part(6,j,m))**2
      part(4,j,m) = dx
      part(5,j,m) = dy
      part(6,j,m) = dz
c new position
      dx = part(1,j,m) + dx*dt
      dy = part(2,j,m) + dy*dt
      dz = part(3,j,m) + dz*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,m)
            part(4,j,m) = -part(4,j,m)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,m)
            part(5,j,m) = -part(5,j,m)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = part(3,j,m)
            part(6,j,m) = -part(6,j,m)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,m)
            part(4,j,m) = -part(4,j,m)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,m)
            part(5,j,m) = -part(5,j,m)
         endif
      endif
c set new position
      part(1,j,m) = dx
      part(2,j,m) = dy
      part(3,j,m) = dz
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSPUSH32(part,fxyz,npp,noff,qbm,dt,ek,nx,ny,nz,idimp,n
     1pmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, for distributed data
c with 2D spatial decomposition
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 238 flops/particle, 87 loads, 6 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
c z(t+dt) = z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (.75-dz**2)*((.75-dy**2)*((.75-dx**2)*fx(n,m,l)+
c (.5*(.5+dx)**2)*fx(n+1,m,l)+(.5*(.5-dx)**2)*fx(n-1,m,l)) +
c (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1,l)+
c (.5*(.5+dx)**2)*fx(n+1,m+1,l)+(.5*(.5-dx)**2)*fx(n-1,m+1,l)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1,l)+
c (.5*(.5+dx)**2)*fx(n+1,m-1,l)+(.5*(.5-dx)**2)*fx(n-1,m-1,l))) +
c (.5*(.5+dz)**2)*((.75-dy**2)*((.75-dx**2)*fx(n,m,l+1)+
c (.5*(.5+dx)**2)*fx(n+1,m,l+1)+(.5*(.5-dx)**2)*fx(n-1,m,l+1)) +
c (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1,l+1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1,l+1)+(.5*(.5-dx)**2)*fx(n-1,m+1,l+1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1,l+1)+
c (.5*(.5+dx)**2)*fx(n+1,m-1,l+1)+(.5*(.5-dx)**2)*fx(n-1,m-1,l+1)))
c (.5*(.5-dz)**2)*((.75-dy**2)*((.75-dx**2)*fx(n,m,l-1)+
c (.5*(.5+dx)**2)*fx(n+1,m,l-1)+(.5*(.5-dx)**2)*fx(n-1,m,l-1)) +
c (.5*(.5+dy)**2)*((.75-dx**2)*fx(n,m+1,l-1)+
c (.5*(.5+dx)**2)*fx(n+1,m+1,l-1)+(.5*(.5-dx)**2)*fx(n-1,m+1,l-1)) +
c (.5*(.5-dy)**2)*((.75-dx**2)*fx(n,m-1,l-1)+
c (.5*(.5+dx)**2)*fx(n+1,m-1,l-1)+(.5*(.5-dx)**2)*fx(n-1,m-1,l-1)))
c where n,m,l = nearest grid points and dx = x-n, dy = y-m, dz = z-l
c and similarly for fy(x,y,z) and fz(x,y,z)
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c fxyz(1,j+1,k,l,m) = x component of force/charge at grid (j,kk,ll)
c fxyz(2,j+1,k,l,m) = y component of force/charge at grid (j,kk,ll)
c fxyz(3,j+1,k,l,m) = z component of force/charge at grid (j,kk,ll)
c in other words, fxyz are the convolutions of the electric field
c over the particle shape,
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
c (vz(t+dt/2)+vz(t-dt/2))**2)
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first virtual dimension of field array, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nxyzp = second dimension of field array, must be >= nxv*nypmx*nzpmx
c idds = dimensionality of domain decomposition
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      double precision sum1
      dimension part(idimp,npmax,mnblok)
      dimension fxyz(3,nxyzp,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      qtm = qbm*dt
      sum1 = 0.0d0
      nxvy = nxv*nypmx
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
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
         edgerx = float(nx-1)
         edgery = float(ny-1)
      endif
      do 20 m = 1, mnblok
      if (npp(m).lt.1) go to 20
      mnoff = noff(1,m)
      lnoff = noff(2,m)
c begin first particle
      nnn = part(1,1,m) + .5
      mmm = part(2,1,m) + .5
      lll = part(3,1,m) + .5
      dxn = part(1,1,m) - float(nnn)
      dyn = part(2,1,m) - float(mmm)
      dzn = part(3,1,m) - float(lll)
      mmm = mmm - mnoff
      lll = lll - lnoff
      nop1 = npp(m) - 1
c find interpolation weights
      do 10 j = 1, nop1
      nn = nnn + 1
      mm = nxv*mmm + nxvy*lll
      nnn = part(1,j+1,m) + .5
      mmm = part(2,j+1,m) + .5
      lll = part(3,j+1,m) + .5
      dx = dxn
      dy = dyn
      dz = dzn
      dxn = part(1,j+1,m) - float(nnn)
      dyn = part(2,j+1,m) - float(mmm)
      dzn = part(3,j+1,m) - float(lll)
      amx = .75 - dx*dx
      dxl = .5*(.5 - dx)**2
      dxp = .5*(.5 + dx)**2
      ml = mm + nn
      dyl = .5*(.5 - dy)**2
      amy = .75 - dy*dy
      dyp = .5*(.5 + dy)**2
      mn = ml + nxv
      dx1 = dxl*dyl
      dx2 = amx*dyl
      dyl = dxp*dyl
      dzl = .5*(.5 - dz)**2
      dy1 = dxl*amy
      dy2 = amx*amy
      amy = dxp*amy
      amz = .75 - dz*dz
      dxl = dxl*dyp
      amx = amx*dyp
      dyp = dxp*dyp
      dzp = .5*(.5 + dz)**2
      mp = mn + nxv
      mmm = mmm - mnoff
      lll = lll - lnoff
c find acceleration
      dx = dzl*(dx1*fxyz(1,ml,m) + dx2*fxyz(1,ml+1,m) + dyl*fxyz(1,ml+2,
     1m) + dy1*fxyz(1,mn,m) + dy2*fxyz(1,mn+1,m) + amy*fxyz(1,mn+2,m) + 
     2dxl*fxyz(1,mp,m) + amx*fxyz(1,mp+1,m) + dyp*fxyz(1,mp+2,m))
      dy = dzl*(dx1*fxyz(2,ml,m) + dx2*fxyz(2,ml+1,m) + dyl*fxyz(2,ml+2,
     1m) + dy1*fxyz(2,mn,m) + dy2*fxyz(2,mn+1,m) + amy*fxyz(2,mn+2,m) + 
     2dxl*fxyz(2,mp,m) + amx*fxyz(2,mp+1,m) + dyp*fxyz(2,mp+2,m))
      dz = dzl*(dx1*fxyz(3,ml,m) + dx2*fxyz(3,ml+1,m) + dyl*fxyz(3,ml+2,
     1m) + dy1*fxyz(3,mn,m) + dy2*fxyz(3,mn+1,m) + amy*fxyz(3,mn+2,m) + 
     2dxl*fxyz(3,mp,m) + amx*fxyz(3,mp+1,m) + dyp*fxyz(3,mp+2,m))
      ml = ml + nxvy
      mn = mn + nxvy
      mp = mp + nxvy
      dx = dx + amz*(dx1*fxyz(1,ml,m) + dx2*fxyz(1,ml+1,m) + dyl*fxyz(1,
     1ml+2,m) + dy1*fxyz(1,mn,m) + dy2*fxyz(1,mn+1,m) + amy*fxyz(1,mn+2,
     2m) + dxl*fxyz(1,mp,m) + amx*fxyz(1,mp+1,m) + dyp*fxyz(1,mp+2,m))
      dy = dy + amz*(dx1*fxyz(2,ml,m) + dx2*fxyz(2,ml+1,m) + dyl*fxyz(2,
     1ml+2,m) + dy1*fxyz(2,mn,m) + dy2*fxyz(2,mn+1,m) + amy*fxyz(2,mn+2,
     2m) + dxl*fxyz(2,mp,m) + amx*fxyz(2,mp+1,m) + dyp*fxyz(2,mp+2,m))
      dz = dz + amz*(dx1*fxyz(3,ml,m) + dx2*fxyz(3,ml+1,m) + dyl*fxyz(3,
     1ml+2,m) + dy1*fxyz(3,mn,m) + dy2*fxyz(3,mn+1,m) + amy*fxyz(3,mn+2,
     2m) + dxl*fxyz(3,mp,m) + amx*fxyz(3,mp+1,m) + dyp*fxyz(3,mp+2,m))
      ml = ml + nxvy
      mn = mn + nxvy
      mp = mp + nxvy
      dx = dx + dzp*(dx1*fxyz(1,ml,m) + dx2*fxyz(1,ml+1,m) + dyl*fxyz(1,
     1ml+2,m) + dy1*fxyz(1,mn,m) + dy2*fxyz(1,mn+1,m) + amy*fxyz(1,mn+2,
     2m) + dxl*fxyz(1,mp,m) + amx*fxyz(1,mp+1,m) + dyp*fxyz(1,mp+2,m))
      dy = dy + dzp*(dx1*fxyz(2,ml,m) + dx2*fxyz(2,ml+1,m) + dyl*fxyz(2,
     1ml+2,m) + dy1*fxyz(2,mn,m) + dy2*fxyz(2,mn+1,m) + amy*fxyz(2,mn+2,
     2m) + dxl*fxyz(2,mp,m) + amx*fxyz(2,mp+1,m) + dyp*fxyz(2,mp+2,m))
      dz = dz + dzp*(dx1*fxyz(3,ml,m) + dx2*fxyz(3,ml+1,m) + dyl*fxyz(3,
     1ml+2,m) + dy1*fxyz(3,mn,m) + dy2*fxyz(3,mn+1,m) + amy*fxyz(3,mn+2,
     2m) + dxl*fxyz(3,mp,m) + amx*fxyz(3,mp+1,m) + dyp*fxyz(3,mp+2,m))
c new velocity
      dx = part(4,j,m) + qtm*dx
      dy = part(5,j,m) + qtm*dy
      dz = part(6,j,m) + qtm*dz
c average kinetic energy
      sum1 = sum1 + (dx + part(4,j,m))**2 + (dy + part(5,j,m))**2 + (dz 
     1+ part(6,j,m))**2
      part(4,j,m) = dx
      part(5,j,m) = dy
      part(6,j,m) = dz
c new position
      dx = part(1,j,m) + dx*dt
      dy = part(2,j,m) + dy*dt
      dz = part(3,j,m) + dz*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,m)
            part(4,j,m) = -part(4,j,m)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,m)
            part(5,j,m) = -part(5,j,m)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = part(3,j,m)
            part(6,j,m) = -part(6,j,m)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,m)
            part(4,j,m) = -part(4,j,m)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,m)
            part(5,j,m) = -part(5,j,m)
         endif
      endif
c set new position
      part(1,j,m) = dx
      part(2,j,m) = dy
      part(3,j,m) = dz
   10 continue
      nop = npp(m)
c push last particle
      nn = nnn + 1
      mm = nxv*mmm + nxvy*lll
      amx = .75 - dxn*dxn
      dxl = .5*(.5 - dxn)**2
      dxp = .5*(.5 + dxn)**2
      ml = mm + nn
      dyl = .5*(.5 - dyn)**2
      amy = .75 - dyn*dyn
      dyp = .5*(.5 + dyn)**2
      mn = ml + nxv
      dx1 = dxl*dyl
      dx2 = amx*dyl
      dyl = dxp*dyl
      dzl = .5*(.5 - dzn)**2
      dy1 = dxl*amy
      dy2 = amx*amy
      amy = dxp*amy
      amz = .75 - dzn*dzn
      dxl = dxl*dyp
      amx = amx*dyp
      dyp = dxp*dyp
      dzp = .5*(.5 + dzn)**2
      mp = mn + nxv
c find acceleration
      dx = dzl*(dx1*fxyz(1,ml,m) + dx2*fxyz(1,ml+1,m) + dyl*fxyz(1,ml+2,
     1m) + dy1*fxyz(1,mn,m) + dy2*fxyz(1,mn+1,m) + amy*fxyz(1,mn+2,m) + 
     2dxl*fxyz(1,mp,m) + amx*fxyz(1,mp+1,m) + dyp*fxyz(1,mp+2,m))
      dy = dzl*(dx1*fxyz(2,ml,m) + dx2*fxyz(2,ml+1,m) + dyl*fxyz(2,ml+2,
     1m) + dy1*fxyz(2,mn,m) + dy2*fxyz(2,mn+1,m) + amy*fxyz(2,mn+2,m) + 
     2dxl*fxyz(2,mp,m) + amx*fxyz(2,mp+1,m) + dyp*fxyz(2,mp+2,m))
      dz = dzl*(dx1*fxyz(3,ml,m) + dx2*fxyz(3,ml+1,m) + dyl*fxyz(3,ml+2,
     1m) + dy1*fxyz(3,mn,m) + dy2*fxyz(3,mn+1,m) + amy*fxyz(3,mn+2,m) + 
     2dxl*fxyz(3,mp,m) + amx*fxyz(3,mp+1,m) + dyp*fxyz(3,mp+2,m))
      ml = ml + nxvy
      mn = mn + nxvy
      mp = mp + nxvy
      dx = dx + amz*(dx1*fxyz(1,ml,m) + dx2*fxyz(1,ml+1,m) + dyl*fxyz(1,
     1ml+2,m) + dy1*fxyz(1,mn,m) + dy2*fxyz(1,mn+1,m) + amy*fxyz(1,mn+2,
     2m) + dxl*fxyz(1,mp,m) + amx*fxyz(1,mp+1,m) + dyp*fxyz(1,mp+2,m))
      dy = dy + amz*(dx1*fxyz(2,ml,m) + dx2*fxyz(2,ml+1,m) + dyl*fxyz(2,
     1ml+2,m) + dy1*fxyz(2,mn,m) + dy2*fxyz(2,mn+1,m) + amy*fxyz(2,mn+2,
     2m) + dxl*fxyz(2,mp,m) + amx*fxyz(2,mp+1,m) + dyp*fxyz(2,mp+2,m))
      dz = dz + amz*(dx1*fxyz(3,ml,m) + dx2*fxyz(3,ml+1,m) + dyl*fxyz(3,
     1ml+2,m) + dy1*fxyz(3,mn,m) + dy2*fxyz(3,mn+1,m) + amy*fxyz(3,mn+2,
     2m) + dxl*fxyz(3,mp,m) + amx*fxyz(3,mp+1,m) + dyp*fxyz(3,mp+2,m))
      ml = ml + nxvy
      mn = mn + nxvy
      mp = mp + nxvy
      dx = dx + dzp*(dx1*fxyz(1,ml,m) + dx2*fxyz(1,ml+1,m) + dyl*fxyz(1,
     1ml+2,m) + dy1*fxyz(1,mn,m) + dy2*fxyz(1,mn+1,m) + amy*fxyz(1,mn+2,
     2m) + dxl*fxyz(1,mp,m) + amx*fxyz(1,mp+1,m) + dyp*fxyz(1,mp+2,m))
      dy = dy + dzp*(dx1*fxyz(2,ml,m) + dx2*fxyz(2,ml+1,m) + dyl*fxyz(2,
     1ml+2,m) + dy1*fxyz(2,mn,m) + dy2*fxyz(2,mn+1,m) + amy*fxyz(2,mn+2,
     2m) + dxl*fxyz(2,mp,m) + amx*fxyz(2,mp+1,m) + dyp*fxyz(2,mp+2,m))
      dz = dz + dzp*(dx1*fxyz(3,ml,m) + dx2*fxyz(3,ml+1,m) + dyl*fxyz(3,
     1ml+2,m) + dy1*fxyz(3,mn,m) + dy2*fxyz(3,mn+1,m) + amy*fxyz(3,mn+2,
     2m) + dxl*fxyz(3,mp,m) + amx*fxyz(3,mp+1,m) + dyp*fxyz(3,mp+2,m))
c new velocity
      dx = part(4,nop,m) + qtm*dx
      dy = part(5,nop,m) + qtm*dy
      dz = part(6,nop,m) + qtm*dz
c average kinetic energy
      sum1 = sum1 + (dx + part(4,nop,m))**2 + (dy + part(5,nop,m))**2 + 
     1(dz + part(6,nop,m))**2
      part(4,nop,m) = dx
      part(5,nop,m) = dy
      part(6,nop,m) = dz
c new position
      dx = part(1,nop,m) + dx*dt
      dy = part(2,nop,m) + dy*dt
      dz = part(3,nop,m) + dz*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,m)
            part(4,nop,m) = -part(4,nop,m)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop,m)
            part(5,nop,m) = -part(5,nop,m)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = part(3,nop,m)
            part(6,nop,m) = -part(6,nop,m)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,m)
            part(4,nop,m) = -part(4,nop,m)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop,m)
            part(5,nop,m) = -part(5,nop,m)
         endif
      endif
c set new position
      part(1,nop,m) = dx
      part(2,nop,m) = dy
      part(3,nop,m) = dz
   20 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PPUSH32L(part,fx,fy,fz,npp,noff,qbm,dt,ek,nx,idimp,npma
     1x,mnblok,nxv,nypmx,nzpmx,idds)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with periodic boundary conditions,
c for distributed data with 2D spatial decomposition
c baseline scalar distributed version
c 98 flops/particle, 30 loads, 6 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
c z(t+dt) = z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
c                + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
c                + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
c fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
c                + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
c                + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c fx(j,k,l,m) = x component of force/charge at grid (j,kk,ll)
c fy(j,k,l,m) = y component of force/charge at grid (j,kk,ll)
c fz(j,k,l,m) = z component of force/charge at grid (j,kk,ll)
c in other words, fx/fy/fz are the convolutions of the electric field
c over the particle shape,
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
c (vz(t+dt/2)+vz(t-dt/2))**2)
c nx = system length in x direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first dimension of field arrays, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
      double precision sum1
      dimension part(idimp,npmax,mnblok)
      dimension fx(nxv,nypmx,nzpmx,mnblok), fy(nxv,nypmx,nzpmx,mnblok)
      dimension fz(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      zero = 0.
      anx = float(nx)
      qtm = qbm*dt
      sum1 = 0.0d0
      do 20 m = 1, mnblok
      mnoff = noff(1,m) - 1
      lnoff = noff(2,m) - 1
c find interpolation weights
      do 10 j = 1, npp(m)
      nn = part(1,j,m)
      mm = part(2,j,m)
      ll = part(3,j,m)
      dxp = part(1,j,m) - float(nn)
      dyp = part(2,j,m) - float(mm)
      dzp = part(3,j,m) - float(ll)
      nn = nn + 1
      amx = 1. - dxp
      np = nn + 1
      if (np.gt.nx) np = np - nx
      mm = mm - mnoff
      amy = 1. - dyp
      mp = mm + 1
      ll = ll - lnoff
      amz = 1. - dzp
      lp = ll + 1
c find acceleration
      dx = amz*(amy*(amx*fx(nn,mm,ll,m) + dxp*fx(np,mm,ll,m)) + dyp*(amx
     1*fx(nn,mp,ll,m) + dxp*fx(np,mp,ll,m))) + dzp*(amy*(amx*fx(nn,mm,lp
     2,m) + dxp*fx(np,mm,lp,m)) + dyp*(amx*fx(nn,mp,lp,m) + dxp*fx(np,mp
     3,lp,m)))
      dy = amz*(amy*(amx*fy(nn,mm,ll,m) + dxp*fy(np,mm,ll,m)) + dyp*(amx
     1*fy(nn,mp,ll,m) + dxp*fy(np,mp,ll,m))) + dzp*(amy*(amx*fy(nn,mm,lp
     2,m) + dxp*fy(np,mm,lp,m)) + dyp*(amx*fy(nn,mp,lp,m) + dxp*fy(np,mp
     3,lp,m)))
      dz = amz*(amy*(amx*fz(nn,mm,ll,m) + dxp*fz(np,mm,ll,m)) + dyp*(amx
     1*fz(nn,mp,ll,m) + dxp*fz(np,mp,ll,m))) + dzp*(amy*(amx*fz(nn,mm,lp
     2,m) + dxp*fz(np,mm,lp,m)) + dyp*(amx*fz(nn,mp,lp,m) + dxp*fz(np,mp
     3,lp,m)))
c new velocity
      dx = part(4,j,m) + qtm*dx
      dy = part(5,j,m) + qtm*dy
      dz = part(6,j,m) + qtm*dz
c average kinetic energy
      sum1 = sum1 + (dx + part(4,j,m))**2 + (dy + part(5,j,m))**2 + (dz
     1+ part(6,j,m))**2
      part(4,j,m) = dx
      part(5,j,m) = dy
      part(6,j,m) = dz
c new position
      dx = part(1,j,m) + dx*dt
      dy = part(2,j,m) + dy*dt
      dz = part(3,j,m) + dz*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j,m) = dx
      part(2,j,m) = dy
      part(3,j,m) = dz
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGPUSH32L(part,fxyz,npp,noff,qbm,dt,ek,nx,ny,nz,idimp,n
     1pmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for distributed data
c with 2D spatial decomposition
c scalar version using guard cells
c 90 flops/particle, 30 loads, 6 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
c z(t+dt) = z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
c                + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
c                + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
c fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
c                + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
c                + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c fxyz(1,j,k,l,m) = x component of force/charge at grid (j,kk,ll)
c fxyz(2,j,k,l,m) = y component of force/charge at grid (j,kk,ll)
c fxyz(3,j,k,l,m) = z component of force/charge at grid (j,kk,ll)
c in other words, fxyz are the convolutions of the electric field
c over the particle shape,
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
c (vz(t+dt/2)+vz(t-dt/2))**2)
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first dimension of field array, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      double precision sum1
      dimension part(idimp,npmax,mnblok)
      dimension fxyz(3,nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      qtm = qbm*dt
      sum1 = 0.0d0
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
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
         edgerx = float(nx-1)
         edgery = float(ny-1)
      endif
      do 20 m = 1, mnblok
      mnoff = noff(1,m) - 1
      lnoff = noff(2,m) - 1
c find interpolation weights
      do 10 j = 1, npp(m)
      nn = part(1,j,m)
      mm = part(2,j,m)
      ll = part(3,j,m)
      dxp = part(1,j,m) - float(nn)
      dyp = part(2,j,m) - float(mm)
      dzp = part(3,j,m) - float(ll)
      nn = nn + 1
      amx = 1. - dxp
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
c find acceleration
      dx = amz*(amx*fxyz(1,nn,mm,ll,m) + amy*fxyz(1,np,mm,ll,m) + dyp*fx
     1yz(1,nn,mp,ll,m) + dx1*fxyz(1,np,mp,ll,m)) + dzp*(amx*fxyz(1,nn,mm
     2,lp,m) + amy*fxyz(1,np,mm,lp,m) + dyp*fxyz(1,nn,mp,lp,m) + dx1*fxy
     3z(1,np,mp,lp,m))
      dy = amz*(amx*fxyz(2,nn,mm,ll,m) + amy*fxyz(2,np,mm,ll,m) + dyp*fx
     1yz(2,nn,mp,ll,m) + dx1*fxyz(2,np,mp,ll,m)) + dzp*(amx*fxyz(2,nn,mm
     2,lp,m) + amy*fxyz(2,np,mm,lp,m) + dyp*fxyz(2,nn,mp,lp,m) + dx1*fxy
     3z(2,np,mp,lp,m))
      dz = amz*(amx*fxyz(3,nn,mm,ll,m) + amy*fxyz(3,np,mm,ll,m) + dyp*fx
     1yz(3,nn,mp,ll,m) + dx1*fxyz(3,np,mp,ll,m)) + dzp*(amx*fxyz(3,nn,mm
     2,lp,m) + amy*fxyz(3,np,mm,lp,m) + dyp*fxyz(3,nn,mp,lp,m) + dx1*fxy
     3z(3,np,mp,lp,m))
c new velocity
      dx = part(4,j,m) + qtm*dx
      dy = part(5,j,m) + qtm*dy
      dz = part(6,j,m) + qtm*dz
c average kinetic energy
      sum1 = sum1 + (dx + part(4,j,m))**2 + (dy + part(5,j,m))**2 + (dz
     1+ part(6,j,m))**2
      part(4,j,m) = dx
      part(5,j,m) = dy
      part(6,j,m) = dz
c new position
      dx = part(1,j,m) + dx*dt
      dy = part(2,j,m) + dy*dt
      dz = part(3,j,m) + dz*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,m)
            part(4,j,m) = -part(4,j,m)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,m)
            part(5,j,m) = -part(5,j,m)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = part(3,j,m)
            part(6,j,m) = -part(6,j,m)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,m)
            part(4,j,m) = -part(4,j,m)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,m)
            part(5,j,m) = -part(5,j,m)
         endif
      endif
c set new position
      part(1,j,m) = dx
      part(2,j,m) = dy
      part(3,j,m) = dz
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSPUSH32L(part,fxyz,npp,noff,qbm,dt,ek,nx,ny,nz,idimp,
     1npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for distributed data
c with 2D spatial decomposition
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 90 flops/particle, 30 loads, 6 stores
c input: all, output: part, ek
c equations used are:
c vx(t+dt/2) = vx(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c vy(t+dt/2) = vy(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c vz(t+dt/2) = vz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + vx(t+dt/2)*dt, y(t+dt) = y(t) + vy(t+dt/2)*dt,
c z(t+dt) = z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c fy(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fy(n,m,l)+dx*fy(n+1,m,l))
c                + dy*((1-dx)*fy(n,m+1,l) + dx*fy(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fy(n,m,l+1)+dx*fy(n+1,m,l+1))
c                + dy*((1-dx)*fy(n,m+1,l+1) + dx*fy(n+1,m+1,l+1)))
c fz(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fz(n,m,l)+dx*fz(n+1,m,l))
c                + dy*((1-dx)*fz(n,m+1,l) + dx*fz(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fz(n,m,l+1)+dx*fz(n+1,m,l+1))
c                + dy*((1-dx)*fz(n,m+1,l+1) + dx*fz(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c fxyz(1,j,k,l,m) = x component of force/charge at grid (j,kk,ll)
c fxyz(2,j,k,l,m) = y component of force/charge at grid (j,kk,ll)
c fxyz(3,j,k,l,m) = z component of force/charge at grid (j,kk,ll)
c in other words, fxyz are the convolutions of the electric field
c over the particle shape,
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .125*sum((vx(t+dt/2)+vx(t-dt/2))**2+(vy(t+dt/2)+vy(t-dt/2))**2+
c (vz(t+dt/2)+vz(t-dt/2))**2)
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first virtual dimension of field array, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nxyzp = dimension of field array, must be >= nxv*nypmx*nzpmx
c idds = dimensionality of domain decomposition
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      double precision sum1
      dimension part(idimp,npmax,mnblok)
      dimension fxyz(3,nxyzp,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      qtm = qbm*dt
      sum1 = 0.0d0
      nxvy = nxv*nypmx
c set boundary values
      if (ipbc.eq.1) then
         edgelx = 0.
         edgerx = float(nx)
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
         edgerx = float(nx-1)
         edgery = float(ny-1)
      endif
      do 20 m = 1, mnblok
      if (npp(m).lt.1) go to 20
      mnoff = noff(1,m)
      lnoff = noff(2,m)
c begin first particle
      nnn = part(1,1,m)
      mmm = part(2,1,m)
      lll = part(3,1,m)
      dxn = part(1,1,m) - float(nnn)
      dyn = part(2,1,m) - float(mmm)
      dzn = part(3,1,m) - float(lll)
      mmm = mmm - mnoff
      lll = lll - lnoff
      nop1 = npp(m) - 1
      do 10 j = 1, nop1
c find interpolation weights
      nn = nnn + 1
      mm = nxv*mmm + nxvy*lll
      nnn = part(1,j+1,m)
      mmm = part(2,j+1,m)
      lll = part(3,j+1,m)
      dxp = dxn
      dyp = dyn
      dzp = dzn
      dxn = part(1,j+1,m) - float(nnn)
      dyn = part(2,j+1,m) - float(mmm)
      dzn = part(3,j+1,m) - float(lll)
      amx = 1. - dxp
      amy = 1. - dyp
      mm = mm + nn
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + nxv
      amx = amx*amy
      amz = 1. - dzp
      ll = mm + nxvy
      amy = dxp*amy
      lp = mp + nxvy
      mmm = mmm - mnoff
      lll = lll - lnoff
c find acceleration
      dx = amz*(amx*fxyz(1,mm,m) + amy*fxyz(1,mm+1,m) + dyp*fxyz(1,mp,m)
     1 + dx1*fxyz(1,mp+1,m)) + dzp*(amx*fxyz(1,ll,m) + amy*fxyz(1,ll+1,m
     2) + dyp*fxyz(1,lp,m) + dx1*fxyz(1,lp+1,m))
      dy = amz*(amx*fxyz(2,mm,m) + amy*fxyz(2,mm+1,m) + dyp*fxyz(2,mp,m)
     1 + dx1*fxyz(2,mp+1,m)) + dzp*(amx*fxyz(2,ll,m) + amy*fxyz(2,ll+1,m
     2) + dyp*fxyz(2,lp,m) + dx1*fxyz(2,lp+1,m))
      dz = amz*(amx*fxyz(3,mm,m) + amy*fxyz(3,mm+1,m) + dyp*fxyz(3,mp,m)
     1 + dx1*fxyz(3,mp+1,m)) + dzp*(amx*fxyz(3,ll,m) + amy*fxyz(3,ll+1,m
     2) + dyp*fxyz(3,lp,m) + dx1*fxyz(3,lp+1,m))
c new velocity
      dx = part(4,j,m) + qtm*dx
      dy = part(5,j,m) + qtm*dy
      dz = part(6,j,m) + qtm*dz
c average kinetic energy
      sum1 = sum1 + (dx + part(4,j,m))**2 + (dy + part(5,j,m))**2 + (dz
     1+ part(6,j,m))**2
      part(4,j,m) = dx
      part(5,j,m) = dy
      part(6,j,m) = dz
c new position
      dx = part(1,j,m) + dx*dt
      dy = part(2,j,m) + dy*dt
      dz = part(3,j,m) + dz*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,m)
            part(4,j,m) = -part(4,j,m)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,m)
            part(5,j,m) = -part(5,j,m)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = part(3,j,m)
            part(6,j,m) = -part(6,j,m)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j,m)
            part(4,j,m) = -part(4,j,m)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j,m)
            part(5,j,m) = -part(5,j,m)
         endif
      endif
c set new position
      part(1,j,m) = dx
      part(2,j,m) = dy
      part(3,j,m) = dz
   10 continue
      nop = npp(m)
c push last particle
      nn = nnn + 1
      mm = nxv*mmm + nxvy*lll
      amx = 1. - dxn
      amy = 1. - dyn
      mm = mm + nn
      dx1 = dxn*dyn
      dyp = amx*dyn
      mp = mm + nxv
      amx = amx*amy
      amz = 1. - dzn
      ll = mm + nxvy
      amy = dxn*amy
      lp = mp + nxvy
c find acceleration
      dx = amz*(amx*fxyz(1,mm,m) + amy*fxyz(1,mm+1,m) + dyp*fxyz(1,mp,m)
     1 + dx1*fxyz(1,mp+1,m)) + dzn*(amx*fxyz(1,ll,m) + amy*fxyz(1,ll+1,m
     2) + dyp*fxyz(1,lp,m) + dx1*fxyz(1,lp+1,m))
      dy = amz*(amx*fxyz(2,mm,m) + amy*fxyz(2,mm+1,m) + dyp*fxyz(2,mp,m)
     1 + dx1*fxyz(2,mp+1,m)) + dzn*(amx*fxyz(2,ll,m) + amy*fxyz(2,ll+1,m
     2) + dyp*fxyz(2,lp,m) + dx1*fxyz(2,lp+1,m))
      dz = amz*(amx*fxyz(3,mm,m) + amy*fxyz(3,mm+1,m) + dyp*fxyz(3,mp,m)
     1 + dx1*fxyz(3,mp+1,m)) + dzn*(amx*fxyz(3,ll,m) + amy*fxyz(3,ll+1,m
     2) + dyp*fxyz(3,lp,m) + dx1*fxyz(3,lp+1,m))
c new velocity
      dx = part(4,nop,m) + qtm*dx
      dy = part(5,nop,m) + qtm*dy
      dz = part(6,nop,m) + qtm*dz
c average kinetic energy
      sum1 = sum1 + (dx + part(4,nop,m))**2 + (dy + part(5,nop,m))**2 + 
     1(dz+ part(6,nop,m))**2
      part(4,nop,m) = dx
      part(5,nop,m) = dy
      part(6,nop,m) = dz
c new position
      dx = part(1,nop,m) + dx*dt
      dy = part(2,nop,m) + dy*dt
      dz = part(3,nop,m) + dz*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,m)
            part(4,nop,m) = -part(4,nop,m)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop,m)
            part(5,nop,m) = -part(5,nop,m)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = part(3,nop,m)
            part(6,nop,m) = -part(6,nop,m)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,nop,m)
            part(4,nop,m) = -part(4,nop,m)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,nop,m)
            part(5,nop,m) = -part(5,nop,m)
         endif
      endif
c set new position
      part(1,nop,m) = dx
      part(2,nop,m) = dy
      part(3,nop,m) = dz
   20 continue
c normalize kinetic energy
      ek = ek + .125*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PSORTP32YZ(part,pt,ip,npic,npp,noff,nyzp,idimp,npmax,mn
     1blok,nyzpm1,idds)
c this subroutine sorts particles by y,z grid
c quadratic interpolation, spatial decomposition in y and z direction
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c pt = scratch array for reordering particles
c ip = index array for reordering particles
c npic = address offset for reordering particles
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions
c nyzpm1 = max(nyzp(1,m)+1)*max(nyzp(2,m)+1)
c idds = dimensionality of domain decomposition
      dimension part(idimp,npmax,mnblok), pt(npmax,mnblok)
      dimension ip(npmax,mnblok), npic(nyzpm1,mnblok)
      dimension npp(mnblok), noff(idds,mnblok), nyzp(idds,mnblok)
      do 80 m = 1, mnblok
      nnoff = noff(1,m)
      lnoff = noff(2,m)
      nyp1 = nyzp(1,m) + 1
      nyzp1 = nyp1*(nyzp(2,m) + 1)
c clear counter array
      do 10 k = 1, nyzp1
      npic(k,m) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, npp(m)
      n = part(2,j,m) + 1.5
      l = part(3,j,m) + .5
      l = n - nnoff + nyp1*(l - lnoff)
      npic(l,m) = npic(l,m) + 1
      ip(j,m) = l
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nyzp1
      ist = npic(k,m)
      npic(k,m) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid
      do 40 j = 1, npp(m)
      l = ip(j,m)
      npic(l,m) = npic(l,m) + 1
      ip(j,m) = npic(l,m)
   40 continue
c reorder particles by copying each component to temporary
      do 70 i = 1, idimp
      do 50 j = 1, npp(m)
      pt(ip(j,m),m) = part(i,j,m)
   50 continue
      do 60 j = 1, npp(m)
      part(i,j,m) = pt(j,m)
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSORTP32YZL(part,pt,ip,npic,npp,noff,nyzp,idimp,npmax,m
     1nblok,nyzpm1,idds)
c this subroutine sorts particles by y,z grid
c linear interpolation, spatial decomposition in y and z direction
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c pt = scratch array for reordering particles
c ip = index array for reordering particles
c npic = address offset for reordering particles
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions
c nyzpm1 = max(nyzp(1,m)+1)*max(nyzp(2,m)+1)
c idds = dimensionality of domain decomposition
      dimension part(idimp,npmax,mnblok), pt(npmax,mnblok)
      dimension ip(npmax,mnblok), npic(nyzpm1,mnblok)
      dimension npp(mnblok), noff(idds,mnblok), nyzp(idds,mnblok)
      do 80 m = 1, mnblok
      nnoff = noff(1,m)
      lnoff = noff(2,m)
      nyp1 = nyzp(1,m) + 1
      nyzp1 = nyp1*(nyzp(2,m) + 1)
c clear counter array
      do 10 k = 1, nyzp1
      npic(k,m) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, npp(m)
      n = part(2,j,m) + 1.0
      l = part(3,j,m)
      l = n - nnoff + nyp1*(l - lnoff)
      npic(l,m) = npic(l,m) + 1
      ip(j,m) = l
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nyzp1
      ist = npic(k,m)
      npic(k,m) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid
      do 40 j = 1, npp(m)
      l = ip(j,m)
      npic(l,m) = npic(l,m) + 1
      ip(j,m) = npic(l,m)
   40 continue
c reorder particles by copying each component to temporary
      do 70 i = 1, idimp
      do 50 j = 1, npp(m)
      pt(ip(j,m),m) = part(i,j,m)
   50 continue
      do 60 j = 1, npp(m)
      part(i,j,m) = pt(j,m)
   60 continue
   70 continue
   80 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDSORTP32YZ(parta,partb,npic,npp,noff,nyzp,idimp,npmax,
     ,mnblok,nyzpm1,idds)
c this subroutine sorts particles by y,z grid
c quadratic interpolation, spatial decomposition in y and z direction
c parta/partb = input/output particle array
c parta(2,n,m) = position y of particle n in partition m
c parta(3,n,m) = position z of particle n in partition m
c npic = address offset for reordering particles
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions
c nyzpm1 = max(nyzp(1,m)+1)*max(nyzp(2,m)+1)
c idds = dimensionality of domain decomposition
      implicit none
      integer idimp, npmax, mnblok, nyzpm1, idds
      integer npic, npp, noff, nyzp
      real parta, partb
      dimension parta(idimp,npmax,mnblok), partb(idimp,npmax,mnblok)
      dimension npic(nyzpm1,mnblok)
      dimension npp(mnblok), noff(idds,mnblok), nyzp(idds,mnblok)
c local data
      integer i, j, k, l, m, n, nnoff, lnoff, nyp1, nyzp1, isum, ist, ip
      do 60 m = 1, mnblok
      nnoff = noff(1,m)
      lnoff = noff(2,m)
      nyp1 = nyzp(1,m) + 1
      nyzp1 = nyp1*(nyzp(2,m) + 1)
c clear counter array
      do 10 k = 1, nyzp1
      npic(k,m) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, npp(m)
      n = parta(2,j,m) + 1.5
      l = parta(3,j,m) + .5
      l = n - nnoff + nyp1*(l - lnoff)
      npic(l,m) = npic(l,m) + 1
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nyzp1
      ist = npic(k,m)
      npic(k,m) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid and reorder particles
      do 50 j = 1, npp(m)
      n = parta(2,j,m) + 1.5
      l = parta(3,j,m) + .5
      l = n - nnoff + nyp1*(l - lnoff)
      ip = npic(l,m) + 1
      do 40 i = 1, idimp
      partb(i,ip,m) = parta(i,j,m)
   40 continue
      npic(l,m) = ip
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDSORTP32YZL(parta,partb,npic,npp,noff,nyzp,idimp,npmax
     1,mnblok,nyzpm1,idds)
c this subroutine sorts particles by y,z grid
c linear interpolation, spatial decomposition in y and z direction
c parta/partb = input/output particle array
c parta(2,n,m) = position y of particle n in partition m
c parta(3,n,m) = position z of particle n in partition m
c npic = address offset for reordering particles
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions
c nyzpm1 = max(nyzp(1,m)+1)*max(nyzp(2,m)+1)
c idds = dimensionality of domain decomposition
      implicit none
      integer idimp, npmax, mnblok, nyzpm1, idds
      integer npic, npp, noff, nyzp
      real parta, partb
      dimension parta(idimp,npmax,mnblok), partb(idimp,npmax,mnblok)
      dimension npic(nyzpm1,mnblok)
      dimension npp(mnblok), noff(idds,mnblok), nyzp(idds,mnblok)
c local data
      integer i, j, k, l, m, n, nnoff, lnoff, nyp1, nyzp1, isum, ist, ip
      do 60 m = 1, mnblok
      nnoff = noff(1,m)
      lnoff = noff(2,m)
      nyp1 = nyzp(1,m) + 1
      nyzp1 = nyp1*(nyzp(2,m) + 1)
c clear counter array
      do 10 k = 1, nyzp1
      npic(k,m) = 0
   10 continue
c find how many particles in each grid
      do 20 j = 1, npp(m)
      n = parta(2,j,m) + 1.0
      l = parta(3,j,m)
      l = n - nnoff + nyp1*(l - lnoff)
      npic(l,m) = npic(l,m) + 1
   20 continue
c find address offset
      isum = 0
      do 30 k = 1, nyzp1
      ist = npic(k,m)
      npic(k,m) = isum
      isum = isum + ist
   30 continue
c find addresses of particles at each grid and reorder particles
      do 50 j = 1, npp(m)
      n = parta(2,j,m) + 1.0
      l = parta(3,j,m)
      l = n - nnoff + nyp1*(l - lnoff)
      ip = npic(l,m) + 1
      do 40 i = 1, idimp
      partb(i,ip,m) = parta(i,j,m)
   40 continue
      npic(l,m) = ip
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCOUNT32YZL(part,isign,npicyz,npp,noff,nyzp,idimp,npmax
     1,mnblok,myzpm1,idds)
c this subroutine counts particles by y,z grid, accumulating into npic.
c These arrays are initialized to zero, unless isign = 0
c part = particle array
c isign = (0,1) = (no,yes) initialize npicy, npicz arrays to zero
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c npicyz = number of particles per grid in y and z direction
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c idimp = size of phase space
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions
c myzpm1 = max(max(nyzp(2,m)+1),max(nyzp(1,m)+1))
c idds = dimensionality of domain decomposition
      implicit none
      real part
      integer npicyz, npp, noff, nyzp
      integer isign, idimp, npmax, mnblok, myzpm1, idds
      dimension part(idimp,npmax,mnblok)
      dimension npicyz(myzpm1,idds,mnblok)
      dimension npp(mnblok), noff(idds,mnblok), nyzp(idds,mnblok)
      integer j, l, m, n, nnoff, lnoff, nyp1, nzp1
      do 40 m = 1, mnblok
      nnoff = noff(1,m) - 1
      lnoff = noff(2,m) - 1
      nyp1 = nyzp(1,m) + 1
      nzp1 = nyzp(2,m) + 1
      if (isign.ne.0) then
         do 10 j = 1, nyp1
         npicyz(j,1,m) = 0
   10    continue
         do 20 j = 1, nzp1
         npicyz(j,2,m) = 0
   20    continue
      endif
c find how many particles in each grid
      do 30 j = 1, npp(m)
      n = part(2,j,m)
      l = part(3,j,m)
      n = n - nnoff
      l = l - lnoff
      npicyz(n,1,m) = npicyz(n,1,m) + 1
      npicyz(l,2,m) = npicyz(l,2,m) + 1
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PNCOUNT32YZL(part,isign,npic,npicy,npp,noff,nyzp,idimp,
     1npmax,mnblok,nyzpm1,nypm1,idds)
c this subroutine counts particles by y,z grid, accumulating into npicy
c and npicz.  These arrays are initialized to zero, unless isign = 0
c part = particle array
c isign = (0,1) = (no,yes) initialize npicy, npicz arrays to zero
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c npic = number of particles per grid in y and z direction
c npicy = number of particles per grid in y direction
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c idimp = size of phase space
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions
c nypm1 = max(nyzp(1,m)+1)
c nyzpm1 = nypm1*max(nyzp(2,m)+1)
c idds = dimensionality of domain decomposition
      implicit none
      real part
      integer npic, npicy, npp, noff, nyzp
      integer isign, idimp, npmax, mnblok, nyzpm1, nypm1, idds
      dimension part(idimp,npmax,mnblok)
      dimension npic(nyzpm1,mnblok), npicy(nypm1,mnblok)
      dimension npp(mnblok), noff(idds,mnblok), nyzp(idds,mnblok)
      integer j, l, m, n, nnoff, lnoff, nyzp1, nyp1
      do 40 m = 1, mnblok
      nnoff = noff(1,m) - 1
      lnoff = noff(2,m)
      nyp1 = nyzp(1,m) + 1
      nyzp1 = nyp1*(nyzp(2,m) + 1)
      if (isign.ne.0) then
         do 10 j = 1, nyzp1
         npic(j,m) = 0
   10    continue
         do 20 j = 1, nyp1
         npicy(j,m) = 0
   20    continue
      endif
c find how many particles in each grid
      do 30 j = 1, npp(m)
      n = part(2,j,m)
      l = part(3,j,m)
      n = n - nnoff
      l = n + nyp1*(l - lnoff)
      npicy(n,m) = npicy(n,m) + 1
      npic(l,m) = npic(l,m) + 1
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PRMOVE32(part,npp,ihole,jss,nx,ny,nz,idimp,npmax,mnblok
     1,idds,ntmax,ipbc)
c this subroutine removes particles which would normally be reflected
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position y of particle n in partition m
c npp(m) = number of particles in partition m
c ihole = location of holes left in particle arrays
c jss(idps,m) = scratch array for particle partition m
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions
c idds = dimensionality of domain decomposition
c ntmax =  size of hole array for particles leaving processors
c ipbc = particle boundary condition = (0,-1,-2,-3) =
c (none,2d periodic,2d reflecting,mixed reflecting/periodic)
      implicit none
      real part
      integer npp, ihole, jss
      integer nx, ny, nz, idimp, npmax, mnblok, idds, ntmax, ipbc
      dimension part(idimp,npmax,mnblok)
      dimension npp(mnblok), ihole(ntmax,mnblok), jss(idds,mnblok)
c local data
      integer i, j, j1, j2, m, nter
      real edgelx, edgely, edgelz, edgerx, edgery, edgerz, dx, dy, dz
c set boundary values
      if (ipbc.eq.(-1)) then
         edgelx = 0.
         edgely = 0.
         edgelz = 0.
         edgerx = float(nx)
         edgery = float(ny)
         edgerz = float(nz)
      else if (ipbc.eq.(-2)) then
         edgelx = 1.
         edgely = 1.
         edgelz = 1.
         edgerx = float(nx-1)
         edgery = float(ny-1)
         edgerz = float(nz-1)
      else if (ipbc.eq.(-3)) then
         edgelx = 1.
         edgely = 1.
         edgelz = 0.
         edgerx = float(nx-1)
         edgery = float(ny-1)
         edgerz = float(nz)
      endif
      nter = 0
c buffer outgoing particles
      do 60 m = 1, mnblok
   10 jss(1,m) = 0
      jss(2,m) = 0
      do 20 j = 1, npp(m)
      dx = part(1,j,m)
      dy = part(2,j,m)
      dz = part(3,j,m)
c periodic boundary conditions
      if (ipbc.eq.(-1)) then
         if (dx.lt.edgelx) part(1,j,m) = dx + edgerx
         if (dx.ge.edgerx) part(1,j,m) = part(1,j,m) - edgerx
         if (dy.lt.edgely) part(2,j,m) = dy + edgery
         if (dy.ge.edgery) part(2,j,m) = part(2,j,m) - edgery
         if (dz.lt.edgelz) part(3,j,m) = dz + edgerz
         if (dz.ge.edgerz) part(3,j,m) = part(3,j,m) - edgerz
c reflecting boundary conditions
      else if (ipbc.eq.(-2)) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            if (jss(1,m).lt.ntmax) then
               jss(1,m) = jss(1,m) + 1
               ihole(jss(1,m),m) = j
            else
               jss(2,m) = 1
               go to 30
            endif
         else if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            if (jss(1,m).lt.ntmax) then
               jss(1,m) = jss(1,m) + 1
               ihole(jss(1,m),m) = j
            else
               jss(2,m) = 1
               go to 30
            endif
         else if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            if (jss(1,m).lt.ntmax) then
               jss(1,m) = jss(1,m) + 1
               ihole(jss(1,m),m) = j
            else
               jss(2,m) = 1
               go to 30
            endif
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.(-3)) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            if (jss(1,m).lt.ntmax) then
               jss(1,m) = jss(1,m) + 1
               ihole(jss(1,m),m) = j
            else
               jss(2,m) = 1
               go to 30
            endif
         else if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            if (jss(1,m).lt.ntmax) then
               jss(1,m) = jss(1,m) + 1
               ihole(jss(1,m),m) = j
            else
               jss(2,m) = 1
               go to 30
            endif
         endif
         if (dz.lt.edgelz) part(3,j,m) = dz + edgerz
         if (dz.ge.edgerz) part(3,j,m) = part(3,j,m) - edgerz
      endif
   20 continue
   30 continue
c fill up holes in particle array with particles from bottom
      do 50 j = 1, jss(1,m)
      j1 = npp(m) - j + 1
      j2 = jss(1,m) - j + 1
      if (j1.gt.ihole(j2,m)) then
c move particle only if it is below current hole
         do 40 i = 1, idimp
         part(i,ihole(j2,m),m) = part(i,j1,m)
   40    continue
      endif
   50 continue
      npp(m) = npp(m) - jss(1,m)
c check if buffer overflowed and more particles remain to be checked
      if (jss(2,m).gt.0) then
         nter = nter + 1
         go to 10
      endif
      jss(2,m) = nter
   60 continue
c information
      nter = 0
      do 70 m = 1, mnblok 
      nter = max0(nter,jss(2,m))
   70 continue
      if (nter.gt.0) then
         write (2,*) 'Info: ', nter, ' buffer overflows, ntmax=', ntmax
      endif
      return
      end
