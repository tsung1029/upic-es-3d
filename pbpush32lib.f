c 3d parallel PIC library for pushing particles with magnetic field
c with 2D domain decomposition
c written by viktor k. decyk, ucla
c copyright 1999, regents of the university of california
c update: december 19, 2003
c-----------------------------------------------------------------------
      subroutine PJDOST32(part,cux,cuy,cuz,npp,noff,qm,dt,nx,idimp,npmax
     1,mnblok,nxv,nypmx,nzpmx,idds)
c for 3d code, this subroutine calculates particle current density
c using second-order spline interpolation, periodic boundaries
c and distributed data with 2D spatial decomposition
c in addition, particle positions are advanced a half time-step
c baseline scalar distributed version
c 207 flops/particle, 87 loads, 84 stores
c input: all, output: part, cux, cuy, cuz
c current density is approximated by values at the nearest grid points
c cui(n,m,l)=qci*(.75-dx**2)*(.75-dy**2)*(.75-dz**2)
c cui(n+1,m,l)=.5*qci*((.5+dx)**2)*(.75-dy**2)*(.75-dz**2)
c cuiq(n-1,m,l)=.5*qci*((.5-dx)**2)*(.75-dy**2)*(.75-dz**2)
c cui(n,m+1,l)=.5*qci*(.75-dx**2)*(.5+dy)**2*(.75-dz**2)
c cui(n+1,m+1,l)=.25*qci*((.5+dx)**2)*((.5+dy)**2)*(.75-dz**2)
c cui(n-1,m+1,l)=.25*qci*((.5-dx)**2)*((.5+dy)**2)*(.75-dz**2)
c cui(n,m-1,l)=.5*qci*(.75-dx**2)*((.5-dy)**2)*(.75-dz**2)
c cui(n+1,m-1,l)=.25*qci*((.5+dx)**2)*((.5-dy)**2)*(.75-dz**2)
c cui(n-1,m-1,l)=.25*qci*((.5-dx)**2)*((.5-dy)**2)*(.75-dz**2)
c cui(n,m,l+1)=.5*qci*(.75-dx**2)*(.75-dy**2)*((.5+dz)**2)
c cui(n+1,m,l+1)=.25*qci*((.5+dx)**2)*(.75-dy**2)*((.5+dz)**2)
c cui(n-1,m,l+1)=.25*qci*((.5-dx)**2)*(.75-dy**2)*((.5+dz)**2)
c cui(n,m+1,l+1)=.25*qci*(.75-dx**2)*(.5+dy)**2*((.5+dz)**2)
c cui(n+1,m+1,l+1)=.125*qci*((.5+dx)**2)*((.5+dy)**2)*((.5+dz)**2)
c cui(n-1,m+1,l+1)=.125*qci*((.5-dx)**2)*((.5+dy)**2)*((.5+dz)**2)
c cui(n,m-1,l+1)=.25*qci*(.75-dx**2)*((.5-dy)**2)*((.5+dz)**2)
c cui(n+1,m-1,l+1)=.125*qci*((.5+dx)**2)*((.5-dy)**2)*((.5+dz)**2)
c cui(n-1,m-1,l+1)=.125*qci*((.5-dx)**2)*((.5-dy)**2)*((.5+dz)**2)
c cui(n,m,l-1)=.5*qci*(.75-dx**2)*(.75-dy**2)*((.5-dz)**2)
c cui(n+1,m,l-1)=.25*qci*((.5+dx)**2)*(.75-dy**2)*((.5-dz)**2)
c cui(n-1,m,l-1)=.25*qci*((.5-dx)**2)*(.75-dy**2)*((.5-dz)**2)
c cui(n,m+1,l-1)=.25*qci*(.75-dx**2)*(.5+dy)**2*((.5-dz)**2)
c cui(n+1,m+1,l-1)=.125*qci*((.5+dx)**2)*((.5+dy)**2)*((.5-dz)**2)
c cui(n-1,m+1,l-1)=.125*qci*((.5-dx)**2)*((.5+dy)**2)*((.5-dz)**2)
c cui(n,m-1,l-1)=.25*qci*(.75-dx**2)*((.5-dy)**2)*((.5-dz)**2)
c cui(n+1,m-1,l-1)=.125*qci*((.5+dx)**2)*((.5-dy)**2)*((.5-dz)**2)
c cui(n-1,m-1,l-1)=.125*qci*((.5-dx)**2)*((.5-dy)**2)*((.5-dz)**2)
c where n,m,l = nearest grid points and dx = x-n, dy = y-m, dz = z-l
c and qci = qm*vi, where i = x,y,z
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = x velocity of particle n in partition m
c part(5,n,m) = y velocity of particle n in partition m
c part(6,n,m) = z velocity of particle n in partition m
c cui(j,k,l,m)=ith component of current density at grid point (j,kk,ll)
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx = system length in x direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first dimension of current array, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
      dimension part(idimp,npmax,mnblok), cux(nxv,nypmx,nzpmx,mnblok)
      dimension cuy(nxv,nypmx,nzpmx,mnblok), cuz(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      zero = 0.
      anx = float(nx)
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
c deposit current
      dx = dxl*dyl*dzl
      dy = amx*dyl*dzl
      dz = dxp*dyl*dzl
      vx = part(4,j,m)
      vy = part(5,j,m)
      vz = part(6,j,m)
      cux(nl,ml,lm,m) = cux(nl,ml,lm,m) + vx*dx
      cux(nn,ml,lm,m) = cux(nn,ml,lm,m) + vx*dy
      cux(np,ml,lm,m) = cux(np,ml,lm,m) + vx*dz
      cuy(nl,ml,lm,m) = cuy(nl,ml,lm,m) + vy*dx
      cuy(nn,ml,lm,m) = cuy(nn,ml,lm,m) + vy*dy
      cuy(np,ml,lm,m) = cuy(np,ml,lm,m) + vy*dz
      cuz(nl,ml,lm,m) = cuz(nl,ml,lm,m) + vz*dx
      cuz(nn,ml,lm,m) = cuz(nn,ml,lm,m) + vz*dy
      cuz(np,ml,lm,m) = cuz(np,ml,lm,m) + vz*dz
      dx = dxl*amy*dzl
      dy = amx*amy*dzl
      dz = dxp*amy*dzl
      cux(nl,mm,lm,m) = cux(nl,mm,lm,m) + vx*dx
      cux(nn,mm,lm,m) = cux(nn,mm,lm,m) + vx*dy
      cux(np,mm,lm,m) = cux(np,mm,lm,m) + vx*dz
      cuy(nl,mm,lm,m) = cuy(nl,mm,lm,m) + vy*dx
      cuy(nn,mm,lm,m) = cuy(nn,mm,lm,m) + vy*dy
      cuy(np,mm,lm,m) = cuy(np,mm,lm,m) + vy*dz
      cuz(nl,mm,lm,m) = cuz(nl,mm,lm,m) + vz*dx
      cuz(nn,mm,lm,m) = cuz(nn,mm,lm,m) + vz*dy
      cuz(np,mm,lm,m) = cuz(np,mm,lm,m) + vz*dz
      dx = dxl*dyp*dzl
      dy = amx*dyp*dzl
      dz = dxp*dyp*dzl
      cux(nl,mp,lm,m) = cux(nl,mp,lm,m) + vx*dx
      cux(nn,mp,lm,m) = cux(nn,mp,lm,m) + vx*dy
      cux(np,mp,lm,m) = cux(np,mp,lm,m) + vx*dz
      cuy(nl,mp,lm,m) = cuy(nl,mp,lm,m) + vy*dx
      cuy(nn,mp,lm,m) = cuy(nn,mp,lm,m) + vy*dy
      cuy(np,mp,lm,m) = cuy(np,mp,lm,m) + vy*dz
      cuz(nl,mp,lm,m) = cuz(nl,mp,lm,m) + vz*dx
      cuz(nn,mp,lm,m) = cuz(nn,mp,lm,m) + vz*dy
      cuz(np,mp,lm,m) = cuz(np,mp,lm,m) + vz*dz
      dx = dxl*dyl*amz
      dy = amx*dyl*amz
      dz = dxp*dyl*amz
      cux(nl,ml,ll,m) = cux(nl,ml,ll,m) + vx*dx
      cux(nn,ml,ll,m) = cux(nn,ml,ll,m) + vx*dy
      cux(np,ml,ll,m) = cux(np,ml,ll,m) + vx*dz
      cuy(nl,ml,ll,m) = cuy(nl,ml,ll,m) + vy*dx
      cuy(nn,ml,ll,m) = cuy(nn,ml,ll,m) + vy*dy
      cuy(np,ml,ll,m) = cuy(np,ml,ll,m) + vy*dz
      cuz(nl,ml,ll,m) = cuz(nl,ml,ll,m) + vz*dx
      cuz(nn,ml,ll,m) = cuz(nn,ml,ll,m) + vz*dy
      cuz(np,ml,ll,m) = cuz(np,ml,ll,m) + vz*dz
      dx = dxl*amy*amz
      dy = amx*amy*amz
      dz = dxp*amy*amz
      cux(nl,mm,ll,m) = cux(nl,mm,ll,m) + vx*dx
      cux(nn,mm,ll,m) = cux(nn,mm,ll,m) + vx*dy
      cux(np,mm,ll,m) = cux(np,mm,ll,m) + vx*dz
      cuy(nl,mm,ll,m) = cuy(nl,mm,ll,m) + vy*dx
      cuy(nn,mm,ll,m) = cuy(nn,mm,ll,m) + vy*dy
      cuy(np,mm,ll,m) = cuy(np,mm,ll,m) + vy*dz
      cuz(nl,mm,ll,m) = cuz(nl,mm,ll,m) + vz*dx
      cuz(nn,mm,ll,m) = cuz(nn,mm,ll,m) + vz*dy
      cuz(np,mm,ll,m) = cuz(np,mm,ll,m) + vz*dz
      dx = dxl*dyp*amz
      dy = amx*dyp*amz
      dz = dxp*dyp*amz
      cux(nl,mp,ll,m) = cux(nl,mp,ll,m) + vx*dx
      cux(nn,mp,ll,m) = cux(nn,mp,ll,m) + vx*dy
      cux(np,mp,ll,m) = cux(np,mp,ll,m) + vx*dz
      cuy(nl,mp,ll,m) = cuy(nl,mp,ll,m) + vy*dx
      cuy(nn,mp,ll,m) = cuy(nn,mp,ll,m) + vy*dy
      cuy(np,mp,ll,m) = cuy(np,mp,ll,m) + vy*dz
      cuz(nl,mp,ll,m) = cuz(nl,mp,ll,m) + vz*dx
      cuz(nn,mp,ll,m) = cuz(nn,mp,ll,m) + vz*dy
      cuz(np,mp,ll,m) = cuz(np,mp,ll,m) + vz*dz
      dx = dxl*dyl*dzp
      dy = amx*dyl*dzp
      dz = dxp*dyl*dzp
      cux(nl,ml,lp,m) = cux(nl,ml,lp,m) + vx*dx
      cux(nn,ml,lp,m) = cux(nn,ml,lp,m) + vx*dy
      cux(np,ml,lp,m) = cux(np,ml,lp,m) + vx*dz
      cuy(nl,ml,lp,m) = cuy(nl,ml,lp,m) + vy*dx
      cuy(nn,ml,lp,m) = cuy(nn,ml,lp,m) + vy*dy
      cuy(np,ml,lp,m) = cuy(np,ml,lp,m) + vy*dz
      cuz(nl,ml,lp,m) = cuz(nl,ml,lp,m) + vz*dx
      cuz(nn,ml,lp,m) = cuz(nn,ml,lp,m) + vz*dy
      cuz(np,ml,lp,m) = cuz(np,ml,lp,m) + vz*dz
      dx = dxl*amy*dzp
      dy = amx*amy*dzp
      dz = dxp*amy*dzp
      cux(nl,mm,lp,m) = cux(nl,mm,lp,m) + vx*dx
      cux(nn,mm,lp,m) = cux(nn,mm,lp,m) + vx*dy
      cux(np,mm,lp,m) = cux(np,mm,lp,m) + vx*dz
      cuy(nl,mm,lp,m) = cuy(nl,mm,lp,m) + vy*dx
      cuy(nn,mm,lp,m) = cuy(nn,mm,lp,m) + vy*dy
      cuy(np,mm,lp,m) = cuy(np,mm,lp,m) + vy*dz
      cuz(nl,mm,lp,m) = cuz(nl,mm,lp,m) + vz*dx
      cuz(nn,mm,lp,m) = cuz(nn,mm,lp,m) + vz*dy
      cuz(np,mm,lp,m) = cuz(np,mm,lp,m) + vz*dz
      dx = dxl*dyp*dzp
      dy = amx*dyp*dzp
      dz = dxp*dyp*dzp
      cux(nl,mp,lp,m) = cux(nl,mp,lp,m) + vx*dx
      cux(nn,mp,lp,m) = cux(nn,mp,lp,m) + vx*dy
      cux(np,mp,lp,m) = cux(np,mp,lp,m) + vx*dz
      cuy(nl,mp,lp,m) = cuy(nl,mp,lp,m) + vy*dx
      cuy(nn,mp,lp,m) = cuy(nn,mp,lp,m) + vy*dy
      cuy(np,mp,lp,m) = cuy(np,mp,lp,m) + vy*dz
      cuz(nl,mp,lp,m) = cuz(nl,mp,lp,m) + vz*dx
      cuz(nn,mp,lp,m) = cuz(nn,mp,lp,m) + vz*dy
      cuz(np,mp,lp,m) = cuz(np,mp,lp,m) + vz*dz
c advance position half a time-step
      dx = part(1,j,m) + vx*dt
      dy = part(2,j,m) + vy*dt
      dz = part(3,j,m) + vz*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j,m) = dx
      part(2,j,m) = dy
      part(3,j,m) = dz
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGJPOST32(part,cu,npp,noff,qm,dt,nx,ny,nz,idimp,npmax,m
     1nblok,nxv,nypmx,nzpmx,idds,ipbc)
c for 3d code, this subroutine calculates particle current density
c using second-order spline interpolation, for distributed data
c with 2D spatial decomposition
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 189 flops/particle, 87 loads, 84 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m,l)=qci*(.75-dx**2)*(.75-dy**2)*(.75-dz**2)
c cu(i,n+1,m,l)=.5*qci*((.5+dx)**2)*(.75-dy**2)*(.75-dz**2)
c cu(i,n-1,m,l)=.5*qci*((.5-dx)**2)*(.75-dy**2)*(.75-dz**2)
c cu(i,n,m+1,l)=.5*qci*(.75-dx**2)*(.5+dy)**2*(.75-dz**2)
c cu(i,n+1,m+1,l)=.25*qci*((.5+dx)**2)*((.5+dy)**2)*(.75-dz**2)
c cu(i,n-1,m+1,l)=.25*qci*((.5-dx)**2)*((.5+dy)**2)*(.75-dz**2)
c cu(i,n,m-1,l)=.5*qci*(.75-dx**2)*((.5-dy)**2)*(.75-dz**2)
c cu(i,n+1,m-1,l)=.25*qci*((.5+dx)**2)*((.5-dy)**2)*(.75-dz**2)
c cu(i,n-1,m-1,l)=.25*qci*((.5-dx)**2)*((.5-dy)**2)*(.75-dz**2)
c cu(i,n,m,l+1)=.5*qci*(.75-dx**2)*(.75-dy**2)*((.5+dz)**2)
c cu(i,n+1,m,l+1)=.25*qci*((.5+dx)**2)*(.75-dy**2)*((.5+dz)**2)
c cu(i,n-1,m,l+1)=.25*qci*((.5-dx)**2)*(.75-dy**2)*((.5+dz)**2)
c cu(i,n,m+1,l+1)=.25*qci*(.75-dx**2)*(.5+dy)**2*((.5+dz)**2)
c cu(i,n+1,m+1,l+1)=.125*qci*((.5+dx)**2)*((.5+dy)**2)*((.5+dz)**2)
c cu(i,n-1,m+1,l+1)=.125*qci*((.5-dx)**2)*((.5+dy)**2)*((.5+dz)**2)
c cu(i,n,m-1,l+1)=.25*qci*(.75-dx**2)*((.5-dy)**2)*((.5+dz)**2)
c cu(i,n+1,m-1,l+1)=.125*qci*((.5+dx)**2)*((.5-dy)**2)*((.5+dz)**2)
c cu(i,n-1,m-1,l+1)=.125*qci*((.5-dx)**2)*((.5-dy)**2)*((.5+dz)**2)
c cu(i,n,m,l-1)=.5*qci*(.75-dx**2)*(.75-dy**2)*((.5-dz)**2)
c cu(i,n+1,m,l-1)=.25*qci*((.5+dx)**2)*(.75-dy**2)*((.5-dz)**2)
c cu(i,n-1,m,l-1)=.25*qci*((.5-dx)**2)*(.75-dy**2)*((.5-dz)**2)
c cu(i,n,m+1,l-1)=.25*qci*(.75-dx**2)*(.5+dy)**2*((.5-dz)**2)
c cu(i,n+1,m+1,l-1)=.125*qci*((.5+dx)**2)*((.5+dy)**2)*((.5-dz)**2)
c cu(i,n-1,m+1,l-1)=.125*qci*((.5-dx)**2)*((.5+dy)**2)*((.5-dz)**2)
c cu(i,n,m-1,l-1)=.25*qci*(.75-dx**2)*((.5-dy)**2)*((.5-dz)**2)
c cu(i,n+1,m-1,l-1)=.125*qci*((.5+dx)**2)*((.5-dy)**2)*((.5-dz)**2)
c cu(i,n-1,m-1,l-1)=.125*qci*((.5-dx)**2)*((.5-dy)**2)*((.5-dz)**2)
c where n,m,l = nearest grid points and dx = x-n, dy = y-m, dz = z-l
c and qci = qm*vi, where i = x,y,z
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = x velocity of particle n in partition m
c part(5,n,m) = y velocity of particle n in partition m
c part(6,n,m) = z velocity of particle n in partition m
c cu(i,j+1,k,l,m) = ith component of current density at grid point 
c (j,kk,ll), where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first dimension of current array, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      dimension part(idimp,npmax,mnblok), cu(3,nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      qmh = .5*qm
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
c deposit current
      dx = dx1*dzl
      dy = dx2*dzl
      dz = dyl*dzl
      vx = part(4,j,m)
      vy = part(5,j,m)
      vz = part(6,j,m)
      cu(1,nl,ml,lm,m) = cu(1,nl,ml,lm,m) + vx*dx
      cu(2,nl,ml,lm,m) = cu(2,nl,ml,lm,m) + vy*dx
      cu(3,nl,ml,lm,m) = cu(3,nl,ml,lm,m) + vz*dx
      dx = dy1*dzl
      cu(1,nn,ml,lm,m) = cu(1,nn,ml,lm,m) + vx*dy
      cu(2,nn,ml,lm,m) = cu(2,nn,ml,lm,m) + vy*dy
      cu(3,nn,ml,lm,m) = cu(3,nn,ml,lm,m) + vz*dy
      dy = dy2*dzl
      cu(1,np,ml,lm,m) = cu(1,np,ml,lm,m) + vx*dz
      cu(2,np,ml,lm,m) = cu(2,np,ml,lm,m) + vy*dz
      cu(3,np,ml,lm,m) = cu(3,np,ml,lm,m) + vz*dz
      dz = amy*dzl
      cu(1,nl,mm,lm,m) = cu(1,nl,mm,lm,m) + vx*dx
      cu(2,nl,mm,lm,m) = cu(2,nl,mm,lm,m) + vy*dx
      cu(3,nl,mm,lm,m) = cu(3,nl,mm,lm,m) + vz*dx
      dx = dxl*dzl
      cu(1,nn,mm,lm,m) = cu(1,nn,mm,lm,m) + vx*dy
      cu(2,nn,mm,lm,m) = cu(2,nn,mm,lm,m) + vy*dy
      cu(3,nn,mm,lm,m) = cu(3,nn,mm,lm,m) + vz*dy
      dy = amx*dzl
      cu(1,np,mm,lm,m) = cu(1,np,mm,lm,m) + vx*dz
      cu(2,np,mm,lm,m) = cu(2,np,mm,lm,m) + vy*dz
      cu(3,np,mm,lm,m) = cu(3,np,mm,lm,m) + vz*dz
      dz = dyp*dzl
      cu(1,nl,mp,lm,m) = cu(1,nl,mp,lm,m) + vx*dx
      cu(2,nl,mp,lm,m) = cu(2,nl,mp,lm,m) + vy*dx
      cu(3,nl,mp,lm,m) = cu(3,nl,mp,lm,m) + vz*dx
      dx = dx1*amz
      cu(1,nn,mp,lm,m) = cu(1,nn,mp,lm,m) + vx*dy
      cu(2,nn,mp,lm,m) = cu(2,nn,mp,lm,m) + vy*dy
      cu(3,nn,mp,lm,m) = cu(3,nn,mp,lm,m) + vz*dy
      dy = dx2*amz
      cu(1,np,mp,lm,m) = cu(1,np,mp,lm,m) + vx*dz
      cu(2,np,mp,lm,m) = cu(2,np,mp,lm,m) + vy*dz
      cu(3,np,mp,lm,m) = cu(3,np,mp,lm,m) + vz*dz
      dz = dyl*amz
      cu(1,nl,ml,ll,m) = cu(1,nl,ml,ll,m) + vx*dx
      cu(2,nl,ml,ll,m) = cu(2,nl,ml,ll,m) + vy*dx
      cu(3,nl,ml,ll,m) = cu(3,nl,ml,ll,m) + vz*dx
      dx = dy1*amz
      cu(1,nn,ml,ll,m) = cu(1,nn,ml,ll,m) + vx*dy
      cu(2,nn,ml,ll,m) = cu(2,nn,ml,ll,m) + vy*dy
      cu(3,nn,ml,ll,m) = cu(3,nn,ml,ll,m) + vz*dy
      dy = dy2*amz
      cu(1,np,ml,ll,m) = cu(1,np,ml,ll,m) + vx*dz
      cu(2,np,ml,ll,m) = cu(2,np,ml,ll,m) + vy*dz
      cu(3,np,ml,ll,m) = cu(3,np,ml,ll,m) + vz*dz
      dz = amy*amz
      cu(1,nl,mm,ll,m) = cu(1,nl,mm,ll,m) + vx*dx
      cu(2,nl,mm,ll,m) = cu(2,nl,mm,ll,m) + vy*dx
      cu(3,nl,mm,ll,m) = cu(3,nl,mm,ll,m) + vz*dx
      dx = dxl*amz
      cu(1,nn,mm,ll,m) = cu(1,nn,mm,ll,m) + vx*dy
      cu(2,nn,mm,ll,m) = cu(2,nn,mm,ll,m) + vy*dy
      cu(3,nn,mm,ll,m) = cu(3,nn,mm,ll,m) + vz*dy
      dy = amx*amz
      cu(1,np,mm,ll,m) = cu(1,np,mm,ll,m) + vx*dz
      cu(2,np,mm,ll,m) = cu(2,np,mm,ll,m) + vy*dz
      cu(3,np,mm,ll,m) = cu(3,np,mm,ll,m) + vz*dz
      dz = dyp*amz
      cu(1,nl,mp,ll,m) = cu(1,nl,mp,ll,m) + vx*dx
      cu(2,nl,mp,ll,m) = cu(2,nl,mp,ll,m) + vy*dx
      cu(3,nl,mp,ll,m) = cu(3,nl,mp,ll,m) + vz*dx
      dx = dx1*dzp
      cu(1,nn,mp,ll,m) = cu(1,nn,mp,ll,m) + vx*dy
      cu(2,nn,mp,ll,m) = cu(2,nn,mp,ll,m) + vy*dy
      cu(3,nn,mp,ll,m) = cu(3,nn,mp,ll,m) + vz*dy
      dy = dx2*dzp
      cu(1,np,mp,ll,m) = cu(1,np,mp,ll,m) + vx*dz
      cu(2,np,mp,ll,m) = cu(2,np,mp,ll,m) + vy*dz
      cu(3,np,mp,ll,m) = cu(3,np,mp,ll,m) + vz*dz
      dz = dyl*dzp
      cu(1,nl,ml,lp,m) = cu(1,nl,ml,lp,m) + vx*dx
      cu(2,nl,ml,lp,m) = cu(2,nl,ml,lp,m) + vy*dx
      cu(3,nl,ml,lp,m) = cu(3,nl,ml,lp,m) + vz*dx
      dx = dy1*dzp
      cu(1,nn,ml,lp,m) = cu(1,nn,ml,lp,m) + vx*dy
      cu(2,nn,ml,lp,m) = cu(2,nn,ml,lp,m) + vy*dy
      cu(3,nn,ml,lp,m) = cu(3,nn,ml,lp,m) + vz*dy
      dy = dy2*dzp
      cu(1,np,ml,lp,m) = cu(1,np,ml,lp,m) + vx*dz
      cu(2,np,ml,lp,m) = cu(2,np,ml,lp,m) + vy*dz
      cu(3,np,ml,lp,m) = cu(3,np,ml,lp,m) + vz*dz
      dz = amy*dzp
      cu(1,nl,mm,lp,m) = cu(1,nl,mm,lp,m) + vx*dx
      cu(2,nl,mm,lp,m) = cu(2,nl,mm,lp,m) + vy*dx
      cu(3,nl,mm,lp,m) = cu(3,nl,mm,lp,m) + vz*dx
      dx = dxl*dzp
      cu(1,nn,mm,lp,m) = cu(1,nn,mm,lp,m) + vx*dy
      cu(2,nn,mm,lp,m) = cu(2,nn,mm,lp,m) + vy*dy
      cu(3,nn,mm,lp,m) = cu(3,nn,mm,lp,m) + vz*dy
      dy = amx*dzp
      cu(1,np,mm,lp,m) = cu(1,np,mm,lp,m) + vx*dz
      cu(2,np,mm,lp,m) = cu(2,np,mm,lp,m) + vy*dz
      cu(3,np,mm,lp,m) = cu(3,np,mm,lp,m) + vz*dz
      dz = dyp*dzp
      cu(1,nl,mp,lp,m) = cu(1,nl,mp,lp,m) + vx*dx
      cu(2,nl,mp,lp,m) = cu(2,nl,mp,lp,m) + vy*dx
      cu(3,nl,mp,lp,m) = cu(3,nl,mp,lp,m) + vz*dx
      cu(1,nn,mp,lp,m) = cu(1,nn,mp,lp,m) + vx*dy
      cu(2,nn,mp,lp,m) = cu(2,nn,mp,lp,m) + vy*dy
      cu(3,nn,mp,lp,m) = cu(3,nn,mp,lp,m) + vz*dy
      cu(1,np,mp,lp,m) = cu(1,np,mp,lp,m) + vx*dz
      cu(2,np,mp,lp,m) = cu(2,np,mp,lp,m) + vy*dz
      cu(3,np,mp,lp,m) = cu(3,np,mp,lp,m) + vz*dz
c advance position half a time-step
      dx = part(1,j,m) + vx*dt
      dy = part(2,j,m) + vy*dt
      dz = part(3,j,m) + vz*dt
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
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSJPOST32(part,cu,npp,noff,qm,dt,nx,ny,nz,idimp,npmax,
     1mnblok,nxv,nypmx,nxyzp,idds,ipbc)
c for 3d code, this subroutine calculates particle current density
c using second-order spline interpolation, for distributed data
c with 2D spatial decomposition
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 189 flops/particle, 87 loads, 84 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m,l)=qci*(.75-dx**2)*(.75-dy**2)*(.75-dz**2)
c cu(i,n+1,m,l)=.5*qci*((.5+dx)**2)*(.75-dy**2)*(.75-dz**2)
c cu(i,n-1,m,l)=.5*qci*((.5-dx)**2)*(.75-dy**2)*(.75-dz**2)
c cu(i,n,m+1,l)=.5*qci*(.75-dx**2)*(.5+dy)**2*(.75-dz**2)
c cu(i,n+1,m+1,l)=.25*qci*((.5+dx)**2)*((.5+dy)**2)*(.75-dz**2)
c cu(i,n-1,m+1,l)=.25*qci*((.5-dx)**2)*((.5+dy)**2)*(.75-dz**2)
c cu(i,n,m-1,l)=.5*qci*(.75-dx**2)*((.5-dy)**2)*(.75-dz**2)
c cu(i,n+1,m-1,l)=.25*qci*((.5+dx)**2)*((.5-dy)**2)*(.75-dz**2)
c cu(i,n-1,m-1,l)=.25*qci*((.5-dx)**2)*((.5-dy)**2)*(.75-dz**2)
c cu(i,n,m,l+1)=.5*qci*(.75-dx**2)*(.75-dy**2)*((.5+dz)**2)
c cu(i,n+1,m,l+1)=.25*qci*((.5+dx)**2)*(.75-dy**2)*((.5+dz)**2)
c cu(i,n-1,m,l+1)=.25*qci*((.5-dx)**2)*(.75-dy**2)*((.5+dz)**2)
c cu(i,n,m+1,l+1)=.25*qci*(.75-dx**2)*(.5+dy)**2*((.5+dz)**2)
c cu(i,n+1,m+1,l+1)=.125*qci*((.5+dx)**2)*((.5+dy)**2)*((.5+dz)**2)
c cu(i,n-1,m+1,l+1)=.125*qci*((.5-dx)**2)*((.5+dy)**2)*((.5+dz)**2)
c cu(i,n,m-1,l+1)=.25*qci*(.75-dx**2)*((.5-dy)**2)*((.5+dz)**2)
c cu(i,n+1,m-1,l+1)=.125*qci*((.5+dx)**2)*((.5-dy)**2)*((.5+dz)**2)
c cu(i,n-1,m-1,l+1)=.125*qci*((.5-dx)**2)*((.5-dy)**2)*((.5+dz)**2)
c cu(i,n,m,l-1)=.5*qci*(.75-dx**2)*(.75-dy**2)*((.5-dz)**2)
c cu(i,n+1,m,l-1)=.25*qci*((.5+dx)**2)*(.75-dy**2)*((.5-dz)**2)
c cu(i,n-1,m,l-1)=.25*qci*((.5-dx)**2)*(.75-dy**2)*((.5-dz)**2)
c cu(i,n,m+1,l-1)=.25*qci*(.75-dx**2)*(.5+dy)**2*((.5-dz)**2)
c cu(i,n+1,m+1,l-1)=.125*qci*((.5+dx)**2)*((.5+dy)**2)*((.5-dz)**2)
c cu(i,n-1,m+1,l-1)=.125*qci*((.5-dx)**2)*((.5+dy)**2)*((.5-dz)**2)
c cu(i,n,m-1,l-1)=.25*qci*(.75-dx**2)*((.5-dy)**2)*((.5-dz)**2)
c cu(i,n+1,m-1,l-1)=.125*qci*((.5+dx)**2)*((.5-dy)**2)*((.5-dz)**2)
c cu(i,n-1,m-1,l-1)=.125*qci*((.5-dx)**2)*((.5-dy)**2)*((.5-dz)**2)
c where n,m,l = nearest grid points and dx = x-n, dy = y-m, dz = z-l
c and qci = qm*vi, where i = x,y,z
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = x velocity of particle n in partition m
c part(5,n,m) = y velocity of particle n in partition m
c part(6,n,m) = z velocity of particle n in partition m
c cu(i,n,m) = ith component of current density at grid point (j,kk,ll)
c where n = j + nxv*kk + nxv*nypmx*ll
c and kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first virtual dimension of current array, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nxyzp = dimension of charge array, must be >= nxv*nypmx*nzpmx
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      dimension part(idimp,npmax,mnblok), cu(3,nxyzp,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      qmh = .5*qm
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
c deposit current
      dx = dx1*dzl
      dy = dx2*dzl
      dz = dyl*dzl
      vx = part(4,j-1,m)
      vy = part(5,j-1,m)
      vz = part(6,j-1,m)
      dx3 = cu(1,ml,m) + vx*dx
      dy3 = cu(2,ml,m) + vy*dx
      dz3 = cu(3,ml,m) + vz*dx
      dxp = cu(1,ml+1,m) + vx*dy
      dy4 = cu(2,ml+1,m) + vy*dy
      dz4 = cu(3,ml+1,m) + vz*dy
      dx = cu(1,ml+2,m) + vx*dz
      dy = cu(2,ml+2,m) + vy*dz
      dz = cu(3,ml+2,m) + vz*dz
      cu(1,ml,m) = dx3
      cu(2,ml,m) = dy3
      cu(3,ml,m) = dz3
      cu(1,ml+1,m) = dxp
      cu(2,ml+1,m) = dy4
      cu(3,ml+1,m) = dz4
      cu(1,ml+2,m) = dx
      cu(2,ml+2,m) = dy
      cu(3,ml+2,m) = dz
      dx = dy1*dzl
      dy = dy2*dzl
      dz = amy*dzl
      dx3 = cu(1,mn,m) + vx*dx
      dy3 = cu(2,mn,m) + vy*dx
      dz3 = cu(3,mn,m) + vz*dx
      dxp = cu(1,mn+1,m) + vx*dy
      dy4 = cu(2,mn+1,m) + vy*dy
      dz4 = cu(3,mn+1,m) + vz*dy
      dx = cu(1,mn+2,m) + vx*dz
      dy = cu(2,mn+2,m) + vy*dz
      dz = cu(3,mn+2,m) + vz*dz
      cu(1,mn,m) = dx3
      cu(2,mn,m) = dy3
      cu(3,mn,m) = dz3
      cu(1,mn+1,m) = dxp
      cu(2,mn+1,m) = dy4
      cu(3,mn+1,m) = dz4
      cu(1,mn+2,m) = dx
      cu(2,mn+2,m) = dy
      cu(3,mn+2,m) = dz
      dx = dxl*dzl
      dy = amx*dzl
      dz = dyp*dzl
      dx3 = cu(1,mp,m) + vx*dx
      dy3 = cu(2,mp,m) + vy*dx
      dz3 = cu(3,mp,m) + vz*dx
      dxp = cu(1,mp+1,m) + vx*dy
      dy4 = cu(2,mp+1,m) + vy*dy
      dzl = cu(3,mp+1,m) + vz*dy
      dx = cu(1,mp+2,m) + vx*dz
      dy = cu(2,mp+2,m) + vy*dz
      dz = cu(3,mp+2,m) + vz*dz
      cu(1,mp,m) = dx3
      cu(2,mp,m) = dy3
      cu(3,mp,m) = dz3
      cu(1,mp+1,m) = dxp
      cu(2,mp+1,m) = dy4
      cu(3,mp+1,m) = dzl
      cu(1,mp+2,m) = dx 
      cu(2,mp+2,m) = dy
      cu(3,mp+2,m) = dz
      ml = ml + nxvy
      mn = mn + nxvy
      mp = mp + nxvy
      dx = dx1*amz
      dy = dx2*amz
      dz = dyl*amz
      dx3 = cu(1,ml,m) + vx*dx
      dy3 = cu(2,ml,m) + vy*dx
      dz3 = cu(3,ml,m) + vz*dx
      dxp = cu(1,ml+1,m) + vx*dy
      dy4 = cu(2,ml+1,m) + vy*dy
      dzl = cu(3,ml+1,m) + vz*dy
      dx = cu(1,ml+2,m) + vx*dz
      dy = cu(2,ml+2,m) + vy*dz
      dz = cu(3,ml+2,m) + vz*dz
      cu(1,ml,m) = dx3
      cu(2,ml,m) = dy3
      cu(3,ml,m) = dz3
      cu(1,ml+1,m) = dxp
      cu(2,ml+1,m) = dy4
      cu(3,ml+1,m) = dzl
      cu(1,ml+2,m) = dx
      cu(2,ml+2,m) = dy
      cu(3,ml+2,m) = dz
      dx = dy1*amz
      dy = dy2*amz
      dz = amy*amz
      dx3 = cu(1,mn,m) + vx*dx
      dy3 = cu(2,mn,m) + vy*dx
      dz3 = cu(3,mn,m) + vz*dx
      dxp = cu(1,mn+1,m) + vx*dy
      dy4 = cu(2,mn+1,m) + vy*dy
      dzl = cu(3,mn+1,m) + vz*dy
      dx = cu(1,mn+2,m) + vx*dz
      dy = cu(2,mn+2,m) + vy*dz
      dz = cu(3,mn+2,m) + vz*dz
      cu(1,mn,m) = dx3
      cu(2,mn,m) = dy3
      cu(3,mn,m) = dz3
      cu(1,mn+1,m) = dxp
      cu(2,mn+1,m) = dy4
      cu(3,mn+1,m) = dzl
      cu(1,mn+2,m) = dx
      cu(2,mn+2,m) = dy
      cu(3,mn+2,m) = dz
      dx = dxl*amz
      dy = amx*amz
      dz = dyp*amz
      dx3 = cu(1,mp,m) + vx*dx
      dy3 = cu(2,mp,m) + vy*dx
      dz3 = cu(3,mp,m) + vz*dx
      dxp = cu(1,mp+1,m) + vx*dy
      dzl = cu(2,mp+1,m) + vy*dy
      amz = cu(3,mp+1,m) + vz*dy
      dx = cu(1,mp+2,m) + vx*dz
      dy = cu(2,mp+2,m) + vy*dz
      dz = cu(3,mp+2,m) + vz*dz
      cu(1,mp,m) = dx3
      cu(2,mp,m) = dy3
      cu(3,mp,m) = dz3
      cu(1,mp+1,m) = dxp
      cu(2,mp+1,m) = dzl
      cu(3,mp+1,m) = amz
      cu(1,mp+2,m) = dx
      cu(2,mp+2,m) = dy
      cu(3,mp+2,m) = dz
      ml = ml + nxvy
      mn = mn + nxvy
      mp = mp + nxvy
      dx = dx1*dzp
      dy = dx2*dzp
      dz = dyl*dzp
      dx3 = cu(1,ml,m) + vx*dx
      dy3 = cu(2,ml,m) + vy*dx
      dz3 = cu(3,ml,m) + vz*dx
      dxp = cu(1,ml+1,m) + vx*dy
      dzl = cu(2,ml+1,m) + vy*dy
      amz = cu(3,ml+1,m) + vz*dy
      dx = cu(1,ml+2,m) + vx*dz
      dy = cu(2,ml+2,m) + vy*dz
      dz = cu(3,ml+2,m) + vz*dz
      cu(1,ml,m) = dx3
      cu(2,ml,m) = dy3
      cu(3,ml,m) = dz3
      cu(1,ml+1,m) = dxp
      cu(2,ml+1,m) = dzl
      cu(3,ml+1,m) = amz
      cu(1,ml+2,m) = dx
      cu(2,ml+2,m) = dy
      cu(3,ml+2,m) = dz
      dx = dy1*dzp
      dy = dy2*dzp
      dz = amy*dzp
      dx3 = cu(1,mn,m) + vx*dx
      dy3 = cu(2,mn,m) + vy*dx
      dz3 = cu(3,mn,m) + vz*dx
      dxp = cu(1,mn+1,m) + vx*dy
      dzl = cu(2,mn+1,m) + vy*dy
      amz = cu(3,mn+1,m) + vz*dy
      dx = cu(1,mn+2,m) + vx*dz
      dy = cu(2,mn+2,m) + vy*dz
      dz = cu(3,mn+2,m) + vz*dz
      cu(1,mn,m) = dx3
      cu(2,mn,m) = dy3
      cu(3,mn,m) = dz3
      cu(1,mn+1,m) = dxp
      cu(2,mn+1,m) = dzl
      cu(3,mn+1,m) = amz
      cu(1,mn+2,m) = dx
      cu(2,mn+2,m) = dy
      cu(3,mn+2,m) = dz
      dx = dxl*dzp
      dy = amx*dzp
      dz = dyp*dzp
      dx3 = cu(1,mp,m) + vx*dx
      dy3 = cu(2,mp,m) + vy*dx
      dz3 = cu(3,mp,m) + vz*dx
      dxp = cu(1,mp+1,m) + vx*dy
      dzl = cu(2,mp+1,m) + vy*dy
      amz = cu(3,mp+1,m) + vz*dy
      dx = cu(1,mp+2,m) + vx*dz
      dy = cu(2,mp+2,m) + vy*dz
      dz = cu(3,mp+2,m) + vz*dz
      cu(1,mp,m) = dx3
      cu(2,mp,m) = dy3
      cu(3,mp,m) = dz3
      cu(1,mp+1,m) = dxp
      cu(2,mp+1,m) = dzl
      cu(3,mp+1,m) = amz
      cu(1,mp+2,m) = dx
      cu(2,mp+2,m) = dy
      cu(3,mp+2,m) = dz
c advance position half a time-step
      dx = part(1,j-1,m) + vx*dt
      dy = part(2,j-1,m) + vy*dt
      dz = part(3,j-1,m) + vz*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1,m)
            part(4,j-1,m) = -part(4,j-1,m)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j-1,m)
            part(5,j-1,m) = -part(5,j-1,m)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = part(3,j-1,m)
            part(6,j-1,m) = -part(6,j-1,m)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1,m)
            part(4,j-1,m) = -part(4,j-1,m)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j-1,m)
            part(5,j-1,m) = -part(5,j-1,m)
         endif
      endif
c set new position
      part(1,j-1,m) = dx
      part(2,j-1,m) = dy
      part(3,j-1,m) = dz
   10 continue
      nop = npp(m)
c deposit current for last particle
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
c deposit current
      dx = dx1*dzl
      dy = dx2*dzl
      dz = dyl*dzl
      vx = part(4,nop,m)
      vy = part(5,nop,m)
      vz = part(6,nop,m)
      cu(1,ml,m) = cu(1,ml,m) + vx*dx
      cu(2,ml,m) = cu(2,ml,m) + vy*dx
      cu(3,ml,m) = cu(3,ml,m) + vz*dx
      dx = dy1*dzl
      cu(1,ml+1,m) = cu(1,ml+1,m) + vx*dy
      cu(2,ml+1,m) = cu(2,ml+1,m) + vy*dy
      cu(3,ml+1,m) = cu(3,ml+1,m) + vz*dy
      dy = dy2*dzl
      cu(1,ml+2,m) = cu(1,ml+2,m) + vx*dz
      cu(2,ml+2,m) = cu(2,ml+2,m) + vy*dz
      cu(3,ml+2,m) = cu(3,ml+2,m) + vz*dz
      dz = amy*dzl
      cu(1,mn,m) = cu(1,mn,m) + vx*dx
      cu(2,mn,m) = cu(2,mn,m) + vy*dx
      cu(3,mn,m) = cu(3,mn,m) + vz*dx
      dx = dxl*dzl
      cu(1,mn+1,m) = cu(1,mn+1,m) + vx*dy
      cu(2,mn+1,m) = cu(2,mn+1,m) + vy*dy
      cu(3,mn+1,m) = cu(3,mn+1,m) + vz*dy
      dy = amx*dzl
      cu(1,mn+2,m) = cu(1,mn+2,m) + vx*dz
      cu(2,mn+2,m) = cu(2,mn+2,m) + vy*dz
      cu(3,mn+2,m) = cu(3,mn+2,m) + vz*dz
      dz = dyp*dzl
      cu(1,mp,m) = cu(1,mp,m) + vx*dx
      cu(2,mp,m) = cu(2,mp,m) + vy*dx
      cu(3,mp,m) = cu(3,mp,m) + vz*dx
      dx = dx1*amz
      cu(1,mp+1,m) = cu(1,mp+1,m) + vx*dy
      cu(2,mp+1,m) = cu(2,mp+1,m) + vy*dy
      cu(3,mp+1,m) = cu(3,mp+1,m) + vz*dy
      dy = dx2*amz
      cu(1,mp+2,m) = cu(1,mp+2,m) + vx*dz
      cu(2,mp+2,m) = cu(2,mp+2,m) + vy*dz
      cu(3,mp+2,m) = cu(3,mp+2,m) + vz*dz
      ml = ml + nxvy
      mn = mn + nxvy
      mp = mp + nxvy
      dz = dyl*amz
      cu(1,ml,m) = cu(1,ml,m) + vx*dx
      cu(2,ml,m) = cu(2,ml,m) + vy*dx
      cu(3,ml,m) = cu(3,ml,m) + vz*dx
      dx = dy1*amz
      cu(1,ml+1,m) = cu(1,ml+1,m) + vx*dy
      cu(2,ml+1,m) = cu(2,ml+1,m) + vy*dy
      cu(3,ml+1,m) = cu(3,ml+1,m) + vz*dy
      dy = dy2*amz
      cu(1,ml+2,m) = cu(1,ml+2,m) + vx*dz
      cu(2,ml+2,m) = cu(2,ml+2,m) + vy*dz
      cu(3,ml+2,m) = cu(3,ml+2,m) + vz*dz
      dz = amy*amz
      cu(1,mn,m) = cu(1,mn,m) + vx*dx
      cu(2,mn,m) = cu(2,mn,m) + vy*dx
      cu(3,mn,m) = cu(3,mn,m) + vz*dx
      dx = dxl*amz
      cu(1,mn+1,m) = cu(1,mn+1,m) + vx*dy
      cu(2,mn+1,m) = cu(2,mn+1,m) + vy*dy
      cu(3,mn+1,m) = cu(3,mn+1,m) + vz*dy
      dy = amx*amz
      cu(1,mn+2,m) = cu(1,mn+2,m) + vx*dz
      cu(2,mn+2,m) = cu(2,mn+2,m) + vy*dz
      cu(3,mn+2,m) = cu(3,mn+2,m) + vz*dz
      dz = dyp*amz
      cu(1,mp,m) = cu(1,mp,m) + vx*dx
      cu(2,mp,m) = cu(2,mp,m) + vy*dx
      cu(3,mp,m) = cu(3,mp,m) + vz*dx
      dx = dx1*dzp
      cu(1,mp+1,m) = cu(1,mp+1,m) + vx*dy
      cu(2,mp+1,m) = cu(2,mp+1,m) + vy*dy
      cu(3,mp+1,m) = cu(3,mp+1,m) + vz*dy
      dy = dx2*dzp
      cu(1,mp+2,m) = cu(1,mp+2,m) + vx*dz
      cu(2,mp+2,m) = cu(2,mp+2,m) + vy*dz
      cu(3,mp+2,m) = cu(3,mp+2,m) + vz*dz
      ml = ml + nxvy
      mn = mn + nxvy
      mp = mp + nxvy
      dz = dyl*dzp
      cu(1,ml,m) = cu(1,ml,m) + vx*dx
      cu(2,ml,m) = cu(2,ml,m) + vy*dx
      cu(3,ml,m) = cu(3,ml,m) + vz*dx
      dx = dy1*dzp
      cu(1,ml+1,m) = cu(1,ml+1,m) + vx*dy
      cu(2,ml+1,m) = cu(2,ml+1,m) + vy*dy
      cu(3,ml+1,m) = cu(3,ml+1,m) + vz*dy
      dy = dy2*dzp
      cu(1,ml+2,m) = cu(1,ml+2,m) + vx*dz
      cu(2,ml+2,m) = cu(2,ml+2,m) + vy*dz
      cu(3,ml+2,m) = cu(3,ml+2,m) + vz*dz
      dz = amy*dzp
      cu(1,mn,m) = cu(1,mn,m) + vx*dx
      cu(2,mn,m) = cu(2,mn,m) + vy*dx
      cu(3,mn,m) = cu(3,mn,m) + vz*dx
      dx = dxl*dzp
      cu(1,mn+1,m) = cu(1,mn+1,m) + vx*dy
      cu(2,mn+1,m) = cu(2,mn+1,m) + vy*dy
      cu(3,mn+1,m) = cu(3,mn+1,m) + vz*dy
      dy = amx*dzp
      cu(1,mn+2,m) = cu(1,mn+2,m) + vx*dz
      cu(2,mn+2,m) = cu(2,mn+2,m) + vy*dz
      cu(3,mn+2,m) = cu(3,mn+2,m) + vz*dz
      dz = dyp*dzp
      cu(1,mp,m) = cu(1,mp,m) + vx*dx
      cu(2,mp,m) = cu(2,mp,m) + vy*dx
      cu(3,mp,m) = cu(3,mp,m) + vz*dx
      cu(1,mp+1,m) = cu(1,mp+1,m) + vx*dy
      cu(2,mp+1,m) = cu(2,mp+1,m) + vy*dy
      cu(3,mp+1,m) = cu(3,mp+1,m) + vz*dy
      cu(1,mp+2,m) = cu(1,mp+2,m) + vx*dz
      cu(2,mp+2,m) = cu(2,mp+2,m) + vy*dz
      cu(3,mp+2,m) = cu(3,mp+2,m) + vz*dz
c advance position half a time-step
      dx = part(1,nop,m) + vx*dt
      dy = part(2,nop,m) + vy*dt
      dz = part(3,nop,m) + vz*dt
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
      return
      end
c-----------------------------------------------------------------------
      subroutine PSJOST32X(part,cu,npp,noff,nn,amxyz,qm,dt,nx,idimp,npma
     1x,mnblok,nxv,nypmx,nxvyzp,idds,npd,n81)
c for 3d code, this subroutine calculates particle current density
c using second-order spline interpolation, periodic boundaries,
c with short vectors over independent weights, and distributed data,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized distributed version with 1d addressing
c for 2D spatial decomposition
c 207 flops/particle, 249 loads, 244 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(n,m,l,i)=qci*(.75-dx**2)*(.75-dy**2)*(.75-dz**2)
c cu(n+1,m,l,i)=.5*qci*((.5+dx)**2)*(.75-dy**2)*(.75-dz**2)
c cu(n-1,m,l,i)=.5*qci*((.5-dx)**2)*(.75-dy**2)*(.75-dz**2)
c cu(n,m+1,l,i)=.5*qci*(.75-dx**2)*(.5+dy)**2*(.75-dz**2)
c cu(n+1,m+1,l,i)=.25*qci*((.5+dx)**2)*((.5+dy)**2)*(.75-dz**2)
c cu(n-1,m+1,l,i)=.25*qci*((.5-dx)**2)*((.5+dy)**2)*(.75-dz**2)
c cu(n,m-1,l,i)=.5*qci*(.75-dx**2)*((.5-dy)**2)*(.75-dz**2)
c cu(n+1,m-1,l,i)=.25*qci*((.5+dx)**2)*((.5-dy)**2)*(.75-dz**2)
c cu(n-1,m-1,l,i)=.25*qci*((.5-dx)**2)*((.5-dy)**2)*(.75-dz**2)
c cu(n,m,l+1,i)=.5*qci*(.75-dx**2)*(.75-dy**2)*((.5+dz)**2)
c cu(n+1,m,l+1,i)=.25*qci*((.5+dx)**2)*(.75-dy**2)*((.5+dz)**2)
c cu(n-1,m,l+1,i)=.25*qci*((.5-dx)**2)*(.75-dy**2)*((.5+dz)**2)
c cu(n,m+1,l+1,i)=.25*qci*(.75-dx**2)*(.5+dy)**2*((.5+dz)**2)
c cu(n+1,m+1,l+1,i)=.125*qci*((.5+dx)**2)*((.5+dy)**2)*((.5+dz)**2)
c cu(n-1,m+1,l+1,i=.125*qci*((.5-dx)**2)*((.5+dy)**2)*((.5+dz)**2)
c cu(n,m-1,l+1,i)=.25*qci*(.75-dx**2)*((.5-dy)**2)*((.5+dz)**2)
c cu(n+1,m-1,l+1,i)=.125*qci*((.5+dx)**2)*((.5-dy)**2)*((.5+dz)**2)
c cu(n-1,m-1,l+1,i)=.125*qci*((.5-dx)**2)*((.5-dy)**2)*((.5+dz)**2)
c cu(n,m,l-1,i)=.5*qci*(.75-dx**2)*(.75-dy**2)*((.5-dz)**2)
c cu(n+1,m,l-1,i)=.25*qci*((.5+dx)**2)*(.75-dy**2)*((.5-dz)**2)
c cu(n-1,m,l-1,i)=.25*qci*((.5-dx)**2)*(.75-dy**2)*((.5-dz)**2)
c cu(n,m+1,l-1,i)=.25*qci*(.75-dx**2)*(.5+dy)**2*((.5-dz)**2)
c cu(n+1,m+1,l-1,i)=.125*qci*((.5+dx)**2)*((.5+dy)**2)*((.5-dz)**2)
c cu(n-1,m+1,l-1,i)=.125*qci*((.5-dx)**2)*((.5+dy)**2)*((.5-dz)**2)
c cu(n,m-1,l-1,i)=.25*qci*(.75-dx**2)*((.5-dy)**2)*((.5-dz)**2)
c cu(n+1,m-1,l-1,i)=.125*qci*((.5+dx)**2)*((.5-dy)**2)*((.5-dz)**2)
c cu(n-1,m-1,l-1,i)=.125*qci*((.5-dx)**2)*((.5-dy)**2)*((.5-dz)**2)
c where n,m,l = nearest grid points and dx = x-n, dy = y-m, dz = z-l
c and qci = qm*vi, where i = x,y,z
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = x velocity of particle n in partition m
c part(5,n,m) = y velocity of particle n in partition m
c part(6,n,m) = z velocity of particle n in partition m
c cu(n,i,m) = ith component of current density at grid point (j,kk,ll),
c where n = j + nxv*(kk - 1) + nxv*nyv*(ll - 1)
c and kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nn = scratch address array for vectorized current deposition
c amxyz = scratch weight array for vectorized current deposition
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx = system length in x direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first virtual dimension of current array, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nxvyzp = nxv*nypmx*nzpmx, first dimension of charge array
c idds = dimensionality of domain decomposition
c npd = size of scratch buffers for vectorized push/current deposition
c n81 = number of independent weights
c version with spatial decomposition and 1d addressing scheme
      dimension part(idimp,npmax,mnblok), cu(3*nxvyzp,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      dimension nn(n81,npd,mnblok), amxyz(n81,npd,mnblok)
      zero = 0.
      anx = float(nx)
      nxvy = nxv*nypmx
      nxvyzp2 = nxvyzp + nxvyzp
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
      in = n + mm + l
      nn(1,i,m) = in
      nn(1+27,i,m) = in + nxvyzp
      nn(1+54,i,m) = in + nxvyzp2
      in = np + mm + l
      nn(2,i,m) = in
      nn(2+27,i,m) = in + nxvyzp
      nn(2+54,i,m) = in + nxvyzp2
      in = nl + mm + l
      nn(3,i,m) = in
      nn(3+27,i,m) = in + nxvyzp
      nn(3+54,i,m) = in + nxvyzp2
      in = n + mp + l
      nn(4,i,m) = in
      nn(4+27,i,m) = in + nxvyzp
      nn(4+54,i,m) = in + nxvyzp2
      in = np + mp + l
      nn(5,i,m) = in
      nn(5+27,i,m) = in + nxvyzp
      nn(5+54,i,m) = in + nxvyzp2
      in = nl + mp + l
      nn(6,i,m) = in
      nn(6+27,i,m) = in + nxvyzp
      nn(6+54,i,m) = in + nxvyzp2
      in = n + ml + l
      nn(7,i,m) = in
      nn(7+27,i,m) = in + nxvyzp
      nn(7+54,i,m) = in + nxvyzp2
      in = np + ml + l
      nn(8,i,m) = in
      nn(8+27,i,m) = in + nxvyzp
      nn(8+54,i,m) = in + nxvyzp2
      in = nl + ml + l
      nn(9,i,m) = in
      nn(9+27,i,m) = in + nxvyzp
      nn(9+54,i,m) = in + nxvyzp2
      in = n + mm + lp
      nn(10,i,m) = in
      nn(10+27,i,m) = in + nxvyzp
      nn(10+54,i,m) = in + nxvyzp2
      in = np + mm + lp
      nn(11,i,m) = in
      nn(11+27,i,m) = in + nxvyzp
      nn(11+54,i,m) = in + nxvyzp2
      in = nl + mm + lp
      nn(12,i,m) = in
      nn(12+27,i,m) = in + nxvyzp
      nn(12+54,i,m) = in + nxvyzp2
      in = n + mp + lp
      nn(13,i,m) = in
      nn(13+27,i,m) = in + nxvyzp
      nn(13+54,i,m) = in + nxvyzp2
      in = np + mp + lp
      nn(14,i,m) = in
      nn(14+27,i,m) = in + nxvyzp
      nn(14+54,i,m) = in + nxvyzp2
      in = nl + mp + lp
      nn(15,i,m) = in
      nn(15+27,i,m) = in + nxvyzp
      nn(15+54,i,m) = in + nxvyzp2
      in = n + ml + lp
      nn(16,i,m) = in
      nn(16+27,i,m) = in + nxvyzp
      nn(16+54,i,m) = in + nxvyzp2
      in = np + ml + lp
      nn(17,i,m) = in
      nn(17+27,i,m) = in + nxvyzp
      nn(17+54,i,m) = in + nxvyzp2
      in = nl + ml + lp
      nn(18,i,m) = in
      nn(18+27,i,m) = in + nxvyzp
      nn(18+54,i,m) = in + nxvyzp2
      in = n + mm + lm
      nn(19,i,m) = in
      nn(19+27,i,m) = in + nxvyzp
      nn(19+54,i,m) = in + nxvyzp2
      in = np + mm + lm
      nn(20,i,m) = in
      nn(20+27,i,m) = in + nxvyzp
      nn(20+54,i,m) = in + nxvyzp2
      in = nl + mm + lm
      nn(21,i,m) = in
      nn(21+27,i,m) = in + nxvyzp
      nn(21+54,i,m) = in + nxvyzp2
      in = n + mp + lm
      nn(22,i,m) = in
      nn(22+27,i,m) = in + nxvyzp
      nn(22+54,i,m) = in + nxvyzp2
      in = np + mp + lm
      nn(23,i,m) = in
      nn(23+27,i,m) = in + nxvyzp
      nn(23+54,i,m) = in + nxvyzp2
      in = nl + mp + lm
      nn(24,i,m) = in
      nn(24+27,i,m) = in + nxvyzp
      nn(24+54,i,m) = in + nxvyzp2
      in = n + ml + lm
      nn(25,i,m) = in
      nn(25+27,i,m) = in + nxvyzp
      nn(25+54,i,m) = in + nxvyzp2
      in = np + ml + lm
      nn(26,i,m) = in
      nn(26+27,i,m) = in + nxvyzp
      nn(26+54,i,m) = in + nxvyzp2
      in = nl + ml + lm
      nn(27,i,m) = in
      nn(27+27,i,m) = in + nxvyzp
      nn(27+54,i,m) = in + nxvyzp2
      dx = amx*amy*amz
      dy = dxp*amy*amz
      dz = dxl*amy*amz
      vx = part(4,i+jb,m)
      vy = part(5,i+jb,m)
      vz = part(6,i+jb,m)
      amxyz(1,i,m) = vx*dx
      amxyz(2,i,m) = vx*dy
      amxyz(3,i,m) = vx*dz
      amxyz(1+27,i,m) = vy*dx
      amxyz(2+27,i,m) = vy*dy
      amxyz(3+27,i,m) = vy*dz
      amxyz(1+54,i,m) = vz*dx
      amxyz(2+54,i,m) = vz*dy
      amxyz(3+54,i,m) = vz*dz
      dx = amx*dyp*amz
      dy = dxp*dyp*amz
      dz = dxl*dyp*amz
      amxyz(4,i,m) = vx*dx
      amxyz(5,i,m) = vx*dy
      amxyz(6,i,m) = vx*dz
      amxyz(4+27,i,m) = vy*dx
      amxyz(5+27,i,m) = vy*dy
      amxyz(6+27,i,m) = vy*dz
      amxyz(4+54,i,m) = vz*dx
      amxyz(5+54,i,m) = vz*dy
      amxyz(6+54,i,m) = vz*dz
      dx = amx*dyl*amz
      dy = dxp*dyl*amz
      dz = dxl*dyl*amz
      amxyz(7,i,m) = vx*dx
      amxyz(8,i,m) = vx*dy
      amxyz(9,i,m) = vx*dz
      amxyz(7+27,i,m) = vy*dx
      amxyz(8+27,i,m) = vy*dy
      amxyz(9+27,i,m) = vy*dz
      amxyz(7+54,i,m) = vz*dx
      amxyz(8+54,i,m) = vz*dy
      amxyz(9+54,i,m) = vz*dz
      dx = amx*amy*dzp
      dy = dxp*amy*dzp
      dz = dxl*amy*dzp
      amxyz(10,i,m) = vx*dx
      amxyz(11,i,m) = vx*dy
      amxyz(12,i,m) = vx*dz
      amxyz(10+27,i,m) = vy*dx
      amxyz(11+27,i,m) = vy*dy
      amxyz(12+27,i,m) = vy*dz
      amxyz(10+54,i,m) = vz*dx
      amxyz(11+54,i,m) = vz*dy
      amxyz(12+54,i,m) = vz*dz
      dx = amx*dyp*dzp
      dy = dxp*dyp*dzp
      dz = dxl*dyp*dzp
      amxyz(13,i,m) = vx*dx
      amxyz(14,i,m) = vx*dy
      amxyz(15,i,m) = vx*dz
      amxyz(13+27,i,m) = vy*dx
      amxyz(14+27,i,m) = vy*dy
      amxyz(15+27,i,m) = vy*dz
      amxyz(13+54,i,m) = vz*dx
      amxyz(14+54,i,m) = vz*dy
      amxyz(15+54,i,m) = vz*dz
      dx = amx*dyl*dzp
      dy = dxp*dyl*dzp
      dz = dxl*dyl*dzp
      amxyz(16,i,m) = vx*dx
      amxyz(17,i,m) = vx*dy
      amxyz(18,i,m) = vx*dz
      amxyz(16+27,i,m) = vy*dx
      amxyz(17+27,i,m) = vy*dy
      amxyz(18+27,i,m) = vy*dz
      amxyz(16+54,i,m) = vz*dx
      amxyz(17+54,i,m) = vz*dy
      amxyz(18+54,i,m) = vz*dz
      dx = amx*amy*dzl
      dy = dxp*amy*dzl
      dz = dxl*amy*dzl
      amxyz(19,i,m) = vx*dx
      amxyz(20,i,m) = vx*dy
      amxyz(21,i,m) = vx*dz
      amxyz(19+27,i,m) = vy*dx
      amxyz(20+27,i,m) = vy*dy
      amxyz(21+27,i,m) = vy*dz
      amxyz(19+54,i,m) = vz*dx
      amxyz(20+54,i,m) = vz*dy
      amxyz(21+54,i,m) = vz*dz
      dx = amx*dyp*dzl
      dy = dxp*dyp*dzl
      dz = dxl*dyp*dzl
      amxyz(22,i,m) = vx*dx
      amxyz(23,i,m) = vx*dy
      amxyz(24,i,m) = vx*dz
      amxyz(22+27,i,m) = vy*dx
      amxyz(23+27,i,m) = vy*dy
      amxyz(24+27,i,m) = vy*dz
      amxyz(22+54,i,m) = vz*dx
      amxyz(23+54,i,m) = vz*dy
      amxyz(24+54,i,m) = vz*dz
      dx = amx*dyl*dzl
      dy = dxp*dyl*dzl
      dz = dxl*dyl*dzl
      amxyz(25,i,m) = vx*dx
      amxyz(26,i,m) = vx*dy
      amxyz(27,i,m) = vx*dz
      amxyz(25+27,i,m) = vy*dx
      amxyz(26+27,i,m) = vy*dy
      amxyz(27+27,i,m) = vy*dz
      amxyz(25+54,i,m) = vz*dx
      amxyz(26+54,i,m) = vz*dy
      amxyz(27+54,i,m) = vz*dz
c advance position half a time-step
      dx = part(1,i+jb,m) + vx*dt
      dy = part(2,i+jb,m) + vy*dt
      dz = part(3,i+jb,m) + vz*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,i+jb,m) = dx
      part(2,i+jb,m) = dy
      part(3,i+jb,m) = dz
   10 continue
c deposit current
      do 30 i = 1, npb
cdir$ ivdep
      do 20 k = 1, 81
      cu(nn(k,i,m),m) = cu(nn(k,i,m),m) + amxyz(k,i,m)
   20 continue
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSJOST32X(part,cu,npp,noff,nn,amxyz,qm,dt,nx,ny,nz,idi
     1mp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n81,ipbc)
c for 3d code, this subroutine calculates particle current density
c using second-order spline interpolation, for distributed data
c with 2D spatial decomposition
c with short vectors over independent weights,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized version with guard cells and 1d addressing
c 189 flops/particle, 249 loads, 244 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m,l)=qci*(.75-dx**2)*(.75-dy**2)*(.75-dz**2)
c cu(i,n+1,m,l)=.5*qci*((.5+dx)**2)*(.75-dy**2)*(.75-dz**2)
c cu(i,n-1,m,l)=.5*qci*((.5-dx)**2)*(.75-dy**2)*(.75-dz**2)
c cu(i,n,m+1,l)=.5*qci*(.75-dx**2)*(.5+dy)**2*(.75-dz**2)
c cu(i,n+1,m+1,l)=.25*qci*((.5+dx)**2)*((.5+dy)**2)*(.75-dz**2)
c cu(i,n-1,m+1,l)=.25*qci*((.5-dx)**2)*((.5+dy)**2)*(.75-dz**2)
c cu(i,n,m-1,l)=.5*qci*(.75-dx**2)*((.5-dy)**2)*(.75-dz**2)
c cu(i,n+1,m-1,l)=.25*qci*((.5+dx)**2)*((.5-dy)**2)*(.75-dz**2)
c cu(i,n-1,m-1,l)=.25*qci*((.5-dx)**2)*((.5-dy)**2)*(.75-dz**2)
c cu(i,n,m,l+1)=.5*qci*(.75-dx**2)*(.75-dy**2)*((.5+dz)**2)
c cu(i,n+1,m,l+1)=.25*qci*((.5+dx)**2)*(.75-dy**2)*((.5+dz)**2)
c cu(i,n-1,m,l+1)=.25*qci*((.5-dx)**2)*(.75-dy**2)*((.5+dz)**2)
c cu(i,n,m+1,l+1)=.25*qci*(.75-dx**2)*(.5+dy)**2*((.5+dz)**2)
c cu(i,n+1,m+1,l+1)=.125*qci*((.5+dx)**2)*((.5+dy)**2)*((.5+dz)**2)
c cu(i,n-1,m+1,l+1)=.125*qci*((.5-dx)**2)*((.5+dy)**2)*((.5+dz)**2)
c cu(i,n,m-1,l+1)=.25*qci*(.75-dx**2)*((.5-dy)**2)*((.5+dz)**2)
c cu(i,n+1,m-1,l+1)=.125*qci*((.5+dx)**2)*((.5-dy)**2)*((.5+dz)**2)
c cu(i,n-1,m-1,l+1)=.125*qci*((.5-dx)**2)*((.5-dy)**2)*((.5+dz)**2)
c cu(i,n,m,l-1)=.5*qci*(.75-dx**2)*(.75-dy**2)*((.5-dz)**2)
c cu(i,n+1,m,l-1)=.25*qci*((.5+dx)**2)*(.75-dy**2)*((.5-dz)**2)
c cu(i,n-1,m,l-1)=.25*qci*((.5-dx)**2)*(.75-dy**2)*((.5-dz)**2)
c cu(i,n,m+1,l-1)=.25*qci*(.75-dx**2)*(.5+dy)**2*((.5-dz)**2)
c cu(i,n+1,m+1,l-1)=.125*qci*((.5+dx)**2)*((.5+dy)**2)*((.5-dz)**2)
c cu(i,n-1,m+1,l-1)=.125*qci*((.5-dx)**2)*((.5+dy)**2)*((.5-dz)**2)
c cu(i,n,m-1,l-1)=.25*qci*(.75-dx**2)*((.5-dy)**2)*((.5-dz)**2)
c cu(i,n+1,m-1,l-1)=.125*qci*((.5+dx)**2)*((.5-dy)**2)*((.5-dz)**2)
c cu(i,n-1,m-1,l-1)=.125*qci*((.5-dx)**2)*((.5-dy)**2)*((.5-dz)**2)
c where n,m,l = nearest grid points and dx = x-n, dy = y-m, dz = z-l
c and qci = qm*vi, where i = x,y,z
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = x velocity of particle n in partition m
c part(5,n,m) = y velocity of particle n in partition m
c part(6,n,m) = z velocity of particle n in partition m
c cu(i,n,m) = ith component of current density at grid point (j,kk,ll),
c where n = j + nxv*k + nxv*nyv*ll
c and kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nn = scratch address array for vectorized current deposition
c amxyz = scratch weight array for vectorized current deposition
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first virtual dimension of current array, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nxvyzp = nxv*nypmx*nzpmx, first dimension of charge array
c idds = dimensionality of domain decomposition
c npd = size of scratch buffers for vectorized push/current deposition
c n81 = number of independent weights
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      dimension part(idimp,npmax,mnblok), cu(3*nxvyzp,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      dimension nn(n81,npd,mnblok), amxyz(n81,npd,mnblok)
      nxv3 = 3*nxv
      nxvy3 = nxv3*nypmx
      qmh = .5*qm
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
      n3 = 3*n + 1
      mm = nxv3*(mm - mnoff) + nxvy3*(ll - lnoff)
      amx = qm*(.75 - dx*dx)
      dxl = qmh*(.5 - dx)**2
      dxp = qmh*(.5 + dx)**2
      ml = mm + n3
      dyl = .5*(.5 - dy)**2
      amy = .75 - dy*dy
      dyp = .5*(.5 + dy)**2
      mn = ml + nxv3
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
      mp = mn + nxv3
      nn(1,i,m) = ml
      nn(2,i,m) = ml + 1
      nn(3,i,m) = ml + 2
      nn(4,i,m) = ml + 3
      nn(5,i,m) = ml + 4
      nn(6,i,m) = ml + 5
      nn(7,i,m) = ml + 6
      nn(8,i,m) = ml + 7
      nn(9,i,m) = ml + 8
      nn(10,i,m) = mn
      nn(11,i,m) = mn + 1
      nn(12,i,m) = mn + 2
      nn(13,i,m) = mn + 3
      nn(14,i,m) = mn + 4
      nn(15,i,m) = mn + 5
      nn(16,i,m) = mn + 6
      nn(17,i,m) = mn + 7
      nn(18,i,m) = mn + 8
      nn(19,i,m) = mp
      nn(20,i,m) = mp + 1
      nn(21,i,m) = mp + 2
      nn(22,i,m) = mp + 3
      nn(23,i,m) = mp + 4
      nn(24,i,m) = mp + 5
      nn(25,i,m) = mp + 6
      nn(26,i,m) = mp + 7
      nn(27,i,m) = mp + 8
      ml = ml + nxvy3
      mn = mn + nxvy3
      mp = mp + nxvy3
      nn(28,i,m) = ml
      nn(29,i,m) = ml + 1
      nn(30,i,m) = ml + 2
      nn(31,i,m) = ml + 3
      nn(32,i,m) = ml + 4
      nn(33,i,m) = ml + 5
      nn(34,i,m) = ml + 6
      nn(35,i,m) = ml + 7
      nn(36,i,m) = ml + 8
      nn(37,i,m) = mn
      nn(38,i,m) = mn + 1
      nn(39,i,m) = mn + 2
      nn(40,i,m) = mn + 3
      nn(41,i,m) = mn + 4
      nn(42,i,m) = mn + 5
      nn(43,i,m) = mn + 6
      nn(44,i,m) = mn + 7
      nn(45,i,m) = mn + 8
      nn(46,i,m) = mp
      nn(47,i,m) = mp + 1
      nn(48,i,m) = mp + 2
      nn(49,i,m) = mp + 3
      nn(50,i,m) = mp + 4
      nn(51,i,m) = mp + 5
      nn(52,i,m) = mp + 6
      nn(53,i,m) = mp + 7
      nn(54,i,m) = mp + 8
      ml = ml + nxvy3
      mn = mn + nxvy3
      mp = mp + nxvy3
      nn(55,i,m) = ml
      nn(56,i,m) = ml + 1
      nn(57,i,m) = ml + 2
      nn(58,i,m) = ml + 3
      nn(59,i,m) = ml + 4
      nn(60,i,m) = ml + 5
      nn(61,i,m) = ml + 6
      nn(62,i,m) = ml + 7
      nn(63,i,m) = ml + 8
      nn(64,i,m) = mn
      nn(65,i,m) = mn + 1
      nn(66,i,m) = mn + 2
      nn(67,i,m) = mn + 3
      nn(68,i,m) = mn + 4
      nn(69,i,m) = mn + 5
      nn(70,i,m) = mn + 6
      nn(71,i,m) = mn + 7
      nn(72,i,m) = mn + 8
      nn(73,i,m) = mp
      nn(74,i,m) = mp + 1
      nn(75,i,m) = mp + 2
      nn(76,i,m) = mp + 3
      nn(77,i,m) = mp + 4
      nn(78,i,m) = mp + 5
      nn(79,i,m) = mp + 6
      nn(80,i,m) = mp + 7
      nn(81,i,m) = mp + 8
c deposit current
      dx = dx1*dzl
      dy = dx2*dzl
      dz = dyl*dzl
      vx = part(4,i+jb,m)
      vy = part(5,i+jb,m)
      vz = part(6,i+jb,m)
      amxyz(1,i,m) = vx*dx
      amxyz(2,i,m) = vy*dx
      amxyz(3,i,m) = vz*dx
      dx = dy1*dzl
      amxyz(4,i,m) = vx*dy
      amxyz(5,i,m) = vy*dy
      amxyz(6,i,m) = vz*dy
      dy = dy2*dzl
      amxyz(7,i,m) = vx*dz
      amxyz(8,i,m) = vy*dz
      amxyz(9,i,m) = vz*dz
      dz = amy*dzl
      amxyz(10,i,m) = vx*dx
      amxyz(11,i,m) = vy*dx
      amxyz(12,i,m) = vz*dx
      dx = dxl*dzl
      amxyz(13,i,m) = vx*dy
      amxyz(14,i,m) = vy*dy
      amxyz(15,i,m) = vz*dy
      dy = amx*dzl
      amxyz(16,i,m) = vx*dz
      amxyz(17,i,m) = vy*dz
      amxyz(18,i,m) = vz*dz
      dz = dyp*dzl
      amxyz(19,i,m) = vx*dx
      amxyz(20,i,m) = vy*dx
      amxyz(21,i,m) = vz*dx
      dx = dx1*amz
      amxyz(22,i,m) = vx*dy
      amxyz(23,i,m) = vy*dy
      amxyz(24,i,m) = vz*dy
      dy = dx2*amz
      amxyz(25,i,m) = vx*dz
      amxyz(26,i,m) = vy*dz
      amxyz(27,i,m) = vz*dz
      dz = dyl*amz
      amxyz(28,i,m) = vx*dx
      amxyz(29,i,m) = vy*dx
      amxyz(30,i,m) = vz*dx
      dx = dy1*amz
      amxyz(31,i,m) = vx*dy
      amxyz(32,i,m) = vy*dy
      amxyz(33,i,m) = vz*dy
      dy = dy2*amz
      amxyz(34,i,m) = vx*dz
      amxyz(35,i,m) = vy*dz
      amxyz(36,i,m) = vz*dz
      dz = amy*amz
      amxyz(37,i,m) = vx*dx
      amxyz(38,i,m) = vy*dx
      amxyz(39,i,m) = vz*dx
      dx = dxl*amz
      amxyz(40,i,m) = vx*dy
      amxyz(41,i,m) = vy*dy
      amxyz(42,i,m) = vz*dy
      dy = amx*amz
      amxyz(43,i,m) = vx*dz
      amxyz(44,i,m) = vy*dz
      amxyz(45,i,m) = vz*dz
      dz = dyp*amz
      amxyz(46,i,m) = vx*dx
      amxyz(47,i,m) = vy*dx
      amxyz(48,i,m) = vz*dx
      dx = dx1*dzp
      amxyz(49,i,m) = vx*dy
      amxyz(50,i,m) = vy*dy
      amxyz(51,i,m) = vz*dy
      dy = dx2*dzp
      amxyz(52,i,m) = vx*dz
      amxyz(53,i,m) = vy*dz
      amxyz(54,i,m) = vz*dz
      dz = dyl*dzp
      amxyz(55,i,m) = vx*dx
      amxyz(56,i,m) = vy*dx
      amxyz(57,i,m) = vz*dx
      dx = dy1*dzp
      amxyz(58,i,m) = vx*dy
      amxyz(59,i,m) = vy*dy
      amxyz(60,i,m) = vz*dy
      dy = dy2*dzp
      amxyz(61,i,m) = vx*dz
      amxyz(62,i,m) = vy*dz
      amxyz(63,i,m) = vz*dz
      dz = amy*dzp
      amxyz(64,i,m) = vx*dx
      amxyz(65,i,m) = vy*dx
      amxyz(66,i,m) = vz*dx
      dx = dxl*dzp
      amxyz(67,i,m) = vx*dy
      amxyz(68,i,m) = vy*dy
      amxyz(69,i,m) = vz*dy
      dy = amx*dzp
      amxyz(70,i,m) = vx*dz
      amxyz(71,i,m) = vy*dz
      amxyz(72,i,m) = vz*dz
      dz = dyp*dzp
      amxyz(73,i,m) = vx*dx
      amxyz(74,i,m) = vy*dx
      amxyz(75,i,m) = vz*dx
      amxyz(76,i,m) = vx*dy
      amxyz(77,i,m) = vy*dy
      amxyz(78,i,m) = vz*dy
      amxyz(79,i,m) = vx*dz
      amxyz(80,i,m) = vy*dz
      amxyz(81,i,m) = vz*dz
c advance position half a time-step
      dx = part(1,i+jb,m) + vx*dt
      dy = part(2,i+jb,m) + vy*dt
      dz = part(3,i+jb,m) + vz*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i+jb,m)
            part(4,i+jb,m) = -part(4,i+jb,m)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,i+jb,m)
            part(5,i+jb,m) = -part(5,i+jb,m)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = part(3,i+jb,m)
            part(6,i+jb,m) = -part(6,i+jb,m)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i+jb,m)
            part(4,i+jb,m) = -part(4,i+jb,m)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,i+jb,m)
            part(5,i+jb,m) = -part(5,i+jb,m)
         endif
      endif
c set new position
      part(1,i+jb,m) = dx
      part(2,i+jb,m) = dy
      part(3,i+jb,m) = dz
   10 continue
c deposit current
      do 30 i = 1, npb
cdir$ ivdep
      do 20 k = 1, 81
      cu(nn(k,i,m),m) = cu(nn(k,i,m),m) + amxyz(k,i,m)
   20 continue
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PJDOST32L(part,cux,cuy,cuz,npp,noff,qm,dt,nx,idimp,npma
     1x,mnblok,nxv,nypmx,nzpmx,idds)
c for 3d code, this subroutine calculates particle current density
c using first-order linear interpolation, periodic boundaries
c and distributed data with 2D spatial decomposition
c in addition, particle positions are advanced a half time-step
c baseline scalar distributed version
c 69 flops/particle, 30 loads, 27 stores
c input: all, output: part, cux, cuy, cuz
c current density is approximated by values at the nearest grid points
c cui(n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
c cui(n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
c cui(n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
c cui(n+1,m+1,l)=qci*dx*dy*(1.-dz)
c cui(n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
c cui(n+1,m,l+1)=qci*dx*(1.-dy)*dz
c cui(n,m+1,l+1)=qci*(1.-dx)*dy*dz
c cui(n+1,m+1,l+1)=qci*dx*dy*dz
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c and qci = qm*vi, where i = x,y,z
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = x velocity of particle n in partition m
c part(5,n,m) = y velocity of particle n in partition m
c part(6,n,m) = z velocity of particle n in partition m
c cui(j,k,l,m)=ith component of current density at grid point (j,kk,ll)
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx = system length in x direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first dimension of current array, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
      dimension part(idimp,npmax,mnblok), cux(nxv,nypmx,nzpmx,mnblok)
      dimension cuy(nxv,nypmx,nzpmx,mnblok), cuz(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      zero = 0.
      anx = float(nx)
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
c deposit current
      dx = amx*amy*amz
      dy = dxp*amy*amz
      vx = part(4,j,m)
      vy = part(5,j,m)
      vz = part(6,j,m)
      cux(nn,mm,ll,m) = cux(nn,mm,ll,m) + vx*dx
      cux(np,mm,ll,m) = cux(np,mm,ll,m) + vx*dy
      cuy(nn,mm,ll,m) = cuy(nn,mm,ll,m) + vy*dx
      cuy(np,mm,ll,m) = cuy(np,mm,ll,m) + vy*dy
      cuz(nn,mm,ll,m) = cuz(nn,mm,ll,m) + vz*dx
      cuz(np,mm,ll,m) = cuz(np,mm,ll,m) + vz*dy
      dx = amx*dyp*amz
      dy = dxp*dyp*amz
      cux(nn,mp,ll,m) = cux(nn,mp,ll,m) + vx*dx
      cux(np,mp,ll,m) = cux(np,mp,ll,m) + vx*dy
      cuy(nn,mp,ll,m) = cuy(nn,mp,ll,m) + vy*dx
      cuy(np,mp,ll,m) = cuy(np,mp,ll,m) + vy*dy
      cuz(nn,mp,ll,m) = cuz(nn,mp,ll,m) + vz*dx
      cuz(np,mp,ll,m) = cuz(np,mp,ll,m) + vz*dy
      dx = amx*amy*dzp
      dy = dxp*amy*dzp
      cux(nn,mm,lp,m) = cux(nn,mm,lp,m) + vx*dx
      cux(np,mm,lp,m) = cux(np,mm,lp,m) + vx*dy
      cuy(nn,mm,lp,m) = cuy(nn,mm,lp,m) + vy*dx
      cuy(np,mm,lp,m) = cuy(np,mm,lp,m) + vy*dy
      cuz(nn,mm,lp,m) = cuz(nn,mm,lp,m) + vz*dx
      cuz(np,mm,lp,m) = cuz(np,mm,lp,m) + vz*dy
      dx = amx*dyp*dzp
      dy = dxp*dyp*dzp
      cux(nn,mp,lp,m) = cux(nn,mp,lp,m) + vx*dx
      cux(np,mp,lp,m) = cux(np,mp,lp,m) + vx*dy
      cuy(nn,mp,lp,m) = cuy(nn,mp,lp,m) + vy*dx
      cuy(np,mp,lp,m) = cuy(np,mp,lp,m) + vy*dy
      cuz(nn,mp,lp,m) = cuz(nn,mp,lp,m) + vz*dx
      cuz(np,mp,lp,m) = cuz(np,mp,lp,m) + vz*dy
c advance position half a time-step
      dx = part(1,j,m) + vx*dt
      dy = part(2,j,m) + vy*dt
      dz = part(3,j,m) + vz*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j,m) = dx
      part(2,j,m) = dy
      part(3,j,m) = dz
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGJPOST32L(part,cu,npp,noff,qm,dt,nx,ny,nz,idimp,npmax,
     1mnblok,nxv,nypmx,nzpmx,idds,ipbc)
c for 3d code, this subroutine calculates particle current density
c using first-order linear interpolation, for distributed data
c with 2D spatial decomposition
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 65 flops/particle, 30 loads, 27 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
c cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
c cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
c cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
c cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
c cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
c cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
c cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
c where n,m,l = nearest grid points and dx = x-n, dy = y-m, dz = z-l
c and qci = qm*vi, where i = x,y,z
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = x velocity of particle n in partition m
c part(5,n,m) = y velocity of particle n in partition m
c part(6,n,m) = z velocity of particle n in partition m
c cu(i,j,k,l,m) = ith component of current density at grid point
c (j,kk,ll), where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first dimension of current array, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      dimension part(idimp,npmax,mnblok), cu(3,nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
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
c deposit current
      dx = amx*amz
      dy = amy*amz
      vx = part(4,j,m)
      vy = part(5,j,m)
      vz = part(6,j,m)
      cu(1,nn,mm,ll,m) = cu(1,nn,mm,ll,m) + vx*dx
      cu(2,nn,mm,ll,m) = cu(2,nn,mm,ll,m) + vy*dx
      cu(3,nn,mm,ll,m) = cu(3,nn,mm,ll,m) + vz*dx
      dx = dyp*amz
      cu(1,np,mm,ll,m) = cu(1,np,mm,ll,m) + vx*dy
      cu(2,np,mm,ll,m) = cu(2,np,mm,ll,m) + vy*dy
      cu(3,np,mm,ll,m) = cu(3,np,mm,ll,m) + vz*dy
      dy = dx1*amz
      cu(1,nn,mp,ll,m) = cu(1,nn,mp,ll,m) + vx*dx
      cu(2,nn,mp,ll,m) = cu(2,nn,mp,ll,m) + vy*dx
      cu(3,nn,mp,ll,m) = cu(3,nn,mp,ll,m) + vz*dx
      dx = amx*dzp
      cu(1,np,mp,ll,m) = cu(1,np,mp,ll,m) + vx*dy
      cu(2,np,mp,ll,m) = cu(2,np,mp,ll,m) + vy*dy
      cu(3,np,mp,ll,m) = cu(3,np,mp,ll,m) + vz*dy
      dy = amy*dzp
      cu(1,nn,mm,lp,m) = cu(1,nn,mm,lp,m) + vx*dx
      cu(2,nn,mm,lp,m) = cu(2,nn,mm,lp,m) + vy*dx
      cu(3,nn,mm,lp,m) = cu(3,nn,mm,lp,m) + vz*dx
      dx = dyp*dzp
      cu(1,np,mm,lp,m) = cu(1,np,mm,lp,m) + vx*dy
      cu(2,np,mm,lp,m) = cu(2,np,mm,lp,m) + vy*dy
      cu(3,np,mm,lp,m) = cu(3,np,mm,lp,m) + vz*dy
      dy = dx1*dzp
      cu(1,nn,mp,lp,m) = cu(1,nn,mp,lp,m) + vx*dx
      cu(2,nn,mp,lp,m) = cu(2,nn,mp,lp,m) + vy*dx
      cu(3,nn,mp,lp,m) = cu(3,nn,mp,lp,m) + vz*dx
      cu(1,np,mp,lp,m) = cu(1,np,mp,lp,m) + vx*dy
      cu(2,np,mp,lp,m) = cu(2,np,mp,lp,m) + vy*dy
      cu(3,np,mp,lp,m) = cu(3,np,mp,lp,m) + vz*dy
c advance position half a time-step
      dx = part(1,j,m) + vx*dt
      dy = part(2,j,m) + vy*dt
      dz = part(3,j,m) + vz*dt
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
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSJPOST32L(part,cu,npp,noff,qm,dt,nx,ny,nz,idimp,npmax
     1,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
c for 3d code, this subroutine calculates particle current density
c using first-order linear interpolation, for distributed data
c with 2D spatial decomposition
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 65 flops/particle, 30 loads, 27 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
c cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
c cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
c cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
c cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
c cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
c cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
c cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
c where n,m,l = nearest grid points and dx = x-n, dy = y-m, dz = z-l
c and qci = qm*vi, where i = x,y,z
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = x velocity of particle n in partition m
c part(5,n,m) = y velocity of particle n in partition m
c part(6,n,m) = z velocity of particle n in partition m
c cu(i,n,m) = ith component of current density at grid point (j,kk,ll),
c where n = j + nxv*(k-1) + nxv*nyv*(ll-1)
c and kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first virtual dimension of current array, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nxyzp = dimension of charge array, must be >= nxv*nypmx*nzpmx
c idds = dimensionality of domain decomposition
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      dimension part(idimp,npmax,mnblok), cu(3,nxyzp,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
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
c deposit current
      dx = amx*amz
      dz = amy*amz
      vx = part(4,j-1,m)
      vy = part(5,j-1,m)
      vz = part(6,j-1,m)
      dxp = cu(1,mm,m) + vx*dx
      dx2 = cu(2,mm,m) + vy*dx
      dx3 = cu(3,mm,m) + vz*dx
      dx = cu(1,mm+1,m) + vx*dz
      dy = cu(2,mm+1,m) + vy*dz
      dz = cu(3,mm+1,m) + vz*dz
      cu(1,mm,m) = dxp
      cu(2,mm,m) = dx2
      cu(3,mm,m) = dx3
      cu(1,mm+1,m) = dx
      cu(2,mm+1,m) = dy
      cu(3,mm+1,m) = dz
      dx = dyp*amz
      dz = dx1*amz
      dxp = cu(1,mp,m) + vx*dx
      dx2 = cu(2,mp,m) + vy*dx
      dx3 = cu(3,mp,m) + vz*dx
      dx = cu(1,mp+1,m) + vx*dz
      dy = cu(2,mp+1,m) + vy*dz
      dz = cu(3,mp+1,m) + vz*dz
      cu(1,mp,m) = dxp
      cu(2,mp,m) = dx2
      cu(3,mp,m) = dx3
      cu(1,mp+1,m) = dx
      cu(2,mp+1,m) = dy
      cu(3,mp+1,m) = dz
      dx = amx*dzp
      dz = amy*dzp
      dxp = cu(1,ll,m) + vx*dx
      dx2 = cu(2,ll,m) + vy*dx
      dx3 = cu(3,ll,m) + vz*dx
      dx = cu(1,ll+1,m) + vx*dz
      dy = cu(2,ll+1,m) + vy*dz
      dz = cu(3,ll+1,m) + vz*dz
      cu(1,ll,m) = dxp
      cu(2,ll,m) = dx2
      cu(3,ll,m) = dx3
      cu(1,ll+1,m) = dx
      cu(2,ll+1,m) = dy
      cu(3,ll+1,m) = dz
      dx = dyp*dzp
      dz = dx1*dzp
      dxp = cu(1,lp,m) + vx*dx
      dx2 = cu(2,lp,m) + vy*dx
      dx3 = cu(3,lp,m) + vz*dx
      dx = cu(1,lp+1,m) + vx*dz
      dy = cu(2,lp+1,m) + vy*dz
      dz = cu(3,lp+1,m) + vz*dz
      cu(1,lp,m) = dxp
      cu(2,lp,m) = dx2
      cu(3,lp,m) = dx3
      cu(1,lp+1,m) = dx
      cu(2,lp+1,m) = dy
      cu(3,lp+1,m) = dz
c advance position half a time-step
      dx = part(1,j-1,m) + vx*dt
      dy = part(2,j-1,m) + vy*dt
      dz = part(3,j-1,m) + vz*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1,m)
            part(4,j-1,m) = -part(4,j-1,m)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j-1,m)
            part(5,j-1,m) = -part(5,j-1,m)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = part(3,j-1,m)
            part(6,j-1,m) = -part(6,j-1,m)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,j-1,m)
            part(4,j-1,m) = -part(4,j-1,m)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,j-1,m)
            part(5,j-1,m) = -part(5,j-1,m)
         endif
      endif
c set new position
      part(1,j-1,m) = dx
      part(2,j-1,m) = dy
      part(3,j-1,m) = dz
   10 continue
      nop = npp(m)
c deposit current for last particle
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
c deposit current
      dx = amx*amz
      dy = amy*amz
      vx = part(4,nop,m)
      vy = part(5,nop,m)
      vz = part(6,nop,m)
      cu(1,mm,m) = cu(1,mm,m) + vx*dx
      cu(2,mm,m) = cu(2,mm,m) + vy*dx
      cu(3,mm,m) = cu(3,mm,m) + vz*dx
      cu(1,mm+1,m) = cu(1,mm+1,m) + vx*dy
      cu(2,mm+1,m) = cu(2,mm+1,m) + vy*dy
      cu(3,mm+1,m) = cu(3,mm+1,m) + vz*dy
      dx = dyp*amz
      dy = dx1*amz
      cu(1,mp,m) = cu(1,mp,m) + vx*dx
      cu(2,mp,m) = cu(2,mp,m) + vy*dx
      cu(3,mp,m) = cu(3,mp,m) + vz*dx
      cu(1,mp+1,m) = cu(1,mp+1,m) + vx*dy
      cu(2,mp+1,m) = cu(2,mp+1,m) + vy*dy
      cu(3,mp+1,m) = cu(3,mp+1,m) + vz*dy
      dx = amx*dzn
      dy = amy*dzn
      cu(1,ll,m) = cu(1,ll,m) + vx*dx
      cu(2,ll,m) = cu(2,ll,m) + vy*dx
      cu(3,ll,m) = cu(3,ll,m) + vz*dx
      cu(1,ll+1,m) = cu(1,ll+1,m) + vx*dy
      cu(2,ll+1,m) = cu(2,ll+1,m) + vy*dy
      cu(3,ll+1,m) = cu(3,ll+1,m) + vz*dy
      dx = dyp*dzn
      dy = dx1*dzn
      cu(1,lp,m) = cu(1,lp,m) + vx*dx
      cu(2,lp,m) = cu(2,lp,m) + vy*dx
      cu(3,lp,m) = cu(3,lp,m) + vz*dx
      cu(1,lp+1,m) = cu(1,lp+1,m) + vx*dy
      cu(2,lp+1,m) = cu(2,lp+1,m) + vy*dy
      cu(3,lp+1,m) = cu(3,lp+1,m) + vz*dy
c advance position half a time-step
      dx = part(1,nop,m) + vx*dt
      dy = part(2,nop,m) + vy*dt
      dz = part(3,nop,m) + vz*dt
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
      return
      end
c-----------------------------------------------------------------------
      subroutine PSJOST32XL(part,cu,npp,noff,nn,amxyz,qm,dt,nx,idimp,npm
     1ax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n24)
c for 3d code, this subroutine calculates particle current density
c using first-order linear interpolation, periodic boundaries,
c with short vectors over independent weights, and distributed data,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized distributed version with 1d addressing
c for 2D spatial decomposition
c 67 flops/particle, 78 loads, 75 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(n,m,l,i)=qci*(1.-dx)*(1.-dy)*(1.-dz)
c cu(n+1,m,l,i)=qci*dx*(1.-dy)*(1.-dz)
c cu(n,m+1,l,i)=qci*(1.-dx)*dy*(1.-dz)
c cu(n+1,m+1,l,i)=qci*dx*dy*(1.-dz)
c cu(n,m,l+1,i)=qci*(1.-dx)*(1.-dy)*dz
c cu(n+1,m,l+1,i)=qci*dx*(1.-dy)*dz
c cu(n,m+1,l+1,i)=qci*(1.-dx)*dy*dz
c cu(n+1,m+1,l+1,i)=qci*dx*dy*dz
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c and qci = qm*vi, where i = x,y,z
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = x velocity of particle n in partition m
c part(5,n,m) = y velocity of particle n in partition m
c part(6,n,m) = z velocity of particle n in partition m
c cu(n,i,m) = ith component of current density at grid point (j,kk,ll),
c where n = j + nxv*(k - 1) + nxv*nyv*(ll - 1)
c and kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nn = scratch address array for vectorized current deposition
c amxyz = scratch weight array for vectorized current deposition
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx = system length in x direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first virtual dimension of current array, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nxvyzp = nxv*nypmx*nzpmx, first dimension of charge array
c idds = dimensionality of domain decomposition
c npd = size of scratch buffers for vectorized push/current deposition
c n24 = number of independent weights
      dimension part(idimp,npmax,mnblok), cu(3*nxvyzp,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      dimension nn(n24,npd,mnblok), amxyz(n24,npd,mnblok)
      nxvy = nxv*nypmx
      nxvyzp2 = nxvyzp + nxvyzp
      zero = 0.
      anx = float(nx)
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
      in = n + mm + l
      nn(1,i,m) = in
      nn(1+8,i,m) = in + nxvyzp
      nn(1+16,i,m) = in + nxvyzp2
      in = np + mm + l
      nn(2,i,m) = in
      nn(2+8,i,m) = in + nxvyzp
      nn(2+16,i,m) = in + nxvyzp2
      in = n + mp + l
      nn(3,i,m) = in
      nn(3+8,i,m) = in + nxvyzp
      nn(3+16,i,m) = in + nxvyzp2
      in = np + mp + l
      nn(4,i,m) = in
      nn(4+8,i,m) = in + nxvyzp
      nn(4+16,i,m) = in + nxvyzp2
      in = n + mm + lp
      nn(5,i,m) = in
      nn(5+8,i,m) = in + nxvyzp
      nn(5+16,i,m) = in + nxvyzp2
      in = np + mm + lp
      nn(6,i,m) = in
      nn(6+8,i,m) = in + nxvyzp
      nn(6+16,i,m) = in + nxvyzp2
      in = n + mp + lp
      nn(7,i,m) = in
      nn(7+8,i,m) = in + nxvyzp
      nn(7+16,i,m) = in + nxvyzp2
      in = np + mp + lp
      nn(8,i,m) = in
      nn(8+8,i,m) = in + nxvyzp
      nn(8+16,i,m) = in + nxvyzp2
      dx = amx*amy*amz
      dy = dxp*amy*amz
      vx = part(4,i+jb,m)
      vy = part(5,i+jb,m)
      vz = part(6,i+jb,m)
      amxyz(1,i,m) = vx*dx
      amxyz(2,i,m) = vx*dy
      amxyz(1+8,i,m) = vy*dx
      amxyz(2+8,i,m) = vy*dy
      amxyz(1+16,i,m) = vz*dx
      amxyz(2+16,i,m) = vz*dy
      dx = amx*dyp*amz
      dy = dxp*dyp*amz
      amxyz(3,i,m) = vx*dx
      amxyz(4,i,m) = vx*dy
      amxyz(3+8,i,m) = vy*dx
      amxyz(4+8,i,m) = vy*dy
      amxyz(3+16,i,m) = vz*dx
      amxyz(4+16,i,m) = vz*dy
      dx = amx*amy*dzp
      dy = dxp*amy*dzp
      amxyz(5,i,m) = vx*dx
      amxyz(6,i,m) = vx*dy
      amxyz(5+8,i,m) = vy*dx
      amxyz(6+8,i,m) = vy*dy
      amxyz(5+16,i,m) = vz*dx
      amxyz(6+16,i,m) = vz*dy
      dx = amx*dyp*dzp
      dy = dxp*dyp*dzp
      amxyz(7,i,m) = vx*dx
      amxyz(8,i,m) = vx*dy
      amxyz(7+8,i,m) = vy*dx
      amxyz(8+8,i,m) = vy*dy
      amxyz(7+16,i,m) = vz*dx
      amxyz(8+16,i,m) = vz*dy   
c advance position half a time-step
      dx = part(1,i+jb,m) + vx*dt
      dy = part(2,i+jb,m) + vy*dt
      dz = part(3,i+jb,m) + vz*dt
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,i+jb,m) = dx
      part(2,i+jb,m) = dy
      part(3,i+jb,m) = dz
   10 continue
c deposit current
      do 30 i = 1, npb
cdir$ ivdep
      do 20 k = 1, 24
      cu(nn(k,i,m),m) = cu(nn(k,i,m),m) + amxyz(k,i,m)
   20 continue
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSJOST32XL(part,cu,npp,noff,nn,amxyz,qm,dt,nx,ny,nz,id
     1imp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n24,ipbc)
c for 3d code, this subroutine calculates particle current density
c using first-order linear interpolation, for distributed data
c with 2D spatial decomposition
c with short vectors over independent weights, 
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized distributed version with guard cells and 1d addressing
c 65 flops/particle, 78 loads, 75 stores
c input: all, output: part, cu
c current density is approximated by values at the nearest grid points
c cu(i,n,m,l)=qci*(1.-dx)*(1.-dy)*(1.-dz)
c cu(i,n+1,m,l)=qci*dx*(1.-dy)*(1.-dz)
c cu(i,n,m+1,l)=qci*(1.-dx)*dy*(1.-dz)
c cu(i,n+1,m+1,l)=qci*dx*dy*(1.-dz)
c cu(i,n,m,l+1)=qci*(1.-dx)*(1.-dy)*dz
c cu(i,n+1,m,l+1)=qci*dx*(1.-dy)*dz
c cu(i,n,m+1,l+1)=qci*(1.-dx)*dy*dz
c cu(i,n+1,m+1,l+1)=qci*dx*dy*dz
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c and qci = qm*vi, where i = x,y,z
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = x velocity of particle n in partition m
c part(5,n,m) = y velocity of particle n in partition m
c part(6,n,m) = z velocity of particle n in partition m
c cu(i,n,m) = ith component of current density at grid point (j,kk,ll),
c where n = j + nxv*(k - 1) + nxv*nyv*(ll - 1)
c and kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c nn = scratch address array for vectorized current deposition
c amxyz = scratch weight array for vectorized current deposition
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first virtual dimension of current array, must be >= nx+1
c nypmx = maximum size of particle partition in y, including guard cells
c nxvyzp = nxv*nypmx*nzpmx, first dimension of charge array
c idds = dimensionality of domain decomposition
c npd = size of scratch buffers for vectorized push/current deposition
c n24 = number of independent weights
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      dimension part(idimp,npmax,mnblok), cu(3*nxvyzp,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      dimension nn(n24,npd,mnblok), amxyz(n24,npd,mnblok)
      nxv3 = 3*nxv
      nxvy3 = nxv3*nypmx
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
      n3 = 3*n + 1
      mm = nxv3*(mm - mnoff) + nxvy3*(ll - lnoff)
      amx = qm - dxp
      amy = 1. - dyp
      mm = mm + n3
      dx1 = dxp*dyp
      dyp = amx*dyp
      mp = mm + nxv3
      amx = amx*amy
      amz = 1. - dzp
      ll = mm + nxvy3
      amy = dxp*amy
      lp = mp + nxvy3
      nn(1,i,m) = mm
      nn(2,i,m) = mm + 1
      nn(3,i,m) = mm + 2
      nn(4,i,m) = mm + 3
      nn(5,i,m) = mm + 4
      nn(6,i,m) = mm + 5
      nn(7,i,m) = mp
      nn(8,i,m) = mp + 1
      nn(9,i,m) = mp + 2
      nn(10,i,m) = mp + 3
      nn(11,i,m) = mp + 4
      nn(12,i,m) = mp + 5
      nn(13,i,m) = ll
      nn(14,i,m) = ll + 1
      nn(15,i,m) = ll + 2
      nn(16,i,m) = ll + 3
      nn(17,i,m) = ll + 4
      nn(18,i,m) = ll + 5
      nn(19,i,m) = lp
      nn(20,i,m) = lp + 1
      nn(21,i,m) = lp + 2
      nn(22,i,m) = lp + 3
      nn(23,i,m) = lp + 4
      nn(24,i,m) = lp + 5
      dx = amx*amz
      dy = amy*amz
      vx = part(4,i+jb,m)
      vy = part(5,i+jb,m)
      vz = part(6,i+jb,m)
      amxyz(1,i,m) = vx*dx
      amxyz(2,i,m) = vy*dx
      amxyz(3,i,m) = vz*dx
      dx = dyp*amz
      amxyz(4,i,m) = vx*dy
      amxyz(5,i,m) = vy*dy
      amxyz(6,i,m) = vz*dy
      dy = dx1*amz
      amxyz(7,i,m) = vx*dx
      amxyz(8,i,m) = vy*dx
      amxyz(9,i,m) = vz*dx
      dx = amx*dzp
      amxyz(10,i,m) = vx*dy
      amxyz(11,i,m) = vy*dy
      amxyz(12,i,m) = vz*dy
      dy = amy*dzp
      amxyz(13,i,m) = vx*dx
      amxyz(14,i,m) = vy*dx
      amxyz(15,i,m) = vz*dx
      dx = dyp*dzp
      amxyz(16,i,m) = vx*dy
      amxyz(17,i,m) = vy*dy
      amxyz(18,i,m) = vz*dy
      dy = dx1*dzp
      amxyz(19,i,m) = vx*dx
      amxyz(20,i,m) = vy*dx
      amxyz(21,i,m) = vz*dx
      amxyz(22,i,m) = vx*dy
      amxyz(23,i,m) = vy*dy
      amxyz(24,i,m) = vz*dy
c advance position half a time-step
      dx = part(1,i+jb,m) + vx*dt
      dy = part(2,i+jb,m) + vy*dt
      dz = part(3,i+jb,m) + vz*dt
c periodic boundary conditions
      if (ipbc.eq.1) then
         if (dx.lt.edgelx) dx = dx + edgerx
         if (dx.ge.edgerx) dx = dx - edgerx
c reflecting boundary conditions
      else if (ipbc.eq.2) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i+jb,m)
            part(4,i+jb,m) = -part(4,i+jb,m)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,i+jb,m)
            part(5,i+jb,m) = -part(5,i+jb,m)
         endif
         if ((dz.lt.edgelz).or.(dz.ge.edgerz)) then
            dz = part(3,i+jb,m)
            part(6,i+jb,m) = -part(6,i+jb,m)
         endif
c mixed reflecting/periodic boundary conditions
      else if (ipbc.eq.3) then
         if ((dx.lt.edgelx).or.(dx.ge.edgerx)) then
            dx = part(1,i+jb,m)
            part(4,i+jb,m) = -part(4,i+jb,m)
         endif
         if ((dy.lt.edgely).or.(dy.ge.edgery)) then
            dy = part(2,i+jb,m)
            part(5,i+jb,m) = -part(5,i+jb,m)
         endif
      endif
c set new position
      part(1,i+jb,m) = dx
      part(2,i+jb,m) = dy
      part(3,i+jb,m) = dz
   10 continue
c deposit current
      do 30 i = 1, npb
cdir$ ivdep
      do 20 k = 1, 24
      cu(nn(k,i,m),m) = cu(nn(k,i,m),m) + amxyz(k,i,m)
   20 continue
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPUSH32(part,fx,fy,fz,bx,by,bz,npp,noff,qbm,dt,ek,nx,i
     1dimp,npmax,mnblok,nxv,nypmx,nzpmx,idds)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field and periodic boundary
c conditions. Using the Boris Mover,
c for distributed data with 2D spatial decomposition
c baseline scalar distributed version
c 494 flops/particle, 1 divide, 168 loads, 6 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
c omz = (q/m)*bz(x(t),y(t),z(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c z(t+dt)=z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
c bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
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
c similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
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
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
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
      dimension fz(nxv,nypmx,nzpmx,mnblok), bx(nxv,nypmx,nzpmx,mnblok)
      dimension by(nxv,nypmx,nzpmx,mnblok), bz(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      zero = 0.
      anx = float(nx)
      qtmh = .5*qbm*dt
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
c find electric field
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
c find magnetic field
      ox = dzl*(dyl*(dxl*bx(nl,ml,lm,m) + amx*bx(nn,ml,lm,m) + dxp*bx(np
     1,ml,lm,m)) + amy*(dxl*bx(nl,mm,lm,m) + amx*bx(nn,mm,lm,m) + dxp*bx
     2(np,mm,lm,m)) + dyp*(dxl*bx(nl,mp,lm,m) + amx*bx(nn,mp,lm,m) + dxp
     3*bx(np,mp,lm,m))) + amz*(dyl*(dxl*bx(nl,ml,ll,m) + amx*bx(nn,ml,ll
     4,m) + dxp*bx(np,ml,ll,m)) + amy*(dxl*bx(nl,mm,ll,m) + amx*bx(nn,mm
     5,ll,m) + dxp*bx(np,mm,ll,m)) + dyp*(dxl*bx(nl,mp,ll,m) + amx*bx(nn
     6,mp,ll,m) + dxp*bx(np,mp,ll,m))) + dzp*(dyl*(dxl*bx(nl,ml,lp,m) + 
     7amx*bx(nn,ml,lp,m) + dxp*bx(np,ml,lp,m)) + amy*(dxl*bx(nl,mm,lp,m) 
     8 + amx*bx(nn,mm,lp,m) + dxp*bx(np,mm,lp,m)) + dyp*(dxl*bx(nl,mp,lp
     9,m) + amx*bx(nn,mp,lp,m) + dxp*bx(np,mp,lp,m)))
      oy = dzl*(dyl*(dxl*by(nl,ml,lm,m) + amx*by(nn,ml,lm,m) + dxp*by(np
     1,ml,lm,m)) + amy*(dxl*by(nl,mm,lm,m) + amx*by(nn,mm,lm,m) + dxp*by
     2(np,mm,lm,m)) + dyp*(dxl*by(nl,mp,lm,m) + amx*by(nn,mp,lm,m) + dxp
     3*by(np,mp,lm,m))) + amz*(dyl*(dxl*by(nl,ml,ll,m) + amx*by(nn,ml,ll
     4,m) + dxp*by(np,ml,ll,m)) + amy*(dxl*by(nl,mm,ll,m) + amx*by(nn,mm
     5,ll,m) + dxp*by(np,mm,ll,m)) + dyp*(dxl*by(nl,mp,ll,m) + amx*by(nn
     6,mp,ll,m) + dxp*by(np,mp,ll,m))) + dzp*(dyl*(dxl*by(nl,ml,lp,m) + 
     7amx*by(nn,ml,lp,m) + dxp*by(np,ml,lp,m)) + amy*(dxl*by(nl,mm,lp,m) 
     8 + amx*by(nn,mm,lp,m) + dxp*by(np,mm,lp,m)) + dyp*(dxl*by(nl,mp,lp
     9,m) + amx*by(nn,mp,lp,m) + dxp*by(np,mp,lp,m)))
      oz = dzl*(dyl*(dxl*bz(nl,ml,lm,m) + amx*bz(nn,ml,lm,m) + dxp*bz(np
     1,ml,lm,m)) + amy*(dxl*bz(nl,mm,lm,m) + amx*bz(nn,mm,lm,m) + dxp*bz
     2(np,mm,lm,m)) + dyp*(dxl*bz(nl,mp,lm,m) + amx*bz(nn,mp,lm,m) + dxp
     3*bz(np,mp,lm,m))) + amz*(dyl*(dxl*bz(nl,ml,ll,m) + amx*bz(nn,ml,ll
     4,m) + dxp*bz(np,ml,ll,m)) + amy*(dxl*bz(nl,mm,ll,m) + amx*bz(nn,mm
     5,ll,m) + dxp*bz(np,mm,ll,m)) + dyp*(dxl*bz(nl,mp,ll,m) + amx*bz(nn
     6,mp,ll,m) + dxp*bz(np,mp,ll,m))) + dzp*(dyl*(dxl*bz(nl,ml,lp,m) + 
     7amx*bz(nn,ml,lp,m) + dxp*bz(np,ml,lp,m)) + amy*(dxl*bz(nl,mm,lp,m) 
     8 + amx*bz(nn,mm,lp,m) + dxp*bz(np,mm,lp,m)) + dyp*(dxl*bz(nl,mp,lp
     9,m) + amx*bz(nn,mp,lp,m) + dxp*bz(np,mp,lp,m)))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j,m) + dx
      acy = part(5,j,m) + dy
      acz = part(6,j,m) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
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
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGBPUSH32(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ek,nx,ny,n
     1z,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field, for distributed data
c with 2D spatial decomposition
c Using the Boris Mover.
c scalar version using guard cells, 
c 451 flops/particle, 1 divide, 168 loads, 6 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
c omz = (q/m)*bz(x(t),y(t),z(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c z(t+dt)=z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
c bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
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
c similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c fxyz(1,j+1,k+1,l,m) = x component of force/charge at grid (j,kk,ll)
c fxyz(2,j+1,k+1,l,m) = y component of force/charge at grid (j,kk,ll)
c fxyz(3,j+1,k+1,l,m) = z component of force/charge at grid (j,kk,ll)
c in other words, fxyz are the convolutions of the electric field
c over the particle shape,
c bxyz(1,j+1,k+1,l,m) = x component of magnetic field at grid (j,kk,ll)
c bxyz(2,j+1,k+1,l,m) = y component of magnetic field at grid (j,kk,ll)
c bxyz(3,j+1,k+1,l,m) = z component of magnetic field at grid (j,kk,ll)
c in other words, bxyz is the convolution of the magnetic field
c over the particle shape,
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
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
      dimension bxyz(3,nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      qtmh = .5*qbm*dt
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
c find electric field
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
c find magnetic field
      ox = dzl*(dx1*bxyz(1,nl,ml,lm,m) + dx2*bxyz(1,nn,ml,lm,m) + dyl*bx
     1yz(1,np,ml,lm,m) + dy1*bxyz(1,nl,mm,lm,m) + dy2*bxyz(1,nn,mm,lm,m)
     2 + amy*bxyz(1,np,mm,lm,m) + dxl*bxyz(1,nl,mp,lm,m) + amx*bxyz(1,nn
     3,mp,lm,m) + dyp*bxyz(1,np,mp,lm,m))
      oy = dzl*(dx1*bxyz(2,nl,ml,lm,m) + dx2*bxyz(2,nn,ml,lm,m) + dyl*bx
     1yz(2,np,ml,lm,m) + dy1*bxyz(2,nl,mm,lm,m) + dy2*bxyz(2,nn,mm,lm,m)
     2 + amy*bxyz(2,np,mm,lm,m) + dxl*bxyz(2,nl,mp,lm,m) + amx*bxyz(2,nn
     3,mp,lm,m) + dyp*bxyz(2,np,mp,lm,m))
      oz = dzl*(dx1*bxyz(3,nl,ml,lm,m) + dx2*bxyz(3,nn,ml,lm,m) + dyl*bx
     1yz(3,np,ml,lm,m) + dy1*bxyz(3,nl,mm,lm,m) + dy2*bxyz(3,nn,mm,lm,m)
     2 + amy*bxyz(3,np,mm,lm,m) + dxl*bxyz(3,nl,mp,lm,m) + amx*bxyz(3,nn
     3,mp,lm,m) + dyp*bxyz(3,np,mp,lm,m))
      ox = ox + amz*(dx1*bxyz(1,nl,ml,ll,m) + dx2*bxyz(1,nn,ml,ll,m) + d
     1yl*bxyz(1,np,ml,ll,m) + dy1*bxyz(1,nl,mm,ll,m) + dy2*bxyz(1,nn,mm,
     2ll,m) + amy*bxyz(1,np,mm,ll,m) + dxl*bxyz(1,nl,mp,ll,m) + amx*bxyz
     3(1,nn,mp,ll,m) + dyp*bxyz(1,np,mp,ll,m))
      oy = oy + amz*(dx1*bxyz(2,nl,ml,ll,m) + dx2*bxyz(2,nn,ml,ll,m) + d
     1yl*bxyz(2,np,ml,ll,m) + dy1*bxyz(2,nl,mm,ll,m) + dy2*bxyz(2,nn,mm,
     2ll,m) + amy*bxyz(2,np,mm,ll,m) + dxl*bxyz(2,nl,mp,ll,m) + amx*bxyz
     3(2,nn,mp,ll,m) + dyp*bxyz(2,np,mp,ll,m))
      oz = oz + amz*(dx1*bxyz(3,nl,ml,ll,m) + dx2*bxyz(3,nn,ml,ll,m) + d
     1yl*bxyz(3,np,ml,ll,m) + dy1*bxyz(3,nl,mm,ll,m) + dy2*bxyz(3,nn,mm,
     2ll,m) + amy*bxyz(3,np,mm,ll,m) + dxl*bxyz(3,nl,mp,ll,m) + amx*bxyz
     3(3,nn,mp,ll,m) + dyp*bxyz(3,np,mp,ll,m))
      ox = ox + dzp*(dx1*bxyz(1,nl,ml,lp,m) + dx2*bxyz(1,nn,ml,lp,m) + d
     1yl*bxyz(1,np,ml,lp,m) + dy1*bxyz(1,nl,mm,lp,m) + dy2*bxyz(1,nn,mm,
     2lp,m) + amy*bxyz(1,np,mm,lp,m) + dxl*bxyz(1,nl,mp,lp,m) + amx*bxyz
     3(1,nn,mp,lp,m) + dyp*bxyz(1,np,mp,lp,m))
      oy = oy + dzp*(dx1*bxyz(2,nl,ml,lp,m) + dx2*bxyz(2,nn,ml,lp,m) + d
     1yl*bxyz(2,np,ml,lp,m) + dy1*bxyz(2,nl,mm,lp,m) + dy2*bxyz(2,nn,mm,
     2lp,m) + amy*bxyz(2,np,mm,lp,m) + dxl*bxyz(2,nl,mp,lp,m) + amx*bxyz
     3(2,nn,mp,lp,m) + dyp*bxyz(2,np,mp,lp,m))
      oz = oz + dzp*(dx1*bxyz(3,nl,ml,lp,m) + dx2*bxyz(3,nn,ml,lp,m) + d
     1yl*bxyz(3,np,ml,lp,m) + dy1*bxyz(3,nl,mm,lp,m) + dy2*bxyz(3,nn,mm,
     2lp,m) + amy*bxyz(3,np,mm,lp,m) + dxl*bxyz(3,nl,mp,lp,m) + amx*bxyz
     3(3,nn,mp,lp,m) + dyp*bxyz(3,np,mp,lp,m))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j,m) + dx
      acy = part(5,j,m) + dy
      acz = part(6,j,m) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(4,j,m) = dx
      part(5,j,m) = dy
      part(6,j,m) = dz
c new position
      dx = part(1,j,m) + dx*dtc
      dy = part(2,j,m) + dy*dtc
      dz = part(3,j,m) + dz*dtc
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
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSBPUSH32(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ek,nx,ny,
     1nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field, for distributed data
c with 2D spatial decomposition
c Using the Boris Mover,
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c 451 flops/particle, 1 divide, 168 loads, 6 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
c omz = (q/m)*bz(x(t),y(t),z(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c z(t+dt)=z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
c bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
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
c similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c fxyz(1,j+1,k+1,l,m) = x component of force/charge at grid (j,kk,ll)
c fxyz(2,j+1,k+1,l,m) = y component of force/charge at grid (j,kk,ll)
c fxyz(3,j+1,k+1,l,m) = z component of force/charge at grid (j,kk,ll)
c in other words, fxyz are the convolutions of the electric field
c over the particle shape,
c bxyz(1,j+1,k+1,l,m) = x component of magnetic field at grid (j,kk,ll)
c bxyz(2,j+1,k+1,l,m) = y component of magnetic field at grid (j,kk,ll)
c bxyz(3,j+1,k+1,l,m) = z component of magnetic field at grid (j,kk,ll)
c over the particle shape,
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)s
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
      dimension fxyz(3,nxyzp,mnblok), bxyz(3,nxyzp,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
      nxyv = nxv*nypmx
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
      mm = nxv*mmm + nxyv*lll
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
c find electric field
      dx = dzl*(dx1*fxyz(1,ml,m) + dx2*fxyz(1,ml+1,m) + dyl*fxyz(1,ml+2,
     1m) + dy1*fxyz(1,mn,m) + dy2*fxyz(1,mn+1,m) + amy*fxyz(1,mn+2,m) + 
     2dxl*fxyz(1,mp,m) + amx*fxyz(1,mp+1,m) + dyp*fxyz(1,mp+2,m))
      dy = dzl*(dx1*fxyz(2,ml,m) + dx2*fxyz(2,ml+1,m) + dyl*fxyz(2,ml+2,
     1m) + dy1*fxyz(2,mn,m) + dy2*fxyz(2,mn+1,m) + amy*fxyz(2,mn+2,m) + 
     2dxl*fxyz(2,mp,m) + amx*fxyz(2,mp+1,m) + dyp*fxyz(2,mp+2,m))
      dz = dzl*(dx1*fxyz(3,ml,m) + dx2*fxyz(3,ml+1,m) + dyl*fxyz(3,ml+2,
     1m) + dy1*fxyz(3,mn,m) + dy2*fxyz(3,mn+1,m) + amy*fxyz(3,mn+2,m) + 
     2dxl*fxyz(3,mp,m) + amx*fxyz(3,mp+1,m) + dyp*fxyz(3,mp+2,m))
      ml = ml + nxyv
      mn = mn + nxyv
      mp = mp + nxyv
      dx = dx + amz*(dx1*fxyz(1,ml,m) + dx2*fxyz(1,ml+1,m) + dyl*fxyz(1,
     1ml+2,m) + dy1*fxyz(1,mn,m) + dy2*fxyz(1,mn+1,m) + amy*fxyz(1,mn+2,
     2m) + dxl*fxyz(1,mp,m) + amx*fxyz(1,mp+1,m) + dyp*fxyz(1,mp+2,m))
      dy = dy + amz*(dx1*fxyz(2,ml,m) + dx2*fxyz(2,ml+1,m) + dyl*fxyz(2,
     1ml+2,m) + dy1*fxyz(2,mn,m) + dy2*fxyz(2,mn+1,m) + amy*fxyz(2,mn+2,
     2m) + dxl*fxyz(2,mp,m) + amx*fxyz(2,mp+1,m) + dyp*fxyz(2,mp+2,m))
      dz = dz + amz*(dx1*fxyz(3,ml,m) + dx2*fxyz(3,ml+1,m) + dyl*fxyz(3,
     1ml+2,m) + dy1*fxyz(3,mn,m) + dy2*fxyz(3,mn+1,m) + amy*fxyz(3,mn+2,
     2m) + dxl*fxyz(3,mp,m) + amx*fxyz(3,mp+1,m) + dyp*fxyz(3,mp+2,m))
      ml = ml + nxyv
      mn = mn + nxyv
      mp = mp + nxyv
      dx = dx + dzp*(dx1*fxyz(1,ml,m) + dx2*fxyz(1,ml+1,m) + dyl*fxyz(1,
     1ml+2,m) + dy1*fxyz(1,mn,m) + dy2*fxyz(1,mn+1,m) + amy*fxyz(1,mn+2,
     2m) + dxl*fxyz(1,mp,m) + amx*fxyz(1,mp+1,m) + dyp*fxyz(1,mp+2,m))
      dy = dy + dzp*(dx1*fxyz(2,ml,m) + dx2*fxyz(2,ml+1,m) + dyl*fxyz(2,
     1ml+2,m) + dy1*fxyz(2,mn,m) + dy2*fxyz(2,mn+1,m) + amy*fxyz(2,mn+2,
     2m) + dxl*fxyz(2,mp,m) + amx*fxyz(2,mp+1,m) + dyp*fxyz(2,mp+2,m))
      dz = dz + dzp*(dx1*fxyz(3,ml,m) + dx2*fxyz(3,ml+1,m) + dyl*fxyz(3,
     1ml+2,m) + dy1*fxyz(3,mn,m) + dy2*fxyz(3,mn+1,m) + amy*fxyz(3,mn+2,
     2m) + dxl*fxyz(3,mp,m) + amx*fxyz(3,mp+1,m) + dyp*fxyz(3,mp+2,m))
c find magnetic field
      ml = mm + nn
      mn = ml + nxv
      mp = mn + nxv
      ox = dzl*(dx1*bxyz(1,ml,m) + dx2*bxyz(1,ml+1,m) + dyl*bxyz(1,ml+2,
     1m) + dy1*bxyz(1,mn,m) + dy2*bxyz(1,mn+1,m) + amy*bxyz(1,mn+2,m) + 
     2dxl*bxyz(1,mp,m) + amx*bxyz(1,mp+1,m) + dyp*bxyz(1,mp+2,m))
      oy = dzl*(dx1*bxyz(2,ml,m) + dx2*bxyz(2,ml+1,m) + dyl*bxyz(2,ml+2,
     1m) + dy1*bxyz(2,mn,m) + dy2*bxyz(2,mn+1,m) + amy*bxyz(2,mn+2,m) + 
     2dxl*bxyz(2,mp,m) + amx*bxyz(2,mp+1,m) + dyp*bxyz(2,mp+2,m))
      oz = dzl*(dx1*bxyz(3,ml,m) + dx2*bxyz(3,ml+1,m) + dyl*bxyz(3,ml+2,
     1m) + dy1*bxyz(3,mn,m) + dy2*bxyz(3,mn+1,m) + amy*bxyz(3,mn+2,m) + 
     2dxl*bxyz(3,mp,m) + amx*bxyz(3,mp+1,m) + dyp*bxyz(3,mp+2,m))
      ml = ml + nxyv
      mn = mn + nxyv
      mp = mp + nxyv
      ox = ox + amz*(dx1*bxyz(1,ml,m) + dx2*bxyz(1,ml+1,m) + dyl*bxyz(1,
     1ml+2,m) + dy1*bxyz(1,mn,m) + dy2*bxyz(1,mn+1,m) + amy*bxyz(1,mn+2,
     2m) + dxl*bxyz(1,mp,m) + amx*bxyz(1,mp+1,m) + dyp*bxyz(1,mp+2,m))
      oy = oy + amz*(dx1*bxyz(2,ml,m) + dx2*bxyz(2,ml+1,m) + dyl*bxyz(2,
     1ml+2,m) + dy1*bxyz(2,mn,m) + dy2*bxyz(2,mn+1,m) + amy*bxyz(2,mn+2,
     2m) + dxl*bxyz(2,mp,m) + amx*bxyz(2,mp+1,m) + dyp*bxyz(2,mp+2,m))
      oz = oz + amz*(dx1*bxyz(3,ml,m) + dx2*bxyz(3,ml+1,m) + dyl*bxyz(3,
     1ml+2,m) + dy1*bxyz(3,mn,m) + dy2*bxyz(3,mn+1,m) + amy*bxyz(3,mn+2,
     2m) + dxl*bxyz(3,mp,m) + amx*bxyz(3,mp+1,m) + dyp*bxyz(3,mp+2,m))
      ml = ml + nxyv
      mn = mn + nxyv
      mp = mp + nxyv
      ox = ox + dzp*(dx1*bxyz(1,ml,m) + dx2*bxyz(1,ml+1,m) + dyl*bxyz(1,
     1ml+2,m) + dy1*bxyz(1,mn,m) + dy2*bxyz(1,mn+1,m) + amy*bxyz(1,mn+2,
     2m) + dxl*bxyz(1,mp,m) + amx*bxyz(1,mp+1,m) + dyp*bxyz(1,mp+2,m))
      oy = oy + dzp*(dx1*bxyz(2,ml,m) + dx2*bxyz(2,ml+1,m) + dyl*bxyz(2,
     1ml+2,m) + dy1*bxyz(2,mn,m) + dy2*bxyz(2,mn+1,m) + amy*bxyz(2,mn+2,
     2m) + dxl*bxyz(2,mp,m) + amx*bxyz(2,mp+1,m) + dyp*bxyz(2,mp+2,m))
      oz = oz + dzp*(dx1*bxyz(3,ml,m) + dx2*bxyz(3,ml+1,m) + dyl*bxyz(3,
     1ml+2,m) + dy1*bxyz(3,mn,m) + dy2*bxyz(3,mn+1,m) + amy*bxyz(3,mn+2,
     2m) + dxl*bxyz(3,mp,m) + amx*bxyz(3,mp+1,m) + dyp*bxyz(3,mp+2,m))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j,m) + dx
      acy = part(5,j,m) + dy
      acz = part(6,j,m) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(4,j,m) = dx
      part(5,j,m) = dy
      part(6,j,m) = dz
c new position
      dx = part(1,j,m) + dx*dtc
      dy = part(2,j,m) + dy*dtc
      dz = part(3,j,m) + dz*dtc
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
      mm = nxv*mmm + nxyv*lll
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
c find electric field
      dx = dzl*(dx1*fxyz(1,ml,m) + dx2*fxyz(1,ml+1,m) + dyl*fxyz(1,ml+2,
     1m) + dy1*fxyz(1,mn,m) + dy2*fxyz(1,mn+1,m) + amy*fxyz(1,mn+2,m) + 
     2dxl*fxyz(1,mp,m) + amx*fxyz(1,mp+1,m) + dyp*fxyz(1,mp+2,m))
      dy = dzl*(dx1*fxyz(2,ml,m) + dx2*fxyz(2,ml+1,m) + dyl*fxyz(2,ml+2,
     1m) + dy1*fxyz(2,mn,m) + dy2*fxyz(2,mn+1,m) + amy*fxyz(2,mn+2,m) + 
     2dxl*fxyz(2,mp,m) + amx*fxyz(2,mp+1,m) + dyp*fxyz(2,mp+2,m))
      dz = dzl*(dx1*fxyz(3,ml,m) + dx2*fxyz(3,ml+1,m) + dyl*fxyz(3,ml+2,
     1m) + dy1*fxyz(3,mn,m) + dy2*fxyz(3,mn+1,m) + amy*fxyz(3,mn+2,m) + 
     2dxl*fxyz(3,mp,m) + amx*fxyz(3,mp+1,m) + dyp*fxyz(3,mp+2,m))
      ml = ml + nxyv
      mn = mn + nxyv
      mp = mp + nxyv
      dx = dx + amz*(dx1*fxyz(1,ml,m) + dx2*fxyz(1,ml+1,m) + dyl*fxyz(1,
     1ml+2,m) + dy1*fxyz(1,mn,m) + dy2*fxyz(1,mn+1,m) + amy*fxyz(1,mn+2,
     2m) + dxl*fxyz(1,mp,m) + amx*fxyz(1,mp+1,m) + dyp*fxyz(1,mp+2,m))
      dy = dy + amz*(dx1*fxyz(2,ml,m) + dx2*fxyz(2,ml+1,m) + dyl*fxyz(2,
     1ml+2,m) + dy1*fxyz(2,mn,m) + dy2*fxyz(2,mn+1,m) + amy*fxyz(2,mn+2,
     2m) + dxl*fxyz(2,mp,m) + amx*fxyz(2,mp+1,m) + dyp*fxyz(2,mp+2,m))
      dz = dz + amz*(dx1*fxyz(3,ml,m) + dx2*fxyz(3,ml+1,m) + dyl*fxyz(3,
     1ml+2,m) + dy1*fxyz(3,mn,m) + dy2*fxyz(3,mn+1,m) + amy*fxyz(3,mn+2,
     2m) + dxl*fxyz(3,mp,m) + amx*fxyz(3,mp+1,m) + dyp*fxyz(3,mp+2,m))
      ml = ml + nxyv
      mn = mn + nxyv
      mp = mp + nxyv
      dx = dx + dzp*(dx1*fxyz(1,ml,m) + dx2*fxyz(1,ml+1,m) + dyl*fxyz(1,
     1ml+2,m) + dy1*fxyz(1,mn,m) + dy2*fxyz(1,mn+1,m) + amy*fxyz(1,mn+2,
     2m) + dxl*fxyz(1,mp,m) + amx*fxyz(1,mp+1,m) + dyp*fxyz(1,mp+2,m))
      dy = dy + dzp*(dx1*fxyz(2,ml,m) + dx2*fxyz(2,ml+1,m) + dyl*fxyz(2,
     1ml+2,m) + dy1*fxyz(2,mn,m) + dy2*fxyz(2,mn+1,m) + amy*fxyz(2,mn+2,
     2m) + dxl*fxyz(2,mp,m) + amx*fxyz(2,mp+1,m) + dyp*fxyz(2,mp+2,m))
      dz = dz + dzp*(dx1*fxyz(3,ml,m) + dx2*fxyz(3,ml+1,m) + dyl*fxyz(3,
     1ml+2,m) + dy1*fxyz(3,mn,m) + dy2*fxyz(3,mn+1,m) + amy*fxyz(3,mn+2,
     2m) + dxl*fxyz(3,mp,m) + amx*fxyz(3,mp+1,m) + dyp*fxyz(3,mp+2,m))
c find magnetic field
      ml = mm + nn
      mn = ml + nxv
      mp = mn + nxv
      ox = dzl*(dx1*bxyz(1,ml,m) + dx2*bxyz(1,ml+1,m) + dyl*bxyz(1,ml+2,
     1m) + dy1*bxyz(1,mn,m) + dy2*bxyz(1,mn+1,m) + amy*bxyz(1,mn+2,m) + 
     2dxl*bxyz(1,mp,m) + amx*bxyz(1,mp+1,m) + dyp*bxyz(1,mp+2,m))
      oy = dzl*(dx1*bxyz(2,ml,m) + dx2*bxyz(2,ml+1,m) + dyl*bxyz(2,ml+2,
     1m) + dy1*bxyz(2,mn,m) + dy2*bxyz(2,mn+1,m) + amy*bxyz(2,mn+2,m) + 
     2dxl*bxyz(2,mp,m) + amx*bxyz(2,mp+1,m) + dyp*bxyz(2,mp+2,m))
      oz = dzl*(dx1*bxyz(3,ml,m) + dx2*bxyz(3,ml+1,m) + dyl*bxyz(3,ml+2,
     1m) + dy1*bxyz(3,mn,m) + dy2*bxyz(3,mn+1,m) + amy*bxyz(3,mn+2,m) + 
     2dxl*bxyz(3,mp,m) + amx*bxyz(3,mp+1,m) + dyp*bxyz(3,mp+2,m))
      ml = ml + nxyv
      mn = mn + nxyv
      mp = mp + nxyv
      ox = ox + amz*(dx1*bxyz(1,ml,m) + dx2*bxyz(1,ml+1,m) + dyl*bxyz(1,
     1ml+2,m) + dy1*bxyz(1,mn,m) + dy2*bxyz(1,mn+1,m) + amy*bxyz(1,mn+2,
     2m) + dxl*bxyz(1,mp,m) + amx*bxyz(1,mp+1,m) + dyp*bxyz(1,mp+2,m))
      oy = oy + amz*(dx1*bxyz(2,ml,m) + dx2*bxyz(2,ml+1,m) + dyl*bxyz(2,
     1ml+2,m) + dy1*bxyz(2,mn,m) + dy2*bxyz(2,mn+1,m) + amy*bxyz(2,mn+2,
     2m) + dxl*bxyz(2,mp,m) + amx*bxyz(2,mp+1,m) + dyp*bxyz(2,mp+2,m))
      oz = oz + amz*(dx1*bxyz(3,ml,m) + dx2*bxyz(3,ml+1,m) + dyl*bxyz(3,
     1ml+2,m) + dy1*bxyz(3,mn,m) + dy2*bxyz(3,mn+1,m) + amy*bxyz(3,mn+2,
     2m) + dxl*bxyz(3,mp,m) + amx*bxyz(3,mp+1,m) + dyp*bxyz(3,mp+2,m))
      ml = ml + nxyv
      mn = mn + nxyv
      mp = mp + nxyv
      ox = ox + dzp*(dx1*bxyz(1,ml,m) + dx2*bxyz(1,ml+1,m) + dyl*bxyz(1,
     1ml+2,m) + dy1*bxyz(1,mn,m) + dy2*bxyz(1,mn+1,m) + amy*bxyz(1,mn+2,
     2m) + dxl*bxyz(1,mp,m) + amx*bxyz(1,mp+1,m) + dyp*bxyz(1,mp+2,m))
      oy = oy + dzp*(dx1*bxyz(2,ml,m) + dx2*bxyz(2,ml+1,m) + dyl*bxyz(2,
     1ml+2,m) + dy1*bxyz(2,mn,m) + dy2*bxyz(2,mn+1,m) + amy*bxyz(2,mn+2,
     2m) + dxl*bxyz(2,mp,m) + amx*bxyz(2,mp+1,m) + dyp*bxyz(2,mp+2,m))
      oz = oz + dzp*(dx1*bxyz(3,ml,m) + dx2*bxyz(3,ml+1,m) + dyl*bxyz(3,
     1ml+2,m) + dy1*bxyz(3,mn,m) + dy2*bxyz(3,mn+1,m) + amy*bxyz(3,mn+2,
     2m) + dxl*bxyz(3,mp,m) + amx*bxyz(3,mp+1,m) + dyp*bxyz(3,mp+2,m))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,nop,m) + dx
      acy = part(5,nop,m) + dy
      acz = part(6,nop,m) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(4,nop,m) = dx
      part(5,nop,m) = dy
      part(6,nop,m) = dz
c new position
      dx = part(1,nop,m) + dx*dtc
      dy = part(2,nop,m) + dy*dtc
      dz = part(3,nop,m) + dz*dtc
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
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPUSH32L(part,fx,fy,fz,bx,by,bz,npp,noff,qbm,dt,ek,nx,
     1idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field and periodic boundary
c conditions. Using the Boris Mover,
c for distributed data with 2D spatial decomposition
c baseline scalar distributed version
c 206 flops/particle, 1 divide, 54 loads, 6 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
c omz = (q/m)*bz(x(t),y(t),z(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c z(t+dt)=z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
c bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
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
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
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
      dimension fz(nxv,nypmx,nzpmx,mnblok), bx(nxv,nypmx,nzpmx,mnblok)
      dimension by(nxv,nypmx,nzpmx,mnblok), bz(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      zero = 0.
      anx = float(nx)
      qtmh = .5*qbm*dt
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
c find electric field
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
c find magnetic field
      ox = amz*(amy*(amx*bx(nn,mm,ll,m) + dxp*bx(np,mm,ll,m)) + dyp*(amx
     1*bx(nn,mp,ll,m) + dxp*bx(np,mp,ll,m))) + dzp*(amy*(amx*bx(nn,mm,lp
     2,m) + dxp*bx(np,mm,lp,m)) + dyp*(amx*bx(nn,mp,lp,m) + dxp*bx(np,mp
     3,lp,m)))
      oy = amz*(amy*(amx*by(nn,mm,ll,m) + dxp*by(np,mm,ll,m)) + dyp*(amx
     1*by(nn,mp,ll,m) + dxp*by(np,mp,ll,m))) + dzp*(amy*(amx*by(nn,mm,lp
     2,m) + dxp*by(np,mm,lp,m)) + dyp*(amx*by(nn,mp,lp,m) + dxp*by(np,mp
     3,lp,m)))
      oz = amz*(amy*(amx*bz(nn,mm,ll,m) + dxp*bz(np,mm,ll,m)) + dyp*(amx
     1*bz(nn,mp,ll,m) + dxp*bz(np,mp,ll,m))) + dzp*(amy*(amx*bz(nn,mm,lp
     2,m) + dxp*bz(np,mm,lp,m)) + dyp*(amx*bz(nn,mp,lp,m) + dxp*bz(np,mp
     3,lp,m)))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j,m) + dx
      acy = part(5,j,m) + dy
      acz = part(6,j,m) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
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
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGBPUSH32L(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ek,nx,ny,
     1nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field, for distributed data
c with 2D spatial decomposition
c Using the Boris Mover.
c scalar version using guard cells
c 188 flops/particle, 1 divide, 54 loads, 6 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
c omz = (q/m)*bz(x(t),y(t),z(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c z(t+dt)=z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
c bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
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
c bxyz(1,j,k,l,m) = x component of magnetic field at grid (j,kk,ll)
c bxyz(2,j,k,l,m) = y component of magnetic field at grid (j,kk,ll)
c bxyz(3,j,k,l,m) = z component of magnetic field at grid (j,kk,ll)
c in other words, bxyz is the convolution of the magnetic field
c over the particle shape,
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
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
      dimension bxyz(3,nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      qtmh = .5*qbm*dt
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
c find electric field
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
c find magnetic field
      ox = amz*(amx*bxyz(1,nn,mm,ll,m) + amy*bxyz(1,np,mm,ll,m) + dyp*bx
     1yz(1,nn,mp,ll,m) + dx1*bxyz(1,np,mp,ll,m)) + dzp*(amx*bxyz(1,nn,mm
     2,lp,m) + amy*bxyz(1,np,mm,lp,m) + dyp*bxyz(1,nn,mp,lp,m) + dx1*bxy
     3z(1,np,mp,lp,m))
      oy = amz*(amx*bxyz(2,nn,mm,ll,m) + amy*bxyz(2,np,mm,ll,m) + dyp*bx
     1yz(2,nn,mp,ll,m) + dx1*bxyz(2,np,mp,ll,m)) + dzp*(amx*bxyz(2,nn,mm
     2,lp,m) + amy*bxyz(2,np,mm,lp,m) + dyp*bxyz(2,nn,mp,lp,m) + dx1*bxy
     3z(2,np,mp,lp,m))
      oz = amz*(amx*bxyz(3,nn,mm,ll,m) + amy*bxyz(3,np,mm,ll,m) + dyp*bx
     1yz(3,nn,mp,ll,m) + dx1*bxyz(3,np,mp,ll,m)) + dzp*(amx*bxyz(3,nn,mm
     2,lp,m) + amy*bxyz(3,np,mm,lp,m) + dyp*bxyz(3,nn,mp,lp,m) + dx1*bxy
     3z(3,np,mp,lp,m))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j,m) + dx
      acy = part(5,j,m) + dy
      acz = part(6,j,m) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(4,j,m) = dx
      part(5,j,m) = dy
      part(6,j,m) = dz
c new position
      dx = part(1,j,m) + dx*dtc
      dy = part(2,j,m) + dy*dtc
      dz = part(3,j,m) + dz*dtc
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
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSBPUSH32L(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ek,nx,ny
     1,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field, for distributed data
c with 2D spatial decomposition
c Using the Boris Mover,
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 188 flops/particle, 1 divide, 54 loads, 6 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
c omz = (q/m)*bz(x(t),y(t),z(t)).
c position equations used are:
c x(t+dt)=x(t) + vx(t+dt/2)*dt
c y(t+dt)=y(t) + vy(t+dt/2)*dt
c z(t+dt)=z(t) + vz(t+dt/2)*dt
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
c bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
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
c bxyz(1,j,k,l,m) = x component of magnetic field at grid (j,kk,ll)
c bxyz(2,j,k,l,m) = y component of magnetic field at grid (j,kk,ll)
c bxyz(3,j,k,l,m) = z component of magnetic field at grid (j,kk,ll)
c in other words, bxyz is the convolution of the magnetic field
c over the particle shape,
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c dtc = time interval between successive co-ordinate calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
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
      dimension fxyz(3,nxyzp,mnblok), bxyz(3,nxyzp,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      qtmh = .5*qbm*dt
      sum1 = 0.0d0
      nxyv = nxv*nypmx
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
c find interpolation weights
      do 10 j = 1, nop1
      nn = nnn + 1
      mm = nxv*mmm + nxyv*lll
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
      ll = mm + nxyv
      amy = dxp*amy
      lp = mp + nxyv
      mmm = mmm - mnoff
      lll = lll - lnoff
c find electric field
      dx = amz*(amx*fxyz(1,mm,m) + amy*fxyz(1,mm+1,m) + dyp*fxyz(1,mp,m)
     1 + dx1*fxyz(1,mp+1,m)) + dzp*(amx*fxyz(1,ll,m) + amy*fxyz(1,ll+1,m
     2) + dyp*fxyz(1,lp,m) + dx1*fxyz(1,lp+1,m))
      dy = amz*(amx*fxyz(2,mm,m) + amy*fxyz(2,mm+1,m) + dyp*fxyz(2,mp,m)
     1 + dx1*fxyz(2,mp+1,m)) + dzp*(amx*fxyz(2,ll,m) + amy*fxyz(2,ll+1,m
     2) + dyp*fxyz(2,lp,m) + dx1*fxyz(2,lp+1,m))
      dz = amz*(amx*fxyz(3,mm,m) + amy*fxyz(3,mm+1,m) + dyp*fxyz(3,mp,m)
     1 + dx1*fxyz(3,mp+1,m)) + dzp*(amx*fxyz(3,ll,m) + amy*fxyz(3,ll+1,m
     2) + dyp*fxyz(3,lp,m) + dx1*fxyz(3,lp+1,m))
c find magnetic field
      ox = amz*(amx*bxyz(1,mm,m) + amy*bxyz(1,mm+1,m) + dyp*bxyz(1,mp,m)
     1 + dx1*bxyz(1,mp+1,m)) + dzp*(amx*bxyz(1,ll,m) + amy*bxyz(1,ll+1,m
     2) + dyp*bxyz(1,lp,m) + dx1*bxyz(1,lp+1,m))
      oy = amz*(amx*bxyz(2,mm,m) + amy*bxyz(2,mm+1,m) + dyp*bxyz(2,mp,m)
     1 + dx1*bxyz(2,mp+1,m)) + dzp*(amx*bxyz(2,ll,m) + amy*bxyz(2,ll+1,m
     2) + dyp*bxyz(2,lp,m) + dx1*bxyz(2,lp+1,m))
      oz = amz*(amx*bxyz(3,mm,m) + amy*bxyz(3,mm+1,m) + dyp*bxyz(3,mp,m)
     1 + dx1*bxyz(3,mp+1,m)) + dzp*(amx*bxyz(3,ll,m) + amy*bxyz(3,ll+1,m
     2) + dyp*bxyz(3,lp,m) + dx1*bxyz(3,lp+1,m))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j,m) + dx
      acy = part(5,j,m) + dy
      acz = part(6,j,m) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(4,j,m) = dx
      part(5,j,m) = dy
      part(6,j,m) = dz
c new position
      dx = part(1,j,m) + dx*dtc
      dy = part(2,j,m) + dy*dtc
      dz = part(3,j,m) + dz*dtc
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
      mm = nxv*mmm + nxyv*lll
      amx = 1. - dxn
      amy = 1. - dyn
      mm = mm + nn
      dx1 = dxn*dyn
      dyp = amx*dyn
      mp = mm + nxv
      amx = amx*amy
      amz = 1. - dzn
      ll = mm + nxyv
      amy = dxn*amy
      lp = mp + nxyv
c find electric field
      dx = amz*(amx*fxyz(1,mm,m) + amy*fxyz(1,mm+1,m) + dyp*fxyz(1,mp,m)
     1 + dx1*fxyz(1,mp+1,m)) + dzn*(amx*fxyz(1,ll,m) + amy*fxyz(1,ll+1,m
     2) + dyp*fxyz(1,lp,m) + dx1*fxyz(1,lp+1,m))
      dy = amz*(amx*fxyz(2,mm,m) + amy*fxyz(2,mm+1,m) + dyp*fxyz(2,mp,m)
     1 + dx1*fxyz(2,mp+1,m)) + dzn*(amx*fxyz(2,ll,m) + amy*fxyz(2,ll+1,m
     2) + dyp*fxyz(2,lp,m) + dx1*fxyz(2,lp+1,m))
      dz = amz*(amx*fxyz(3,mm,m) + amy*fxyz(3,mm+1,m) + dyp*fxyz(3,mp,m)
     1 + dx1*fxyz(3,mp+1,m)) + dzn*(amx*fxyz(3,ll,m) + amy*fxyz(3,ll+1,m
     2) + dyp*fxyz(3,lp,m) + dx1*fxyz(3,lp+1,m))
c find magnetic field
      ox = amz*(amx*bxyz(1,mm,m) + amy*bxyz(1,mm+1,m) + dyp*bxyz(1,mp,m)
     1 + dx1*bxyz(1,mp+1,m)) + dzn*(amx*bxyz(1,ll,m) + amy*bxyz(1,ll+1,m
     2) + dyp*bxyz(1,lp,m) + dx1*bxyz(1,lp+1,m))
      oy = amz*(amx*bxyz(2,mm,m) + amy*bxyz(2,mm+1,m) + dyp*bxyz(2,mp,m)
     1 + dx1*bxyz(2,mp+1,m)) + dzn*(amx*bxyz(2,ll,m) + amy*bxyz(2,ll+1,m
     2) + dyp*bxyz(2,lp,m) + dx1*bxyz(2,lp+1,m))
      oz = amz*(amx*bxyz(3,mm,m) + amy*bxyz(3,mm+1,m) + dyp*bxyz(3,mp,m)
     1 + dx1*bxyz(3,mp+1,m)) + dzn*(amx*bxyz(3,ll,m) + amy*bxyz(3,ll+1,m
     2) + dyp*bxyz(3,lp,m) + dx1*bxyz(3,lp+1,m))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,nop,m) + dx
      acy = part(5,nop,m) + dy
      acz = part(6,nop,m) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      anorm = 2./(1. + omt)
      omt = .5*(1. - omt)
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(4,nop,m) = dx
      part(5,nop,m) = dy
      part(6,nop,m) = dz
c new position
      dx = part(1,nop,m) + dx*dtc
      dy = part(2,nop,m) + dy*dtc
      dz = part(3,nop,m) + dz*dtc
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
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPUSH32C(part,fx,fy,fz,bx,by,bz,npp,noff,qbm,dt,ek,nx,
     1idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field and periodic boundary
c conditions.  Using the Boris Mover with correction,
c for distributed data with 2D spatial decomposition
c baseline scalar distributed version
c 512 flops/particle, 2 divides, 2 sqrts, 168 loads, 6 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
c omz = (q/m)*bz(x(t),y(t),z(t)).
c position equations used are:
c x(t+dt)=x(t) + (vx(t+dt/2) - vplx)*dt/sqrt(1 + (om*dt/2)**2) + vplx*dt
c y(t+dt)=y(t) + (vy(t+dt/2) - vply)*dt/sqrt(1 + (om*dt/2)**2) + vply*dt
c z(t+dt)=z(t) + (vz(t+dt/2) - vplz)*dt/sqrt(1 + (om*dt/2)**2) + vplz*dt
c where vplx = omx*vpl/om, vply = omy*vpl/om, and vplz = omz*vpl/om,
c and vpl = (omx*vx(t+dt/2) + omy*vy(t+dt/2) + omz*vz(t+dt/2))/om
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
c bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
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
c similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
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
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx = system length in x direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first dimension of charge array, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
      double precision sum1
      dimension part(idimp,npmax,mnblok)
      dimension fx(nxv,nypmx,nzpmx,mnblok), fy(nxv,nypmx,nzpmx,mnblok)
      dimension fz(nxv,nypmx,nzpmx,mnblok), bx(nxv,nypmx,nzpmx,mnblok)
      dimension by(nxv,nypmx,nzpmx,mnblok), bz(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      zero = 0.
      anx = float(nx)
      qtmh = .5*qbm*dt
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
c find electric field
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
c find magnetic field
      ox = dzl*(dyl*(dxl*bx(nl,ml,lm,m) + amx*bx(nn,ml,lm,m) + dxp*bx(np
     1,ml,lm,m)) + amy*(dxl*bx(nl,mm,lm,m) + amx*bx(nn,mm,lm,m) + dxp*bx
     2(np,mm,lm,m)) + dyp*(dxl*bx(nl,mp,lm,m) + amx*bx(nn,mp,lm,m) + dxp
     3*bx(np,mp,lm,m))) + amz*(dyl*(dxl*bx(nl,ml,ll,m) + amx*bx(nn,ml,ll
     4,m) + dxp*bx(np,ml,ll,m)) + amy*(dxl*bx(nl,mm,ll,m) + amx*bx(nn,mm
     5,ll,m) + dxp*bx(np,mm,ll,m)) + dyp*(dxl*bx(nl,mp,ll,m) + amx*bx(nn
     6,mp,ll,m) + dxp*bx(np,mp,ll,m))) + dzp*(dyl*(dxl*bx(nl,ml,lp,m) + 
     7amx*bx(nn,ml,lp,m) + dxp*bx(np,ml,lp,m)) + amy*(dxl*bx(nl,mm,lp,m) 
     8 + amx*bx(nn,mm,lp,m) + dxp*bx(np,mm,lp,m)) + dyp*(dxl*bx(nl,mp,lp
     9,m) + amx*bx(nn,mp,lp,m) + dxp*bx(np,mp,lp,m)))
      oy = dzl*(dyl*(dxl*by(nl,ml,lm,m) + amx*by(nn,ml,lm,m) + dxp*by(np
     1,ml,lm,m)) + amy*(dxl*by(nl,mm,lm,m) + amx*by(nn,mm,lm,m) + dxp*by
     2(np,mm,lm,m)) + dyp*(dxl*by(nl,mp,lm,m) + amx*by(nn,mp,lm,m) + dxp
     3*by(np,mp,lm,m))) + amz*(dyl*(dxl*by(nl,ml,ll,m) + amx*by(nn,ml,ll
     4,m) + dxp*by(np,ml,ll,m)) + amy*(dxl*by(nl,mm,ll,m) + amx*by(nn,mm
     5,ll,m) + dxp*by(np,mm,ll,m)) + dyp*(dxl*by(nl,mp,ll,m) + amx*by(nn
     6,mp,ll,m) + dxp*by(np,mp,ll,m))) + dzp*(dyl*(dxl*by(nl,ml,lp,m) + 
     7amx*by(nn,ml,lp,m) + dxp*by(np,ml,lp,m)) + amy*(dxl*by(nl,mm,lp,m) 
     8 + amx*by(nn,mm,lp,m) + dxp*by(np,mm,lp,m)) + dyp*(dxl*by(nl,mp,lp
     9,m) + amx*by(nn,mp,lp,m) + dxp*by(np,mp,lp,m)))
      oz = dzl*(dyl*(dxl*bz(nl,ml,lm,m) + amx*bz(nn,ml,lm,m) + dxp*bz(np
     1,ml,lm,m)) + amy*(dxl*bz(nl,mm,lm,m) + amx*bz(nn,mm,lm,m) + dxp*bz
     2(np,mm,lm,m)) + dyp*(dxl*bz(nl,mp,lm,m) + amx*bz(nn,mp,lm,m) + dxp
     3*bz(np,mp,lm,m))) + amz*(dyl*(dxl*bz(nl,ml,ll,m) + amx*bz(nn,ml,ll
     4,m) + dxp*bz(np,ml,ll,m)) + amy*(dxl*bz(nl,mm,ll,m) + amx*bz(nn,mm
     5,ll,m) + dxp*bz(np,mm,ll,m)) + dyp*(dxl*bz(nl,mp,ll,m) + amx*bz(nn
     6,mp,ll,m) + dxp*bz(np,mp,ll,m))) + dzp*(dyl*(dxl*bz(nl,ml,lp,m) + 
     7amx*bz(nn,ml,lp,m) + dxp*bz(np,ml,lp,m)) + amy*(dxl*bz(nl,mm,lp,m) 
     8 + amx*bz(nn,mm,lp,m) + dxp*bz(np,mm,lp,m)) + dyp*(dxl*bz(nl,mp,lp
     9,m) + amx*bz(nn,mp,lp,m) + dxp*bz(np,mp,lp,m)))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j,m) + dx
      acy = part(5,j,m) + dy
      acz = part(6,j,m) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      amt = 0.
      if (omt.ne.0.) amt = 1./sqrt(omt)
      anorm = 1./(1. + omt)
      dtt = dt*sqrt(anorm)
      omt = .5*(1. - omt)
      anorm = anorm + anorm
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c direction cosines
      omxt = omxt*amt
      omyt = omyt*amt
      omzt = omzt*amt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(4,j,m) = dx
      part(5,j,m) = dy
      part(6,j,m) = dz
c parallel component of velocity
      omt = (dt - dtt)*(omxt*dx + omyt*dy + omzt*dz)
c new position
      dx = part(1,j,m) + (dx*dtt + omxt*omt)
      dy = part(2,j,m) + (dy*dtt + omyt*omt)
      dz = part(3,j,m) + (dz*dtt + omzt*omt)
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j,m) = dx
      part(2,j,m) = dy
      part(3,j,m) = dz
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGBPUSH32C(part,fxyz,bxyz,npp,noff,qbm,dt,ek,nx,ny,nz,i
     1dimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field, for distributed data 
c with 2D spatial decomposition
c Using the Boris Mover with correction.
c scalar version using guard cells,
c 469 flops/particle, 2 divides, 2 sqrts, 168 loads, 6 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
c omz = (q/m)*bz(x(t),y(t),z(t)).
c position equations used are:
c x(t+dt)=x(t) + (vx(t+dt/2) - vplx)*dt/sqrt(1 + (om*dt/2)**2) + vplx*dt
c y(t+dt)=y(t) + (vy(t+dt/2) - vply)*dt/sqrt(1 + (om*dt/2)**2) + vply*dt
c z(t+dt)=z(t) + (vz(t+dt/2) - vplz)*dt/sqrt(1 + (om*dt/2)**2) + vplz*dt
c where vplx = omx*vpl/om, vply = omy*vpl/om, and vplz = omz*vpl/om,
c and vpl = (omx*vx(t+dt/2) + omy*vy(t+dt/2) + omz*vz(t+dt/2))/om
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
c bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
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
c similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c fxyz(1,j+1,k+1,l,m) = x component of force/charge at grid (j,kk,ll)
c fxyz(2,j+1,k+1,l,m) = y component of force/charge at grid (j,kk,ll)
c fxyz(3,j+1,k+1,l,m) = z component of force/charge at grid (j,kk,ll)
c in other words, fxyz are the convolutions of the electric field
c over the particle shape, where ll = l + noff(m) - 1
c bxyz(1,j+1,k+1,l,m) = x component of magnetic field at grid (j,kk,ll)
c bxyz(2,j+1,k+1,l,m) = y component of magnetic field at grid (j,kk,ll)
c bxyz(3,j+1,k+1,l,m) = z component of magnetic field at grid (j,kk,ll)
c in other words, bxyz is the convolution of the magnetic field
c over the particle shape,
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m.
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c nxv = first dimension of charge array, must be >= nx+3
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      double precision sum1
      dimension part(idimp,npmax,mnblok)
      dimension fxyz(3,nxv,nypmx,nzpmx,mnblok)
      dimension bxyz(3,nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      qtmh = .5*qbm*dt
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
c find electric field
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
c find magnetic field
      ox = dzl*(dx1*bxyz(1,nl,ml,lm,m) + dx2*bxyz(1,nn,ml,lm,m) + dyl*bx
     1yz(1,np,ml,lm,m) + dy1*bxyz(1,nl,mm,lm,m) + dy2*bxyz(1,nn,mm,lm,m)
     2 + amy*bxyz(1,np,mm,lm,m) + dxl*bxyz(1,nl,mp,lm,m) + amx*bxyz(1,nn
     3,mp,lm,m) + dyp*bxyz(1,np,mp,lm,m))
      oy = dzl*(dx1*bxyz(2,nl,ml,lm,m) + dx2*bxyz(2,nn,ml,lm,m) + dyl*bx
     1yz(2,np,ml,lm,m) + dy1*bxyz(2,nl,mm,lm,m) + dy2*bxyz(2,nn,mm,lm,m)
     2 + amy*bxyz(2,np,mm,lm,m) + dxl*bxyz(2,nl,mp,lm,m) + amx*bxyz(2,nn
     3,mp,lm,m) + dyp*bxyz(2,np,mp,lm,m))
      oz = dzl*(dx1*bxyz(3,nl,ml,lm,m) + dx2*bxyz(3,nn,ml,lm,m) + dyl*bx
     1yz(3,np,ml,lm,m) + dy1*bxyz(3,nl,mm,lm,m) + dy2*bxyz(3,nn,mm,lm,m)
     2 + amy*bxyz(3,np,mm,lm,m) + dxl*bxyz(3,nl,mp,lm,m) + amx*bxyz(3,nn
     3,mp,lm,m) + dyp*bxyz(3,np,mp,lm,m))
      ox = ox + amz*(dx1*bxyz(1,nl,ml,ll,m) + dx2*bxyz(1,nn,ml,ll,m) + d
     1yl*bxyz(1,np,ml,ll,m) + dy1*bxyz(1,nl,mm,ll,m) + dy2*bxyz(1,nn,mm,
     2ll,m) + amy*bxyz(1,np,mm,ll,m) + dxl*bxyz(1,nl,mp,ll,m) + amx*bxyz
     3(1,nn,mp,ll,m) + dyp*bxyz(1,np,mp,ll,m))
      oy = oy + amz*(dx1*bxyz(2,nl,ml,ll,m) + dx2*bxyz(2,nn,ml,ll,m) + d
     1yl*bxyz(2,np,ml,ll,m) + dy1*bxyz(2,nl,mm,ll,m) + dy2*bxyz(2,nn,mm,
     2ll,m) + amy*bxyz(2,np,mm,ll,m) + dxl*bxyz(2,nl,mp,ll,m) + amx*bxyz
     3(2,nn,mp,ll,m) + dyp*bxyz(2,np,mp,ll,m))
      oz = oz + amz*(dx1*bxyz(3,nl,ml,ll,m) + dx2*bxyz(3,nn,ml,ll,m) + d
     1yl*bxyz(3,np,ml,ll,m) + dy1*bxyz(3,nl,mm,ll,m) + dy2*bxyz(3,nn,mm,
     2ll,m) + amy*bxyz(3,np,mm,ll,m) + dxl*bxyz(3,nl,mp,ll,m) + amx*bxyz
     3(3,nn,mp,ll,m) + dyp*bxyz(3,np,mp,ll,m))
      ox = ox + dzp*(dx1*bxyz(1,nl,ml,lp,m) + dx2*bxyz(1,nn,ml,lp,m) + d
     1yl*bxyz(1,np,ml,lp,m) + dy1*bxyz(1,nl,mm,lp,m) + dy2*bxyz(1,nn,mm,
     2lp,m) + amy*bxyz(1,np,mm,lp,m) + dxl*bxyz(1,nl,mp,lp,m) + amx*bxyz
     3(1,nn,mp,lp,m) + dyp*bxyz(1,np,mp,lp,m))
      oy = oy + dzp*(dx1*bxyz(2,nl,ml,lp,m) + dx2*bxyz(2,nn,ml,lp,m) + d
     1yl*bxyz(2,np,ml,lp,m) + dy1*bxyz(2,nl,mm,lp,m) + dy2*bxyz(2,nn,mm,
     2lp,m) + amy*bxyz(2,np,mm,lp,m) + dxl*bxyz(2,nl,mp,lp,m) + amx*bxyz
     3(2,nn,mp,lp,m) + dyp*bxyz(2,np,mp,lp,m))
      oz = oz + dzp*(dx1*bxyz(3,nl,ml,lp,m) + dx2*bxyz(3,nn,ml,lp,m) + d
     1yl*bxyz(3,np,ml,lp,m) + dy1*bxyz(3,nl,mm,lp,m) + dy2*bxyz(3,nn,mm,
     2lp,m) + amy*bxyz(3,np,mm,lp,m) + dxl*bxyz(3,nl,mp,lp,m) + amx*bxyz
     3(3,nn,mp,lp,m) + dyp*bxyz(3,np,mp,lp,m))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j,m) + dx
      acy = part(5,j,m) + dy
      acz = part(6,j,m) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      amt = 0.
      if (omt.ne.0.) amt = 1./sqrt(omt)
      anorm = 1./(1. + omt)
      dtt = dt*sqrt(anorm)
      omt = .5*(1. - omt)
      anorm = anorm + anorm
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c direction cosines
      omxt = omxt*amt
      omyt = omyt*amt
      omzt = omzt*amt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(4,j,m) = dx
      part(5,j,m) = dy
      part(6,j,m) = dz
c parallel component of velocity
      omt = (dt - dtt)*(omxt*dx + omyt*dy + omzt*dz)
c new position
      dx = part(1,j,m) + (dx*dtt + omxt*omt)
      dy = part(2,j,m) + (dy*dtt + omyt*omt)
      dz = part(3,j,m) + (dz*dtt + omzt*omt)
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
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPUSH32CL(part,fx,fy,fz,bx,by,bz,npp,noff,qbm,dt,ek,nx
     1,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field and periodic boundary
c conditions. Using the Boris Mover with correction,
c for distributed data with 2D spatial decomposition
c baseline scalar distributed version
c 224 flops/particle, 2 divides, 2 sqrts, 54 loads, 6 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
c omz = (q/m)*bz(x(t),y(t),z(t)).
c position equations used are:
c x(t+dt)=x(t) + (vx(t+dt/2) - vplx)*dt/sqrt(1 + (om*dt/2)**2) + vplx*dt
c y(t+dt)=y(t) + (vy(t+dt/2) - vply)*dt/sqrt(1 + (om*dt/2)**2) + vply*dt
c z(t+dt)=z(t) + (vz(t+dt/2) - vplz)*dt/sqrt(1 + (om*dt/2)**2) + vplz*dt
c where vplx = omx*vpl/om, vply = omy*vpl/om, and vplz = omz*vpl/om,
c and vpl = (omx*vx(t+dt/2) + omy*vy(t+dt/2) + omz*vz(t+dt/2))/om
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
c bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
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
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
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
      dimension fz(nxv,nypmx,nzpmx,mnblok), bx(nxv,nypmx,nzpmx,mnblok)
      dimension by(nxv,nypmx,nzpmx,mnblok), bz(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      zero = 0.
      anx = float(nx)
      qtmh = .5*qbm*dt
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
c find electric field
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
c find magnetic field
      ox = amz*(amy*(amx*bx(nn,mm,ll,m) + dxp*bx(np,mm,ll,m)) + dyp*(amx
     1*bx(nn,mp,ll,m) + dxp*bx(np,mp,ll,m))) + dzp*(amy*(amx*bx(nn,mm,lp
     2,m) + dxp*bx(np,mm,lp,m)) + dyp*(amx*bx(nn,mp,lp,m) + dxp*bx(np,mp
     3,lp,m)))
      oy = amz*(amy*(amx*by(nn,mm,ll,m) + dxp*by(np,mm,ll,m)) + dyp*(amx
     1*by(nn,mp,ll,m) + dxp*by(np,mp,ll,m))) + dzp*(amy*(amx*by(nn,mm,lp
     2,m) + dxp*by(np,mm,lp,m)) + dyp*(amx*by(nn,mp,lp,m) + dxp*by(np,mp
     3,lp,m)))
      oz = amz*(amy*(amx*bz(nn,mm,ll,m) + dxp*bz(np,mm,ll,m)) + dyp*(amx
     1*bz(nn,mp,ll,m) + dxp*bz(np,mp,ll,m))) + dzp*(amy*(amx*bz(nn,mm,lp
     2,m) + dxp*bz(np,mm,lp,m)) + dyp*(amx*bz(nn,mp,lp,m) + dxp*bz(np,mp
     3,lp,m)))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j,m) + dx
      acy = part(5,j,m) + dy
      acz = part(6,j,m) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      amt = 0.
      if (omt.ne.0.) amt = 1./sqrt(omt)
      anorm = 1./(1. + omt)
      dtt = dt*sqrt(anorm)
      omt = .5*(1. - omt)
      anorm = anorm + anorm
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c direction cosines
      omxt = omxt*amt
      omyt = omyt*amt
      omzt = omzt*amt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(4,j,m) = dx
      part(5,j,m) = dy
      part(6,j,m) = dz
c parallel component of velocity
      omt = (dt - dtt)*(omxt*dx + omyt*dy + omzt*dz)
c new position
      dx = part(1,j,m) + (dx*dtt + omxt*omt)
      dy = part(2,j,m) + (dy*dtt + omyt*omt)
      dz = part(3,j,m) + (dz*dtt + omzt*omt)
c periodic boundary conditions
      if (dx.lt.zero) dx = dx + anx
      if (dx.ge.anx) dx = dx - anx
      part(1,j,m) = dx
      part(2,j,m) = dy
      part(3,j,m) = dz
   10 continue
   20 continue
c normalize kinetic energy
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGBPUSH32CL(part,fxyz,bxyz,npp,noff,qbm,dt,ek,nx,ny,nz,
     1idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field, for distributed data
c with 2D spatial decomposition
c Using the Boris Mover with correction.
c scalar version using guard cells,
c 206 flops/particle, 2 divides, 2 sqrts, 54 loads, 6 stores
c input: all, output: part, ek
c velocity equations used are:
c vx(t+dt/2) = rot(1)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(2)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(3)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
c vy(t+dt/2) = rot(4)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(5)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(6)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
c vz(t+dt/2) = rot(7)*(vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(8)*(vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(9)*(vz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fz(x(t),y(t),z(t))*dt)
c where q/m is charge/mass, and the rotation matrix is given by:
c    rot(1) = (1 - (om*dt/2)**2 + 2*(omx*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(2) = 2*(omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(3) = 2*(-omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(4) = 2*(-omz*dt/2 + (omx*dt/2)*(omy*dt/2))/(1 + (om*dt/2)**2)
c    rot(5) = (1 - (om*dt/2)**2 + 2*(omy*dt/2)**2)/(1 + (om*dt/2)**2)
c    rot(6) = 2*(omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(7) = 2*(omy*dt/2 + (omx*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(8) = 2*(-omx*dt/2 + (omy*dt/2)*(omz*dt/2))/(1 + (om*dt/2)**2)
c    rot(9) = (1 - (om*dt/2)**2 + 2*(omz*dt/2)**2)/(1 + (om*dt/2)**2)
c and om**2 = omx**2 + omy**2 + omz**2
c the rotation matrix is determined by:
c omx = (q/m)*bx(x(t),y(t),z(t)), omy = (q/m)*by(x(t),y(t),z(t)), and
c omz = (q/m)*bz(x(t),y(t),z(t)).
c position equations used are:
c x(t+dt)=x(t) + (vx(t+dt/2) - vplx)*dt/sqrt(1 + (om*dt/2)**2) + vplx*dt
c y(t+dt)=y(t) + (vy(t+dt/2) - vply)*dt/sqrt(1 + (om*dt/2)**2) + vply*dt
c z(t+dt)=z(t) + (vz(t+dt/2) - vplz)*dt/sqrt(1 + (om*dt/2)**2) + vplz*dt
c where vplx = omx*vpl/om, vply = omy*vpl/om, and vplz = omz*vpl/om,
c and vpl = (omx*vx(t+dt/2) + omy*vy(t+dt/2) + omz*vz(t+dt/2))/om
c fx(x(t),y(t),z(t)), fy(x(t),y(t),z(t)), and fz(x(t),y(t),z(t)),
c bx(x(t),y(t),z(t)), by(x(t),y(t),z(t)), and bz(x(t),y(t),z(t))
c are approximated by interpolation from the nearest grid points:
c fx(x,y,z) = (1-dz)*((1-dy)*((1-dx)*fx(n,m,l)+dx*fx(n+1,m,l))
c                + dy*((1-dx)*fx(n,m+1,l) + dx*fx(n+1,m+1,l)))
c           + dz*((1-dy)*((1-dx)*fx(n,m,l+1)+dx*fx(n+1,m,l+1))
c                + dy*((1-dx)*fx(n,m+1,l+1) + dx*fx(n+1,m+1,l+1)))
c where n,m,l = leftmost grid points and dx = x-n, dy = y-m, dz = z-l
c similarly for fy(x,y,z), fz(x,y,z), bx(x,y,z), by(x,y,z), bz(x,y,z)
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
c bxyz(1,j,k,l,m) = x component of magnetic field at grid (j,kk,ll)
c bxyz(2,j,k,l,m) = y component of magnetic field at grid (j,kk,ll)
c bxyz(3,j,k,l,m) = z component of magnetic field at grid (j,kk,ll)
c in other words, bxyz is the convolution of the magnetic field
c over the particle shape,
c where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qbm = particle charge/mass ratio
c dt = time interval between successive calculations
c kinetic energy/mass at time t is also calculated, using
c ek = .5*sum((vx(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (vy(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 + 
c      .25*(vz(t+dt/2) + vz(t-dt/2))**2)
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
      dimension bxyz(3,nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), noff(idds,mnblok)
      qtmh = .5*qbm*dt
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
c find electric field
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
c find magnetic field
      ox = amz*(amx*bxyz(1,nn,mm,ll,m) + amy*bxyz(1,np,mm,ll,m) + dyp*bx
     1yz(1,nn,mp,ll,m) + dx1*bxyz(1,np,mp,ll,m)) + dzp*(amx*bxyz(1,nn,mm
     2,lp,m) + amy*bxyz(1,np,mm,lp,m) + dyp*bxyz(1,nn,mp,lp,m) + dx1*bxy
     3z(1,np,mp,lp,m))
      oy = amz*(amx*bxyz(2,nn,mm,ll,m) + amy*bxyz(2,np,mm,ll,m) + dyp*bx
     1yz(2,nn,mp,ll,m) + dx1*bxyz(2,np,mp,ll,m)) + dzp*(amx*bxyz(2,nn,mm
     2,lp,m) + amy*bxyz(2,np,mm,lp,m) + dyp*bxyz(2,nn,mp,lp,m) + dx1*bxy
     3z(2,np,mp,lp,m))
      oz = amz*(amx*bxyz(3,nn,mm,ll,m) + amy*bxyz(3,np,mm,ll,m) + dyp*bx
     1yz(3,nn,mp,ll,m) + dx1*bxyz(3,np,mp,ll,m)) + dzp*(amx*bxyz(3,nn,mm
     2,lp,m) + amy*bxyz(3,np,mm,lp,m) + dyp*bxyz(3,nn,mp,lp,m) + dx1*bxy
     3z(3,np,mp,lp,m))
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j,m) + dx
      acy = part(5,j,m) + dy
      acz = part(6,j,m) + dz
c time-centered kinetic energy
      sum1 = sum1 + (acx*acx + acy*acy + acz*acz)
c calculate cyclotron frequency
      omxt = qtmh*ox
      omyt = qtmh*oy
      omzt = qtmh*oz
c calculate rotation matrix
      omt = omxt*omxt + omyt*omyt + omzt*omzt
      amt = 0.
      if (omt.ne.0.) amt = 1./sqrt(omt)
      anorm = 1./(1. + omt)
      dtt = dt*sqrt(anorm)
      omt = .5*(1. - omt)
      anorm = anorm + anorm
      rot4 = omxt*omyt
      rot7 = omxt*omzt
      rot8 = omyt*omzt
      rot1 = omt + omxt*omxt
      rot5 = omt + omyt*omyt
      rot9 = omt + omzt*omzt
      rot2 = omzt + rot4
      rot4 = -omzt + rot4
      rot3 = -omyt + rot7
      rot7 = omyt + rot7
      rot6 = omxt + rot8
      rot8 = -omxt + rot8
c direction cosines
      omxt = omxt*amt
      omyt = omyt*amt
      omzt = omzt*amt
c new velocity
      dx = (rot1*acx + rot2*acy + rot3*acz)*anorm + dx
      dy = (rot4*acx + rot5*acy + rot6*acz)*anorm + dy
      dz = (rot7*acx + rot8*acy + rot9*acz)*anorm + dz
      part(4,j,m) = dx
      part(5,j,m) = dy
      part(6,j,m) = dz
c parallel component of velocity
      omt = (dt - dtt)*(omxt*dx + omyt*dy + omzt*dz)
c new position
      dx = part(1,j,m) + (dx*dtt + omxt*omt)
      dy = part(2,j,m) + (dy*dtt + omyt*omt)
      dz = part(3,j,m) + (dz*dtt + omzt*omt)
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
      ek = ek + .5*sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PRETARD32(part,npp,dtc,nx,ny,nz,idimp,npmax,mnblok,ipbc
     1)
c for 3d code, particle positions are retarded a half time-step
c for distributed data with 2D spatial decomposition.
c input: all, output: part
c equations used are:
c x(t+dt) = x(t) - vx(t+dt/2)*dtc, y(t+dt) = y(t) - vy(t+dt/2)*dtc,
c z(t+dt) = z(t) - vz(t+dt/2)*dtc
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = velocity vx of particle n in partition m
c part(5,n,m) = velocity vy of particle n in partition m
c part(6,n,m) = velocity vz of particle n in partition m
c npp(m) = number of particles in partition m
c dtc = time interval between successive co-ordinate calculations
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      dimension part(idimp,npmax,mnblok)
      dimension npp(mnblok)
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
c retard position half a time-step for current deposit
      do 10 j = 1, npp(m)
      dx = part(1,j,m) - part(4,j,m)*dtc
      dy = part(2,j,m) - part(5,j,m)*dtc
      dz = part(3,j,m) - part(6,j,m)*dtc
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
      return
      end
