c 3d parallel PIC library for pushing relativistic particles with
c magnetic field and depositing current
c with 2D domain decomposition
c written by viktor k. decyk, ucla
c copyright 1999, regents of the university of california
c update: december 20, 2003
c-----------------------------------------------------------------------
      subroutine PGRJPOST32(part,cu,npp,noff,qm,dt,ci,nx,ny,nz,idimp,npm
     1ax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
c for 3d code, this subroutine calculates particle current density
c using second-order spline interpolation, for relativistic particles
c and distributed data with 2D spatial decomposition
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 199 flops/particle, 1 divide, 1 sqrt, 87 loads, 84 stores
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
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = x momentum of particle n in partition m
c part(5,n,m) = y momentum of particle n in partition m
c part(6,n,m) = z momentum of particle n in partition m
c cu(i,j+1,k,l,m) = ith component of current density at grid point 
c (j,kk,ll), where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprical of velocity of light
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
      ci2 = ci*ci
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
c find inverse gamma
      vx = part(4,j,m)
      vy = part(5,j,m)
      vz = part(6,j,m)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
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
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
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
      subroutine PGSRJPOST32(part,cu,npp,noff,qm,dt,ci,nx,ny,nz,idimp,np
     1max,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
c for 3d code, this subroutine calculates particle current density
c using second-order spline interpolation, for relativistic particles
c and distributed data with 2D spatial decomposition
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 199 flops/particle, 1 divide, 1 sqrt, 87 loads, 84 stores
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
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = x momentum of particle n in partition m
c part(5,n,m) = y momentum of particle n in partition m
c part(6,n,m) = z momentum of particle n in partition m
c cu(i,n,m) = ith component of current density at grid point (j,kk,ll)
c where n = j + nxv*kk + nxv*nypmx*ll
c and kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprical of velocity of light
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
      ci2 = ci*ci
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
      vxn = part(4,1,m)
      vyn = part(5,1,m)
      vzn = part(6,1,m)
      p2 = vxn*vxn + vyn*vyn + vzn*vzn
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
c calculate weights
      dx = dx1*dzl
      dy = dx2*dzl
      dz = dyl*dzl
      vx = vxn*gami
      vy = vyn*gami
      vz = vzn*gami
c get momentum for next particle
      vxn = part(4,j,m)
      vyn = part(5,j,m)
      vzn = part(6,j,m)
      p2 = vxn*vxn + vyn*vyn + vzn*vzn
c deposit current
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
c find inverse gamma for next particle
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
      vx = vxn*gami
      vy = vyn*gami
      vz = vzn*gami
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
      subroutine PGSRJOST32X(part,cu,npp,noff,nn,amxyz,qm,dt,ci,nx,ny,nz
     1,idimp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n81,ipbc)
c for 3d code, this subroutine calculates particle current density
c using second-order spline interpolation, for relativistic particles
c and distributed data with 2D spatial decomposition
c with short vectors over independent weights,
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized version with guard cells and 1d addressing
c 199 flops/particle, 1 divide, 1 sqrt, 249 loads, 244 stores
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
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = x momentum of particle n in partition m
c part(5,n,m) = y momentum of particle n in partition m
c part(6,n,m) = z momentum of particle n in partition m
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
c ci = reciprical of velocity of light
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
      ci2 = ci*ci
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
c find inverse gamma
      vx = part(4,i+jb,m)
      vy = part(5,i+jb,m)
      vz = part(6,i+jb,m)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
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
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
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
      subroutine PGRJPOST32L(part,cu,npp,noff,qm,dt,ci,nx,ny,nz,idimp,np
     1max,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
c for 3d code, this subroutine calculates particle current density
c using first-order linear interpolation, for relativistic particles
c and distributed data with 2D spatial decomposition
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells
c 75 flops/particle, 1 divide, 1 sqrt, 30 loads, 27 stores
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
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = x momentum of particle n in partition m
c part(5,n,m) = y momentum of particle n in partition m
c part(6,n,m) = z momentum of particle n in partition m
c cu(i,j,k,l,m) = ith component of current density at grid point
c (j,kk,ll), where kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprical of velocity of light
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
      ci2 = ci*ci
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
c find inverse gamma
      vx = part(4,j,m)
      vy = part(5,j,m)
      vz = part(6,j,m)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
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
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
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
      subroutine PGSRJPOST32L(part,cu,npp,noff,qm,dt,ci,nx,ny,nz,idimp,n
     1pmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
c for 3d code, this subroutine calculates particle current density
c using first-order linear interpolation, for relativistic particles
c and distributed data with 2D spatial decomposition
c in addition, particle positions are advanced a half time-step
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 75 flops/particle, 1 divide, 1 sqrt, 30 loads, 27 stores
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
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = x momentum of particle n in partition m
c part(5,n,m) = y momentum of particle n in partition m
c part(6,n,m) = z momentum of particle n in partition m
c cu(i,n,m) = ith component of current density at grid point (j,kk,ll),
c where n = j + nxv*(k-1) + nxv*nyv*(ll-1)
c and kk = k + noff(1,m) - 1, and ll = l + noff(2,m) - 1
c npp(m) = number of particles in partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c qm = charge on particle, in units of e
c dt = time interval between successive calculations
c ci = reciprical of velocity of light
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
      ci2 = ci*ci
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
c find inverse gamma
      vxn = part(4,1,m)
      vyn = part(5,1,m)
      vzn = part(6,1,m)
      p2 = vxn*vxn + vyn*vyn + vzn*vzn
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
c calculate weights
      dx = amx*amz
      dz = amy*amz
      vx = vxn*gami
      vy = vyn*gami
      vz = vzn*gami
c get momentum for next particle
      vxn = part(4,j,m)
      vyn = part(5,j,m)
      vzn = part(6,j,m)
      p2 = vxn*vxn + vyn*vyn + vzn*vzn
c deposit current
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
c find inverse gamma for next particle
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
      vx = vxn*gami
      vy = vyn*gami
      vz = vzn*gami
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
      subroutine PGSRJOST32XL(part,cu,npp,noff,nn,amxyz,qm,dt,ci,nx,ny,n
     1z,idimp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n24,ipbc)
c for 3d code, this subroutine calculates particle current density
c using first-order linear interpolation, for relativistic particles
c and distributed data with 2D spatial decomposition
c with short vectors over independent weights, 
c as in j. schwartzmeier and t. hewitt, proc. 12th conf. on numerical
c simulation of plasmas, san francisco, ca, 1987.
c in addition, particle positions are advanced a half time-step
c vectorized distributed version with guard cells and 1d addressing
c 75 flops/particle, 1 divide, 1 sqrt, 78 loads, 75 stores
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
c and qci = qm*pi*gami, where i = x,y,z
c where gami = 1./sqrt(1.+sum(pi**2)*ci*ci)
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = x momentum of particle n in partition m
c part(5,n,m) = y momentum of particle n in partition m
c part(6,n,m) = z momentum of particle n in partition m
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
c ci = reciprical of velocity of light
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
      ci2 = ci*ci
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
c find inverse gamma
      vx = part(4,i+jb,m)
      vy = part(5,i+jb,m)
      vz = part(6,i+jb,m)
      p2 = vx*vx + vy*vy + vz*vz
      gami = 1.0/sqrt(1.0 + p2*ci2)
c calculate weights
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
      vx = vx*gami
      vy = vy*gami
      vz = vz*gami
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
      subroutine PGRPUSH32(part,fxyz,npp,noff,qbm,dt,ci,ek,nx,ny,nz,idim
     1p,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, for relativistic particles
c and distributed data with 2D spatial decomposition
c scalar version using guard cells,
c 248 flops/particle, 1 divide, 1 sqrt, 87 loads, 6 stores
c input: all, output: part, ek
c equations used are:
c px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c pz(t+dt/2) = pz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c z(t+dt) = z(t) + pz(t+dt/2)*dtg, where
c dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
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
c part(4,n,m) = momentum px of particle n in partition m
c part(5,n,m) = momentum py of particle n in partition m
c part(6,n,m) = momentum pz of particle n in partition m
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
c ci = reciprical of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gamma)
c where gamma = sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
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
      qtmh = .5*qbm*dt
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j,m) + dx
      acy = part(5,j,m) + dy
      acz = part(6,j,m) + dz
c time-centered kinetic energy
      p2 = acx*acx + acy*acy + acz*acz
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
c new velocity
      dx = acx + dx
      dy = acy + dy
      dz = acz + dz
      part(4,j,m) = dx
      part(5,j,m) = dy
      part(6,j,m) = dz
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dt/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j,m) + dx*dtg
      dy = part(2,j,m) + dy*dtg
      dz = part(3,j,m) + dz*dtg
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
      ek = ek + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSRPUSH32(part,fxyz,npp,noff,qbm,dt,ci,ek,nx,ny,nz,idi
     1mp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, for relativistic particles
c and distributed data with 2D spatial decomposition
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 248 flops/particle, 1 divide, 1 sqrt, 87 loads, 6 stores
c input: all, output: part, ek
c equations used are:
c px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c pz(t+dt/2) = pz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c z(t+dt) = z(t) + pz(t+dt/2)*dtg, where
c dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
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
c part(4,n,m) = momentum px of particle n in partition m
c part(5,n,m) = momentum py of particle n in partition m
c part(6,n,m) = momentum pz of particle n in partition m
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
c ci = reciprical of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gamma)
c where gamma = sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
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
      qtmh = .5*qbm*dt
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j,m) + dx
      acy = part(5,j,m) + dy
      acz = part(6,j,m) + dz
c time-centered kinetic energy
      p2 = acx*acx + acy*acy + acz*acz
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
c new velocity
      dx = acx + dx
      dy = acy + dy
      dz = acz + dz
      part(4,j,m) = dx
      part(5,j,m) = dy
      part(6,j,m) = dz
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dt/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j,m) + dx*dtg
      dy = part(2,j,m) + dy*dtg
      dz = part(3,j,m) + dz*dtg
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,nop,m) + dx
      acy = part(5,nop,m) + dy
      acz = part(6,nop,m) + dz
c time-centered kinetic energy
      p2 = acx*acx + acy*acy + acz*acz
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
c new velocity
      dx = acx + dx
      dy = acy + dy
      dz = acz + dz
      part(4,nop,m) = dx
      part(5,nop,m) = dy
      part(6,nop,m) = dz
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dt/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,nop,m) + dx*dtg
      dy = part(2,nop,m) + dy*dtg
      dz = part(3,nop,m) + dz*dtg
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
      ek = ek + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSRPUSH32F(part,fxyz,npp,noff,qbm,dt,ci,ek,nx,ny,nz,id
     1imp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, for relativistic particles
c and distributed data with 2D spatial decomposition
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c cases 9-10 in v.k.decyk et al, computers in physics 10, 290 (1996).
c 248 flops/particle, 1 divide, 1 sqrt, 87 loads, 6 stores
c input: all, output: part, ek
c equations used are:
c px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c pz(t+dt/2) = pz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c z(t+dt) = z(t) + pz(t+dt/2)*dtg, where
c dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
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
c part(4,n,m) = momentum px of particle n in partition m
c part(5,n,m) = momentum py of particle n in partition m
c part(6,n,m) = momentum pz of particle n in partition m
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
c ci = reciprical of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gamma)
c where gamma = sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space + tag storage + forces storage = 10
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
      qtmh = .5*qbm*dt
      ci2 = ci*ci
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

c save forces
      part(idimp-2,j,m) = dx*qbm
      part(idimp-1,j,m) = dy*qbm
      part(idimp,j,m) = dz*qbm

c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j,m) + dx
      acy = part(5,j,m) + dy
      acz = part(6,j,m) + dz
c time-centered kinetic energy
      p2 = acx*acx + acy*acy + acz*acz
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
c new velocity
      dx = acx + dx
      dy = acy + dy
      dz = acz + dz
      part(4,j,m) = dx
      part(5,j,m) = dy
      part(6,j,m) = dz
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dt/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j,m) + dx*dtg
      dy = part(2,j,m) + dy*dtg
      dz = part(3,j,m) + dz*dtg
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

c save forces
      part(idimp-2,nop,m) = dx*qbm
      part(idimp-1,nop,m) = dy*qbm
      part(idimp,nop,m) = dz*qbm

c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,nop,m) + dx
      acy = part(5,nop,m) + dy
      acz = part(6,nop,m) + dz
c time-centered kinetic energy
      p2 = acx*acx + acy*acy + acz*acz
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
c new velocity
      dx = acx + dx
      dy = acy + dy
      dz = acz + dz
      part(4,nop,m) = dx
      part(5,nop,m) = dy
      part(6,nop,m) = dz
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dt/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,nop,m) + dx*dtg
      dy = part(2,nop,m) + dy*dtg
      dz = part(3,nop,m) + dz*dtg
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
      ek = ek + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGRPUSH32L(part,fxyz,npp,noff,qbm,dt,ci,ek,nx,ny,nz,idi
     1mp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles
c and distributed data with 2D spatial decomposition
c scalar version using guard cells
c 100 flops/particle, 1 divide, 1 sqrt, 30 loads, 6 stores
c input: all, output: part, ek
c equations used are:
c px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c pz(t+dt/2) = pz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c z(t+dt) = z(t) + pz(t+dt/2)*dtg, where
c dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
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
c part(4,n,m) = momentum px of particle n in partition m
c part(5,n,m) = momentum py of particle n in partition m
c part(6,n,m) = momentum pz of particle n in partition m
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
c ci = reciprical of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gamma)
c where gamma = sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
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
      qtmh = .5*qbm*dt
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j,m) + dx
      acy = part(5,j,m) + dy
      acz = part(6,j,m) + dz
c time-centered kinetic energy
      p2 = acx*acx + acy*acy + acz*acz
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
c new velocity
      dx = acx + dx
      dy = acy + dy
      dz = acz + dz
      part(4,j,m) = dx
      part(5,j,m) = dy
      part(6,j,m) = dz
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dt/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j,m) + dx*dtg
      dy = part(2,j,m) + dy*dtg
      dz = part(3,j,m) + dz*dtg
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
      ek = ek + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSRPUSH32L(part,fxyz,npp,noff,qbm,dt,ci,ek,nx,ny,nz,id
     1imp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, for relativistic particles
c and distributed data with 2D spatial decomposition
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 100 flops/particle, 1 divide, 1 sqrt, 30 loads, 6 stores
c input: all, output: part, ek
c equations used are:
c px(t+dt/2) = px(t-dt/2) + (q/m)*fx(x(t),y(t),z(t))*dt,
c py(t+dt/2) = py(t-dt/2) + (q/m)*fy(x(t),y(t),z(t))*dt,
c pz(t+dt/2) = pz(t-dt/2) + (q/m)*fz(x(t),y(t),z(t))*dt,
c where q/m is charge/mass, and
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c z(t+dt) = z(t) + pz(t+dt/2)*dtg, where
c dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
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
c part(4,n,m) = momentum px of particle n in partition m
c part(5,n,m) = momentum py of particle n in partition m
c part(6,n,m) = momentum pz of particle n in partition m
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
c ci = reciprical of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gamma)
c where gamma = sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
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
      qtmh = .5*qbm*dt
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j,m) + dx
      acy = part(5,j,m) + dy
      acz = part(6,j,m) + dz
c time-centered kinetic energy
      p2 = acx*acx + acy*acy + acz*acz
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
c new velocity
      dx = acx + dx
      dy = acy + dy
      dz = acz + dz
      part(4,j,m) = dx
      part(5,j,m) = dy
      part(6,j,m) = dz
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dt/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j,m) + dx*dtg
      dy = part(2,j,m) + dy*dtg
      dz = part(3,j,m) + dz*dtg
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,nop,m) + dx
      acy = part(5,nop,m) + dy
      acz = part(6,nop,m) + dz
c time-centered kinetic energy
      p2 = acx*acx + acy*acy + acz*acz
      sum1 = sum1 + p2/(1.0 + sqrt(1.0 + p2*ci2))
c new velocity
      dx = acx + dx
      dy = acy + dy
      dz = acz + dz
      part(4,nop,m) = dx
      part(5,nop,m) = dy
      part(6,nop,m) = dz
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dt/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,nop,m) + dx*dtg
      dy = part(2,nop,m) + dy*dtg
      dz = part(3,nop,m) + dz*dtg
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
      ek = ek + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGRBPUSH32(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ci,ek,nx,
     1ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field for relativistic particles
c and distributed data with 2D spatial decomposition
c Using the Boris Mover.
c scalar version using guard cells, 
c 461 flops/particle, 4 divides, 2 sqrts, 168 loads, 6 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
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
c omx = (q/m)*bx(x(t),y(t),z(t))*gami, 
c omy = (q/m)*by(x(t),y(t),z(t))*gami,
c omz = (q/m)*bz(x(t),y(t),z(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c z(t+dt) = z(t) + pz(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
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
c part(4,n,m) = momentum px of particle n in partition m
c part(5,n,m) = momentum py of particle n in partition m
c part(6,n,m) = momentum pz of particle n in partition m
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
c ci = reciprical of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
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
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j,m) + dx
      acy = part(5,j,m) + dy
      acz = part(6,j,m) + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
c renormalize magnetic field
      qtmg = qtmh*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j,m) + dx*dtg
      dy = part(2,j,m) + dy*dtg
      dz = part(3,j,m) + dz*dtg
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
      ek = ek + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSRBPUSH32(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ci,ek,nx
     1,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and second-order spline
c interpolation in space, with magnetic field for relativistic particles
c and distributed data with 2D spatial decomposition
c Using the Boris Mover,
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing, for distributed data
c 461 flops/particle, 4 divides, 2 sqrts, 168 loads, 6 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
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
c omx = (q/m)*bx(x(t),y(t),z(t))*gami, 
c omy = (q/m)*by(x(t),y(t),z(t))*gami,
c omz = (q/m)*bz(x(t),y(t),z(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c z(t+dt) = z(t) + pz(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
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
c part(4,n,m) = momentum px of particle n in partition m
c part(5,n,m) = momentum py of particle n in partition m
c part(6,n,m) = momentum pz of particle n in partition m
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
c ci = reciprical of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
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
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j,m) + dx
      acy = part(5,j,m) + dy
      acz = part(6,j,m) + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
c renormalize magnetic field
      qtmg = qtmh*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j,m) + dx*dtg
      dy = part(2,j,m) + dy*dtg
      dz = part(3,j,m) + dz*dtg
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,nop,m) + dx
      acy = part(5,nop,m) + dy
      acz = part(6,nop,m) + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
c renormalize magnetic field
      qtmg = qtmh*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,nop,m) + dx*dtg
      dy = part(2,nop,m) + dy*dtg
      dz = part(3,nop,m) + dz*dtg
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
      ek = ek + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGRBPUSH32L(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ci,ek,nx
     1,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field for relativistic particles
c and distributed data with 2D spatial decomposition
c Using the Boris Mover.
c scalar version using guard cells
c 198 flops/particle, 4 divides, 2 sqrts, 54 loads, 6 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
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
c omx = (q/m)*bx(x(t),y(t),z(t))*gami, 
c omy = (q/m)*by(x(t),y(t),z(t))*gami,
c omz = (q/m)*bz(x(t),y(t),z(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c z(t+dt) = z(t) + pz(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
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
c part(4,n,m) = momentum px of particle n in partition m
c part(5,n,m) = momentum py of particle n in partition m
c part(6,n,m) = momentum pz of particle n in partition m
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
c ci = reciprical of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
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
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j,m) + dx
      acy = part(5,j,m) + dy
      acz = part(6,j,m) + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
c renormalize magnetic field
      qtmg = qtmh*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j,m) + dx*dtg
      dy = part(2,j,m) + dy*dtg
      dz = part(3,j,m) + dz*dtg
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
      ek = ek + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PGSRBPUSH32L(part,fxyz,bxyz,npp,noff,qbm,dt,dtc,ci,ek,n
     1x,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
c for 3d code, this subroutine updates particle co-ordinates and
c velocities using leap-frog scheme in time and first-order linear
c interpolation in space, with magnetic field for relativistic particles
c and distributed data with 2D spatial decomposition
c Using the Boris Mover,
c scalar version using guard cells, integer conversion precalculation,
c and 1d addressing
c 198 flops/particle, 4 divides, 2 sqrts, 54 loads, 6 stores
c input: all, output: part, ek
c momentum equations used are:
c px(t+dt/2) = rot(1)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(2)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(3)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fx(x(t),y(t),z(t))*dt)
c py(t+dt/2) = rot(4)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(5)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(6)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
c    .5*(q/m)*fy(x(t),y(t),z(t))*dt)
c pz(t+dt/2) = rot(7)*(px(t-dt/2) + .5*(q/m)*fx(x(t),y(t),z(t))*dt) +
c    rot(8)*(py(t-dt/2) + .5*(q/m)*fy(x(t),y(t),z(t))*dt) +
c    rot(9)*(pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t),z(t))*dt) +
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
c omx = (q/m)*bx(x(t),y(t),z(t))*gami, 
c omy = (q/m)*by(x(t),y(t),z(t))*gami,
c omz = (q/m)*bz(x(t),y(t),z(t))*gami,
c where gami = 1./sqrt(1.+(px(t)*px(t)+py(t)*py(t)+pz(t)*pz(t))*ci*ci)
c position equations used are:
c x(t+dt) = x(t) + px(t+dt/2)*dtg
c y(t+dt) = y(t) + py(t+dt/2)*dtg
c z(t+dt) = z(t) + pz(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
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
c part(4,n,m) = momentum px of particle n in partition m
c part(5,n,m) = momentum py of particle n in partition m
c part(6,n,m) = momentum pz of particle n in partition m
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
c ci = reciprical of velocity of light
c kinetic energy/mass at time t is also calculated, using
c ek = gami*sum((px(t-dt/2) + .5*(q/m)*fx(x(t),y(t))*dt)**2 +
c      (py(t-dt/2) + .5*(q/m)*fy(x(t),y(t))*dt)**2 +
c      (pz(t-dt/2) + .5*(q/m)*fz(x(t),y(t))*dt)**2)/(1. + gami)
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
      ci2 = ci*ci
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,j,m) + dx
      acy = part(5,j,m) + dy
      acz = part(6,j,m) + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
c renormalize magnetic field
      qtmg = qtmh*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,j,m) + dx*dtg
      dy = part(2,j,m) + dy*dtg
      dz = part(3,j,m) + dz*dtg
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
c calculate half impulse
      dx = qtmh*dx
      dy = qtmh*dy
      dz = qtmh*dz
c half acceleration
      acx = part(4,nop,m) + dx
      acy = part(5,nop,m) + dy
      acz = part(6,nop,m) + dz
c find inverse gamma
      p2 = acx*acx + acy*acy + acz*acz
      gami = 1.0/sqrt(1.0 + p2*ci2)
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
c renormalize magnetic field
      qtmg = qtmh*gami
c time-centered kinetic energy
      sum1 = sum1 + gami*p2/(1.0 + gami)
c calculate cyclotron frequency
      omxt = qtmg*ox
      omyt = qtmg*oy
      omzt = qtmg*oz
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
c update inverse gamma
      p2 = dx*dx + dy*dy + dz*dz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c new position
      dx = part(1,nop,m) + dx*dtg
      dy = part(2,nop,m) + dy*dtg
      dz = part(3,nop,m) + dz*dtg
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
      ek = ek + sum1
      return
      end
c-----------------------------------------------------------------------
      subroutine PRRETARD32(part,npp,dtc,ci,nx,ny,nz,idimp,npmax,mnblok,
     1ipbc)
c for 3d code, particle positions are retarded a half time-step
c for distributed data with 2D spatial decomposition.
c input: all, output: part
c equations used are:
c x(t+dt) = x(t) - px(t+dt/2)*dtg
c y(t+dt) = y(t) - py(t+dt/2)*dtg
c z(t+dt) = z(t) - pz(t+dt/2)*dtg
c where dtg = dtc/sqrt(1.+(px(t+dt/2)*px(t+dt/2)+py(t+dt/2)*py(t+dt/2)+
c pz(t+dt/2)*pz(t+dt/2))*ci*ci)
c part(1,n,m) = position x of particle n in partition m
c part(2,n,m) = position y of particle n in partition m
c part(3,n,m) = position z of particle n in partition m
c part(4,n,m) = momentum px of particle n in partition m
c part(5,n,m) = momentum py of particle n in partition m
c part(6,n,m) = momentum pz of particle n in partition m
c npp(m) = number of particles in partition m
c dtc = time interval between successive co-ordinate calculations
c ci = reciprical of velocity of light
c nx/ny/nz = system length in x/y/z direction
c idimp = size of phase space = 6
c npmax = maximum number of particles in each partition
c mnblok = number of particle partitions.
c ipbc = particle boundary condition = (0,1,2,3) =
c (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      dimension part(idimp,npmax,mnblok)
      dimension npp(mnblok)
      ci2 = ci*ci
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
      do 10 j = 1, npp(m)
c find inverse gamma
      vx = part(4,j,m)
      vy = part(5,j,m)
      vz = part(6,j,m)
      p2 = vx*vx + vy*vy + vz*vz
      dtg = dtc/sqrt(1.0 + p2*ci2)
c retard position half a time-step for current deposit
      dx = part(1,j,m) - part(4,j,m)*dtg
      dy = part(2,j,m) - part(5,j,m)*dtg
      dz = part(3,j,m) - part(6,j,m)*dtg
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
