c 3d parallel PIC library for solving field equations
c with 2D domain decomposition
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: april 19, 2008
c-----------------------------------------------------------------------
      subroutine PCGUARD32X(fxyz,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
c replicate extended field
c quadratic interpolation, for distributed data with 2D decomposition
      implicit none
      real fxyz
      integer nyzp, nx, nxe, nypmx, nzpmx, idds, mnblok
      dimension fxyz(3,nxe,nypmx,nzpmx,mnblok), nyzp(idds,mnblok)
      integer i, k, l, m, nyp3, nzp3
      do 40 m = 1, mnblok
      nyp3 = nyzp(1,m) + 3
      nzp3 = nyzp(2,m) + 3
      do 30 l = 1, nzp3
      do 20 k = 1, nyp3
      do 10 i = 1, 3
      fxyz(i,1,k,l,m) = fxyz(i,nx+1,k,l,m)
      fxyz(i,nx+2,k,l,m) = fxyz(i,2,k,l,m)
      fxyz(i,nx+3,k,l,m) = fxyz(i,3,k,l,m)
   10 continue
   20 continue
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDGUARD32X(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
c replicate extended scalar field
c quadratic interpolation, for distributed data with 2D decomposition
      implicit none
      real q
      integer nyzp, nx, nxe, nypmx, nzpmx, idds, mnblok
      dimension q(nxe,nypmx,nzpmx,mnblok), nyzp(idds,mnblok)
      integer  k, l, m, nyp3, nzp3
      do 30 m = 1, mnblok
      nyp3 = nyzp(1,m) + 3
      nzp3 = nyzp(2,m) + 3
      do 20 l = 1, nzp3
      do 10 k = 1, nyp3
      q(1,k,l,m) = q(nx+1,k,l,m)
      q(nx+2,k,l,m) = q(2,k,l,m)
      q(nx+3,k,l,m) = q(3,k,l,m)
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSCGUARD32(cu,nyzp,xj0,yj0,zj0,nx,nxe,nypmx,nzpmx,idds,
     1mnblok)
c initialize extended periodic field
c quadratic interpolation, for distributed data with 2D decomposition
      implicit none
      real cu, xj0, yj0, zj0
      integer nyzp, nx, nxe, nypmx, nzpmx, idds, mnblok
      dimension cu(3,nxe,nypmx,nzpmx,mnblok), nyzp(idds,mnblok)
      integer i, j, k, l, m, nyp3, nzp3, nx3
c initialize extended field, with zero in the edges
      nx3 = nx + 3
      do 100 m = 1, mnblok
      nyp3 = nyzp(1,m) + 3
      nzp3 = nyzp(2,m) + 3
      do 60 l = 1, nyzp(2,m)
      do 30 k = 1, nyzp(1,m)
      do 10 j = 1, nx
      cu(1,j+1,k+1,l+1,m) = xj0
      cu(2,j+1,k+1,l+1,m) = yj0
      cu(3,j+1,k+1,l+1,m) = zj0
   10 continue
      do 20 i = 1, 3
      cu(i,1,k+1,l+1,m) = 0.
      cu(i,nx+2,k+1,l+1,m) = 0.
      cu(i,nx+3,k+1,l+1,m) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx3
      do 40 i = 1, 3
      cu(i,j,1,l+1,m) = 0.
      cu(i,j,nyp3-1,l+1,m) = 0.
      cu(i,j,nyp3,l+1,m) = 0.
   40 continue
   50 continue
   60 continue
      do 90 k = 1, nyp3
      do 80 j = 1, nx3
      do 70 i = 1, 3
      cu(i,j,k,1,m) = 0.
      cu(i,j,k,nzp3-1,m) = 0.
      cu(i,j,k,nzp3,m) = 0.
   70 continue
   80 continue
   90 continue
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSGUARD32(q,nyzp,qi0,nx,nxe,nypmx,nzpmx,idds,mnblok)
c initialize extended peridoc scalar field
c quadratic interpolation, for distributed data with 2D decomposition
      implicit none
      real q, qi0
      integer nyzp, nx, nxe, nypmx, nzpmx, idds, mnblok
      dimension q(nxe,nypmx,nzpmx,mnblok), nyzp(idds,mnblok)
      integer j, k, l, m, nyp3, nzp3, nx3
c initialize extended field, with zero in the edges
      nx3 = nx + 3
      do 70 m = 1, mnblok
      nyp3 = nyzp(1,m) + 3
      nzp3 = nyzp(2,m) + 3
      do 40 l = 1, nyzp(2,m)
      do 20 k = 1, nyzp(1,m)
      do 10 j = 1, nx
      q(j+1,k+1,l+1,m) = qi0
   10 continue
      q(1,k+1,l+1,m) = 0.
      q(nx+2,k+1,l+1,m) = 0.
      q(nx+3,k+1,l+1,m) = 0.
   20 continue
      do 30 j = 1, nx3
      q(j,1,l+1,m) = 0.
      q(j,nyp3-1,l+1,m) = 0.
      q(j,nyp3,l+1,m) = 0. 
   30 continue
   40 continue
      do 60 k = 1, nyp3
      do 50 j = 1, nx3
      q(j,k,1,m) = 0.
      q(j,k,nzp3-1,m) = 0.
      q(j,k,nzp3,m) = 0.
   50 continue
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PACGUARD32X(cu,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
c accumulate extended periodic vector field
c quadratic interpolation, for distributed data with 2D decomposition
      implicit none
      real cu
      integer nyzp, nx, nxe, nypmx, nzpmx, idds, mnblok
      dimension cu(3,nxe,nypmx,nzpmx,mnblok), nyzp(idds,mnblok)
      integer i, k, l, m, nyp3, nzp3
c accumulate edges of extended field
      do 40 m = 1, mnblok
      nyp3 = nyzp(1,m) + 3
      nzp3 = nyzp(2,m) + 3
      do 30 l = 1, nzp3
      do 20 k = 1, nyp3
      do 10 i = 1, 3
      cu(i,2,k,l,m) = cu(i,2,k,l,m) + cu(i,nx+2,k,l,m)
      cu(i,3,k,l,m) = cu(i,3,k,l,m) + cu(i,nx+3,k,l,m)
      cu(i,nx+1,k,l,m) = cu(i,nx+1,k,l,m) + cu(i,1,k,l,m)
      cu(i,1,k,l,m) = 0.
      cu(i,nx+2,k,l,m) = 0.
      cu(i,nx+3,k,l,m) = 0.
   10 continue
   20 continue
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PAGUARD32X(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
c accumulate extended periodic scalar field
c quadratic interpolation, for distributed data with 2D decomposition
      implicit none
      real q
      integer nyzp, nx, nxe, nypmx, nzpmx, idds, mnblok
      dimension q(nxe,nypmx,nzpmx,mnblok), nyzp(idds,mnblok)
      integer k, l, m, nyp3, nzp3
c accumulate edges of extended field
      do 30 m = 1, mnblok
      nyp3 = nyzp(1,m) + 3
      nzp3 = nyzp(2,m) + 3
      do 20 l = 1, nzp3
      do 10 k = 1, nyp3
      q(2,k,l,m) = q(2,k,l,m) + q(nx+2,k,l,m)
      q(3,k,l,m) = q(3,k,l,m) + q(nx+3,k,l,m)
      q(nx+1,k,l,m) = q(nx+1,k,l,m) + q(1,k,l,m)
      q(1,k,l,m) = 0.
      q(nx+2,k,l,m) = 0.
      q(nx+3,k,l,m) = 0.
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PZCGUARD32(cu,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
c zero out guard cells in extended periodic vector field
c quadratic interpolation, for distributed data with 2D decomposition
      implicit none
      real cu
      integer nyzp, nx, nxe, nypmx, nzpmx, idds, mnblok
      dimension cu(3,nxe,nypmx,nzpmx,mnblok), nyzp(idds,mnblok)
      integer i, j, k, l, m, nx3, nyp3
      nx3 = nx + 3
      do 90 m = 1, mnblok
      nyp3 = nyzp(1,m) + 3
c zero out guard cells in x
      do 50 l = 1, nyzp(2,m)
      do 20 k = 1, nyzp(1,m)
      do 10 i = 1, 3
      cu(i,1,k+1,l+1,m) = 0.
      cu(i,nx+2,k+1,l+1,m) = 0.
      cu(i,nx+3,k+1,l+1,m) = 0.
   10 continue
   20 continue
c zero out guard cells in y
      do 40 j = 1, nx3
      do 30 i = 1, 3
      cu(i,j,1,l+1,m) = 0.
      cu(i,j,nyzp(1,m)+2,l+1,m) = 0.
      cu(i,j,nyzp(1,m)+3,l+1,m) = 0.
   30 continue
   40 continue
   50 continue
c zero out guard cells in z
      do 80 k = 1, nyp3
      do 70 j = 1, nx3
      do 60 i = 1, 3
      cu(i,j,k,1,m) = 0.
      cu(i,j,k,nyzp(2,m)+2,m) = 0.
      cu(i,j,k,nyzp(2,m)+3,m) = 0.
   60 continue
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PZGUARD32(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
c zero out guard cells in extended periodic scalar field
c quadratic interpolation, for distributed data with 2D decomposition
      implicit none
      real q
      integer nyzp, nx, nxe, nypmx, nzpmx, idds, mnblok
      dimension q(nxe,nypmx,nzpmx,mnblok), nyzp(idds,mnblok)
      integer j, k, l, m, nx3, nyp3
      nx3 = nx + 3
      do 60 m = 1, mnblok
      nyp3 = nyzp(1,m) + 3
c zero out guard cells in x
      do 30 l = 1, nyzp(2,m)
      do 10 k = 1, nyzp(1,m)
      q(1,k+1,l+1,m) = 0.
      q(nx+2,k+1,l+1,m) = 0.
      q(nx+3,k+1,l+1,m) = 0.
   10 continue
c zero out guard cells in y
      do 20 j = 1, nx3
      q(j,1,l+1,m) = 0.
      q(j,nyzp(1,m)+2,l+1,m) = 0.
      q(j,nyzp(1,m)+3,l+1,m) = 0.
   20 continue
   30 continue
c zero out guard cells in z
      do 50 k = 1, nyp3
      do 40 j = 1, nx3
      q(j,k,1,m) = 0.
      q(j,k,nyzp(2,m)+2,m) = 0.
      q(j,k,nyzp(2,m)+3,m) = 0.
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCGUARD32XL(fxyz,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
c replicate extended field
c linear interpolation, for distributed data with 2D decomposition
      implicit none
      real fxyz
      integer nyzp, nx, nxe, nypmx, nzpmx, idds, mnblok
      dimension fxyz(3,nxe,nypmx,nzpmx,mnblok), nyzp(idds,mnblok)
      integer i,  k, l, m, nyp1, nzp1
      do 40 m = 1, mnblok
      nyp1 = nyzp(1,m) + 1
      nzp1 = nyzp(2,m) + 1
      do 30 l = 1, nzp1
      do 20 k = 1, nyp1
      do 10 i = 1, 3
      fxyz(i,nx+1,k,l,m) = fxyz(i,1,k,l,m)
   10 continue
   20 continue
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDGUARD32XL(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
c replicate extended scalar field
c linear interpolation, for distributed data with 2D decomposition
      implicit none
      real q
      integer nyzp, nx, nxe, nypmx, nzpmx, idds, mnblok
      dimension q(nxe,nypmx,nzpmx,mnblok), nyzp(idds,mnblok)
      integer k, l, m, nyp1, nzp1
      do 30 m = 1, mnblok
      nyp1 = nyzp(1,m) + 1
      nzp1 = nyzp(2,m) + 1
      do 20 l = 1, nzp1
      do 10 k = 1, nyp1
      q(nx+1,k,l,m) = q(1,k,l,m)
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSCGUARD32L(cu,nyzp,xj0,yj0,zj0,nx,nxe,nypmx,nzpmx,idds
     1,mnblok)
c initialize extended periodic field
c linear interpolation, for distributed data with 2D decomposition
      implicit none
      real cu, xj0, yj0, zj0
      integer nyzp, nx, nxe, nypmx, nzpmx, idds, mnblok
      dimension cu(3,nxe,nypmx,nzpmx,mnblok), nyzp(idds,mnblok)
      integer i, j, k, l, m, nyp1, nzp1, nx1
c initialize extended field, with zero in the edges
      nx1 = nx + 1
      do 100 m = 1, mnblok
      nyp1 = nyzp(1,m) + 1
      nzp1 = nyzp(2,m) + 1
      do 60 l = 1, nyzp(2,m)
      do 30 k = 1, nyzp(1,m)
      do 10 j = 1, nx
      cu(1,j,k,l,m) = xj0
      cu(2,j,k,l,m) = yj0
      cu(3,j,k,l,m) = zj0
   10 continue
      do 20 i = 1, 3
      cu(i,nx+1,k,l,m) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx1
      do 40 i = 1, 3
      cu(i,j,nyp1,l,m) = 0.
   40 continue
   50 continue
   60 continue
      do 90 k = 1, nyp1
      do 80 j = 1, nx1
      do 70 i = 1, 3
      cu(i,j,k,nzp1,m) = 0.
   70 continue
   80 continue
   90 continue
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PSGUARD32L(q,nyzp,qi0,nx,nxe,nypmx,nzpmx,idds,mnblok)
c initialize extended periodic scalar field
c linear interpolation, for distributed data with 2D decomposition
      implicit none
      real q, qi0
      integer nyzp, nx, nxe, nypmx, nzpmx, idds, mnblok
      dimension q(nxe,nypmx,nzpmx,mnblok), nyzp(idds,mnblok)
      integer j, k, l, m, nyp1, nzp1, nx1
c initialize extended field, with zero in the edges
      nx1 = nx + 1
      do 70 m = 1, mnblok
      nyp1 = nyzp(1,m) + 1
      nzp1 = nyzp(2,m) + 1
      do 40 l = 1, nyzp(2,m)
      do 20 k = 1, nyzp(1,m)
      do 10 j = 1, nx
      q(j,k,l,m) = qi0
   10 continue
      q(nx+1,k,l,m) = 0.
   20 continue
      do 30 j = 1, nx1
      q(j,nyp1,l,m) = 0.
   30 continue
   40 continue
      do 60 k = 1, nyp1
      do 50 j = 1, nx1
      q(j,k,nzp1,m) = 0.
   50 continue
   60 continue
   70 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PACGUARD32XL(cu,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
c accumulate extended periodic vector field
c linear interpolation, for distributed data with 2D decomposition
      implicit none
      real cu
      integer nyzp, nx, nxe, nypmx, nzpmx, idds, mnblok
      dimension cu(3,nxe,nypmx,nzpmx,mnblok), nyzp(idds,mnblok)
      integer i, k, l, m, nyp1, nzp1
c accumulate edges of extended field
      do 40 m = 1, mnblok
      nyp1 = nyzp(1,m) + 1
      nzp1 = nyzp(2,m) + 1
      do 30 l = 1, nzp1
      do 20 k = 1, nyp1
      do 10 i = 1, 3
      cu(i,1,k,l,m) = cu(i,1,k,l,m) + cu(i,nx+1,k,l,m)
      cu(i,nx+1,k,l,m) = 0.
   10 continue
   20 continue
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PAGUARD32XL(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
c accumulate extended periodic scalar field
c linear interpolation, for distributed data with 2D decomposition
      implicit none
      real q
      integer nyzp, nx, nxe, nypmx, nzpmx, idds, mnblok
      dimension q(nxe,nypmx,nzpmx,mnblok), nyzp(idds,mnblok)
      integer k, l, m, nyp1, nzp1
c accumulate edges of extended field
      do 30 m = 1, mnblok
      nyp1 = nyzp(1,m) + 1
      nzp1 = nyzp(2,m) + 1
      do 20 l = 1, nzp1
      do 10 k = 1, nyp1
      q(1,k,l,m) = q(1,k,l,m) + q(nx+1,k,l,m)
      q(nx+1,k,l,m) = 0.
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PZCGUARD32L(cu,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
c zero out guard cells in extended periodic vector field
c linear interpolation, for distributed data with 2D decomposition
      implicit none
      real cu
      integer nyzp, nx, nxe, nypmx, nzpmx, idds, mnblok
      dimension cu(3,nxe,nypmx,nzpmx,mnblok), nyzp(idds,mnblok)
      integer i, j, k, l, m, nx1, nyp1
      nx1 = nx + 1
c accumulate edges of extended field
      do 90 m = 1, mnblok
      nyp1 = nyzp(1,m) + 1
c zero out guard cells in x
      do 50 l = 1, nyzp(2,m)
      do 20 k = 1, nyzp(1,m)
      do 10 i = 1, 3
      cu(i,nx+1,k,l,m) = 0.
   10 continue
   20 continue
c zero out guard cells in y
      do 40 j = 1, nx1
      do 30 i = 1, 3
      cu(i,j,nyzp(1,m)+1,l,m) = 0.
   30 continue
   40 continue
   50 continue
c zero out guard cells in z
      do 80 k = 1, nyp1
      do 70 j = 1, nx1
      do 60 i = 1, 3
      cu(i,j,k,nyzp(2,m)+1,m) = 0.
   60 continue
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PZGUARD32L(q,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
c zero out guard cells in extended periodic scalar field
c linear interpolation, for distributed data with 2D decomposition
      implicit none
      real q
      integer nyzp, nx, nxe, nypmx, nzpmx, idds, mnblok
      dimension q(nxe,nypmx,nzpmx,mnblok), nyzp(idds,mnblok)
      integer j, k, l, m, nx1, nyp1
      nx1 = nx + 1
c accumulate edges of extended field
      do 60 m = 1, mnblok
      nyp1 = nyzp(1,m) + 1
c zero out guard cells in x
      do 30 l = 1, nyzp(2,m)
      do 10 k = 1, nyzp(1,m)
      q(nx+1,k,l,m) = 0.
   10 continue
c zero out guard cells in y
      do 20 j = 1, nx1
      q(j,nyzp(1,m)+1,l,m) = 0.
   20 continue
   30 continue
c zero out guard cells in z
      do 50 k = 1, nyp1
      do 40 j = 1, nx1
      q(j,k,nyzp(2,m)+1,m) = 0.
   40 continue
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOIS32INIT(ffc,ax,ay,az,affp,nx,ny,nz,kstrt,kxyp,kyzp,
     1jblok,mblok,nzhd)
c this subroutine calculates form factor array need to solve 3d
c poisson's equation in fourier space  with periodic boundary conditions
c for distributed data, with 2D spatial decomposition
c input: isign,ax,ay,az,affp,nx,ny,nz,kstrt,kxyp,kyzp,jblok,mblok,nzhd
c output: ffc
c aimag(ffc(l,j,k,m)) = finite-size particle shape factor s
c real(ffc(l,j,k,m)) = potential green's function g
c for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*(mx - 1) and
c kk = k + kyzp*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c ax/ay/az = half-width of particle in x/y/z direction
c affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c nx/ny/nz = system length in x/y/z direction
c jblok/mblok = number of field partitions in x/y
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c nzhd = first dimension of form factor array, must be >= nzh
      implicit none
      integer nx, ny, nz, kstrt, kxyp, kyzp, jblok, mblok, nzhd
      real ax, ay, az, affp
      complex ffc
      dimension ffc(nzhd,kxyp,kyzp,jblok*mblok)
      integer j, k, l, mx, my, m, joff, koff, moff, nxh, nyh, nzh
      integer kxb, kyzb, js, ks, k1
      real dnx, dny, dnz, dkx, dky, dkz, at1, at2, at3, at4, at5, at6
      complex zero
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      dnz = 6.28318530717959/float(nz)
      zero = cmplx(0.,0.)
      if (kstrt.gt.(kxb*kyzb)) return
c prepare form factor array
      do 50 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 40 mx = 1, jblok
      joff = kxyp*(mx + js) - 1
      m = mx + moff
      do 30 k = 1, kyzp
      k1 = k + koff
      if (k1.gt.nyh) k1 = k1 - ny
      dky = dny*float(k1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 20 j = 1, kxyp
      dkx = dnx*float(j + joff)
      at3 = dkx*dkx + at1
      at4 = (dkx*ax)**2 + at2
      do 10 l = 1, nzh
      dkz = dnz*float(l - 1)
      at5 = dkz*dkz + at3
      at6 = exp(-.5*((dkz*az)**2 + at4))
      if (at5.eq.0.) then
         ffc(l,j,k,m) = cmplx(affp,1.0)
      else
         ffc(l,j,k,m) = cmplx(affp*at6/at5,at6)
      endif
   10 continue
   20 continue
   30 continue
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISP32(q,fx,fy,fz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz
     1,kstrt,nzv,kxyp,kyzp,jblok,mblok,nzhd)
c this subroutine solves 3d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides smoothing a smoothing function,
c with periodic boundary conditions for distributed data,
c with 2D spatial decomposition
c for isign = 0, input: isign,ax,ay,az,affp,nx,ny,nz,kstrt,kxyp,kyzp,
c                       jblok,mblok,nzhd
c output: ffc
c for isign = -1, input: q,ffc,isign,nx,ny,nz,kstrt,nzv,kxyp,kyzp,
c                        jblok,mblok,nzhd
c output: fx,fy,fz,we
c approximate flop count is:
c 62*nxc*nyc*nzc + 33*(nxc*nyc + nxc*nzc + nyc*nzc)
c for isign = 1, input: q,ffc,isign,nx,ny,nz,kstrt,nzv,kxyp,kyzp,
c                        jblok,mblok,nzhd
c output: fx,we
c approximate flop count is:
c 41*nxc*nyc*nzc + 21*(nxc*nyc + nxc*nzc + nyc*nzc)
c for isign = 2, input: q,ffc,isign,nx,ny,nz,kstrt,nzv,kxyp,kyzp,
c                        jblok,mblok,nzhd
c output: fy
c approximate flop count is:
c 8*nxc*nyc*nzc + 4*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
c nvpy/nvpz = number of procs in y/z
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky,kz) = -sqrt(-1)*kx*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c fy(kx,ky,kz) = -sqrt(-1)*ky*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c fz(kx,ky,kz) = -sqrt(-1)*kz*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fz(kx=pi) = 0,
c fx(ky=pi) = fy(ky=pi) = fx(ky=pi) = 0,
c fx(kz=pi) = fy(kz=pi) = fz(kz=pi) = 0,
c fx(kx=0,ky=0,kz=0) = fy(kx=0,ky=0,kz=0) = fz(kx=0,ky=0,kz=0) = 0.
c if isign = 1, potential is calculated using the equation:
c fx(kx,ky,kz) = g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz)
c if isign = 2, smoothing is calculated using the equation:
c fy(kx,ky,kz) = q(kx,ky,kz)*s(kx,ky,kz)
c q(l,j,k,m) = complex charge density for fourier mode jj-1,kk-1,l-1
c fx(l,j,k,m) = x component of force/charge
c fy(l,j,k,m) = y component of force/charge
c fz(l,j,k,m) = z component of force/charge
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*(mx - 1) and
c kk = k + kyzp*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c jblok/mblok = number of field partitions in x/y
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c if isign = 0, form factor array is prepared
c aimag(ffc(l,j,k,m)) = finite-size particle shape factor s
c real(ffc(l,j,k,m)) = potential green's function g
c for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*(mx - 1) and
c kk = k + kyzp*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c ax/ay/az = half-width of particle in x/y/z direction
c affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*
c    |q(kx,ky,kz)*s(kx,ky,kz)|**2)
c nx/ny/nz = system length in x/y/z direction
c nzv = first dimension of field arrays, must be >= nz
c nzhd = first dimension of form factor array, must be >= nzh
      double precision wp
      complex q, fx, fy, fz, ffc, zero
      dimension q(nzv,kxyp,kyzp,jblok*mblok)
      dimension fx(nzv,kxyp,kyzp,jblok*mblok)
      dimension fy(nzv,kxyp,kyzp,jblok*mblok)
      dimension fz(nzv,kxyp,kyzp,jblok*mblok)
      dimension ffc(nzhd,kxyp,kyzp,jblok*mblok)
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2
      nz2 = nz + 2
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      dnz = 6.28318530717959/float(nz)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 60
      if (kstrt.gt.(kxb*kyzb)) return
c prepare form factor array
      do 50 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 40 mx = 1, jblok
      joff = kxyp*(mx + js) - 1
      m = mx + moff
      do 30 k = 1, kyzp
      k1 = k + koff
      if (k1.gt.nyh) k1 = k1 - ny
      dky = dny*float(k1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 20 j = 1, kxyp
      dkx = dnx*float(j + joff)
      at3 = dkx*dkx + at1
      at4 = (dkx*ax)**2 + at2
      do 10 l = 1, nzh
      dkz = dnz*float(l - 1)
      at5 = dkz*dkz + at3
      at6 = exp(-.5*((dkz*az)**2 + at4))
      if (at5.eq.0.) then
         ffc(l,j,k,m) = cmplx(affp,1.)
      else
         ffc(l,j,k,m) = cmplx(affp*at6/at5,at6)
      endif
   10 continue
   20 continue
   30 continue
   40 continue
   50 continue
      return
   60 if (isign.gt.0) go to 210
c calculate force/charge and sum field energy
      wp = 0.0d0
      if (kstrt.gt.(kxb*kyzb)) go to 200
      do 190 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 180 mx = 1, jblok
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      joff = kxyp*(mx + js) - 1
      m = mx + moff
      do 110 k = 1, kyzp
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*float(k1)
         do 80 j = 1, kxyp
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 70 l = 2, nzh
            l1 = nz2 - l
            at1 = real(ffc(l,j,k,m))*aimag(ffc(l,j,k,m))
            at2 = dkx*at1
            at3 = dky*at1
            at4 = dnz*float(l - 1)*at1
            fx(l,j,k,m) = at2*cmplx(aimag(q(l,j,k,m)),-real(q(l,j,k,m)))
            fx(l1,j,k,m) = at2*cmplx(aimag(q(l1,j,k,m)),-real(q(l1,j,k,m
     1)))
            fy(l,j,k,m) = at3*cmplx(aimag(q(l,j,k,m)),-real(q(l,j,k,m)))
            fy(l1,j,k,m) = at3*cmplx(aimag(q(l1,j,k,m)),-real(q(l1,j,k,m
     1)))
            fz(l,j,k,m) = at4*cmplx(aimag(q(l,j,k,m)),-real(q(l,j,k,m)))
            fz(l1,j,k,m) = at4*cmplx(-aimag(q(l1,j,k,m)),real(q(l1,j,k,m
     1)))
            wp = wp + at1*(q(l,j,k,m)*conjg(q(l,j,k,m)) + q(l1,j,k,m)*co
     1njg(q(l1,j,k,m)))
   70       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = real(ffc(1,j,k,m))*aimag(ffc(1,j,k,m))
            at2 = dkx*at1
            at3 = dky*at1
            fx(1,j,k,m) = at2*cmplx(aimag(q(1,j,k,m)),-real(q(1,j,k,m)))
            fx(l1,j,k,m) = zero
            fy(1,j,k,m) = at3*cmplx(aimag(q(1,j,k,m)),-real(q(1,j,k,m)))
            fy(l1,j,k,m) = zero
            fz(1,j,k,m) = zero
            fz(l1,j,k,m) = zero
            wp = wp + at1*(q(1,j,k,m)*conjg(q(1,j,k,m)))
         endif
   80    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 90 l = 2, nzh
               l1 = nz2 - l
               at1 = real(ffc(l,1,k,m))*aimag(ffc(l,1,k,m))
               at3 = dky*at1
               at4 = dnz*float(l - 1)*at1
               fx(l,1,k,m) = zero
               fx(l1,1,k,m) = zero
               fy(l,1,k,m) = at3*cmplx(aimag(q(l,1,k,m)),-real(q(l,1,k,m
     1)))
               fy(l1,1,k,m) = at3*cmplx(aimag(q(l1,1,k,m)),-real(q(l1,1,
     1k,m)))
               fz(l,1,k,m) = at4*cmplx(aimag(q(l,1,k,m)),-real(q(l,1,k,m
     1)))
               fz(l1,1,k,m) = at4*cmplx(-aimag(q(l1,1,k,m)),real(q(l1,1,
     1k,m)))
               wp = wp + at1*(q(l,1,k,m)*conjg(q(l,1,k,m)) + q(l1,1,k,m)
     1*conjg(q(l1,1,k,m)))
   90          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at1 = real(ffc(1,1,k,m))*aimag(ffc(1,1,k,m))
               at3 = dky*at1
               fx(1,1,k,m) = zero
               fx(l1,1,k,m) = zero
               fy(1,1,k,m) = at3*cmplx(aimag(q(1,1,k,m)),-real(q(1,1,k,m
     1)))
               fy(l1,1,k,m) = zero
               fz(1,1,k,m) = zero
               fz(l1,1,k,m) = zero
               wp = wp + at1*(q(1,1,k,m)*conjg(q(1,1,k,m)))
c throw away kx = nx/2
            else
               do 100 l = 1, nz
               fx(l,1,k,m) = zero
               fy(l,1,k,m) = zero
               fz(l,1,k,m) = zero
  100          continue
            endif
         endif
      endif
  110 continue
c mode numbers ky = 0, ny/2
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 130 j = 1, kxyp
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 120 l = 2, nzh
            l1 = nz2 - l
            at1 = real(ffc(l,j,1,m))*aimag(ffc(l,j,1,m))
            at2 = dkx*at1
            at4 = dnz*float(l - 1)*at1
            fx(l,j,1,m) = at2*cmplx(aimag(q(l,j,1,m)),-real(q(l,j,1,m)))
            fx(l1,j,1,m) = at2*cmplx(aimag(q(l1,j,1,m)),-real(q(l1,j,1,m
     1)))
            fy(l,j,1,m) = zero
            fy(l1,j,1,m) = zero
            fz(l,j,1,m) = at4*cmplx(aimag(q(l,j,1,m)),-real(q(l,j,1,m)))
            fz(l1,j,1,m) = at4*cmplx(-aimag(q(l1,j,1,m)),real(q(l1,j,1,m
     1)))
            wp = wp + at1*(q(l,j,1,m)*conjg(q(l,j,1,m)) + q(l1,j,1,m)*co
     1njg(q(l1,j,1,m)))
  120       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = real(ffc(1,j,1,m))*aimag(ffc(1,j,1,m))
            at2 = dkx*at1
            fx(1,j,1,m) = at2*cmplx(aimag(q(1,j,1,m)),-real(q(1,j,1,m)))
            fx(l1,j,1,m) = zero
            fy(1,j,1,m) = zero
            fy(l1,j,1,m) = zero
            fz(1,j,1,m) = zero
            fz(l1,j,1,m) = zero
            wp = wp + at1*(q(1,j,1,m)*conjg(q(1,j,1,m)))
         endif
  130    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 140 l = 2, nzh
            l1 = nz2 - l
            at1 = real(ffc(l,1,1,m))*aimag(ffc(l,1,1,m))
            at4 = dnz*float(l - 1)*at1
            fx(l,1,1,m) = zero
            fx(l1,1,1,m) = zero
            fy(l,1,1,m) = zero
            fy(l1,1,1,m) = zero
            fz(l,1,1,m) = at4*cmplx(aimag(q(l,1,1,m)),-real(q(l,1,1,m)))
            fz(l1,1,1,m) = zero
            wp = wp + at1*(q(l,1,1,m)*conjg(q(l,1,1,m)))
  140       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            fx(1,1,1,m) = zero
            fx(l1,1,1,m) = zero
            fy(1,1,1,m) = zero
            fy(l1,1,1,m) = zero
            fz(1,1,1,m) = zero
            fz(l1,1,1,m) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = (nyh/kyzp)*kyzp
      if (n2.eq.k1) then
         k1 = nyh - k1 + 1
         do 160 j = 1, kxyp
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 150 l = 1, nz
            fx(l,j,k1,m) = zero
            fy(l,j,k1,m) = zero
            fz(l,j,k1,m) = zero
  150       continue
         endif
  160    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 170 l = 1, nz
            fx(l,1,k1,m) = zero
            fy(l,1,k1,m) = zero
            fz(l,1,k1,m) = zero
  170       continue
         endif
      endif
  180 continue
  190 continue
  200 continue
      we = float(nx*ny*nz)*wp
      return
c calculate potential and sum field energy
  210 if (isign.gt.1) go to 360
c calculate force/charge and sum field energy
      wp = 0.0d0
      if (kstrt.gt.(kxb*kyzb)) go to 350
      do 340 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 330 mx = 1, jblok
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      joff = kxyp*(mx + js) - 1
      m = mx + moff
      do 260 k = 1, kyzp
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         do 230 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 220 l = 2, nzh
            l1 = nz2 - l
            at2 = real(ffc(l,j,k,m))
            at1 = at2*aimag(ffc(l,j,k,m))
            fx(l,j,k,m) = at2*q(l,j,k,m)
            fx(l1,j,k,m) = at2*q(l1,j,k,m)
            wp = wp + at1*(q(l,j,k,m)*conjg(q(l,j,k,m)) + q(l1,j,k,m)*co
     1njg(q(l1,j,k,m)))
  220       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at2 = real(ffc(1,j,k,m))
            at1 = at2*aimag(ffc(1,j,k,m))
            fx(1,j,k,m) = at2*q(1,j,k,m)
            fx(l1,j,k,m) = zero
            wp = wp + at1*(q(1,j,k,m)*conjg(q(1,j,k,m)))
         endif
  230    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 240 l = 2, nzh
               l1 = nz2 - l
               at2 = real(ffc(l,1,k,m))
               at1 = at2*aimag(ffc(l,1,k,m))
               fx(l,1,k,m) = at2*q(l,1,k,m)
               fx(l1,1,k,m) = at2*q(l1,1,k,m)
               wp = wp + at1*(q(l,1,k,m)*conjg(q(l,1,k,m)) + q(l1,1,k,m)
     1*conjg(q(l1,1,k,m)))
  240          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at2 = real(ffc(1,1,k,m))
               at1 = at2*aimag(ffc(1,1,k,m))
               fx(1,1,k,m) = at2*q(1,1,k,m)
               fx(l1,1,k,m) = zero
               wp = wp + at1*(q(1,1,k,m)*conjg(q(1,1,k,m)))
c throw away kx = nx/2
            else
               do 250 l = 1, nz
               fx(l,1,k,m) = zero
  250          continue
            endif
         endif
      endif
  260 continue
c mode numbers ky = 0, ny/2
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 280 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 270 l = 2, nzh
            l1 = nz2 - l
            at2 = real(ffc(l,j,1,m))
            at1 = at2*aimag(ffc(l,j,1,m))
            fx(l,j,1,m) = at2*q(l,j,1,m)
            fx(l1,j,1,m) = at2*q(l1,j,1,m)
            wp = wp + at1*(q(l,j,1,m)*conjg(q(l,j,1,m)) + q(l1,j,1,m)*co
     1njg(q(l1,j,1,m)))
  270       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at2 = real(ffc(1,j,1,m))
            at1 = at2*aimag(ffc(1,j,1,m))
            fx(1,j,1,m) = at2*q(1,j,1,m)
            fx(l1,j,1,m) = zero
            wp = wp + at1*(q(1,j,1,m)*conjg(q(1,j,1,m)))
         endif
  280    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 290 l = 2, nzh
            l1 = nz2 - l
            at2 = real(ffc(l,1,1,m))
            at1 = at2*aimag(ffc(l,1,1,m))
            fx(l,1,1,m) = at2*q(l,1,1,m)
            fx(l1,1,1,m) = zero
            wp = wp + at1*(q(l,1,1,m)*conjg(q(l,1,1,m)))
  290       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            fx(1,1,1,m) = zero
            fx(l1,1,1,m) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = (nyh/kyzp)*kyzp
      if (n2.eq.k1) then
         k1 = nyh - k1 + 1
         do 310 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 300 l = 1, nz
            fx(l,j,k1,m) = zero
  300       continue
         endif
  310    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 320 l = 1, nz
            fx(l,1,k1,m) = zero
  320       continue
         endif
      endif
  330 continue
  340 continue
  350 continue
      we = float(nx*ny*nz)*wp
      return
c calculate smoothing
  360 if (kstrt.gt.(kxb*kyzb)) go to 500
      do 490 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 480 mx = 1, jblok
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      joff = kxyp*(mx + js) - 1
      m = mx + moff
      do 410 k = 1, kyzp
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         do 380 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 370 l = 2, nzh
            l1 = nz2 - l
            at1 = aimag(ffc(l,j,k,m))
            fy(l,j,k,m) = at1*q(l,j,k,m)
            fy(l1,j,k,m) = at1*q(l1,j,k,m)
  370       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = aimag(ffc(1,j,k,m))
            fy(1,j,k,m) = at1*q(1,j,k,m)
            fy(l1,j,k,m) = zero
         endif
  380    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 390 l = 2, nzh
               l1 = nz2 - l
               at1 = aimag(ffc(l,1,k,m))
               fy(l,1,k,m) = at1*q(l,1,k,m)
               fy(l1,1,k,m) = at1*q(l1,1,k,m)
  390          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at1 = aimag(ffc(1,1,k,m))
               fy(1,1,k,m) = at1*q(1,1,k,m)
               fy(l1,1,k,m) = zero
c throw away kx = nx/2
            else
               do 400 l = 1, nz
               fy(l,1,k,m) = zero
  400          continue
            endif
         endif
      endif
  410 continue
c mode numbers ky = 0, ny/2
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 430 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 420 l = 2, nzh
            l1 = nz2 - l
            at1 = aimag(ffc(l,j,1,m))
            fy(l,j,1,m) = at1*q(l,j,1,m)
            fy(l1,j,1,m) = at1*q(l1,j,1,m)
  420       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = aimag(ffc(1,j,1,m))
            fy(1,j,1,m) = at1*q(1,j,1,m)
            fy(l1,j,1,m) = zero
         endif
  430    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 440 l = 2, nzh
            l1 = nz2 - l
            at1 = aimag(ffc(l,1,1,m))
            fy(l,1,1,m) = at1*q(l,1,1,m)
            fy(l1,1,1,m) = zero
  440       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = aimag(ffc(1,1,1,m))
            fy(1,1,1,m) = cmplx(at1*real(q(1,1,1,m)),0.)
            fy(l1,1,1,m) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = (nyh/kyzp)*kyzp
      if (n2.eq.k1) then
         k1 = nyh - k1 + 1
         do 460 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 450 l = 1, nz
            fy(l,j,k1,m) = zero
  450       continue
         endif
  460    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 470 l = 1, nz
            fy(l,1,k1,m) = zero
  470       continue
         endif
      endif
  480 continue
  490 continue
  500 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISP332(q,fxyz,isign,ffc,ax,ay,az,affp,we,nx,ny,nz,ks
     1trt,nzv,kxyp,kyzp,jblok,mblok,nzhd)
c this subroutine solves 3d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with periodic boundary conditions for distributed data,
c with 2D spatial decomposition
c for isign = 0, input: isign,ax,ay,az,affp,nx,ny,nz,kstrt,kxyp,kyzp,
c                       jblok,mblok,nzhd
c output: ffc
c for isign =/ 0, input: q,ffc,isign,nx,ny,nz,kstrt,nzv,kxyp,kyzp,
c                        jblok,mblok,nzhd
c output: fxyz,we
c approximate flop count is:
c 62*nxc*nyc*nzc + 33*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
c nvpy/nvpz = number of procs in y/z
c the equation used is:
c fx(kx,ky,kz) = -sqrt(-1)*kx*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c fy(kx,ky,kz) = -sqrt(-1)*ky*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c fz(kx,ky,kz) = -sqrt(-1)*kz*g(kx,ky,kz)*q(kx,ky,kz)*s(kx,ky,kz),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers,
c g(kx,ky,kz) = (affp/(kx**2+ky**2+kz**2))*s(kx,ky,kz),
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
c fx(kx=pi) = fy(kx=pi) = fz(kx=pi) = 0,
c fx(ky=pi) = fy(ky=pi) = fx(ky=pi) = 0,
c fx(kz=pi) = fy(kz=pi) = fz(kz=pi) = 0,
c fx(kx=0,ky=0,kz=0) = fy(kx=0,ky=0,kz=0) = fz(kx=0,ky=0,kz=0) = 0.
c q(l,j,k,m) = complex charge density for fourier mode jj-1,kk-1,l-1
c fxyz(1,l,j,k,m) = x component of force/charge
c fxyz(2,l,j,k,m) = y component of force/charge
c fxyz(3,l,j,k,m) = z component of force/charge
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*(mx - 1) and
c kk = k + kyzp*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c jblok/mblok = number of field partitions in x/y
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c if isign = 0, form factor array is prepared
c aimag(ffc(l,j,k,m)) = finite-size particle shape factor s
c real(ffc(l,j,k,m)) = potential green's function g
c for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*(mx - 1) and
c kk = k + kyzp*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c ax/ay/az = half-width of particle in x/y/z direction
c affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c electric field energy is also calculated, using
c we = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*
c    |q(kx,ky,kz)*s(kx,ky,kz)|**2)
c nx/ny/nz = system length in x/y/z direction
c nzv = first dimension of field arrays, must be >= nz
c nzhd = first dimension of form factor array, must be >= nzh
      double precision wp
      complex q, fxyz, ffc, zero, zt1, zt2
      dimension q(nzv,kxyp,kyzp,jblok*mblok)
      dimension fxyz(3,nzv,kxyp,kyzp,jblok*mblok)
      dimension ffc(nzhd,kxyp,kyzp,jblok*mblok)
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2
      nz2 = nz + 2
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      dnz = 6.28318530717959/float(nz)
      zero = cmplx(0.,0.)
      if (isign.ne.0) go to 60
      if (kstrt.gt.(kxb*kyzb)) return
c prepare form factor array
      do 50 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 40 mx = 1, jblok
      joff = kxyp*(mx + js) - 1
      m = mx + moff
      do 30 k = 1, kyzp
      k1 = k + koff
      if (k1.gt.nyh) k1 = k1 - ny
      dky = dny*float(k1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 20 j = 1, kxyp
      dkx = dnx*float(j + joff)
      at3 = dkx*dkx + at1
      at4 = (dkx*ax)**2 + at2
      do 10 l = 1, nzh
      dkz = dnz*float(l - 1)
      at5 = dkz*dkz + at3
      at6 = exp(-.5*((dkz*az)**2 + at4))
      if (at5.eq.0.) then
         ffc(l,j,k,m) = cmplx(0.,1.)
      else
         ffc(l,j,k,m) = cmplx(affp*at6/at5,at6)
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
      do 190 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 180 mx = 1, jblok
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      joff = kxyp*(mx + js) - 1
      m = mx + moff
      do 110 k = 1, kyzp
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*float(k1)
         do 80 j = 1, kxyp
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 70 l = 2, nzh
            l1 = nz2 - l
            at1 = real(ffc(l,j,k,m))*aimag(ffc(l,j,k,m))
            at2 = dkx*at1
            at3 = dky*at1
            at4 = dnz*float(l - 1)*at1
            zt1 = cmplx(aimag(q(l,j,k,m)),-real(q(l,j,k,m)))
            zt2 = cmplx(aimag(q(l1,j,k,m)),-real(q(l1,j,k,m)))
            fxyz(1,l,j,k,m) = at2*zt1
            fxyz(2,l,j,k,m) = at3*zt1
            fxyz(3,l,j,k,m) = at4*zt1
            fxyz(1,l1,j,k,m) = at2*zt2
            fxyz(2,l1,j,k,m) = at3*zt2
            fxyz(3,l1,j,k,m) = -at4*zt2
            wp = wp + at1*(q(l,j,k,m)*conjg(q(l,j,k,m)) + q(l1,j,k,m)*co
     1njg(q(l1,j,k,m)))
   70       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = real(ffc(1,j,k,m))*aimag(ffc(1,j,k,m))
            at2 = dkx*at1
            at3 = dky*at1
            zt1 = cmplx(aimag(q(1,j,k,m)),-real(q(1,j,k,m)))
            fxyz(1,1,j,k,m) = at2*zt1
            fxyz(2,1,j,k,m) = at3*zt1
            fxyz(3,1,j,k,m) = zero
            fxyz(1,l1,j,k,m) = zero
            fxyz(2,l1,j,k,m) = zero
            fxyz(3,l1,j,k,m) = zero
            wp = wp + at1*(q(1,j,k,m)*conjg(q(1,j,k,m)))
         endif
   80    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 90 l = 2, nzh
               l1 = nz2 - l
               at1 = real(ffc(l,1,k,m))*aimag(ffc(l,1,k,m))
               at3 = dky*at1
               at4 = dnz*float(l - 1)*at1
               zt1 = cmplx(aimag(q(l,1,k,m)),-real(q(l,1,k,m)))
               zt2 = cmplx(aimag(q(l1,1,k,m)),-real(q(l1,1,k,m)))
               fxyz(1,l,1,k,m) = zero
               fxyz(2,l,1,k,m) = at3*zt1
               fxyz(3,l,1,k,m) = at4*zt1
               fxyz(1,l1,1,k,m) = zero
               fxyz(2,l1,1,k,m) = at3*zt2
               fxyz(3,l1,1,k,m) = -at4*zt2
               wp = wp + at1*(q(l,1,k,m)*conjg(q(l,1,k,m)) + q(l1,1,k,m)
     1*conjg(q(l1,1,k,m)))
   90          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at1 = real(ffc(1,1,k,m))*aimag(ffc(1,1,k,m))
               at3 = dky*at1
               zt1 = cmplx(aimag(q(1,1,k,m)),-real(q(1,1,k,m)))
               fxyz(1,1,1,k,m) = zero
               fxyz(2,1,1,k,m) = at3*zt1
               fxyz(3,1,1,k,m) = zero
               fxyz(1,l1,1,k,m) = zero
               fxyz(2,l1,1,k,m) = zero
               fxyz(3,l1,1,k,m) = zero
               wp = wp + at1*(q(1,1,k,m)*conjg(q(1,1,k,m)))
c throw away kx = nx/2
            else
               do 100 l = 1, nz
               fxyz(1,l,1,k,m) = zero
               fxyz(2,l,1,k,m) = zero
               fxyz(3,l,1,k,m) = zero
  100          continue
            endif
         endif
      endif
  110 continue
c mode numbers ky = 0, ny/2
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 130 j = 1, kxyp
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 120 l = 2, nzh
            l1 = nz2 - l
            at1 = real(ffc(l,j,1,m))*aimag(ffc(l,j,1,m))
            at2 = dkx*at1
            at4 = dnz*float(l - 1)*at1
            zt1 = cmplx(aimag(q(l,j,1,m)),-real(q(l,j,1,m)))
            zt2 = cmplx(aimag(q(l1,j,1,m)),-real(q(l1,j,1,m)))
            fxyz(1,l,j,1,m) = at2*zt1
            fxyz(2,l,j,1,m) = zero
            fxyz(3,l,j,1,m) = at4*zt1
            fxyz(1,l1,j,1,m) = at2*zt2
            fxyz(2,l1,j,1,m) = zero
            fxyz(3,l1,j,1,m) = -at4*zt2
            wp = wp + at1*(q(l,j,1,m)*conjg(q(l,j,1,m)) + q(l1,j,1,m)*co
     1njg(q(l1,j,1,m)))
  120       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = real(ffc(1,j,1,m))*aimag(ffc(1,j,1,m))
            at2 = dkx*at1
            zt1 = cmplx(aimag(q(1,j,1,m)),-real(q(1,j,1,m)))
            fxyz(1,1,j,1,m) = at2*zt1
            fxyz(2,1,j,1,m) = zero
            fxyz(3,1,j,1,m) = zero
            fxyz(1,l1,j,1,m) = zero
            fxyz(2,l1,j,1,m) = zero
            fxyz(3,l1,j,1,m) = zero
            wp = wp + at1*(q(1,j,1,m)*conjg(q(1,j,1,m)))
         endif
  130    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 140 l = 2, nzh
            l1 = nz2 - l
            at1 = real(ffc(l,1,1,m))*aimag(ffc(l,1,1,m))
            at4 = dnz*float(l - 1)*at1
            zt1 = cmplx(aimag(q(l,1,1,m)),-real(q(l,1,1,m)))
            fxyz(1,l,1,1,m) = zero
            fxyz(2,l,1,1,m) = zero
            fxyz(3,l,1,1,m) = at4*zt1
            fxyz(1,l1,1,1,m) = zero
            fxyz(2,l1,1,1,m) = zero
            fxyz(3,l1,1,1,m) = zero
            wp = wp + at1*(q(l,1,1,m)*conjg(q(l,1,1,m)))
  140       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            fxyz(1,1,1,1,m) = zero
            fxyz(2,1,1,1,m) = zero
            fxyz(3,1,1,1,m) = zero
            fxyz(1,l1,1,1,m) = zero
            fxyz(2,l1,1,1,m) = zero
            fxyz(3,l1,1,1,m) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = (nyh/kyzp)*kyzp
      if (n2.eq.k1) then
         k1 = nyh - k1 + 1
         do 160 j = 1, kxyp
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 150 l = 1, nz
            fxyz(1,l,j,k1,m) = zero
            fxyz(2,l,j,k1,m) = zero
            fxyz(3,l,j,k1,m) = zero
  150       continue
         endif
  160    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 170 l = 1, nz
            fxyz(1,l,1,k1,m) = zero
            fxyz(2,l,1,k1,m) = zero
            fxyz(3,l,1,k1,m) = zero
  170       continue
         endif
      endif
  180 continue
  190 continue
  200 continue
      we = float(nx*ny*nz)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PDIVF32(f,df,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mblok)
c this subroutine calculates the divergence in fourier space
c for distributed data with 2D spatial decomposition
c input: all except df, output: df
c approximate flop count is:
c 35*nxc*nyc*nzc + 16*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
c nvpy/nvpz = number of procs in y/z
c the divergence is calculated using the equations:
c df(kx,ky,kz) = sqrt(-1)*(kx*fx(kx,ky,kz)+ky*fy(kx,ky,kz)
c                       +kz*fz(kx,ky,kz))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers, except for
c df(kx=pi) = 0, df(ky=pi) = 0, df(kz=pi) = 0
c and df(kx=0,ky=0,kz=0) = 0.
c nx/ny/nz = system length in x/y/z direction
c jblok/mblok = number of field partitions in x/y
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c nzv = first dimension of field arrays, must be >= nz
      complex f, df, zero, zt1
      dimension f(3,nzv,kxyp,kyzp,jblok*mblok)
      dimension df(nzv,kxyp,kyzp,jblok*mblok)
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2
      nz2 = nz + 2
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      dnz = 6.28318530717959/float(nz)
      zero = cmplx(0.,0.)
c calculate the divergence
      if (kstrt.gt.(kxb*kyzb)) return
      do 130 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 120 mx = 1, jblok
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      joff = kxyp*(mx + js) - 1
      m = mx + moff
      do 50 k = 1, kyzp
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*float(k1)
         do 20 j = 1, kxyp
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            zt1 = dkx*f(1,l,j,k,m) + dky*f(2,l,j,k,m) + dkz*f(3,l,j,k,m)
            df(l,j,k,m) = cmplx(-aimag(zt1),real(zt1))
            zt1 = dkx*f(1,l1,j,k,m) + dky*f(2,l1,j,k,m) - dkz*f(3,l1,j,k
     1,m)
            df(l1,j,k,m) = cmplx(-aimag(zt1),real(zt1))
   10       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            zt1 = dkx*f(1,1,j,k,m) + dky*f(2,1,j,k,m)
            df(1,j,k,m) = cmplx(-aimag(zt1),real(zt1))
            df(l1,j,k,m) = zero
         endif
   20    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               dkz = dnz*float(l - 1)
               zt1 = dky*f(2,l,1,k,m) + dkz*f(3,l,1,k,m)
               df(l,1,k,m) = cmplx(-aimag(zt1),real(zt1))
               zt1 = dky*f(2,l1,1,k,m) - dkz*f(3,l1,1,k,m)
               df(l1,1,k,m) = cmplx(-aimag(zt1),real(zt1))
   30          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               zt1 = dky*f(2,1,1,k,m)
               df(1,1,k,m) = cmplx(-aimag(zt1),real(zt1))
               df(l1,1,k,m) = zero
c throw away kx = nx/2
            else
               do 40 l = 1, nz
               df(l,1,k,m) = zero
   40          continue
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny/2
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 70 j = 1, kxyp
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            zt1 = dkx*f(1,l,j,1,m) + dkz*f(3,l,j,1,m)
            df(l,j,1,m) = cmplx(-aimag(zt1),real(zt1))
            zt1 = dkx*f(1,l1,j,1,m) - dkz*f(3,l1,j,1,m)
            df(l1,j,1,m) = cmplx(-aimag(zt1),real(zt1))
   60       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            zt1 = dkx*f(1,1,j,1,m)
            df(1,j,1,m) = cmplx(-aimag(zt1),real(zt1))
            df(l1,j,1,m) = zero
         endif
   70    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            zt1 = dkz*f(3,l,1,1,m)
            df(l,1,1,m) = cmplx(-aimag(zt1),real(zt1))
            df(l1,1,1,m) = zero
   80       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            df(1,1,1,m) = zero
            df(l1,1,1,m) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = (nyh/kyzp)*kyzp
      if (n2.eq.k1) then
         k1 = nyh - k1 + 1
         do 100 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            df(l,j,k1,m) = zero
   90       continue
         endif
  100    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 110 l = 1, nz
            df(l,1,k1,m) = zero
  110       continue
         endif
      endif
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGRADF32(df,f,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mblok)
c this subroutine calculates the gradient in fourier space
c for distributed data with 2D spatial decomposition
c input: all except f, output: f
c approximate flop count is:
c 30*nxc*nyc*nzc + 12*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
c nvpy/nvpz = number of procs in y/z
c the gradient is calculated using the equations:
c fx(kx,ky,kz) = sqrt(-1)*kx*df(kx,ky,kz)
c fy(kx,ky,kz) = sqrt(-1)*ky*df(kx,ky,kz)
c fz(kx,ky,kz) = sqrt(-1)*kz*df(kx,ky,kz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers, except for
c fx(kx=pi) = fy(kx=pi) = fz(kx=pi) = 0,
c fx(ky=pi) = fy(ky=pi) = fx(ky=pi) = 0,
c fx(kz=pi) = fy(kz=pi) = fz(kz=pi) = 0,
c fx(kx=0,ky=0,kz=0) = fy(kx=0,ky=0,kz=0) = fz(kx=0,ky=0,kz=0) = 0.
c nx/ny/nz = system length in x/y/z direction
c jblok/mblok = number of field partitions in x/y
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c nzv = first dimension of field arrays, must be >= nz
      complex df, f, zero, zt1
      dimension df(nzv,kxyp,kyzp,jblok*mblok)
      dimension f(3,nzv,kxyp,kyzp,jblok*mblok)
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2
      nz2 = nz + 2
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      dnz = 6.28318530717959/float(nz)
      zero = cmplx(0.,0.)
c calculate the gradient
      if (kstrt.gt.(kxb*kyzb)) return
      do 130 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 120 mx = 1, jblok
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      joff = kxyp*(mx + js) - 1
      m = mx + moff
      do 50 k = 1, kyzp
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*float(k1)
         do 20 j = 1, kxyp
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            zt1 = cmplx(-aimag(df(l,j,k,m)),real(df(l,j,k,m)))
            f(1,l,j,k,m) = dkx*zt1
            f(2,l,j,k,m) = dky*zt1
            f(3,l,j,k,m) = dkz*zt1
            zt1 = cmplx(-aimag(df(l1,j,k,m)),real(df(l1,j,k,m)))
            f(1,l1,j,k,m) = dkx*zt1
            f(2,l1,j,k,m) = dky*zt1
            f(3,l1,j,k,m) = -dkz*zt1
   10       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            zt1 = cmplx(-aimag(df(1,j,k,m)),real(df(1,j,k,m)))
            f(1,1,j,k,m) = dkx*zt1
            f(2,1,j,k,m) = dky*zt1
            f(3,1,j,k,m) = zero
            f(1,l1,j,k,m) = zero
            f(2,l1,j,k,m) = zero
            f(3,l1,j,k,m) = zero
         endif
   20    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               dkz = dnz*float(l - 1)
               zt1 = cmplx(-aimag(df(l,1,k,m)),real(df(l,1,k,m)))
               f(1,l,1,k,m) = zero
               f(2,l,1,k,m) = dky*zt1
               f(3,l,1,k,m) = dkz*zt1
               zt1 = cmplx(-aimag(df(l1,1,k,m)),real(df(l1,1,k,m)))
               f(1,l1,1,k,m) = zero
               f(2,l1,1,k,m) = dky*zt1
               f(3,l1,1,k,m) = -dkz*zt1
   30          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               zt1 = cmplx(-aimag(df(1,1,k,m)),real(df(1,1,k,m)))
               f(1,1,1,k,m) = zero
               f(2,1,1,k,m) = dky*zt1
               f(3,1,1,k,m) = zero
               f(1,l1,1,k,m) = zero
               f(2,l1,1,k,m) = zero
               f(3,l1,1,k,m) = zero
c throw away kx = nx/2
            else
               do 40 l = 1, nz
               f(1,l,1,k,m) = zero
               f(2,l,1,k,m) = zero
               f(3,l,1,k,m) = zero
   40          continue
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny/2
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 70 j = 1, kxyp
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            zt1 = cmplx(-aimag(df(l,j,1,m)),real(df(l,j,1,m)))
            f(1,l,j,1,m) = dkx*zt1
            f(2,l,j,1,m) = zero
            f(3,l,j,1,m) = dkz*zt1
            zt1 = cmplx(-aimag(df(l1,j,1,m)),real(df(l1,j,1,m)))
            f(1,l1,j,1,m) = dkx*zt1
            f(2,l1,j,1,m) = zero
            f(3,l1,j,1,m) = -dkz*zt1
   60       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            zt1 = cmplx(-aimag(df(1,j,1,m)),real(df(1,j,1,m)))
            f(1,1,j,1,m) = dkx*zt1
            f(2,1,j,1,m) = zero
            f(3,1,j,1,m) = zero
            f(1,l1,j,1,m) = zero
            f(2,l1,j,1,m) = zero
            f(3,l1,j,1,m) = zero
         endif
   70    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            zt1 = cmplx(-aimag(df(l,1,1,m)),real(df(l,1,1,m)))
            f(1,l,1,1,m) = zero
            f(2,l,1,1,m) = zero
            f(3,l,1,1,m) = dkz*zt1
            f(1,l1,1,1,m) = zero
            f(2,l1,1,1,m) = zero
            f(3,l1,1,1,m) = zero
   80       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            f(1,1,1,1,m) = zero
            f(2,1,1,1,m) = zero
            f(3,1,1,1,m) = zero
            f(1,l1,1,1,m) = zero
            f(2,l1,1,1,m) = zero
            f(3,l1,1,1,m) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = (nyh/kyzp)*kyzp
      if (n2.eq.k1) then
         k1 = nyh - k1 + 1
         do 100 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            f(1,l,j,k1,m) = zero
            f(2,l,j,k1,m) = zero
            f(3,l,j,k1,m) = zero
   90       continue
         endif
  100    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 110 l = 1, nz
            f(1,l,1,k1,m) = zero
            f(2,l,1,k1,m) = zero
            f(3,l,1,k1,m) = zero
  110       continue
         endif
      endif
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCURLF32(f,g,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mblok)
c this subroutine calculates the curl in fourier space
c for distributed data with 2D spatial decomposition
c input: all except g, output: g
c approximate flop count is:
c 86*nxc*nyc*nzc + 32*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
c nvpy/nvpz = number of procs in y/z
c the curl is calculated using the equations:
c gx(kx,ky,kz) = sqrt(-1)*(ky*fz(kx,ky,kz)-kz*fy(kx,ky,kz))
c gy(kx,ky,kz) = sqrt(-1)*(kz*fx(kx,ky,kz)-kx*fz(kx,ky,kz))
c gz(kx,ky,kz) = sqrt(-1)*(kx*fy(kx,ky,kz)-ky*fx(kx,ky,kz))
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers, except for
c gx(kx=pi) = gy(kx=pi) = gz(kx=pi) = 0,
c gx(ky=pi) = gy(ky=pi) = gx(ky=pi) = 0,
c gx(kz=pi) = gy(kz=pi) = gz(kz=pi) = 0,
c gx(kx=0,ky=0,kz=0) = gy(kx=0,ky=0,kz=0) = gz(kx=0,ky=0,kz=0) = 0.
c nx/ny/nz = system length in x/y/z direction
c jblok/mblok = number of field partitions in x/y
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c nzv = first dimension of field arrays, must be >= nz
      complex f, g, zero, zt1, zt2, zt3
      dimension f(3,nzv,kxyp,kyzp,jblok*mblok)
      dimension g(3,nzv,kxyp,kyzp,jblok*mblok)
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2
      nz2 = nz + 2
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      dnz = 6.28318530717959/float(nz)
      zero = cmplx(0.,0.)
c calculate the curl
      if (kstrt.gt.(kxb*kyzb)) return
      do 130 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 120 mx = 1, jblok
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      joff = kxyp*(mx + js) - 1
      m = mx + moff
      do 50 k = 1, kyzp
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*float(k1)
         do 20 j = 1, kxyp
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            zt1 = cmplx(-aimag(f(3,l,j,k,m)),real(f(3,l,j,k,m)))
            zt2 = cmplx(-aimag(f(2,l,j,k,m)),real(f(2,l,j,k,m)))
            zt3 = cmplx(-aimag(f(1,l,j,k,m)),real(f(1,l,j,k,m)))
            g(1,l,j,k,m) = dky*zt1 - dkz*zt2
            g(2,l,j,k,m) = dkz*zt3 - dkx*zt1
            g(3,l,j,k,m) = dkx*zt2 - dky*zt3
            zt1 = cmplx(-aimag(f(3,l1,j,k,m)),real(f(3,l1,j,k,m)))
            zt2 = cmplx(-aimag(f(2,l1,j,k,m)),real(f(2,l1,j,k,m)))
            zt3 = cmplx(-aimag(f(1,l1,j,k,m)),real(f(1,l1,j,k,m)))
            g(1,l1,j,k,m) = dky*zt1 + dkz*zt2
            g(2,l1,j,k,m) = -dkz*zt3 - dkx*zt1
            g(3,l1,j,k,m) = dkx*zt2 - dky*zt3
   10       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            zt1 = cmplx(-aimag(f(3,1,j,k,m)),real(f(3,1,j,k,m)))
            zt2 = cmplx(-aimag(f(2,1,j,k,m)),real(f(2,1,j,k,m)))
            zt3 = cmplx(-aimag(f(1,1,j,k,m)),real(f(1,1,j,k,m)))
            g(1,1,j,k,m) = dky*zt1
            g(2,1,j,k,m) = -dkx*zt1
            g(3,1,j,k,m) = dkx*zt2 - dky*zt3
            g(1,l1,j,k,m) = zero
            g(2,l1,j,k,m) = zero
            g(3,l1,j,k,m) = zero
         endif
   20    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               dkz = dnz*float(l - 1)
               zt1 = cmplx(-aimag(f(3,l,1,k,m)),real(f(3,l,1,k,m)))
               zt2 = cmplx(-aimag(f(2,l,1,k,m)),real(f(2,l,1,k,m)))
               zt3 = cmplx(-aimag(f(1,l,1,k,m)),real(f(1,l,1,k,m)))
               g(1,l,1,k,m) = dky*zt1 - dkz*zt2
               g(2,l,1,k,m) = dkz*zt3
               g(3,l,1,k,m) = -dky*zt3
               zt1 = cmplx(-aimag(f(3,l1,1,k,m)),real(f(3,l1,1,k,m)))
               zt2 = cmplx(-aimag(f(2,l1,1,k,m)),real(f(2,l1,1,k,m)))
               zt3 = cmplx(-aimag(f(1,l1,1,k,m)),real(f(1,l1,1,k,m)))
               g(1,l1,1,k,m) = dky*zt1 + dkz*zt2
               g(2,l1,1,k,m) = -dkz*zt3
               g(3,l1,1,k,m) = -dky*zt3
   30          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               zt1 = cmplx(-aimag(f(3,1,1,k,m)),real(f(3,1,1,k,m)))
               zt3 = cmplx(-aimag(f(1,1,1,k,m)),real(f(1,1,1,k,m)))
               g(1,1,1,k,m) = dky*zt1
               g(2,1,1,k,m) = zero
               g(3,1,1,k,m) = -dky*zt3
               g(1,l1,1,k,m) = zero
               g(2,l1,1,k,m) = zero
               g(3,l1,1,k,m) = zero
c throw away kx = nx/2
            else
               do 40 l = 1, nz
               g(1,l,1,k,m) = zero
               g(2,l,1,k,m) = zero
               g(3,l,1,k,m) = zero
   40          continue
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny/2
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 70 j = 1, kxyp
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            zt1 = cmplx(-aimag(f(3,l,j,1,m)),real(f(3,l,j,1,m)))
            zt2 = cmplx(-aimag(f(2,l,j,1,m)),real(f(2,l,j,1,m)))
            zt3 = cmplx(-aimag(f(1,l,j,1,m)),real(f(1,l,j,1,m)))
            g(1,l,j,1,m) = -dkz*zt2
            g(2,l,j,1,m) = dkz*zt3 - dkx*zt1
            g(3,l,j,1,m) = dkx*zt2
            zt1 = cmplx(-aimag(f(3,l1,j,1,m)),real(f(3,l1,j,1,m)))
            zt2 = cmplx(-aimag(f(2,l1,j,1,m)),real(f(2,l1,j,1,m)))
            zt3 = cmplx(-aimag(f(1,l1,j,1,m)),real(f(1,l1,j,1,m)))
            g(1,l1,j,1,m) = dkz*zt2
            g(2,l1,j,1,m) = -dkz*zt3 - dkx*zt1
            g(3,l1,j,1,m) = dkx*zt2
   60       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            zt1 = cmplx(-aimag(f(3,1,j,1,m)),real(f(3,1,j,1,m)))
            zt2 = cmplx(-aimag(f(2,1,j,1,m)),real(f(2,1,j,1,m)))
            g(1,1,j,1,m) = zero
            g(2,1,j,1,m) = -dkx*zt1
            g(3,1,j,1,m) = dkx*zt2
            g(1,l1,j,1,m) = zero
            g(2,l1,j,1,m) = zero
            g(3,l1,j,1,m) = zero
         endif
   70    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            zt2 = cmplx(-aimag(f(2,l,1,1,m)),real(f(2,l,1,1,m)))
            zt3 = cmplx(-aimag(f(1,l,1,1,m)),real(f(1,l,1,1,m)))
            g(1,l,1,1,m) = -dkz*zt2
            g(2,l,1,1,m) = dkz*zt3
            g(3,l,1,1,m) = zero
            g(1,l1,1,1,m) = zero
            g(2,l1,1,1,m) = zero
            g(3,l1,1,1,m) = zero
   80       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            g(1,1,1,1,m) = zero
            g(2,1,1,1,m) = zero
            g(3,1,1,1,m) = zero
            g(1,l1,1,1,m) = zero
            g(2,l1,1,1,m) = zero
            g(3,l1,1,1,m) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = (nyh/kyzp)*kyzp
      if (n2.eq.k1) then
         k1 = nyh - k1 + 1
         do 100 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            g(1,l,j,k1,m) = zero
            g(2,l,j,k1,m) = zero
            g(3,l,j,k1,m) = zero
   90       continue
         endif
  100    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 110 l = 1, nz
            g(1,l,1,k1,m) = zero
            g(2,l,1,k1,m) = zero
            g(3,l,1,k1,m) = zero
  110       continue
         endif
      endif
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PCUPERP32(cu,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mblok)
c this subroutine calculates the transverse current in fourier space
c for distributed data with 2D spatial decomposition
c input: all output: cu
c approximate flop count is:
c 100*nxc*nyc*nzc + 36*(nxc*nyc + nxc*nzc + nyc*nzc)
c and (nx/2)*nyc*nzc divides
c where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
c nvpy/nvpz = number of procs in y/z
c the transverse current is calculated using the equation:
c cux(kx,ky,kz) = cux(kx,ky,kz) - kx*(kx*cux(kx,ky,kz)+ky*cuy(kx,ky,kz)+
c                                 kz*cuz(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c cuy(kx,ky,kz) = cuy(kx,ky,kz) - ky*(kx*cux(kx,ky,kz)+ky*cuy(kx,ky,kz)+
c                                 kz*cuz(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c cuz(kx,ky,kz) = cuz(kx,ky,kz) - kz*(kx*cux(kx,ky,kz)+ky*cuy(kx,ky,kz)+
c                                 kz*cuz(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers, except for
c cux(kx=pi) = cuy(kx=pi) = cuz(kx=pi) = 0,
c cux(ky=pi) = cuy(ky=pi) = cux(ky=pi) = 0,
c cux(kz=pi) = cuy(kz=pi) = cuz(kz=pi) = 0,
c cux(kx=0,ky=0,kz=0) = cuy(kx=0,ky=0,kz=0) = cuz(kx=0,ky=0,kz=0) = 0.
c nx/ny/nz = system length in x/y/z direction
c jblok/mblok = number of field partitions in x/y
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c nzv = first dimension of field arrays, must be >= nz
      complex cu, zero, zt1
      dimension cu(3,nzv,kxyp,kyzp,jblok*mblok)
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2
      nz2 = nz + 2
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      dnz = 6.28318530717959/float(nz)
      zero = cmplx(0.,0.)
c calculate transverse part of current
      if (kstrt.gt.(kxb*kyzb)) return
      do 130 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 120 mx = 1, jblok
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      joff = kxyp*(mx + js) - 1
      m = mx + moff
      do 50 k = 1, kyzp
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*float(k1)
         dky2 = dky*dky
         do 20 j = 1, kxyp
         dkx = dnx*float(j + joff)
         dkxy2 = dkx*dkx + dky2
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            at1 = 1./(dkz*dkz + dkxy2)
            zt1 = at1*(dkx*cu(1,l,j,k,m) + dky*cu(2,l,j,k,m) + dkz*cu(3,
     1l,j,k,m))
            cu(1,l,j,k,m) = cu(1,l,j,k,m) - dkx*zt1
            cu(2,l,j,k,m) = cu(2,l,j,k,m) - dky*zt1
            cu(3,l,j,k,m) = cu(3,l,j,k,m) - dkz*zt1
            zt1 = at1*(dkx*cu(1,l1,j,k,m) + dky*cu(2,l1,j,k,m) - dkz*cu(
     13,l1,j,k,m))
            cu(1,l1,j,k,m) = cu(1,l1,j,k,m) - dkx*zt1
            cu(2,l1,j,k,m) = cu(2,l1,j,k,m) - dky*zt1
            cu(3,l1,j,k,m) = cu(3,l1,j,k,m) + dkz*zt1
   10       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = 1./dkxy2
            zt1 = at1*(dkx*cu(1,1,j,k,m) + dky*cu(2,1,j,k,m))
            cu(1,1,j,k,m) = cu(1,1,j,k,m) - dkx*zt1
            cu(2,1,j,k,m) = cu(2,1,j,k,m) - dky*zt1
            cu(1,l1,j,k,m) = zero
            cu(2,l1,j,k,m) = zero
            cu(3,l1,j,k,m) = zero
         endif
   20    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               dkz = dnz*float(l - 1)
               at1 = 1./(dkz*dkz + dky2)
               zt1 = at1*(dky*cu(2,l,1,k,m) + dkz*cu(3,l,1,k,m))
               cu(2,l,1,k,m) = cu(2,l,1,k,m) - dky*zt1
               cu(3,l,1,k,m) = cu(3,l,1,k,m) - dkz*zt1
               zt1 = at1*(dky*cu(2,l1,1,k,m) - dkz*cu(3,l1,1,k,m))
               cu(2,l1,1,k,m) = cu(2,l1,1,k,m) - dky*zt1
               cu(3,l1,1,k,m) = cu(3,l1,1,k,m) + dkz*zt1
   30          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               cu(2,1,1,k,m) = zero
               cu(1,l1,1,k,m) = zero
               cu(2,l1,1,k,m) = zero
               cu(3,l1,1,k,m) = zero
c throw away kx = nx/2
            else
               do 40 l = 1, nz
               cu(1,l,1,k,m) = zero
               cu(2,l,1,k,m) = zero
               cu(3,l,1,k,m) = zero
   40          continue
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny/2
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 70 j = 1, kxyp
         dkx = dnx*float(j + joff)
         dkx2 = dkx*dkx
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            at1 = 1./(dkz*dkz + dkx2)
            zt1 = at1*(dkx*cu(1,l,j,1,m) + dkz*cu(3,l,j,1,m))
            cu(1,l,j,1,m) = cu(1,l,j,1,m) - dkx*zt1
            cu(3,l,j,1,m) = cu(3,l,j,1,m) - dkz*zt1
            zt1 = at1*(dkx*cu(1,l1,j,1,m) - dkz*cu(3,l1,j,1,m))
            cu(1,l1,j,1,m) = cu(1,l1,j,1,m) - dkx*zt1
            cu(3,l1,j,1,m) = cu(3,l1,j,1,m) + dkz*zt1
   60       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            cu(1,1,j,1,m) = zero
            cu(1,l1,j,1,m) = zero
            cu(2,l1,j,1,m) = zero
            cu(3,l1,j,1,m) = zero
         endif
   70    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            cu(3,l,1,1,m) = zero
            cu(1,l1,1,1,m) = zero
            cu(2,l1,1,1,m) = zero
            cu(3,l1,1,1,m) = zero
   80       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            cu(1,1,1,1,m) = zero
            cu(2,1,1,1,m) = zero
            cu(3,1,1,1,m) = zero
            cu(1,l1,1,1,m) = zero
            cu(2,l1,1,1,m) = zero
            cu(3,l1,1,1,m) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = (nyh/kyzp)*kyzp
      if (n2.eq.k1) then
         k1 = nyh - k1 + 1
         do 100 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            cu(1,l,j,k1,m) = zero
            cu(2,l,j,k1,m) = zero
            cu(3,l,j,k1,m) = zero
   90       continue
         endif
  100    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 110 l = 1, nz
            cu(1,l,1,k1,m) = zero
            cu(2,l,1,k1,m) = zero
            cu(3,l,1,k1,m) = zero
  110       continue
         endif
      endif
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PBPOISP332(cu,bxyz,isign,ffc,ax,ay,az,affp,ci,wm,nx,ny,
     1nz,kstrt,nzv,kxyp,kyzp,jblok,mblok,nzhd)
c this subroutine solves 3d poisson's equation in fourier space for
c magnetic field (or convolution of magnetic field over particle shape)
c or for vector potential, or provides a smoothing function,
c with periodic boundary conditions, for distributed data,
c with 2D spatial decomposition
c for isign = 0, output: ffc
c input: isign,ax,ay,az,affp,nx,ny,nz,kstrt,kxyp,kyzp,jblok,mblok,nzhd
c for isign = -1, output: bxyz,wm
c input: cu,ffc,isign,ci,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mblok,nzhd
c approximate flop count is:
c 193*nxc*nyc*nzc + 84*(nxc*nyc + nxc*nzc + nyc*nzc)
c for isign = 1, output: bxyz,wm
c input: cu,ffc,isign,ci,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mblok,nzhd
c approximate flop count is:
c 128*nxc*nyc*nzc + 66*(nxc*nyc + nxc*nzc + nyc*nzc)
c for isign = 2, output: bxyz
c input: cu,ffc,isign,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mblok,nzhd
c approximate flop count is:
c 24*nxc*nyc*nzc + 12*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
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
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
c bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
c bx(ky=pi) = by(ky=pi) = bx(ky=pi) = 0,
c bx(kz=pi) = by(kz=pi) = bz(kz=pi) = 0,
c bx(kx=0,ky=0,kz=0) = by(kx=0,ky=0,kz=0) = bz(kx=0,ky=0,kz=0) = 0.
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
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*(mx - 1) and
c kk = k + kyzp*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c aimag(ffc(l,j,k,m)) = finite-size particle shape factor s
c real(ffc(l,j,k,m)) = potential green's function g
c ax/ay/az = half-width of particle in x/y/z direction
c affp = normalization constant = nx*ny*nz/np,
c where np=number of particles
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny/nz = system length in x/y/z direction
c jblok/mblok = number of field partitions in x/y
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c nzv = first dimension of field arrays, must be >= nz
c nzhd = first dimension of form factor array, must be >= nzh
      double precision wp
      complex cu, bxyz, ffc, zero, zt1, zt2, zt3
      dimension cu(3,nzv,kxyp,kyzp,jblok*mblok)
      dimension bxyz(3,nzv,kxyp,kyzp,jblok*mblok)
      dimension ffc(nzhd,kxyp,kyzp,jblok*mblok)
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2
      nz2 = nz + 2
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      dnz = 6.28318530717959/float(nz)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
      if (isign.ne.0) go to 60
      if (kstrt.gt.(kxb*kyzb)) return
c prepare form factor array
      do 50 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 40 mx = 1, jblok
      joff = kxyp*(mx + js) - 1
      m = mx + moff
      do 30 k = 1, kyzp
      k1 = k + koff
      if (k1.gt.nyh) k1 = k1 - ny
      dky = dny*float(k1)
      at1 = dky*dky
      at2 = (dky*ay)**2
      do 20 j = 1, kxyp
      dkx = dnx*float(j + joff)
      at3 = dkx*dkx + at1
      at4 = (dkx*ax)**2 + at2
      do 10 l = 1, nzh
      dkz = dnz*float(l - 1)
      at5 = dkz*dkz + at3
      at6 = exp(-.5*((dkz*az)**2 + at4))
      if (at5.eq.0.) then
         ffc(l,j,k,m) = cmplx(affp,1.)
      else
         ffc(l,j,k,m) = cmplx(affp*at6/at5,at6)
      endif
   10 continue
   20 continue
   30 continue
   40 continue
   50 continue
      return
   60 if (isign.gt.0) go to 210
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.(kxb*kyzb)) go to 200
      do 190 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 180 mx = 1, jblok
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      joff = kxyp*(mx + js) - 1
      m = mx + moff
      do 110 k = 1, kyzp
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*float(k1)
         do 80 j = 1, kxyp
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 70 l = 2, nzh
            l1 = nz2 - l
            at1 = ci2*real(ffc(l,j,k,m))*aimag(ffc(l,j,k,m))
            at2 = dkx*at1
            at3 = dky*at1
            at4 = dnz*float(l - 1)*at1
            zt1 = cmplx(-aimag(cu(3,l,j,k,m)),real(cu(3,l,j,k,m)))
            zt2 = cmplx(-aimag(cu(2,l,j,k,m)),real(cu(2,l,j,k,m)))
            zt3 = cmplx(-aimag(cu(1,l,j,k,m)),real(cu(1,l,j,k,m)))
            bxyz(1,l,j,k,m) = at3*zt1 - at4*zt2
            bxyz(2,l,j,k,m) = at4*zt3 - at2*zt1
            bxyz(3,l,j,k,m) = at2*zt2 - at3*zt3
            zt1 = cmplx(-aimag(cu(3,l1,j,k,m)),real(cu(3,l1,j,k,m)))
            zt2 = cmplx(-aimag(cu(2,l1,j,k,m)),real(cu(2,l1,j,k,m)))
            zt3 = cmplx(-aimag(cu(1,l1,j,k,m)),real(cu(1,l1,j,k,m)))
            bxyz(1,l1,j,k,m) = at3*zt1 + at4*zt2
            bxyz(2,l1,j,k,m) = -at4*zt3 - at2*zt1
            bxyz(3,l1,j,k,m) = at2*zt2 - at3*zt3
            wp = wp + at1*(cu(1,l,j,k,m)*conjg(cu(1,l,j,k,m)) + cu(2,l,j
     1,k,m)*conjg(cu(2,l,j,k,m)) + cu(3,l,j,k,m)*conjg(cu(3,l,j,k,m)) + 
     2cu(1,l1,j,k,m)*conjg(cu(1,l1,j,k,m)) + cu(2,l1,j,k,m)*conjg(cu(2,l
     31,j,k,m)) + cu(3,l1,j,k,m)*conjg(cu(3,l1,j,k,m)))
   70       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = ci2*real(ffc(1,j,k,m))*aimag(ffc(1,j,k,m))
            at2 = dkx*at1
            at3 = dky*at1
            zt1 = cmplx(-aimag(cu(3,1,j,k,m)),real(cu(3,1,j,k,m)))
            zt2 = cmplx(-aimag(cu(2,1,j,k,m)),real(cu(2,1,j,k,m)))
            zt3 = cmplx(-aimag(cu(1,1,j,k,m)),real(cu(1,1,j,k,m)))
            bxyz(1,1,j,k,m) = at3*zt1
            bxyz(2,1,j,k,m) = -at2*zt1
            bxyz(3,1,j,k,m) = at2*zt2 - at3*zt3
            bxyz(1,l1,j,k,m) = zero
            bxyz(2,l1,j,k,m) = zero
            bxyz(3,l1,j,k,m) = zero
            wp = wp + at1*(cu(1,1,j,k,m)*conjg(cu(1,1,j,k,m)) + cu(2,1,j
     1,k,m)*conjg(cu(2,1,j,k,m)) + cu(3,1,j,k,m)*conjg(cu(3,1,j,k,m)))
         endif
   80    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 90 l = 2, nzh
               l1 = nz2 - l
               at1 = ci2*real(ffc(l,1,k,m))*aimag(ffc(l,1,k,m))
               at3 = dky*at1
               at4 = dnz*float(l - 1)*at1
               zt1 = cmplx(-aimag(cu(3,l,1,k,m)),real(cu(3,l,1,k,m)))
               zt2 = cmplx(-aimag(cu(2,l,1,k,m)),real(cu(2,l,1,k,m)))
               zt3 = cmplx(-aimag(cu(1,l,1,k,m)),real(cu(1,l,1,k,m)))
               bxyz(1,l,1,k,m) = at3*zt1 - at4*zt2
               bxyz(2,l,1,k,m) = at4*zt3
               bxyz(3,l,1,k,m) = -at3*zt3
               zt1 = cmplx(-aimag(cu(3,l1,1,k,m)),real(cu(3,l1,1,k,m)))
               zt2 = cmplx(-aimag(cu(2,l1,1,k,m)),real(cu(2,l1,1,k,m)))
               zt3 = cmplx(-aimag(cu(1,l1,1,k,m)),real(cu(1,l1,1,k,m)))
               bxyz(1,l1,1,k,m) = at3*zt1 + at4*zt2
               bxyz(2,l1,1,k,m) = -at4*zt3
               bxyz(3,l1,1,k,m) = -at3*zt3
               wp = wp + at1*(cu(1,l,1,k,m)*conjg(cu(1,l,1,k,m)) + cu(2,
     1l,1,k,m)*conjg(cu(2,l,1,k,m)) + cu(3,l,1,k,m)*conjg(cu(3,l,1,k,m))
     2 + cu(1,l1,1,k,m)*conjg(cu(1,l1,1,k,m)) + cu(2,l1,1,k,m)*conjg(cu(
     32,l1,1,k,m)) + cu(3,l1,1,k,m)*conjg(cu(3,l1,1,k,m)))
   90          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at1 = ci2*real(ffc(1,1,k,m))*aimag(ffc(1,1,k,m))
               at3 = dky*at1
               zt1 = cmplx(-aimag(cu(3,1,1,k,m)),real(cu(3,1,1,k,m)))
               zt3 = cmplx(-aimag(cu(1,1,1,k,m)),real(cu(1,1,1,k,m)))
               bxyz(1,1,1,k,m) = at3*zt1
               bxyz(2,1,1,k,m) = zero
               bxyz(3,1,1,k,m) = -at3*zt3
               bxyz(1,l1,1,k,m) = zero
               bxyz(2,l1,1,k,m) = zero
               bxyz(3,l1,1,k,m) = zero
               wp = wp + at1*(cu(1,1,1,k,m)*conjg(cu(1,1,1,k,m)) + cu(2,
     11,1,k,m)*conjg(cu(2,1,1,k,m)) + cu(3,1,1,k,m)*conjg(cu(3,1,1,k,m))
     2)
c throw away kx = nx/2
            else
               do 100 l = 1, nz
               bxyz(1,l,1,k,m) = zero
               bxyz(2,l,1,k,m) = zero
               bxyz(3,l,1,k,m) = zero
  100          continue
            endif
         endif
      endif
  110 continue
c mode numbers ky = 0, ny/2
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 130 j = 1, kxyp
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 120 l = 2, nzh
            l1 = nz2 - l
            at1 = ci2*real(ffc(l,j,1,m))*aimag(ffc(l,j,1,m))
            at2 = dkx*at1
            at4 = dnz*float(l - 1)*at1
            zt1 = cmplx(-aimag(cu(3,l,j,1,m)),real(cu(3,l,j,1,m)))
            zt2 = cmplx(-aimag(cu(2,l,j,1,m)),real(cu(2,l,j,1,m)))
            zt3 = cmplx(-aimag(cu(1,l,j,1,m)),real(cu(1,l,j,1,m)))
            bxyz(1,l,j,1,m) = -at4*zt2
            bxyz(2,l,j,1,m) = at4*zt3 - at2*zt1
            bxyz(3,l,j,1,m) = at2*zt2
            zt1 = cmplx(-aimag(cu(3,l1,j,1,m)),real(cu(3,l1,j,1,m)))
            zt2 = cmplx(-aimag(cu(2,l1,j,1,m)),real(cu(2,l1,j,1,m)))
            zt3 = cmplx(-aimag(cu(1,l1,j,1,m)),real(cu(1,l1,j,1,m)))
            bxyz(1,l1,j,1,m) = at4*zt2
            bxyz(2,l1,j,1,m) = -at4*zt3 - at2*zt1
            bxyz(3,l1,j,1,m) = at2*zt2
            wp = wp + at1*(cu(1,l,j,1,m)*conjg(cu(1,l,j,1,m)) + cu(2,l,j
     1,1,m)*conjg(cu(2,l,j,1,m)) + cu(3,l,j,1,m)*conjg(cu(3,l,j,1,m)) + 
     2cu(1,l1,j,1,m)*conjg(cu(1,l1,j,1,m)) + cu(2,l1,j,1,m)*conjg(cu(2,l
     31,j,1,m)) + cu(3,l1,j,1,m)*conjg(cu(3,l1,j,1,m)))
  120       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = ci2*real(ffc(1,j,1,m))*aimag(ffc(1,j,1,m))
            at2 = dkx*at1
            zt1 = cmplx(-aimag(cu(3,1,j,1,m)),real(cu(3,1,j,1,m)))
            zt2 = cmplx(-aimag(cu(2,1,j,1,m)),real(cu(2,1,j,1,m)))
            bxyz(1,1,j,1,m) = zero
            bxyz(2,1,j,1,m) = -at2*zt1
            bxyz(3,1,j,1,m) = at2*zt2
            bxyz(1,l1,j,1,m) = zero
            bxyz(2,l1,j,1,m) = zero
            bxyz(3,l1,j,1,m) = zero
            wp = wp + at1*(cu(1,1,j,1,m)*conjg(cu(1,1,j,1,m)) + cu(2,1,j
     1,1,m)*conjg(cu(2,1,j,1,m)) + cu(3,1,j,1,m)*conjg(cu(3,1,j,1,m)))
         endif
  130    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 140 l = 2, nzh
            l1 = nz2 - l
            at1 = ci2*real(ffc(l,1,1,m))*aimag(ffc(l,1,1,m))
            at4 = dnz*float(l - 1)*at1
            zt2 = cmplx(-aimag(cu(2,l,1,1,m)),real(cu(2,l,1,1,m)))
            zt3 = cmplx(-aimag(cu(1,l,1,1,m)),real(cu(1,l,1,1,m)))
            bxyz(1,l,1,1,m) = -at4*zt2
            bxyz(2,l,1,1,m) = at4*zt3
            bxyz(3,l,1,1,m) = zero
            bxyz(1,l1,1,1,m) = zero
            bxyz(2,l1,1,1,m) = zero
            bxyz(3,l1,1,1,m) = zero
            wp = wp + at1*(cu(1,l,1,1,m)*conjg(cu(1,l,1,1,m)) + cu(2,l,1
     1,1,m)*conjg(cu(2,l,1,1,m)) + cu(3,l,1,1,m)*conjg(cu(3,l,1,1,m)))
  140       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            bxyz(1,1,1,1,m) = zero
            bxyz(2,1,1,1,m) = zero
            bxyz(3,1,1,1,m) = zero
            bxyz(1,l1,1,1,m) = zero
            bxyz(2,l1,1,1,m) = zero
            bxyz(3,l1,1,1,m) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = (nyh/kyzp)*kyzp
      if (n2.eq.k1) then
         k1 = nyh - k1 + 1
         do 160 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 150 l = 1, nz
            bxyz(1,l,j,k1,m) = zero
            bxyz(2,l,j,k1,m) = zero
            bxyz(3,l,j,k1,m) = zero
  150       continue
         endif
  160    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 170 l = 1, nz
            bxyz(1,l,1,k1,m) = zero
            bxyz(2,l,1,k1,m) = zero
            bxyz(3,l,1,k1,m) = zero
  170       continue
         endif
      endif
  180 continue
  190 continue
  200 continue
      wm = float(nx*ny*nz)*wp
      return
c calculate vector potential and sum field energy
  210 if (isign.gt.1) go to 360
      wp = 0.0d0
      if (kstrt.gt.(kxb*kyzb)) go to 350
      do 340 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 330 mx = 1, jblok
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      joff = kxyp*(mx + js) - 1
      m = mx + moff
      do 260 k = 1, kyzp
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         do 230 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 220 l = 2, nzh
            l1 = nz2 - l
            at2 = ci2*real(ffc(l,j,k,m))
            at1 = at2*aimag(ffc(l,j,k,m))
            bxyz(1,l,j,k,m) = at2*cu(1,l,j,k,m)
            bxyz(2,l,j,k,m) = at2*cu(2,l,j,k,m)
            bxyz(3,l,j,k,m) = at2*cu(3,l,j,k,m)
            bxyz(1,l1,j,k,m) = at2*cu(1,l1,j,k,m)
            bxyz(2,l1,j,k,m) = at2*cu(2,l1,j,k,m)
            bxyz(3,l1,j,k,m) = at2*cu(3,l1,j,k,m)
            wp = wp + at1*(cu(1,l,j,k,m)*conjg(cu(1,l,j,k,m)) + cu(2,l,j
     1,k,m)*conjg(cu(2,l,j,k,m)) + cu(3,l,j,k,m)*conjg(cu(3,l,j,k,m)) + 
     2cu(1,l1,j,k,m)*conjg(cu(1,l1,j,k,m)) + cu(2,l1,j,k,m)*conjg(cu(2,l
     31,j,k,m)) + cu(3,l1,j,k,m)*conjg(cu(3,l1,j,k,m)))
  220       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at2 = ci2*real(ffc(1,j,k,m))
            at1 = at2*aimag(ffc(1,j,k,m))
            bxyz(1,1,j,k,m) = at2*cu(1,1,j,k,m)
            bxyz(2,1,j,k,m) = at2*cu(2,1,j,k,m)
            bxyz(3,1,j,k,m) = at2*cu(3,1,j,k,m)
            bxyz(1,l1,j,k,m) = zero
            bxyz(2,l1,j,k,m) = zero
            bxyz(3,l1,j,k,m) = zero
            wp = wp + at1*(cu(1,1,j,k,m)*conjg(cu(1,1,j,k,m)) + cu(2,1,j
     1,k,m)*conjg(cu(2,1,j,k,m)) + cu(3,1,j,k,m)*conjg(cu(3,1,j,k,m)))
         endif
  230    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 240 l = 2, nzh
               l1 = nz2 - l
               at2 = ci2*real(ffc(l,1,k,m))
               at1 = at2*aimag(ffc(l,1,k,m))
               bxyz(1,l,1,k,m) = at2*cu(1,l,1,k,m)
               bxyz(2,l,1,k,m) = at2*cu(2,l,1,k,m)
               bxyz(3,l,1,k,m) = at2*cu(3,l,1,k,m)
               bxyz(1,l1,1,k,m) = at2*cu(1,l1,1,k,m)
               bxyz(2,l1,1,k,m) = at2*cu(2,l1,1,k,m)
               bxyz(3,l1,1,k,m) = at2*cu(3,l1,1,k,m)
               wp = wp + at1*(cu(1,l,1,k,m)*conjg(cu(1,l,1,k,m)) + cu(2,
     1l,1,k,m)*conjg(cu(2,l,1,k,m)) + cu(3,l,1,k,m)*conjg(cu(3,l,1,k,m))
     2 + cu(1,l1,1,k,m)*conjg(cu(1,l1,1,k,m)) + cu(2,l1,1,k,m)*conjg(cu(
     32,l1,1,k,m)) + cu(3,l1,1,k,m)*conjg(cu(3,l1,1,k,m)))
  240          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at2 = ci2*real(ffc(1,1,k,m))
               at1 = at2*aimag(ffc(1,1,k,m))
               bxyz(1,1,1,k,m) = at2*cu(1,1,1,k,m)
               bxyz(2,1,1,k,m) = at2*cu(2,1,1,k,m)
               bxyz(3,1,1,k,m) = at2*cu(3,1,1,k,m)
               bxyz(1,l1,1,k,m) = zero
               bxyz(2,l1,1,k,m) = zero
               bxyz(3,l1,1,k,m) = zero
               wp = wp + at1*(cu(1,1,1,k,m)*conjg(cu(1,1,1,k,m)) + cu(2,
     11,1,k,m)*conjg(cu(2,1,1,k,m)) + cu(3,1,1,k,m)*conjg(cu(3,1,1,k,m))
     2)
c throw away kx = nx/2
            else
               do 250 l = 1, nz
               bxyz(1,l,1,k,m) = zero
               bxyz(2,l,1,k,m) = zero
               bxyz(3,l,1,k,m) = zero
  250          continue
            endif
         endif
      endif
  260 continue
c mode numbers ky = 0, ny/2
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 280 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 270 l = 2, nzh
            l1 = nz2 - l
            at2 = ci2*real(ffc(l,j,1,m))
            at1 = at2*aimag(ffc(l,j,1,m))
            bxyz(1,l,j,1,m) = at2*cu(1,l,j,1,m)
            bxyz(2,l,j,1,m) = at2*cu(2,l,j,1,m)
            bxyz(3,l,j,1,m) = at2*cu(3,l,j,1,m)
            bxyz(1,l1,j,1,m) = at2*cu(1,l1,j,1,m)
            bxyz(2,l1,j,1,m) = at2*cu(2,l1,j,1,m)
            bxyz(3,l1,j,1,m) = at2*cu(3,l1,j,1,m)
            wp = wp + at1*(cu(1,l,j,1,m)*conjg(cu(1,l,j,1,m)) + cu(2,l,j
     1,1,m)*conjg(cu(2,l,j,1,m)) + cu(3,l,j,1,m)*conjg(cu(3,l,j,1,m)) + 
     2cu(1,l1,j,1,m)*conjg(cu(1,l1,j,1,m)) + cu(2,l1,j,1,m)*conjg(cu(2,l
     31,j,1,m)) + cu(3,l1,j,1,m)*conjg(cu(3,l1,j,1,m)))
  270       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at2 = ci2*real(ffc(1,j,1,m))
            at1 = at2*aimag(ffc(1,j,1,m))
            bxyz(1,1,j,1,m) = at2*cu(1,1,j,1,m)
            bxyz(2,1,j,1,m) = at2*cu(2,1,j,1,m)
            bxyz(3,1,j,1,m) = at2*cu(3,1,j,1,m)
            bxyz(1,l1,j,1,m) = zero
            bxyz(2,l1,j,1,m) = zero
            bxyz(3,l1,j,1,m) = zero
            wp = wp + at1*(cu(1,1,j,1,m)*conjg(cu(1,1,j,1,m)) + cu(2,1,j
     1,1,m)*conjg(cu(2,1,j,1,m)) + cu(3,1,j,1,m)*conjg(cu(3,1,j,1,m)))
         endif
  280    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 290 l = 2, nzh
            l1 = nz2 - l
            at2 = ci2*real(ffc(l,1,1,m))
            at1 = at2*aimag(ffc(l,1,1,m))
            bxyz(1,l,1,1,m) = at2*cu(1,l,1,1,m)
            bxyz(2,l,1,1,m) = at2*cu(2,l,1,1,m)
            bxyz(3,l,1,1,m) = at2*cu(3,l,1,1,m)
            bxyz(1,l1,1,1,m) = zero
            bxyz(2,l1,1,1,m) = zero
            bxyz(3,l1,1,1,m) = zero
            wp = wp + at1*(cu(1,l,1,1,m)*conjg(cu(1,l,1,1,m)) + cu(2,l,1
     1,1,m)*conjg(cu(2,l,1,1,m)) + cu(3,l,1,1,m)*conjg(cu(3,l,1,1,m)))
  290       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            bxyz(1,1,1,1,m) = zero
            bxyz(2,1,1,1,m) = zero
            bxyz(3,1,1,1,m) = zero
            bxyz(1,l1,1,1,m) = zero
            bxyz(2,l1,1,1,m) = zero
            bxyz(3,l1,1,1,m) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = (nyh/kyzp)*kyzp
      if (n2.eq.k1) then
         k1 = nyh - k1 + 1
         do 310 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 300 l = 1, nz
            bxyz(1,l,j,k1,m) = zero
            bxyz(2,l,j,k1,m) = zero
            bxyz(3,l,j,k1,m) = zero
  300       continue
         endif
  310    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 320 l = 1, nz
            bxyz(1,l,1,k1,m) = zero
            bxyz(2,l,1,k1,m) = zero
            bxyz(3,l,1,k1,m) = zero
  320       continue
         endif
      endif
  330 continue
  340 continue
  350 continue
      wm = float(nx*ny*nz)*wp
      return
c calculate smoothing
  360 if (kstrt.gt.(kxb*kyzb)) go to 500
      do 490 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 480 mx = 1, jblok
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      joff = kxyp*(mx + js) - 1
      m = mx + moff
      do 410 k = 1, kyzp
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         do 380 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 370 l = 2, nzh
            l1 = nz2 - l
            at1 = aimag(ffc(l,j,k,m))
            bxyz(1,l,j,k,m) = at1*cu(1,l,j,k,m)
            bxyz(2,l,j,k,m) = at1*cu(2,l,j,k,m)
            bxyz(3,l,j,k,m) = at1*cu(3,l,j,k,m)
            bxyz(1,l1,j,k,m) = at1*cu(1,l1,j,k,m)
            bxyz(2,l1,j,k,m) = at1*cu(2,l1,j,k,m)
            bxyz(3,l1,j,k,m) = at1*cu(3,l1,j,k,m)
  370       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = aimag(ffc(1,j,k,m))
            bxyz(1,1,j,k,m) = at1*cu(1,1,j,k,m)
            bxyz(2,1,j,k,m) = at1*cu(2,1,j,k,m)
            bxyz(3,1,j,k,m) = at1*cu(3,1,j,k,m)
            bxyz(1,l1,j,k,m) = zero
            bxyz(2,l1,j,k,m) = zero
            bxyz(3,l1,j,k,m) = zero
         endif
  380    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 390 l = 2, nzh
               l1 = nz2 - l
               at1 = aimag(ffc(l,1,k,m))
               bxyz(1,l,1,k,m) = at1*cu(1,l,1,k,m)
               bxyz(2,l,1,k,m) = at1*cu(2,l,1,k,m)
               bxyz(3,l,1,k,m) = at1*cu(3,l,1,k,m)
               bxyz(1,l1,1,k,m) = at1*cu(1,l1,1,k,m)
               bxyz(2,l1,1,k,m) = at1*cu(2,l1,1,k,m)
               bxyz(3,l1,1,k,m) = at1*cu(3,l1,1,k,m)
  390          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at1 = aimag(ffc(1,1,k,m))
               bxyz(1,1,1,k,m) = at1*cu(1,1,1,k,m)
               bxyz(2,1,1,k,m) = at1*cu(2,1,1,k,m)
               bxyz(3,1,1,k,m) = at1*cu(3,1,1,k,m)
               bxyz(1,l1,1,k,m) = zero
               bxyz(2,l1,1,k,m) = zero
               bxyz(3,l1,1,k,m) = zero
c throw away kx = nx/2
            else
               do 400 l = 1, nz
               bxyz(1,l,1,k,m) = zero
               bxyz(2,l,1,k,m) = zero
               bxyz(3,l,1,k,m) = zero
  400          continue
            endif
         endif
      endif
  410 continue
c mode numbers ky = 0, ny/2
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 430 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 420 l = 2, nzh
            l1 = nz2 - l
            at1 = aimag(ffc(l,j,1,m))
            bxyz(1,l,j,1,m) = at1*cu(1,l,j,1,m)
            bxyz(2,l,j,1,m) = at1*cu(2,l,j,1,m)
            bxyz(3,l,j,1,m) = at1*cu(3,l,j,1,m)
            bxyz(1,l1,j,1,m) = at1*cu(1,l1,j,1,m)
            bxyz(2,l1,j,1,m) = at1*cu(2,l1,j,1,m)
            bxyz(3,l1,j,1,m) = at1*cu(3,l1,j,1,m)
  420       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = aimag(ffc(1,j,1,m))
            bxyz(1,1,j,1,m) = at1*cu(1,1,j,1,m)
            bxyz(2,1,j,1,m) = at1*cu(2,1,j,1,m)
            bxyz(3,1,j,1,m) = at1*cu(3,1,j,1,m)
            bxyz(1,l1,j,1,m) = zero
            bxyz(2,l1,j,1,m) = zero
            bxyz(3,l1,j,1,m) = zero
         endif
  430    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 440 l = 2, nzh
            l1 = nz2 - l
            at1 = aimag(ffc(l,1,1,m))
            bxyz(1,l,1,1,m) = at1*cu(1,l,1,1,m)
            bxyz(2,l,1,1,m) = at1*cu(2,l,1,1,m)
            bxyz(3,l,1,1,m) = at1*cu(3,l,1,1,m)
            bxyz(1,l1,1,1,m) = zero
            bxyz(2,l1,1,1,m) = zero
            bxyz(3,l1,1,1,m) = zero
  440       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = aimag(ffc(1,1,1,m))
            bxyz(1,1,1,1,m) = cmplx(at1*real(cu(1,1,1,1,m)),0.)
            bxyz(2,1,1,1,m) = cmplx(at1*real(cu(2,1,1,1,m)),0.)
            bxyz(3,1,1,1,m) = cmplx(at1*real(cu(3,1,1,1,m)),0.)
            bxyz(1,l1,1,1,m) = zero
            bxyz(2,l1,1,1,m) = zero
            bxyz(3,l1,1,1,m) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = (nyh/kyzp)*kyzp
      if (n2.eq.k1) then
         k1 = nyh - k1 + 1
         do 460 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 450 l = 1, nz
            bxyz(1,l,j,k1,m) = zero
            bxyz(2,l,j,k1,m) = zero
            bxyz(3,l,j,k1,m) = zero
  450       continue
         endif
  460    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 470 l = 1, nz
            bxyz(1,l,1,k1,m) = zero
            bxyz(2,l,1,k1,m) = zero
            bxyz(3,l,1,k1,m) = zero
  470       continue
         endif
      endif
  480 continue
  490 continue
  500 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine IPBPOISP332(cu,bxyz,ffc,ci,wm,nx,ny,nz,kstrt,nzv,kxyp,k
     1yzp,jblok,mblok,nzhd)
c this subroutine solves 3d poisson's equation in fourier space for
c magnetic field with periodic boundary conditions, for distributed data
c with 2D spatial decomposition
c input: cu,ffc,ci,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mblok,nzhd
c output: bxyz, wm
c approximate flop count is:
c 193*nxc*nyc*nzc + 84*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
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
c s(kx,ky,kz) = exp(-((kx*ax)**2+(ky*ay)**2+(kz*az)**2)/2), except for
c bx(kx=pi) = by(kx=pi) = bz(kx=pi) = 0,
c bx(ky=pi) = by(ky=pi) = bx(ky=pi) = 0,
c bx(kz=pi) = by(kz=pi) = bz(kz=pi) = 0,
c bx(kx=0,ky=0,kz=0) = by(kx=0,ky=0,kz=0) = bz(kx=0,ky=0,kz=0) = 0.
c cu(l,j,k,m) = complex current density for fourier mode jj-1,kk-1,l-1
c bxyz(i,l,j,k,m) = i component of complex magnetic field
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*(mx - 1) and
c kk = k + kyzp*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c aimag(ffc(l,j,k,m)) = finite-size particle shape factor s
c real(ffc(l,j,k,m)) = potential green's function g
c ci = reciprical of velocity of light
c magnetic field energy is also calculated, using
c wm = nx*ny*nz*sum((affp/(kx**2+ky**2+kz**2))*ci*ci
c    |cu(kx,ky,kz)*s(kx,ky,kz)|**2)
c this expression is valid only if the current is divergence-free
c nx/ny/nz = system length in x/y/z direction
c jblok/mblok = number of field partitions in x/y
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c nzv = first dimension of field arrays, must be >= nz
c nzhd = first dimension of form factor array, must be >= nzh
      double precision wp
      complex cu, bxyz, ffc, zero, zt1, zt2, zt3
      dimension cu(3,nzv,kxyp,kyzp,jblok*mblok)
      dimension bxyz(3,nzv,kxyp,kyzp,jblok*mblok)
      dimension ffc(nzhd,kxyp,kyzp,jblok*mblok)
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2
      nz2 = nz + 2
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      dnz = 6.28318530717959/float(nz)
      zero = cmplx(0.,0.)
      ci2 = ci*ci
c calculate magnetic field and sum field energy
      wp = 0.0d0
      if (kstrt.gt.(kxb*kyzb)) go to 140
      do 130 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 120 mx = 1, jblok
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      joff = kxyp*(mx + js) - 1
      m = mx + moff
      do 50 k = 1, kyzp
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*float(k1)
         do 20 j = 1, kxyp
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            at1 = ci2*real(ffc(l,j,k,m))
            at2 = dkx*at1
            at3 = dky*at1
            at4 = dnz*float(l - 1)*at1
            at1 = at1*aimag(ffc(l,j,k,m))
            zt1 = cmplx(-aimag(cu(3,l,j,k,m)),real(cu(3,l,j,k,m)))
            zt2 = cmplx(-aimag(cu(2,l,j,k,m)),real(cu(2,l,j,k,m)))
            zt3 = cmplx(-aimag(cu(1,l,j,k,m)),real(cu(1,l,j,k,m)))
            bxyz(1,l,j,k,m) = at3*zt1 - at4*zt2
            bxyz(2,l,j,k,m) = at4*zt3 - at2*zt1
            bxyz(3,l,j,k,m) = at2*zt2 - at3*zt3
            zt1 = cmplx(-aimag(cu(3,l1,j,k,m)),real(cu(3,l1,j,k,m)))
            zt2 = cmplx(-aimag(cu(2,l1,j,k,m)),real(cu(2,l1,j,k,m)))
            zt3 = cmplx(-aimag(cu(1,l1,j,k,m)),real(cu(1,l1,j,k,m)))
            bxyz(1,l1,j,k,m) = at3*zt1 + at4*zt2
            bxyz(2,l1,j,k,m) = -at4*zt3 - at2*zt1
            bxyz(3,l1,j,k,m) = at2*zt2 - at3*zt3
            wp = wp + at1*(cu(1,l,j,k,m)*conjg(cu(1,l,j,k,m)) + cu(2,l,j
     1,k,m)*conjg(cu(2,l,j,k,m)) + cu(3,l,j,k,m)*conjg(cu(3,l,j,k,m)) + 
     2cu(1,l1,j,k,m)*conjg(cu(1,l1,j,k,m)) + cu(2,l1,j,k,m)*conjg(cu(2,l
     31,j,k,m)) + cu(3,l1,j,k,m)*conjg(cu(3,l1,j,k,m)))
   10       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = ci2*real(ffc(1,j,k,m))
            at2 = dkx*at1
            at3 = dky*at1
            at1 = at1*aimag(ffc(1,j,k,m))
            zt1 = cmplx(-aimag(cu(3,1,j,k,m)),real(cu(3,1,j,k,m)))
            zt2 = cmplx(-aimag(cu(2,1,j,k,m)),real(cu(2,1,j,k,m)))
            zt3 = cmplx(-aimag(cu(1,1,j,k,m)),real(cu(1,1,j,k,m)))
            bxyz(1,1,j,k,m) = at3*zt1
            bxyz(2,1,j,k,m) = -at2*zt1
            bxyz(3,1,j,k,m) = at2*zt2 - at3*zt3
            bxyz(1,l1,j,k,m) = zero
            bxyz(2,l1,j,k,m) = zero
            bxyz(3,l1,j,k,m) = zero
            wp = wp + at1*(cu(1,1,j,k,m)*conjg(cu(1,1,j,k,m)) + cu(2,1,j
     1,k,m)*conjg(cu(2,1,j,k,m)) + cu(3,1,j,k,m)*conjg(cu(3,1,j,k,m)))
         endif
   20    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               at1 = ci2*real(ffc(l,1,k,m))
               at3 = dky*at1
               at4 = dnz*float(l - 1)*at1
               at1 = at1*aimag(ffc(l,1,k,m))
               zt1 = cmplx(-aimag(cu(3,l,1,k,m)),real(cu(3,l,1,k,m)))
               zt2 = cmplx(-aimag(cu(2,l,1,k,m)),real(cu(2,l,1,k,m)))
               zt3 = cmplx(-aimag(cu(1,l,1,k,m)),real(cu(1,l,1,k,m)))
               bxyz(1,l,1,k,m) = at3*zt1 - at4*zt2
               bxyz(2,l,1,k,m) = at4*zt3
               bxyz(3,l,1,k,m) = -at3*zt3
               zt1 = cmplx(-aimag(cu(3,l1,1,k,m)),real(cu(3,l1,1,k,m)))
               zt2 = cmplx(-aimag(cu(2,l1,1,k,m)),real(cu(2,l1,1,k,m)))
               zt3 = cmplx(-aimag(cu(1,l1,1,k,m)),real(cu(1,l1,1,k,m)))
               bxyz(1,l1,1,k,m) = at3*zt1 + at4*zt2
               bxyz(2,l1,1,k,m) = -at4*zt3
               bxyz(3,l1,1,k,m) = -at3*zt3
               wp = wp + at1*(cu(1,l,1,k,m)*conjg(cu(1,l,1,k,m)) + cu(2,
     1l,1,k,m)*conjg(cu(2,l,1,k,m)) + cu(3,l,1,k,m)*conjg(cu(3,l,1,k,m))
     2 + cu(1,l1,1,k,m)*conjg(cu(1,l1,1,k,m)) + cu(2,l1,1,k,m)*conjg(cu(
     32,l1,1,k,m)) + cu(3,l1,1,k,m)*conjg(cu(3,l1,1,k,m)))
   30          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at1 = ci2*real(ffc(1,1,k,m))
               at3 = dky*at1
               at1 = at1*aimag(ffc(1,1,k,m))
               zt1 = cmplx(-aimag(cu(3,1,1,k,m)),real(cu(3,1,1,k,m)))
               zt3 = cmplx(-aimag(cu(1,1,1,k,m)),real(cu(1,1,1,k,m)))
               bxyz(1,1,1,k,m) = at3*zt1
               bxyz(2,1,1,k,m) = zero
               bxyz(3,1,1,k,m) = -at3*zt3
               bxyz(1,l1,1,k,m) = zero
               bxyz(2,l1,1,k,m) = zero
               bxyz(3,l1,1,k,m) = zero
               wp = wp + at1*(cu(1,1,1,k,m)*conjg(cu(1,1,1,k,m)) + cu(2,
     11,1,k,m)*conjg(cu(2,1,1,k,m)) + cu(3,1,1,k,m)*conjg(cu(3,1,1,k,m))
     2)
c throw away kx = nx/2
            else
               do 40 l = 1, nz
               bxyz(1,l,1,k,m) = zero
               bxyz(2,l,1,k,m) = zero
               bxyz(3,l,1,k,m) = zero
   40          continue
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny/2
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 70 j = 1, kxyp
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            at1 = ci2*real(ffc(l,j,1,m))
            at2 = dkx*at1
            at4 = dnz*float(l - 1)*at1
            at1 = at1*aimag(ffc(l,j,1,m))
            zt1 = cmplx(-aimag(cu(3,l,j,1,m)),real(cu(3,l,j,1,m)))
            zt2 = cmplx(-aimag(cu(2,l,j,1,m)),real(cu(2,l,j,1,m)))
            zt3 = cmplx(-aimag(cu(1,l,j,1,m)),real(cu(1,l,j,1,m)))
            bxyz(1,l,j,1,m) = -at4*zt2
            bxyz(2,l,j,1,m) = at4*zt3 - at2*zt1
            bxyz(3,l,j,1,m) = at2*zt2
            zt1 = cmplx(-aimag(cu(3,l1,j,1,m)),real(cu(3,l1,j,1,m)))
            zt2 = cmplx(-aimag(cu(2,l1,j,1,m)),real(cu(2,l1,j,1,m)))
            zt3 = cmplx(-aimag(cu(1,l1,j,1,m)),real(cu(1,l1,j,1,m)))
            bxyz(1,l1,j,1,m) = at4*zt2
            bxyz(2,l1,j,1,m) = -at4*zt3 - at2*zt1
            bxyz(3,l1,j,1,m) = at2*zt2
            wp = wp + at1*(cu(1,l,j,1,m)*conjg(cu(1,l,j,1,m)) + cu(2,l,j
     1,1,m)*conjg(cu(2,l,j,1,m)) + cu(3,l,j,1,m)*conjg(cu(3,l,j,1,m)) + 
     2cu(1,l1,j,1,m)*conjg(cu(1,l1,j,1,m)) + cu(2,l1,j,1,m)*conjg(cu(2,l
     31,j,1,m)) + cu(3,l1,j,1,m)*conjg(cu(3,l1,j,1,m)))
   60       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = ci2*real(ffc(1,j,1,m))
            at2 = dkx*at1
            at1 = at1*aimag(ffc(1,j,1,m))
            zt1 = cmplx(-aimag(cu(3,1,j,1,m)),real(cu(3,1,j,1,m)))
            zt2 = cmplx(-aimag(cu(2,1,j,1,m)),real(cu(2,1,j,1,m)))
            bxyz(1,1,j,1,m) = zero
            bxyz(2,1,j,1,m) = -at2*zt1
            bxyz(3,1,j,1,m) = at2*zt2
            bxyz(1,l1,j,1,m) = zero
            bxyz(2,l1,j,1,m) = zero
            bxyz(3,l1,j,1,m) = zero
            wp = wp + at1*(cu(1,1,j,1,m)*conjg(cu(1,1,j,1,m)) + cu(2,1,j
     1,1,m)*conjg(cu(2,1,j,1,m)) + cu(3,1,j,1,m)*conjg(cu(3,1,j,1,m)))
         endif
   70    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            at1 = ci2*real(ffc(l,1,1,m))
            at4 = dnz*float(l - 1)*at1
            at1 = at1*aimag(ffc(l,1,1,m))
            zt2 = cmplx(-aimag(cu(2,l,1,1,m)),real(cu(2,l,1,1,m)))
            zt3 = cmplx(-aimag(cu(1,l,1,1,m)),real(cu(1,l,1,1,m)))
            bxyz(1,l,1,1,m) = -at4*zt2
            bxyz(2,l,1,1,m) = at4*zt3
            bxyz(3,l,1,1,m) = zero
            bxyz(1,l1,1,1,m) = zero
            bxyz(2,l1,1,1,m) = zero
            bxyz(3,l1,1,1,m) = zero
            wp = wp + at1*(cu(1,l,1,1,m)*conjg(cu(1,l,1,1,m)) + cu(2,l,1
     1,1,m)*conjg(cu(2,l,1,1,m)) + cu(3,l,1,1,m)*conjg(cu(3,l,1,1,m)))
   80       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            bxyz(1,1,1,1,m) = zero
            bxyz(2,1,1,1,m) = zero
            bxyz(3,1,1,1,m) = zero
            bxyz(1,l1,1,1,m) = zero
            bxyz(2,l1,1,1,m) = zero
            bxyz(3,l1,1,1,m) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = (nyh/kyzp)*kyzp
      if (n2.eq.k1) then
         k1 = nyh - k1 + 1
         do 100 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            bxyz(1,l,j,k1,m) = zero
            bxyz(2,l,j,k1,m) = zero
            bxyz(3,l,j,k1,m) = zero
   90       continue
         endif
  100    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 110 l = 1, nz
            bxyz(1,l,1,k1,m) = zero
            bxyz(2,l,1,k1,m) = zero
            bxyz(3,l,1,k1,m) = zero
  110       continue
         endif
      endif
  120 continue
  130 continue
  140 continue
      wm = float(nx*ny*nz)*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PMAXWEL32(exyz,bxyz,cu,ffc,affp,ci,dt,wf,wm,nx,ny,nz,ks
     1trt,nzv,kxyp,kyzp,jblok,mblok,nzhd)
c this subroutine solves 3d maxwell's equation in fourier space for
c transverse electric and magnetic fields with periodic boundary
c conditions, for distributed data with 2D spatial decomposition
c input: all, output: wf, wm, exyz, bxyz
c approximate flop count is:
c 679*nxc*nyc*nzc + 149*(nxc*nyc + nxc*nzc + nyc*nzc)
c plus nxc*nyc*nzc divides
c where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz - 1, and
c nvpy/nvpz = number of procs in y/z
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
c j,k,l = fourier mode numbers, except for
c ex(kx=pi) = ey(kx=pi) = ez(kx=pi) = 0,
c ex(ky=pi) = ey(ky=pi) = ex(ky=pi) = 0,
c ex(kz=pi) = ey(kz=pi) = ez(kz=pi) = 0,
c ex(kx=0,ky=0,kz=0) = ey(kx=0,ky=0,kz=0) = ez(kx=0,ky=0,kz=0) = 0.
c and similarly for bx, by, bz.
c exyz(i,l,j,k,m) = i component of complex transverse electric field
c bxyz(i,l,j,k,m) = i component of complex magnetic field
c cu(i,l,j,k,m) = i component of complex current density
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*(mx - 1) and
c kk = k + kyzp*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c aimag(ffc(l,j,k,m)) = finite-size particle shape factor s,
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
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c nzv = first dimension of field arrays, must be >= nz
c jblok/mblok = number of field partitions in x/y
c nzhd = first dimension of form factor array, must be >= nzh
      double precision wp, ws
      complex exyz, bxyz, cu, ffc
      complex zero, zt1, zt2, zt3, zt4, zt5, zt6, zt7, zt8, zt9
      dimension exyz(3,nzv,kxyp,kyzp,jblok*mblok)
      dimension bxyz(3,nzv,kxyp,kyzp,jblok*mblok)
      dimension cu(3,nzv,kxyp,kyzp,jblok*mblok)
      dimension ffc(nzhd,kxyp,kyzp,jblok*mblok)
      if (ci.le.0.) return
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2
      nz2 = nz + 2
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      dnz = 6.28318530717959/float(nz)
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
      if (kstrt.gt.(kxb*kyzb)) go to 140
      do 130 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 120 mx = 1, jblok
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      joff = kxyp*(mx + js) - 1
      m = mx + moff
      do 50 k = 1, kyzp
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*float(k1)
         do 20 j = 1, kxyp
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            afdt = adt*aimag(ffc(l,j,k,m))
c update magnetic field half time step, ky > 0, kz > 0
            zt1 = cmplx(-aimag(exyz(3,l,j,k,m)),real(exyz(3,l,j,k,m)))
            zt2 = cmplx(-aimag(exyz(2,l,j,k,m)),real(exyz(2,l,j,k,m)))
            zt3 = cmplx(-aimag(exyz(1,l,j,k,m)),real(exyz(1,l,j,k,m)))
            zt4 = bxyz(1,l,j,k,m) - dth*(dky*zt1 - dkz*zt2)
            zt5 = bxyz(2,l,j,k,m) - dth*(dkz*zt3 - dkx*zt1)
            zt6 = bxyz(3,l,j,k,m) - dth*(dkx*zt2 - dky*zt3)
c update electric field whole time step
            zt1 = cmplx(-aimag(zt6),real(zt6))
            zt2 = cmplx(-aimag(zt5),real(zt5))
            zt3 = cmplx(-aimag(zt4),real(zt4))
            zt7 = exyz(1,l,j,k,m) + cdt*(dky*zt1 - dkz*zt2) - afdt*cu(1,
     1l,j,k,m)
            zt8 = exyz(2,l,j,k,m) + cdt*(dkz*zt3 - dkx*zt1) - afdt*cu(2,
     1l,j,k,m)
            zt9 = exyz(3,l,j,k,m) + cdt*(dkx*zt2 - dky*zt3) - afdt*cu(3,
     1l,j,k,m)
c update magnetic field half time step and store electric field
            zt1 = cmplx(-aimag(zt9),real(zt9))
            zt2 = cmplx(-aimag(zt8),real(zt8))
            zt3 = cmplx(-aimag(zt7),real(zt7))
            exyz(1,l,j,k,m) = zt7
            exyz(2,l,j,k,m) = zt8
            exyz(3,l,j,k,m) = zt9
            ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg
     1(zt9))
            zt4 = zt4 - dth*(dky*zt1 - dkz*zt2)
            zt5 = zt5 - dth*(dkz*zt3 - dkx*zt1)
            zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
            bxyz(1,l,j,k,m) = zt4
            bxyz(2,l,j,k,m) = zt5
            bxyz(3,l,j,k,m) = zt6
            wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg
     1(zt6))
c update magnetic field half time step, ky > 0, kz < 0
            zt1 = cmplx(-aimag(exyz(3,l1,j,k,m)),real(exyz(3,l1,j,k,m)))
            zt2 = cmplx(-aimag(exyz(2,l1,j,k,m)),real(exyz(2,l1,j,k,m)))
            zt3 = cmplx(-aimag(exyz(1,l1,j,k,m)),real(exyz(1,l1,j,k,m)))
            zt4 = bxyz(1,l1,j,k,m) - dth*(dky*zt1 + dkz*zt2)
            zt5 = bxyz(2,l1,j,k,m) + dth*(dkz*zt3 + dkx*zt1)
            zt6 = bxyz(3,l1,j,k,m) - dth*(dkx*zt2 - dky*zt3)
c update electric field whole time step
            zt1 = cmplx(-aimag(zt6),real(zt6))
            zt2 = cmplx(-aimag(zt5),real(zt5))
            zt3 = cmplx(-aimag(zt4),real(zt4))
            zt7 = exyz(1,l1,j,k,m) + cdt*(dky*zt1 + dkz*zt2) - afdt*cu(1
     1,l1,j,k,m)
            zt8 = exyz(2,l1,j,k,m) - cdt*(dkz*zt3 + dkx*zt1) - afdt*cu(2
     1,l1,j,k,m)
            zt9 = exyz(3,l1,j,k,m) + cdt*(dkx*zt2 - dky*zt3) - afdt*cu(3
     1,l1,j,k,m)
c update magnetic field half time step and store electric field
            zt1 = cmplx(-aimag(zt9),real(zt9))
            zt2 = cmplx(-aimag(zt8),real(zt8))
            zt3 = cmplx(-aimag(zt7),real(zt7))
            exyz(1,l1,j,k,m) = zt7
            exyz(2,l1,j,k,m) = zt8
            exyz(3,l1,j,k,m) = zt9
            ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg
     1(zt9))
            zt4 = zt4 - dth*(dky*zt1 + dkz*zt2)
            zt5 = zt5 + dth*(dkz*zt3 + dkx*zt1)
            zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
            bxyz(1,l1,j,k,m) = zt4
            bxyz(2,l1,j,k,m) = zt5
            bxyz(3,l1,j,k,m) = zt6
            wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg
     1(zt6))
   10       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            afdt = adt*aimag(ffc(1,j,k,m))
c update magnetic field half time step, ky > 0
            zt1 = cmplx(-aimag(exyz(3,1,j,k,m)),real(exyz(3,1,j,k,m)))
            zt2 = cmplx(-aimag(exyz(2,1,j,k,m)),real(exyz(2,1,j,k,m)))
            zt3 = cmplx(-aimag(exyz(1,1,j,k,m)),real(exyz(1,1,j,k,m)))
            zt4 = bxyz(1,1,j,k,m) - dth*(dky*zt1)
            zt5 = bxyz(2,1,j,k,m) + dth*(dkx*zt1)
            zt6 = bxyz(3,1,j,k,m) - dth*(dkx*zt2 - dky*zt3)
c update electric field whole time step
            zt1 = cmplx(-aimag(zt6),real(zt6))
            zt2 = cmplx(-aimag(zt5),real(zt5))
            zt3 = cmplx(-aimag(zt4),real(zt4))
            zt7 = exyz(1,1,j,k,m) + cdt*(dky*zt1) - afdt*cu(1,1,j,k,m)
            zt8 = exyz(2,1,j,k,m) - cdt*(dkx*zt1) - afdt*cu(2,1,j,k,m)
            zt9 = exyz(3,1,j,k,m) + cdt*(dkx*zt2 - dky*zt3) - afdt*cu(3,
     11,j,k,m)
c update magnetic field half time step and store electric field
            zt1 = cmplx(-aimag(zt9),real(zt9))
            zt2 = cmplx(-aimag(zt8),real(zt8))
            zt3 = cmplx(-aimag(zt7),real(zt7))
            exyz(1,1,j,k,m) = zt7
            exyz(2,1,j,k,m) = zt8
            exyz(3,1,j,k,m) = zt9
            ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg
     1(zt9))
            zt4 = zt4 - dth*(dky*zt1)
            zt5 = zt5 + dth*(dkx*zt1)
            zt6 = zt6 - dth*(dkx*zt2 - dky*zt3)
            bxyz(1,1,j,k,m) = zt4
            bxyz(2,1,j,k,m) = zt5
            bxyz(3,1,j,k,m) = zt6
            wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg
     1(zt6))
            bxyz(1,l1,j,k,m) = zero
            bxyz(2,l1,j,k,m) = zero
            bxyz(3,l1,j,k,m) = zero
            exyz(1,l1,j,k,m) = zero
            exyz(2,l1,j,k,m) = zero
            exyz(3,l1,j,k,m) = zero
         endif
   20    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               dkz = dnz*float(l - 1)
               afdt = adt*aimag(ffc(l,1,k,m))
c update magnetic field half time step, kz > 0
               zt1 = cmplx(-aimag(exyz(3,l,1,k,m)),real(exyz(3,l,1,k,m))
     1)
               zt2 = cmplx(-aimag(exyz(2,l,1,k,m)),real(exyz(2,l,1,k,m))
     1)
               zt3 = cmplx(-aimag(exyz(1,l,1,k,m)),real(exyz(1,l,1,k,m))
     1)
               zt4 = bxyz(1,l,1,k,m) - dth*(dky*zt1 - dkz*zt2)
               zt5 = bxyz(2,l,1,k,m) - dth*(dkz*zt3)
               zt6 = bxyz(3,l,1,k,m) + dth*(dky*zt3)
c update electric field whole time step
               zt1 = cmplx(-aimag(zt6),real(zt6))
               zt2 = cmplx(-aimag(zt5),real(zt5))
               zt3 = cmplx(-aimag(zt4),real(zt4))
               zt7 = exyz(1,l,1,k,m) + cdt*(dky*zt1 - dkz*zt2) - afdt*cu
     1(1,l,1,k,m)
               zt8 = exyz(2,l,1,k,m) + cdt*(dkz*zt3) - afdt*cu(2,l,1,k,m
     1)
               zt9 = exyz(3,l,1,k,m) - cdt*(dky*zt3) - afdt*cu(3,l,1,k,m
     1)
c update magnetic field half time step and store electric field
               zt1 = cmplx(-aimag(zt9),real(zt9))
               zt2 = cmplx(-aimag(zt8),real(zt8))
               zt3 = cmplx(-aimag(zt7),real(zt7))
               exyz(1,l,1,k,m) = zt7
               exyz(2,l,1,k,m) = zt8
               exyz(3,l,1,k,m) = zt9
               ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*co
     1njg(zt9))
               zt4 = zt4 - dth*(dky*zt1 - dkz*zt2)
               zt5 = zt5 - dth*(dkz*zt3)
               zt6 = zt6 + dth*(dky*zt3) 
               bxyz(1,l,1,k,m) = zt4
               bxyz(2,l,1,k,m) = zt5
               bxyz(3,l,1,k,m) = zt6
               wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*co
     1njg(zt6))
c update magnetic field half time step, kz < 0
               zt1 = cmplx(-aimag(exyz(3,l1,1,k,m)),real(exyz(3,l1,1,k,m
     1)))
               zt2 = cmplx(-aimag(exyz(2,l1,1,k,m)),real(exyz(2,l1,1,k,m
     1)))
               zt3 = cmplx(-aimag(exyz(1,l1,1,k,m)),real(exyz(1,l1,1,k,m
     1)))
               zt4 = bxyz(1,l1,1,k,m) - dth*(dky*zt1 + dkz*zt2)
               zt5 = bxyz(2,l1,1,k,m) + dth*(dkz*zt3)
               zt6 = bxyz(3,l1,1,k,m) + dth*(dky*zt3)
c update electric field whole time step
               zt1 = cmplx(-aimag(zt6),real(zt6))
               zt2 = cmplx(-aimag(zt5),real(zt5))
               zt3 = cmplx(-aimag(zt4),real(zt4))
               zt7 = exyz(1,l1,1,k,m) + cdt*(dky*zt1 + dkz*zt2) - afdt*c
     1u(1,l1,1,k,m)
               zt8 = exyz(2,l1,1,k,m) - cdt*(dkz*zt3) - afdt*cu(2,l1,1,k
     1,m)
               zt9 = exyz(3,l1,1,k,m) - cdt*(dky*zt3) - afdt*cu(3,l1,1,k
     1,m)
c update magnetic field half time step and store electric field
               zt1 = cmplx(-aimag(zt9),real(zt9))
               zt2 = cmplx(-aimag(zt8),real(zt8))
               zt3 = cmplx(-aimag(zt7),real(zt7))
               exyz(1,l1,1,k,m) = zt7
               exyz(2,l1,1,k,m) = zt8
               exyz(3,l1,1,k,m) = zt9
               ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*co
     1njg(zt9))
               zt4 = zt4 - dth*(dky*zt1 + dkz*zt2)
               zt5 = zt5 + dth*(dkz*zt3)
               zt6 = zt6 + dth*(dky*zt3)
               bxyz(1,l1,1,k,m) = zt4
               bxyz(2,l1,1,k,m) = zt5
               bxyz(3,l1,1,k,m) = zt6
               wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*co
     1njg(zt6))
   30          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               afdt = adt*aimag(ffc(1,1,k,m))
c update magnetic field half time step
               zt1 = cmplx(-aimag(exyz(3,1,1,k,m)),real(exyz(3,1,1,k,m))
     1)
               zt3 = cmplx(-aimag(exyz(1,1,1,k,m)),real(exyz(1,1,1,k,m))
     1)
               zt4 = bxyz(1,1,1,k,m) - dth*(dky*zt1)
               zt6 = bxyz(3,1,1,k,m) + dth*(dky*zt3)
c update electric field whole time step
               zt1 = cmplx(-aimag(zt6),real(zt6))
               zt3 = cmplx(-aimag(zt4),real(zt4))
               zt7 = exyz(1,1,1,k,m) + cdt*(dky*zt1) - afdt*cu(1,1,1,k,m
     1)
               zt9 = exyz(3,1,1,k,m) - cdt*(dky*zt3) - afdt*cu(3,1,1,k,m
     1)
c update magnetic field half time step and store electric field
               zt1 = cmplx(-aimag(zt9),real(zt9))
               zt3 = cmplx(-aimag(zt7),real(zt7))
               exyz(1,1,1,k,m) = zt7
               exyz(2,1,1,k,m) = zero
               exyz(3,1,1,k,m) = zt9
               ws = ws + anorm*(zt7*conjg(zt7) + zt9*conjg(zt9))
               zt4 = zt4 - dth*(dky*zt1)
               zt6 = zt6 + dth*(dky*zt3)
               bxyz(1,1,1,k,m) = zt4
               bxyz(2,1,1,k,m) = zero
               bxyz(3,1,1,k,m) = zt6
               wp = wp + anorm*(zt4*conjg(zt4) + zt6*conjg(zt6))
               bxyz(1,l1,1,k,m) = zero
               bxyz(2,l1,1,k,m) = zero
               bxyz(3,l1,1,k,m) = zero
               exyz(1,l1,1,k,m) = zero
               exyz(2,l1,1,k,m) = zero
               exyz(3,l1,1,k,m) = zero
c throw away kx = nx/2
            else
               do 40 l = 1, nz
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
c mode numbers ky = 0, ny/2
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 70 j = 1, kxyp
         dkx = dnx*float(j + joff)
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            afdt = adt*aimag(ffc(l,j,1,m))
            zt1 = cmplx(-aimag(exyz(3,l,j,1,m)),real(exyz(3,l,j,1,m)))
            zt2 = cmplx(-aimag(exyz(2,l,j,1,m)),real(exyz(2,l,j,1,m)))
            zt3 = cmplx(-aimag(exyz(1,l,j,1,m)),real(exyz(1,l,j,1,m)))
            zt4 = bxyz(1,l,j,1,m) + dth*(dkz*zt2)
            zt5 = bxyz(2,l,j,1,m) - dth*(dkz*zt3 - dkx*zt1)
            zt6 = bxyz(3,l,j,1,m) - dth*(dkx*zt2)
c update electric field whole time step
            zt1 = cmplx(-aimag(zt6),real(zt6))
            zt2 = cmplx(-aimag(zt5),real(zt5))
            zt3 = cmplx(-aimag(zt4),real(zt4))
            zt7 = exyz(1,l,j,1,m) - cdt*(dkz*zt2) - afdt*cu(1,l,j,1,m)
            zt8 = exyz(2,l,j,1,m) + cdt*(dkz*zt3 - dkx*zt1) - afdt*cu(2,
     1l,j,1,m)
            zt9 = exyz(3,l,j,1,m) + cdt*(dkx*zt2) - afdt*cu(3,l,j,1,m)
c update magnetic field half time step and store electric field
            zt1 = cmplx(-aimag(zt9),real(zt9))
            zt2 = cmplx(-aimag(zt8),real(zt8))
            zt3 = cmplx(-aimag(zt7),real(zt7))
            exyz(1,l,j,1,m) = zt7
            exyz(2,l,j,1,m) = zt8
            exyz(3,l,j,1,m) = zt9
            ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg
     1(zt9))
            zt4 = zt4 + dth*(dkz*zt2)
            zt5 = zt5 - dth*(dkz*zt3 - dkx*zt1)
            zt6 = zt6 - dth*(dkx*zt2)
            bxyz(1,l,j,1,m) = zt4
            bxyz(2,l,j,1,m) = zt5
            bxyz(3,l,j,1,m) = zt6
            wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg
     1(zt6))
c update magnetic field half time step, kz < 0
            zt1 = cmplx(-aimag(exyz(3,l1,j,1,m)),real(exyz(3,l1,j,1,m)))
            zt2 = cmplx(-aimag(exyz(2,l1,j,1,m)),real(exyz(2,l1,j,1,m)))
            zt3 = cmplx(-aimag(exyz(1,l1,j,1,m)),real(exyz(1,l1,j,1,m)))
            zt4 = bxyz(1,l1,j,1,m) - dth*(dkz*zt2)
            zt5 = bxyz(2,l1,j,1,m) + dth*(dkz*zt3 + dkx*zt1)
            zt6 = bxyz(3,l1,j,1,m) - dth*(dkx*zt2)
c update electric field whole time step
            zt1 = cmplx(-aimag(zt6),real(zt6))
            zt2 = cmplx(-aimag(zt5),real(zt5))
            zt3 = cmplx(-aimag(zt4),real(zt4))
            zt7 = exyz(1,l1,j,1,m) + cdt*(dkz*zt2) - afdt*cu(1,l1,j,1,m)
            zt8 = exyz(2,l1,j,1,m) - cdt*(dkz*zt3 + dkx*zt1) - afdt*cu(2
     1,l1,j,1,m)
            zt9 = exyz(3,l1,j,1,m) + cdt*(dkx*zt2) - afdt*cu(3,l1,j,1,m)
c update magnetic field half time step and store electric field
            zt1 = cmplx(-aimag(zt9),real(zt9))
            zt2 = cmplx(-aimag(zt8),real(zt8))
            zt3 = cmplx(-aimag(zt7),real(zt7))
            exyz(1,l1,j,1,m) = zt7
            exyz(2,l1,j,1,m) = zt8
            exyz(3,l1,j,1,m) = zt9
            ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8) + zt9*conjg
     1(zt9))
            zt4 = zt4 - dth*(dkz*zt2)
            zt5 = zt5 + dth*(dkz*zt3 + dkx*zt1)
            zt6 = zt6 - dth*(dkx*zt2)
            bxyz(1,l1,j,1,m) = zt4
            bxyz(2,l1,j,1,m) = zt5
            bxyz(3,l1,j,1,m) = zt6
            wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5) + zt6*conjg
     1(zt6))
   60       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            afdt = adt*aimag(ffc(1,j,1,m))
c update magnetic field half time step
            zt1 = cmplx(-aimag(exyz(3,1,j,1,m)),real(exyz(3,1,j,1,m)))
            zt2 = cmplx(-aimag(exyz(2,1,j,1,m)),real(exyz(2,1,j,1,m)))
            zt5 = bxyz(2,1,j,1,m) + dth*(dkx*zt1)
            zt6 = bxyz(3,1,j,1,m) - dth*(dkx*zt2)
c update electric field whole time step
            zt1 = cmplx(-aimag(zt6),real(zt6))
            zt2 = cmplx(-aimag(zt5),real(zt5))
            zt8 = exyz(2,1,j,1,m) - cdt*(dkx*zt1) - afdt*cu(2,1,j,1,m)
            zt9 = exyz(3,1,j,1,m) + cdt*(dkx*zt2) - afdt*cu(3,1,j,1,m)
c update magnetic field half time step and store electric field
            zt1 = cmplx(-aimag(zt9),real(zt9))
            zt2 = cmplx(-aimag(zt8),real(zt8))
            exyz(1,1,j,1,m) = zero
            exyz(2,1,j,1,m) = zt8
            exyz(3,1,j,1,m) = zt9
            ws = ws + anorm*(zt8*conjg(zt8) + zt9*conjg(zt9))
            zt5 = zt5 + dth*(dkx*zt1)
            zt6 = zt6 - dth*(dkx*zt2)
            bxyz(1,1,j,1,m) = zero
            bxyz(2,1,j,1,m) = zt5
            bxyz(3,1,j,1,m) = zt6
            wp = wp + anorm*(zt5*conjg(zt5) + zt6*conjg(zt6))
            bxyz(1,l1,j,1,m) = zero
            bxyz(2,l1,j,1,m) = zero
            bxyz(3,l1,j,1,m) = zero
            exyz(1,l1,j,1,m) = zero
            exyz(2,l1,j,1,m) = zero
            exyz(3,l1,j,1,m) = zero
         endif
   70    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            afdt = adt*aimag(ffc(l,1,1,m))
c update magnetic field half time step
            zt2 = cmplx(-aimag(exyz(2,l,1,1,m)),real(exyz(2,l,1,1,m)))
            zt3 = cmplx(-aimag(exyz(1,l,1,1,m)),real(exyz(1,l,1,1,m)))
            zt4 = bxyz(1,l,1,1,m) + dth*(dkz*zt2)
            zt5 = bxyz(2,l,1,1,m) - dth*(dkz*zt3)
c update electric field whole time step
            zt2 = cmplx(-aimag(zt5),real(zt5))
            zt3 = cmplx(-aimag(zt4),real(zt4))
            zt7 = exyz(1,l,1,1,m) - cdt*(dkz*zt2) - afdt*cu(1,l,1,1,m)
            zt8 = exyz(2,l,1,1,m) + cdt*(dkz*zt3) - afdt*cu(2,l,1,1,m)
c update magnetic field half time step and store electric field
            zt2 = cmplx(-aimag(zt8),real(zt8))
            zt3 = cmplx(-aimag(zt7),real(zt7))
            exyz(1,l,1,1,m) = zt7
            exyz(2,l,1,1,m) = zt8
            exyz(3,l,1,1,m) = zero
            ws = ws + anorm*(zt7*conjg(zt7) + zt8*conjg(zt8))
            zt4 = zt4 + dth*(dkz*zt2)
            zt5 = zt5 - dth*(dkz*zt3)
            bxyz(1,l,1,1,m) = zt4
            bxyz(2,l,1,1,m) = zt5
            bxyz(3,l,1,1,m) = zero
            wp = wp + anorm*(zt4*conjg(zt4) + zt5*conjg(zt5))
            bxyz(1,l1,1,1,m) = zero
            bxyz(2,l1,1,1,m) = zero
            bxyz(3,l1,1,1,m) = zero
            exyz(1,l1,1,1,m) = zero
            exyz(2,l1,1,1,m) = zero
            exyz(3,l1,1,1,m) = zero
   80       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            bxyz(1,1,1,1,m) = zero
            bxyz(2,1,1,1,m) = zero
            bxyz(3,1,1,1,m) = zero
            exyz(1,1,1,1,m) = zero
            exyz(2,1,1,1,m) = zero
            exyz(3,1,1,1,m) = zero
            bxyz(1,l1,1,1,m) = zero
            bxyz(2,l1,1,1,m) = zero
            bxyz(3,l1,1,1,m) = zero
            exyz(1,l1,1,1,m) = zero
            exyz(2,l1,1,1,m) = zero
            exyz(3,l1,1,1,m) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = (nyh/kyzp)*kyzp
      if (n2.eq.k1) then
         k1 = nyh - k1 + 1
         do 100 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            bxyz(1,l,j,k1,m) = zero
            bxyz(2,l,j,k1,m) = zero
            bxyz(3,l,j,k1,m) = zero
            exyz(1,l,j,k1,m) = zero
            exyz(2,l,j,k1,m) = zero
            exyz(3,l,j,k1,m) = zero
   90       continue
         endif
  100    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 110 l = 1, nz
            bxyz(1,l,1,k1,m) = zero
            bxyz(2,l,1,k1,m) = zero
            bxyz(3,l,1,k1,m) = zero
            exyz(1,l,1,k1,m) = zero
            exyz(2,l,1,k1,m) = zero
            exyz(3,l,1,k1,m) = zero
  110       continue
         endif
      endif
  120 continue
  130 continue
  140 continue      
      wf = float(nx*ny*nz)*ws
      wm = float(nx*ny*nz)*c2*wp
      return
      end
c-----------------------------------------------------------------------
      subroutine PEMFIELD32(fxyz,exyz,ffc,isign,nx,ny,nz,kstrt,nzv,kxyp,
     1kyzp,jblok,mblok,nzhd)
c this subroutine either adds complex vector fields if isign > 0
c or copies complex vector fields if isign < 0
c includes additional smoothing
      implicit none
      integer isign, nx, ny, nz, kstrt, nzv, kxyp, kyzp, jblok, mblok
      integer nzhd
      complex fxyz, exyz, ffc
      dimension fxyz(3,nzv,kxyp,kyzp,jblok*mblok)
      dimension exyz(3,nzv,kxyp,kyzp,jblok*mblok)
      dimension ffc(nzhd,kxyp,kyzp,jblok*mblok)
      integer i, j, k, l, m, nxh, nzh, nz2, kxb, kyzb, l1
      real at1
      nxh = nx/2
      nzh = nz/2
      nz2 = nz + 2
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      if (kstrt.gt.(kxb*kyzb)) return
c add the fields
      if (isign.gt.0) then
         do 60 m = 1, jblok*mblok
c mode numbers 0 < kz < nz/2
         do 50 k = 1, kyzp
         do 40 j = 1, kxyp
         do 20 l = 2, nzh
         l1 = nz2 - l
         at1 = aimag(ffc(l,j,k,m))
         do 10 i = 1, 3
         fxyz(i,l,j,k,m) = fxyz(i,l,j,k,m) + exyz(i,l,j,k,m)*at1
         fxyz(i,l1,j,k,m) = fxyz(i,l1,j,k,m) + exyz(i,l1,j,k,m)*at1
   10    continue
   20    continue
c mode numbers kz = 0, nz/2
         l1 = nzh + 1
         at1 = aimag(ffc(1,j,k,m))
         do 30 i = 1, 3
         fxyz(i,1,j,k,m) = fxyz(i,1,j,k,m) + exyz(i,1,j,k,m)*at1
         fxyz(i,l1,j,k,m) = fxyz(i,l1,j,k,m) + exyz(i,l1,j,k,m)*at1
   30    continue
   40    continue
   50    continue
   60    continue
c copy the fields
      else if (isign.lt.0) then
         do 120 m = 1, jblok*mblok
c mode numbers 0 < kz < nz/2
         do 110 k = 1, kyzp
         do 100 j = 1, kxyp
         do 80 l = 2, nzh
         l1 = nz2 - l
         at1 = aimag(ffc(l,j,k,m))
         do 70 i = 1, 3
         fxyz(i,l,j,k,m) = exyz(i,l,j,k,m)*at1
         fxyz(i,l1,j,k,m) = exyz(i,l1,j,k,m)*at1
   70    continue
   80    continue
c mode numbers kz = 0, nz/2
         l1 = nzh + 1
         at1 = aimag(ffc(1,j,k,m))
         do 90 i = 1, 3
         fxyz(i,1,j,k,m) = exyz(i,1,j,k,m)*at1
         fxyz(i,l1,j,k,m) = exyz(i,l1,j,k,m)*at1
   90    continue
  100    continue
  110    continue
  120    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PAVPOT332(bxyz,axyz,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,
     1mblok)
c this subroutine calculates 3d vector potential from magnetic field
c in fourier space with periodic boundary conditions,
c for distributed data
c with 2D spatial decomposition
c input: bxyz,nx,ny,nz,kstrt,nzv,kxyp,kyzp,jblok,mblok
c output: axyz
c approximate flop count is:
c 99*nxc*nyc*nzc + 84*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx/2-1)/nvpy, nyc = (ny/2-1)/nvpz, nzc = nz/2 - 1, and
c nvpy/nvpz = number of procs in y/z
c the vector potential is calculated using the equations:
c ax(kx,ky,kz) = sqrt(-1)*
c                (ky*bz(kx,ky,kz)-kz*by(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c ay(kx,ky,kz) = sqrt(-1)*
c                (kz*bx(kx,ky,kz)-kx*bz(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c az(kx,ky,kz) = sqrt(-1)*
c                (kx*by(kx,ky,kz)-ky*bx(kx,ky,kz))/(kx*kx+ky*ky+kz*kz)
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz, and
c j,k,l = fourier mode numbers, except for
c ax(kx=pi) = ay(kx=pi) = az(kx=pi) = 0,
c ax(ky=pi) = ay(ky=pi) = ax(ky=pi) = 0,
c ax(kz=pi) = ay(kz=pi) = az(kz=pi) = 0,
c ax(kx=0,ky=0,kz=0) = ay(kx=0,ky=0,kz=0) = az(kx=0,ky=0,kz=0) = 0.
c bxyz(i,l,j,k,m) = i component of complex magnetic field
c axyz(i,l,j,k,m) = i component of complex vector potential
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*(mx - 1) and
c kk = k + kyzp*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c nx/ny/nz = system length in x/y/z direction
c jblok/mblok = number of field partitions in x/y
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c kstrt = starting data block number
c nzv = first dimension of field arrays, must be >= nz
      complex bxyz, axyz, zero, zt1, zt2, zt3
      dimension bxyz(3,nzv,kxyp,kyzp,jblok*mblok)
      dimension axyz(3,nzv,kxyp,kyzp,jblok*mblok)
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2
      nz2 = nz + 2
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      dnx = 6.28318530717959/float(nx)
      dny = 6.28318530717959/float(ny)
      dnz = 6.28318530717959/float(nz)
      zero = cmplx(0.,0.)
c calculate vector potential
      if (kstrt.gt.(kxb*kyzb)) go to 140
      do 130 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 120 mx = 1, jblok
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      joff = kxyp*(mx + js) - 1
      m = mx + moff
      do 50 k = 1, kyzp
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         dky = dny*float(k1)
         dky2 = dky*dky
         do 20 j = 1, kxyp
         dkx = dnx*float(j + joff)
         dkxy2 = dkx*dkx + dky2
         if ((j+joff).gt.0) then
            do 10 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            at1 = 1./(dkz*dkz + dkxy2)
            at2 = dkx*at1
            at3 = dky*at1
            at4 = dkz*at1
            zt1 = cmplx(-aimag(bxyz(3,l,j,k,m)),real(bxyz(3,l,j,k,m)))
            zt2 = cmplx(-aimag(bxyz(2,l,j,k,m)),real(bxyz(2,l,j,k,m)))
            zt3 = cmplx(-aimag(bxyz(1,l,j,k,m)),real(bxyz(1,l,j,k,m)))
            axyz(1,l,j,k,m) = at3*zt1 - at4*zt2
            axyz(2,l,j,k,m) = at4*zt3 - at2*zt1
            axyz(3,l,j,k,m) = at2*zt2 - at3*zt3
            zt1 = cmplx(-aimag(bxyz(3,l1,j,k,m)),real(bxyz(3,l1,j,k,m)))
            zt2 = cmplx(-aimag(bxyz(2,l1,j,k,m)),real(bxyz(2,l1,j,k,m)))
            zt3 = cmplx(-aimag(bxyz(1,l1,j,k,m)),real(bxyz(1,l1,j,k,m)))
            axyz(1,l1,j,k,m) = at3*zt1 + at4*zt2
            axyz(2,l1,j,k,m) = -at4*zt3 - at2*zt1
            axyz(3,l1,j,k,m) = at2*zt2 - at3*zt3
   10       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at1 = 1./dkxy2
            at2 = dkx*at1
            at3 = dky*at1
            zt1 = cmplx(-aimag(bxyz(3,1,j,k,m)),real(bxyz(3,1,j,k,m)))
            zt2 = cmplx(-aimag(bxyz(2,1,j,k,m)),real(bxyz(2,1,j,k,m)))
            zt3 = cmplx(-aimag(bxyz(1,1,j,k,m)),real(bxyz(1,1,j,k,m)))
            axyz(1,1,j,k,m) = at3*zt1
            axyz(2,1,j,k,m) = -at2*zt1
            axyz(3,1,j,k,m) = at2*zt2 - at3*zt3
            axyz(1,l1,j,k,m) = zero
            axyz(2,l1,j,k,m) = zero
            axyz(3,l1,j,k,m) = zero
         endif
   20    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 30 l = 2, nzh
               l1 = nz2 - l
               dkz = dnz*float(l - 1)
               at1 = 1./(dkz*dkz + dky2)
               at3 = dky*at1
               at4 = dkz*at1
               zt1 = cmplx(-aimag(bxyz(3,l,1,k,m)),real(bxyz(3,l,1,k,m))
     1)
               zt2 = cmplx(-aimag(bxyz(2,l,1,k,m)),real(bxyz(2,l,1,k,m))
     1)
               zt3 = cmplx(-aimag(bxyz(1,l,1,k,m)),real(bxyz(1,l,1,k,m))
     1)
               axyz(1,l,1,k,m) = at3*zt1 - at4*zt2
               axyz(2,l,1,k,m) = at4*zt3
               axyz(3,l,1,k,m) = -at3*zt3
               zt1 = cmplx(-aimag(bxyz(3,l1,1,k,m)),real(bxyz(3,l1,1,k,m
     1)))
               zt2 = cmplx(-aimag(bxyz(2,l1,1,k,m)),real(bxyz(2,l1,1,k,m
     1)))
               zt3 = cmplx(-aimag(bxyz(1,l1,1,k,m)),real(bxyz(1,l1,1,k,m
     1)))
               axyz(1,l1,1,k,m) = at3*zt1 + at4*zt2
               axyz(2,l1,1,k,m) = -at4*zt3
               axyz(3,l1,1,k,m) = -at3*zt3

   30          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               at3 = 1.0/dky
               zt1 = cmplx(-aimag(bxyz(3,1,1,k,m)),real(bxyz(3,1,1,k,m))
     1)
               zt3 = cmplx(-aimag(bxyz(1,1,1,k,m)),real(bxyz(1,1,1,k,m))
     1)
               axyz(1,1,1,k,m) = at3*zt1
               axyz(2,1,1,k,m) = zero
               axyz(3,1,1,k,m) = -at3*zt3
               axyz(1,l1,1,k,m) = zero
               axyz(2,l1,1,k,m) = zero
               axyz(3,l1,1,k,m) = zero

c throw away kx = nx/2
            else
               do 40 l = 1, nz
               axyz(1,l,1,k,m) = zero
               axyz(2,l,1,k,m) = zero
               axyz(3,l,1,k,m) = zero
   40          continue
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny/2
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 70 j = 1, kxyp
         dkx = dnx*float(j + joff)
         dkx2 = dkx*dkx
         if ((j+joff).gt.0) then
            do 60 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            at1 = 1./(dkz*dkz + dkx2)
            at2 = dkx*at1
            at4 = dkz*at1
            zt1 = cmplx(-aimag(bxyz(3,l,j,1,m)),real(bxyz(3,l,j,1,m)))
            zt2 = cmplx(-aimag(bxyz(2,l,j,1,m)),real(bxyz(2,l,j,1,m)))
            zt3 = cmplx(-aimag(bxyz(1,l,j,1,m)),real(bxyz(1,l,j,1,m)))
            axyz(1,l,j,1,m) = -at4*zt2
            axyz(2,l,j,1,m) = at4*zt3 - at2*zt1
            axyz(3,l,j,1,m) = at2*zt2
            zt1 = cmplx(-aimag(bxyz(3,l1,j,1,m)),real(bxyz(3,l1,j,1,m)))
            zt2 = cmplx(-aimag(bxyz(2,l1,j,1,m)),real(bxyz(2,l1,j,1,m)))
            zt3 = cmplx(-aimag(bxyz(1,l1,j,1,m)),real(bxyz(1,l1,j,1,m)))
            axyz(1,l1,j,1,m) = at4*zt2
            axyz(2,l1,j,1,m) = -at4*zt3 - at2*zt1
            axyz(3,l1,j,1,m) = at2*zt2
   60       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            at2 = 1.0/dkx
            zt1 = cmplx(-aimag(bxyz(3,1,j,1,m)),real(bxyz(3,1,j,1,m)))
            zt2 = cmplx(-aimag(bxyz(2,1,j,1,m)),real(bxyz(2,1,j,1,m)))
            axyz(1,1,j,1,m) = zero
            axyz(2,1,j,1,m) = -at2*zt1
            axyz(3,1,j,1,m) = at2*zt2
            axyz(1,l1,j,1,m) = zero
            axyz(2,l1,j,1,m) = zero
            axyz(3,l1,j,1,m) = zero
         endif
   70    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 80 l = 2, nzh
            l1 = nz2 - l
            dkz = dnz*float(l - 1)
            at4 = 1.0/dkz
            zt2 = cmplx(-aimag(bxyz(2,l,1,1,m)),real(bxyz(2,l,1,1,m)))
            zt3 = cmplx(-aimag(bxyz(1,l,1,1,m)),real(bxyz(1,l,1,1,m)))
            axyz(1,l,1,1,m) = -at4*zt2
            axyz(2,l,1,1,m) = at4*zt3
            axyz(3,l,1,1,m) = zero
            axyz(1,l1,1,1,m) = zero
            axyz(2,l1,1,1,m) = zero
            axyz(3,l1,1,1,m) = zero
   80       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            axyz(1,1,1,1,m) = zero
            axyz(2,1,1,1,m) = zero
            axyz(3,1,1,1,m) = zero
            axyz(1,l1,1,1,m) = zero
            axyz(2,l1,1,1,m) = zero
            axyz(3,l1,1,1,m) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = (nyh/kyzp)*kyzp
      if (n2.eq.k1) then
         k1 = nyh - k1 + 1
         do 100 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 90 l = 1, nz
            axyz(1,l,j,k1,m) = zero
            axyz(2,l,j,k1,m) = zero
            axyz(3,l,j,k1,m) = zero
   90       continue
         endif
  100    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 110 l = 1, nz
            axyz(1,l,1,k1,m) = zero
            axyz(2,l,1,k1,m) = zero
            axyz(3,l,1,k1,m) = zero
  110       continue
         endif
      endif
  120 continue
  130 continue
  140 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGTMODES32(pot,pott,nx,ny,nz,it,modesx,modesy,modesz,ks
     1trt,nzv,kxyp,kyzp,jblok,mblok,nt,modesxpd,modesypd,modeszd)
c this subroutine extracts lowest order modes from complex array pot
c and stores them into a location in a time history array pott
c modes stored: kx = (kxyp*(idprocx)+(0,1,...kxyp-1)),
c and ky = (kyzp*(idprocy)+(0,1,...kyzp-1)), when idprocy < nvpy/2
c and ky = (kyzp*(idprocy-nvpy+1)-(kyzp-1,...,1,0)-1),
c when idprocy >= nvpy/2
c where idprocy=idproc/nvpy, and idprocx = idprocx-nvpy*idprocy
c and kz=(0,+-1,+-2,...,+-(NZ/2-1),NZ/2)
c except kx = NX/2 is stored at location kxyp+1 when idprocx=0,
c and ky = NY/2 is stored at location 1 when idprocy=nvp/2.
c nx/ny/nz = system length in x/y/z direction
c it = current time
c modesx/modesy/modesz = number of modes to store in x/y/z direction,
c where modesx <= nx/2+1, modesy <= ny/2+1, modesz <= nz/2+1
c kstrt = starting data block number
c nzv = first dimension of field arrays, must be >= nz
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c jblok/mblok = number of field partitions in x/y
c nt = first dimension of output array pott, nt >= it
c modeszd = second dimension of output array pott, modeszd  = 2*modesz
c modesypd = fourth dimension of array pott, modesypd >= min(modesy,kyzp)
c modesxpd = third dimension of array pott,
c modesxpd >= min(modesx,kxyp), unless modesx = nx/2+1,
c in which case modesxpd = kxyp+1
      implicit none
      integer nx, ny, nz, it, modesx, modesy ,modesz, kstrt, nzv
      integer kxyp, kyzp, jblok, mblok, nt, modesxpd, modesypd, modeszd
      complex pot, pott
      dimension pot(nzv,kxyp,kyzp,jblok*mblok)
      dimension pott(nt,modeszd,modesxpd,modesypd,jblok*mblok)
      integer nxh, nyh, nzh, jmax, kmax, lmax, nz2, kxb, kyzb, js, ks
      integer j, k, l, m, mx, my, j1, k1, l1, n1, n2, joff, koff, moff
      integer kmin
      if (it.gt.nt) return
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      if ((modesz.le.0).or.(modesz.gt.(nzh+1))) return
      nz2 = nz + 2
      kmax = min0(modesy,nyh)
      kmin = ny - kmax
      lmax = min0(modesz,nzh)
      j1 = kxyp + 1
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      if (kstrt.gt.(kxb*kyzb)) go to 160
      do 150 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 140 mx = 1, jblok
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      joff = kxyp*(mx + js)
      jmax = modesx - joff
      if (jmax.gt.kxyp) then
         jmax = kxyp
      else if (jmax.le.0) then
         jmax = 0
      endif
      m = mx + moff
      do 50 k = 1, kyzp
      k1 = k + koff
      if (((k1.gt.0).and.(k1.lt.kmax)).or.((k1.gt.kmin).and.(k1.lt.ny)))
     1 then
         if (k1.gt.nyh) k1 = k1 - ny
         do 20 j = 1, jmax
         if ((j+joff).gt.1) then
            do 10 l = 2, lmax
            l1 = nz2 - l
            pott(it,2*l-2,j,k,m) = pot(l,j,k,m)
            pott(it,2*l-1,j,k,m) = pot(l1,j,k,m)
   10       continue
c mode numbers kz = 0, nz/2
            pott(it,1,j,k,m) = pot(1,j,k,m)
            if (modesz.gt.nzh) then
               l1 = nzh + 1
               pott(it,nz,j,k,m) = pot(l1,j,k,m)
            endif
         endif
   20    continue
c mode numbers kx = 0, nx/2
         n1 = joff
         if (n1.eq.0) then
c kx = 0
            if (k1.gt.0) then
               do 30 l = 2, lmax
               l1 = nz2 - l
               pott(it,2*l-2,1,k,m) = pot(l,1,k,m)
               pott(it,2*l-1,1,k,m) = pot(l1,1,k,m)
   30          continue
c mode numbers kz = 0, nz/2
               pott(it,1,1,k,m) = pot(1,1,k,m)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pott(it,nz,1,k,m) = pot(l1,1,k,m)
               endif
c kx = nx/2
            else
               if (modesx.gt.nxh) then
                  do 40 l = 2, lmax
                  l1 = nz2 - l
                  pott(it,2*l-2,1,k,m) = pot(l,1,k,m)
                  pott(it,2*l-1,1,k,m) = pot(l1,1,k,m)
   40             continue
c mode numbers kz = 0, nz/2
                  pott(it,1,1,k,m) = pot(1,1,k,m)
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     pott(it,nz,1,k,m) = pot(l1,1,k,m)
                  endif
               endif
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny/2
      n2 = koff + 1
c ky = 0
      if (n2.eq.0) then
         do 70 j = 1, jmax
         if ((j+joff).gt.1) then
            do 60 l = 2, lmax
            l1 = nz2 - l
            pott(it,2*l-2,j,1,m) = pot(l,j,1,m)
            pott(it,2*l-1,j,1,m) = pot(l1,j,1,m)
   60       continue
c mode numbers kz = 0, nz/2
            pott(it,1,j,1,m) = pot(1,j,1,m)
            if (modesz.gt.nzh) then
               l1 = nzh + 1
               pott(it,nz,j,1,m) = pot(l1,j,1,m)
            endif
         endif
   70    continue
c mode numbers kx = 0, nx/2
         n1 = joff
         if (n1.eq.0) then
c kx = 0
            do 80 l = 2, lmax
            pott(it,2*l-2,1,1,m) = pot(l,1,1,m)
            pott(it,2*l-1,1,1,m) = conjg(pot(l,1,1,m))
   80       continue
c mode numbers kz = 0, nz/2
            pott(it,1,1,1,m) = cmplx(real(pot(1,1,1,m)),0.)
            if (modesz.gt.nzh) then
               l1 = nzh + 1
               pott(it,nz,1,1,m) = cmplx(real(pot(l1,1,1,m)),0.)
            endif
c kx = nx/2
            if (modesx.gt.nxh) then
               do 90 l = 2, lmax
               l1 = nz2 - l
               pott(it,2*l-2,j1,1,m) = conjg(pot(l1,1,1,m))
               pott(it,2*l-1,j1,1,m) = pot(l1,1,1,m)
   90          continue
c mode numbers kz = 0, nz/2
               pott(it,1,j1,1,m) = cmplx(aimag(pot(1,1,1,m)),0.)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pott(it,nz,j1,1,m) = cmplx(aimag(pot(l1,1,1,m)),0.)
               endif
            endif
         endif
      endif
c ky = ny/2
      k1 = (nyh/kyzp)*kyzp
      if (n2.eq.k1) then
         k1 = nyh - k1 + 1
         if (modesy.gt.nyh) then
            do 110 j = 1, jmax
            if ((j+joff).gt.1) then
               do 100 l = 2, lmax
               l1 = nz2 - l
               pott(it,2*l-2,j,k1,m) = pot(l,j,k1,m)
               pott(it,2*l-1,j,k1,m) = pot(l1,j,k1,m)
  100          continue
c mode numbers kz = 0, nz/2
               pott(it,1,j,k1,m) = pot(1,j,k1,m)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pott(it,nz,j,k1,m) = pot(l1,j,k1,m)
               endif
            endif
  110       continue
c mode numbers kx = 0, nx/2
            n1 = joff
            if (n1.eq.0) then
c kx = 0
               do 120 l = 2, lmax
               pott(it,2*l-2,1,k1,m) = pot(l,1,k1,m)
               pott(it,2*l-1,1,k1,m) = conjg(pot(l,1,k1,m))
  120          continue
c mode numbers kz = 0, nz/2
               pott(it,1,1,k1,m) = cmplx(real(pot(1,1,k1,m)),0.)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pott(it,nz,1,k1,m) = cmplx(real(pot(l1,1,k1,m)),0.)
               endif
c kx  = nx/2
               if (modesx.gt.nxh) then
                  do 130 l = 2, lmax
                  l1 = nz2 - l
                  pott(it,2*l-2,j1,k1,m) = conjg(pot(l1,1,k1,m))
                  pott(it,2*l-1,j1,k1,m) = pot(l1,1,k1,m)
  130             continue
c mode numbers kz = 0, nz/2
                  pott(it,1,j1,k1,m) = cmplx(aimag(pot(1,1,k1,m)),0.)
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     pott(it,nz,j1,k1,m) = cmplx(aimag(pot(l1,1,k1,m)),0
     1.)
                  endif
               endif
            endif
         endif
      endif
  140 continue
  150 continue
  160 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPTMODES32(pot,pott,nx,ny,nz,it,modesx,modesy,modesz,ks
     1trt,nzv,kxyp,kyzp,jblok,mblok,nt,modesxpd,modesypd,modeszd)
c this subroutine extracts lowest order modes from a location in a time
c history array pott and stores them into complex array pot
c modes stored: kx = (kxyp*(idprocx)+(0,1,...kxyp-1)),
c and ky = (kyzp*(idprocy)+(0,1,...kyzp-1)), when idprocy < nvpy/2
c and ky = (kyzp*(idprocy-nvpy+1)-(kyzp-1,...,1,0)-1),
c when idprocy >= nvpy/2
c where idprocy=idproc/nvpy, and idprocx = idprocx-nvpy*idprocy
c and kz=(0,+-1,+-2,...,+-(NZ/2-1),NZ/2)
c except kx = NX/2 is stored at location kxyp+1 when idprocx=0,
c and ky = NY/2 is stored at location 1 when idprocy=nvp/2.
c nx/ny/nz = system length in x/y/z direction
c it = current time
c modesx/modesy/modesz = number of modes to store in x/y/z direction,
c where modesx <= nx/2+1, modesy <= ny/2+1, modesz <= nz/2+1
c kstrt = starting data block number
c nzv = first dimension of field arrays, must be >= nz
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c jblok/mblok = number of field partitions in x/y
c nt = first dimension of output array pott, nt >= it
c modeszd = second dimension of output array pott, modeszd  = 2*modesz
c modesypd = fourth dimension of array pott, modesypd >= min(modesy,kyzp)
c modesxpd = third dimension of array pott,
c modesxpd >= min(modesx,kxyp), unless modesx = nx/2+1,
c in which case modesxpd = kxyp+1
      implicit none
      integer nx, ny, nz, it, modesx, modesy ,modesz, kstrt, nzv
      integer kxyp, kyzp, jblok, mblok, nt, modesxpd, modesypd, modeszd
      complex pot, pott
      dimension pot(nzv,kxyp,kyzp,jblok*mblok)
      dimension pott(nt,modeszd,modesxpd,modesypd,jblok*mblok)
      integer nxh, nyh, nzh, jmax, kmax, lmax, nz2, kxb, kyzb, js, ks
      integer j, k, l, m, mx, my, j1, k1, l1, n1, n2, joff, koff, moff
      integer kmin
      complex zero
      if (it.gt.nt) return
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      if ((modesz.le.0).or.(modesz.gt.(nzh+1))) return
      nz2 = nz + 2
      kmax = min0(modesy,nyh)
      kmin = ny - kmax
      lmax = min0(modesz,nzh)
      j1 = kxyp + 1
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      zero = cmplx(0.,0.)
      if (kstrt.gt.(kxb*kyzb)) go to 340
      do 330 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 320 mx = 1, jblok
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      joff = kxyp*(mx + js)
      jmax = modesx - joff
      if (jmax.gt.kxyp) then
         jmax = kxyp
      else if (jmax.le.0) then
         jmax = 0
      endif
      m = mx + moff
      do 100 k = 1, kyzp
      k1 = k + koff
      if (((k1.gt.0).and.(k1.lt.kmax)).or.((k1.gt.kmin).and.(k1.lt.ny)))
     1 then
         if (k1.gt.nyh) k1 = k1 - ny
         do 30 j = 1, jmax
         if ((j+joff).gt.1) then
            do 10 l = 2, lmax
            l1 = nz2 - l
            pot(l,j,k,m) = pott(it,2*l-2,j,k,m)
            pot(l1,j,k,m) = pott(it,2*l-1,j,k,m)
   10       continue
            do 20 l = lmax+1, nzh
            l1 = nz2 - l
            pot(l,j,k,m) = zero
            pot(l1,j,k,m) = zero
   20       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1 
            pot(1,j,k,m) = pott(it,1,j,k,m)
            pot(l1,j,k,m) = zero
            if (modesz.gt.nzh) then
               pot(l1,j,k,m) = pott(it,nz,j,k,m)
            endif
         endif
   30    continue
         do 50 j = jmax+1, kxyp
         if ((j+joff).gt.1) then
            do 40 l = 1, nz
            pot(l,j,k,m) = zero
   40       continue
         endif
   50    continue
c mode numbers kx = 0, nx/2
         n1 = joff
         if (n1.eq.0) then
c kx = 0
            if (k1.gt.0) then
               do 60 l = 2, lmax
               l1 = nz2 - l
               pot(l,1,k,m) = pott(it,2*l-2,1,k,m)
               pot(l1,1,k,m) = pott(it,2*l-1,1,k,m)
   60          continue
               do 70 l = lmax+1, nzh
               l1 = nz2 - l
               pot(l,1,k,m) = zero
               pot(l1,1,k,m) = zero
   70          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               pot(1,1,k,m) = pott(it,1,1,k,m)
               pot(l1,1,k,m) = zero
               if (modesz.gt.nzh) then
                  pot(l1,1,k,m) = pott(it,nz,1,k,m)
               endif
c kx = nx/2
            else
               do 80 l = 1, nz
               pot(l,1,k,m) = zero
   80          continue
               if (modesx.gt.nxh) then
                  do 90 l = 2, lmax
                  l1 = nz2 - l
                  pot(l,1,k,m) = pott(it,2*l-2,1,k,m)
                  pot(l1,1,k,m) = pott(it,2*l-1,1,k,m)
   90             continue
c mode numbers kz = 0, nz/2
                  pot(1,1,k,m) = pott(it,1,1,k,m)
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     pot(l1,1,k,m) = pott(it,nz,1,k,m)
                  endif
               endif
            endif
         endif
      endif
  100 continue
      do 140 k = 1, kyzp
      k1 = k + koff
      if ((k1.ge.kmax).and.(k1.le.kmin).and.(k1.ne.nyh)) then
         do 120 j = 1, kxyp
         if ((j+joff).gt.1) then
            do 110 l = 1, nz
            pot(l,j,k,m) = zero
  110       continue
         endif
  120    continue
c mode numbers kx = 0, nx/2
         n1 = joff
         if (n1.eq.0) then
            do 130 l = 1, nz
            pot(l,1,k,m) = zero
  130       continue
         endif
      endif
  140 continue
c mode numbers ky = 0, ny/2
      n2 = koff + 1
c ky = 0
      if (n2.eq.0) then
         do 170 j = 1, jmax
         if ((j+joff).gt.1) then
            do 150 l = 2, lmax
            l1 = nz2 - l
            pot(l,j,1,m) = pott(it,2*l-2,j,1,m)
            pot(l1,j,1,m) = pott(it,2*l-1,j,1,m)
  150       continue
            do 160 l = lmax+1, nzh
            l1 = nz2 - l
            pot(l,j,1,m) = zero
            pot(l1,j,1,m) = zero
  160       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            pot(1,j,1,m) = pott(it,1,j,1,m)
            pot(l1,j,1,m) = zero
            if (modesz.gt.nzh) then
               pot(l1,j,1,m) = pott(it,nz,j,1,m)
            endif
         endif
  170    continue
         do 190 j = jmax+1, kxyp
         if ((j+joff).gt.1) then
            do 180 l = 1, nz
            pot(l,j,1,m) = zero
  180       continue
         endif
  190    continue
c mode numbers kx = 0, nx/2
         n1 = joff
         if (n1.eq.0) then
c kx = 0
            do 200 l = 2, lmax
            l1 = nz2 - l
            pot(l,1,1,m) = pott(it,2*l-2,1,1,m)
            pot(l1,1,1,m) = zero
  200       continue
            do 210 l = lmax+1, nzh
            l1 = nz2 - l
            pot(l,1,1,m) = zero
            pot(l1,1,1,m) = zero
  210       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            pot(1,1,1,m) = cmplx(real(pott(it,1,1,1,m)),0.)
            pot(l1,1,1,m) = zero
            if (modesz.gt.nzh) then
               pot(l1,1,1,m) = cmplx(real(pott(it,nz,1,1,m)),0.)
            endif
c kx = nx/2
            if (modesx.gt.nxh) then
               do 220 l = 2, lmax
               l1 = nz2 - l
               pot(l1,1,1,m) = conjg(pott(it,2*l-2,j1,1,m))
  220          continue
c mode numbers kz = 0, nz/2
               pot(1,1,1,m) = cmplx(real(pot(1,1,1,m)),real(pott(it,1,j1
     1,1,m)))
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pot(l1,1,1,m) = cmplx(real(pot(l1,1,1,m)),real(pott(it
     1,nz,j1,1,m)))
               endif
            endif
         endif
      endif
c ky = ny/2
      k1 = (nyh/kyzp)*kyzp
      if (n2.eq.k1) then
         k1 = nyh - k1 + 1
         do 240 j = 1, jmax
         if ((j+joff).gt.1) then
            do 230 l = 1, nz
            pot(l,j,k1,m) = zero
  230       continue
         endif
  240    continue
         do 260 j = jmax+1, kxyp
         if ((j+joff).gt.1) then
            do 250 l = 1, nz
            pot(l,j,k1,m) = zero
  250       continue
         endif
  260    continue
c mode numbers kx = 0, nx/2
         n1 = joff
         if (n1.eq.0) then
            do 270 l = 1, nz
            pot(l,1,k1,m) = zero
  270       continue
         endif
         if (modesy.gt.nyh) then
            do 290 j = 1, jmax
            if ((j+joff).gt.1) then
               do 280 l = 2, lmax
               l1 = nz2 - l
               pot(l,j,k1,m) = pott(it,2*l-2,j,k1,m)
               pot(l1,j,k1,m) = pott(it,2*l-1,j,k1,m)
  280          continue
c mode numbers kz = 0, nz/2
               pot(1,j,k1,m) = pott(it,1,j,k1,m)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pot(l1,j,k1,m) = pott(it,nz,j,k1,m)
               endif
            endif
  290       continue
c mode numbers kx = 0, nx/2
            n1 = joff
            if (n1.eq.0) then
c kx = 0
               do 300 l = 2, lmax
               pot(l,1,k1,m) = pott(it,2*l-2,1,k1,m)
  300          continue
c mode numbers kz = 0, nz/2
               pot(1,1,k1,m) = cmplx(real(pott(it,1,1,k1,m)),0.)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  pot(l1,1,k1,m) = cmplx(real(pott(it,nz,1,k1,m)),0.)
               endif
c kx  = nx/2
               if (modesx.gt.nxh) then
                  do 310 l = 2, lmax
                  l1 = nz2 - l
                  pot(l1,1,k1,m) = conjg(pott(it,2*l-2,j1,k1,m))
  310             continue
c mode numbers kz = 0, nz/2
                  pot(1,1,k1,m) = cmplx(real(pot(1,1,k1,m)),real(pott(it
     1,1,j1,k1,m)))

                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     pot(l1,1,k1,m) = cmplx(real(pot(l1,1,k1,m)),real(po
     1tt(it,nz,j1,k1,m)))
                  endif
               endif
            endif
         endif
      endif
  320 continue
  330 continue
  340 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PGTVMODES32(vpot,vpott,nx,ny,nz,it,modesx,modesy,modesz
     1,ndim,kstrt,nzv,kxyp,kyzp,jblok,mblok,nt,modesxpd,modesypd,modeszd
     2)
c this subroutine extracts lowest order modes from complex vector array
c vpot and stores them into a location in a time history array vpott
c modes stored: kx = (kxyp*(idprocx)+(0,1,...kxyp-1)),
c and ky = (kyzp*(idprocy)+(0,1,...kyzp-1)), when idprocy < nvpy/2
c and ky = (kyzp*(idprocy-nvpy+1)-(kyzp-1,...,1,0)-1),
c when idprocy >= nvpy/2
c where idprocy=idproc/nvpy, and idprocx = idprocx-nvpy*idprocy
c and kz=(0,+-1,+-2,...,+-(NZ/2-1),NZ/2)
c except kx = NX/2 is stored at location kxyp+1 when idprocx=0,
c and ky = NY/2 is stored at location 1 when idprocy=nvp/2.
c nx/ny/nz = system length in x/y/z direction
c it = current time
c modesx/modesy/modesz = number of modes to store in x/y/z direction,
c where modesx <= nx/2+1, modesy <= ny/2+1, modesz <= nz/2+1
c ndim = number of field arrays, must be >= 1
c kstrt = starting data block number
c nzv = first dimension of field arrays, must be >= nz
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c jblok/mblok = number of field partitions in x/y
c nt = first dimension of output array vpott, nt >= it
c modeszd = third dimension of output array vpott, modeszd  = 2*modesz
c modesypd = fifth dimension of array vpott,
c modesypd >= min(modesy,kyzp)
c modesxpd = fourth dimension of array vpott,
c modesxpd >= min(modesx,kxyp), unless modesx = nx/2+1,
c in which case modesxpd = kxyp+1
      implicit none
      integer nx, ny, nz, it, modesx, modesy ,modesz, ndim, kstrt, nzv
      integer kxyp, kyzp, jblok, mblok, nt, modesxpd, modesypd, modeszd
      complex vpot, vpott
      dimension vpot(ndim,nzv,kxyp,kyzp,jblok*mblok)
      dimension vpott(nt,ndim,modeszd,modesxpd,modesypd,jblok*mblok)
      integer nxh, nyh, nzh, jmax, kmax, lmax, nz2, kxb, kyzb, js, ks
      integer j, k, l, m, mx, my, j1, k1, l1, n1, n2, joff, koff, moff
      integer i, kmin
      if (it.gt.nt) return
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      if ((modesz.le.0).or.(modesz.gt.(nzh+1))) return
      nz2 = nz + 2
      kmax = min0(modesy,nyh)
      kmin = ny - kmax
      lmax = min0(modesz,nzh)
      j1 = kxyp + 1
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      if (kstrt.gt.(kxb*kyzb)) go to 340
      do 330 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 320 mx = 1, jblok
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      joff = kxyp*(mx + js)
      jmax = modesx - joff
      if (jmax.gt.kxyp) then
         jmax = kxyp
      else if (jmax.le.0) then
         jmax = 0
      endif
      m = mx + moff
      do 110 k = 1, kyzp
      k1 = k + koff
      if (((k1.gt.0).and.(k1.lt.kmax)).or.((k1.gt.kmin).and.(k1.lt.ny)))
     1 then
         if (k1.gt.nyh) k1 = k1 - ny
         do 40 j = 1, jmax
         if ((j+joff).gt.1) then
            do 20 l = 2, lmax
            l1 = nz2 - l
            do 10 i = 1, ndim
            vpott(it,i,2*l-2,j,k,m) = vpot(i,l,j,k,m)
            vpott(it,i,2*l-1,j,k,m) = vpot(i,l1,j,k,m)
   10       continue
   20       continue
c mode numbers kz = 0, nz/2
            do 30 i = 1, ndim
            vpott(it,i,1,j,k,m) = vpot(i,1,j,k,m)
            if (modesz.gt.nzh) then
               l1 = nzh + 1
               vpott(it,i,nz,j,k,m) = vpot(i,l1,j,k,m)
            endif
   30       continue
         endif
   40    continue
c mode numbers kx = 0, nx/2
         n1 = joff
         if (n1.eq.0) then
c kx = 0
            if (k1.gt.0) then
               do 60 l = 2, lmax
               l1 = nz2 - l
               do 50 i = 1, ndim
               vpott(it,i,2*l-2,1,k,m) = vpot(i,l,1,k,m)
               vpott(it,i,2*l-1,1,k,m) = vpot(i,l1,1,k,m)
   50          continue
   60          continue
c mode numbers kz = 0, nz/2
               do 70 i = 1, ndim
               vpott(it,i,1,1,k,m) = vpot(i,1,1,k,m)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpott(it,i,nz,1,k,m) = vpot(i,l1,1,k,m)
               endif
   70          continue
c kx = nx/2
            else
               if (modesx.gt.nxh) then
                  do 90 l = 2, lmax
                  l1 = nz2 - l
                  do 80 i = 1, ndim
                  vpott(it,i,2*l-2,1,k,m) = vpot(i,l,1,k,m)
                  vpott(it,i,2*l-1,1,k,m) = vpot(i,l1,1,k,m)
   80             continue
   90             continue
c mode numbers kz = 0, nz/2
                  do 100 i = 1, ndim
                  vpott(it,i,1,1,k,m) = vpot(i,1,1,k,m)
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     vpott(it,i,nz,1,k,m) = vpot(i,l1,1,k,m)
                  endif
  100             continue
               endif
            endif
         endif
      endif
  110 continue
c mode numbers ky = 0, ny/2
      n2 = koff + 1
c ky = 0
      if (n2.eq.0) then
         do 150 j = 1, jmax
         if ((j+joff).gt.1) then
            do 130 l = 2, lmax
            l1 = nz2 - l
            do 120 i = 1, ndim
            vpott(it,i,2*l-2,j,1,m) = vpot(i,l,j,1,m)
            vpott(it,i,2*l-1,j,1,m) = vpot(i,l1,j,1,m)
  120       continue
  130       continue
c mode numbers kz = 0, nz/2
            do 140 i = 1, ndim
            vpott(it,i,1,j,1,m) = vpot(i,1,j,1,m)
            if (modesz.gt.nzh) then
               l1 = nzh + 1
               vpott(it,i,nz,j,1,m) = vpot(i,l1,j,1,m)
            endif
  140       continue
         endif
  150    continue
c mode numbers kx = 0, nx/2
         n1 = joff
         if (n1.eq.0) then
c kx = 0
            do 170 l = 2, lmax
            do 160 i = 1, ndim
            vpott(it,i,2*l-2,1,1,m) = vpot(i,l,1,1,m)
            vpott(it,i,2*l-1,1,1,m) = conjg(vpot(i,l,1,1,m))
  160       continue
  170       continue
c mode numbers kz = 0, nz/2
            do 180 i = 1, ndim
            vpott(it,i,1,1,1,m) = cmplx(real(vpot(i,1,1,1,m)),0.)
            if (modesz.gt.nzh) then
               l1 = nzh + 1
               vpott(it,i,nz,1,1,m) = cmplx(real(vpot(i,l1,1,1,m)),0.)
            endif
  180       continue
c kx = nx/2
            if (modesx.gt.nxh) then
               do 200 l = 2, lmax
               l1 = nz2 - l
               do 190 i = 1, ndim
               vpott(it,i,2*l-2,j1,1,m) = conjg(vpot(i,l1,1,1,m))
               vpott(it,i,2*l-1,j1,1,m) = vpot(i,l1,1,1,m)
  190          continue
  200          continue
c mode numbers kz = 0, nz/2
               do 210 i = 1, ndim
               vpott(it,i,1,j1,1,m) = cmplx(aimag(vpot(i,1,1,1,m)),0.)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpott(it,i,nz,j1,1,m) = cmplx(aimag(vpot(i,l1,1,1,m)),
     10.)
               endif
  210          continue
            endif
         endif
      endif
c ky = ny/2
      k1 = (nyh/kyzp)*kyzp
      if (n2.eq.k1) then
         k1 = nyh - k1 + 1
         if (modesy.gt.nyh) then
            do 250 j = 1, jmax
            if ((j+joff).gt.1) then
               do 230 l = 2, lmax
               l1 = nz2 - l
               do 220 i = 1, ndim
               vpott(it,i,2*l-2,j,k1,m) = vpot(i,l,j,k1,m)
               vpott(it,i,2*l-1,j,k1,m) = vpot(i,l1,j,k1,m)
  220          continue
  230          continue
c mode numbers kz = 0, nz/2
               do 240 i = 1, ndim
               vpott(it,i,1,j,k1,m) = vpot(i,1,j,k1,m)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpott(it,i,nz,j,k1,m) = vpot(i,l1,j,k1,m)
               endif
  240          continue
            endif
  250       continue
c mode numbers kx = 0, nx/2
            n1 = joff
            if (n1.eq.0) then
c kx = 0
               do 270 l = 2, lmax
               do 260 i = 1, ndim
               vpott(it,i,2*l-2,1,k1,m) = vpot(i,l,1,k1,m)
               vpott(it,i,2*l-1,1,k1,m) = conjg(vpot(i,l,1,k1,m))
  260          continue
  270          continue
c mode numbers kz = 0, nz/2
               do 280 i = 1, ndim
               vpott(it,i,1,1,k1,m) = cmplx(real(vpot(i,1,1,k1,m)),0.)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpott(it,i,nz,1,k1,m) = cmplx(real(vpot(i,l1,1,k1,m)),
     10.)
               endif
  280          continue
c kx  = nx/2
               if (modesx.gt.nxh) then
                  do 300 l = 2, lmax
                  l1 = nz2 - l
                  do 290 i = 1, ndim
                  vpott(it,i,2*l-2,j1,k1,m) = conjg(vpot(i,l1,1,k1,m))
                  vpott(it,i,2*l-1,j1,k1,m) = vpot(i,l1,1,k1,m)
  290             continue
  300             continue
c mode numbers kz = 0, nz/2
                  do 310 i = 1, ndim
                  vpott(it,i,1,j1,k1,m) = cmplx(aimag(vpot(i,1,1,k1,m)),
     10.)
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     vpott(it,i,nz,j1,k1,m) = cmplx(aimag(vpot(i,l1,1,k1
     1,m)),0.)
                  endif
  310             continue
               endif
            endif
         endif
      endif
  320 continue
  330 continue
  340 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPTVMODES32(vpot,vpott,nx,ny,nz,it,modesx,modesy,modesz
     1,ndim,kstrt,nzv,kxyp,kyzp,jblok,mblok,nt,modesxpd,modesypd,modeszd
     2)
c this subroutine extracts lowest order modes from a location in a time
c history array vpott and stores them into complex vector array vpot
c modes stored: kx = (kxyp*(idprocx)+(0,1,...kxyp-1)),
c and ky = (kyzp*(idprocy)+(0,1,...kyzp-1)), when idprocy < nvpy/2
c and ky = (kyzp*(idprocy-nvpy+1)-(kyzp-1,...,1,0)-1),
c when idprocy >= nvpy/2
c where idprocy=idproc/nvpy, and idprocx = idprocx-nvpy*idprocy
c and kz=(0,+-1,+-2,...,+-(NZ/2-1),NZ/2)
c except kx = NX/2 is stored at location kxyp+1 when idprocx=0,
c and ky = NY/2 is stored at location 1 when idprocy=nvp/2.
c nx/ny/nz = system length in x/y/z direction
c it = current time
c modesx/modesy/modesz = number of modes to store in x/y/z direction,
c where modesx <= nx/2+1, modesy <= ny/2+1, modesz <= nz/2+1
c ndim = number of field arrays, must be >= 1
c kstrt = starting data block number
c nzv = first dimension of field arrays, must be >= nz
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c jblok/mblok = number of field partitions in x/y
c nt = first dimension of output array vpott, nt >= it
c modeszd = third dimension of output array vpott, modeszd  = 2*modesz
c modesypd = fifth dimension of array vpott,
c modesypd >= min(modesy,kyzp)
c modesxpd = fourth dimension of array vpott,
c modesxpd >= min(modesx,kxyp), unless modesx = nx/2+1,
c in which case modesxpd = kxyp+1
      implicit none
      integer nx, ny, nz, it, modesx, modesy ,modesz, ndim, kstrt, nzv
      integer kxyp, kyzp, jblok, mblok, nt, modesxpd, modesypd, modeszd
      complex vpot, vpott
      dimension vpot(ndim,nzv,kxyp,kyzp,jblok*mblok)
      dimension vpott(nt,ndim,modeszd,modesxpd,modesypd,jblok*mblok)
      integer nxh, nyh, nzh, jmax, kmax, lmax, nz2, kxb, kyzb, js, ks
      integer j, k, l, m, mx, my, j1, k1, l1, n1, n2, joff, koff, moff
      integer i, kmin
      complex zero
      if (it.gt.nt) return
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2
      if ((modesx.le.0).or.(modesx.gt.(nxh+1))) return
      if ((modesy.le.0).or.(modesy.gt.(nyh+1))) return
      if ((modesz.le.0).or.(modesz.gt.(nzh+1))) return
      nz2 = nz + 2
      kmax = min0(modesy,nyh)
      kmin = ny - kmax
      lmax = min0(modesz,nzh)
      j1 = kxyp + 1
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      zero = cmplx(0.,0.)
      if (kstrt.gt.(kxb*kyzb)) go to 640
      do 630 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 620 mx = 1, jblok
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      joff = kxyp*(mx + js)
      jmax = modesx - joff
      if (jmax.gt.kxyp) then
         jmax = kxyp
      else if (jmax.le.0) then
         jmax = 0
      endif
      m = mx + moff
      do 200 k = 1, kyzp
      k1 = k + koff
      if (((k1.gt.0).and.(k1.lt.kmax)).or.((k1.gt.kmin).and.(k1.lt.ny)))
     1 then
         if (k1.gt.nyh) k1 = k1 - ny
         do 60 j = 1, jmax
         if ((j+joff).gt.1) then
            do 20 l = 2, lmax
            l1 = nz2 - l
            do 10 i = 1, ndim
            vpot(i,l,j,k,m) = vpott(it,i,2*l-2,j,k,m)
            vpot(i,l1,j,k,m) = vpott(it,i,2*l-1,j,k,m)
   10       continue
   20       continue
            do 40 l = lmax+1, nzh
            l1 = nz2 - l
            do 30 i = 1, ndim
            vpot(i,l,j,k,m) = zero
            vpot(i,l1,j,k,m) = zero
   30       continue
   40       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1 
            do 50 i = 1, ndim
            vpot(i,1,j,k,m) = vpott(it,i,1,j,k,m)
            vpot(i,l1,j,k,m) = zero
            if (modesz.gt.nzh) then
               vpot(i,l1,j,k,m) = vpott(it,i,nz,j,k,m)
            endif
   50       continue
         endif
   60    continue
         do 90 j = jmax+1, kxyp
         if ((j+joff).gt.1) then
            do 80 l = 1, nz
            do 70 i = 1, ndim
            vpot(i,l,j,k,m) = zero
   70       continue
   80       continue
         endif
   90    continue
c mode numbers kx = 0, nx/2
         n1 = joff
         if (n1.eq.0) then
c kx = 0
            if (k1.gt.0) then
               do 110 l = 2, lmax
               l1 = nz2 - l
               do 100 i = 1, ndim
               vpot(i,l,1,k,m) = vpott(it,i,2*l-2,1,k,m)
               vpot(i,l1,1,k,m) = vpott(it,i,2*l-1,1,k,m)
  100          continue
  110          continue
               do 130 l = lmax+1, nzh
               l1 = nz2 - l
               do 120 i = 1, ndim
               vpot(i,l,1,k,m) = zero
               vpot(i,l1,1,k,m) = zero
  120          continue
  130          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               do 140 i = 1, ndim
               vpot(i,1,1,k,m) = vpott(it,i,1,1,k,m)
               vpot(i,l1,1,k,m) = zero
               if (modesz.gt.nzh) then
                  vpot(i,l1,1,k,m) = vpott(it,i,nz,1,k,m)
               endif
  140          continue
c kx = nx/2
            else
               do 160 l = 1, nz
               do 150 i = 1, ndim
               vpot(i,l,1,k,m) = zero
  150          continue
  160          continue
               if (modesx.gt.nxh) then
                  do 180 l = 2, lmax
                  l1 = nz2 - l
                  do 170 i = 1, ndim
                  vpot(i,l,1,k,m) = vpott(it,i,2*l-2,1,k,m)
                  vpot(i,l1,1,k,m) = vpott(it,i,2*l-1,1,k,m)
  170             continue
  180             continue
c mode numbers kz = 0, nz/2
                  do 190 i = 1, ndim
                  vpot(i,1,1,k,m) = vpott(it,i,1,1,k,m)
                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     vpot(i,l1,1,k,m) = vpott(it,i,nz,1,k,m)
                  endif
  190             continue
               endif
            endif
         endif
      endif
  200 continue
      do 260 k = 1, kyzp
      k1 = k + koff
      if ((k1.ge.kmax).and.(k1.le.kmin).and.(k1.ne.nyh)) then
         do 230 j = 1, kxyp
         if ((j+joff).gt.1) then
            do 220 l = 1, nz
            do 210 i = 1, ndim
            vpot(i,l,j,k,m) = zero
  210       continue
  220       continue
         endif
  230    continue
c mode numbers kx = 0, nx/2
         n1 = joff
         if (n1.eq.0) then
            do 250 l = 1, nz
            do 240 i = 1, ndim
            vpot(i,l,1,k,m) = zero
  240       continue
  250       continue
         endif
      endif
  260 continue
c mode numbers ky = 0, ny/2
      n2 = koff + 1
c ky = 0
      if (n2.eq.0) then
         do 320 j = 1, jmax
         if ((j+joff).gt.1) then
            do 280 l = 2, lmax
            l1 = nz2 - l
            do 270 i = 1, ndim
            vpot(i,l,j,1,m) = vpott(it,i,2*l-2,j,1,m)
            vpot(i,l1,j,1,m) = vpott(it,i,2*l-1,j,1,m)
  270       continue
  280       continue
            do 300 l = lmax+1, nzh
            l1 = nz2 - l
            do 290 i = 1, ndim
            vpot(i,l,j,1,m) = zero
            vpot(i,l1,j,1,m) = zero
  290       continue
  300       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            do 310 i = 1, ndim
            vpot(i,1,j,1,m) = vpott(it,i,1,j,1,m)
            vpot(i,l1,j,1,m) = zero
            if (modesz.gt.nzh) then
               vpot(i,l1,j,1,m) = vpott(it,i,nz,j,1,m)
            endif
  310       continue
         endif
  320    continue
         do 350 j = jmax+1, kxyp
         if ((j+joff).gt.1) then
            do 340 l = 1, nz
            do 330 i = 1, ndim
            vpot(i,l,j,1,m) = zero
  330       continue
  340       continue
         endif
  350    continue
c mode numbers kx = 0, nx/2
         n1 = joff
         if (n1.eq.0) then
c kx = 0
            do 370 l = 2, lmax
            l1 = nz2 - l
            do 360 i = 1, ndim
            vpot(i,l,1,1,m) = vpott(it,i,2*l-2,1,1,m)
            vpot(i,l1,1,1,m) = zero
  360       continue
  370       continue
            do 390 l = lmax+1, nzh
            l1 = nz2 - l
            do 380 i = 1, ndim
            vpot(i,l,1,1,m) = zero
            vpot(i,l1,1,1,m) = zero
  380       continue
  390       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            do 400 i = 1, ndim
            vpot(i,1,1,1,m) = cmplx(real(vpott(it,i,1,1,1,m)),0.)
            vpot(i,l1,1,1,m) = zero
            if (modesz.gt.nzh) then
               vpot(i,l1,1,1,m) = cmplx(real(vpott(it,i,nz,1,1,m)),0.)
            endif
  400       continue
c kx = nx/2
            if (modesx.gt.nxh) then
               do 420 l = 2, lmax
               l1 = nz2 - l
               do 410 i = 1, ndim
               vpot(i,l1,1,1,m) = conjg(vpott(it,i,2*l-2,j1,1,m))
  410          continue
  420          continue
c mode numbers kz = 0, nz/2
               do 430 i = 1, ndim
               vpot(i,1,1,1,m) = cmplx(real(vpot(i,1,1,1,m)),real(vpott(
     1it,i,1,j1,1,m)))
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpot(i,l1,1,1,m) = cmplx(real(vpot(i,l1,1,1,m)),real(v
     1pott(it,i,nz,j1,1,m)))
               endif
  430          continue
            endif
         endif
      endif
c ky = ny/2
      k1 = (nyh/kyzp)*kyzp
      if (n2.eq.k1) then
         k1 = nyh - k1 + 1
         do 460 j = 1, jmax
         if ((j+joff).gt.1) then
            do 450 l = 1, nz
            do 440 i = 1, ndim
            vpot(i,l,j,k1,m) = zero
  440       continue
  450       continue
         endif
  460    continue
         do 490 j = jmax+1, kxyp
         if ((j+joff).gt.1) then
            do 480 l = 1, nz
            do 470 i = 1, ndim
            vpot(i,l,j,k1,m) = zero
  470       continue
  480       continue
         endif
  490    continue
c mode numbers kx = 0, nx/2
         n1 = joff
         if (n1.eq.0) then
            do 510 l = 1, nz
            do 500 i = 1, ndim
            vpot(i,l,1,k1,m) = zero
  500       continue
  510       continue
         endif
         if (modesy.gt.nyh) then
            do 550 j = 1, jmax
            if ((j+joff).gt.1) then
               do 530 l = 2, lmax
               l1 = nz2 - l
               do 520 i = 1, ndim
               vpot(i,l,j,k1,m) = vpott(it,i,2*l-2,j,k1,m)
               vpot(i,l1,j,k1,m) = vpott(it,i,2*l-1,j,k1,m)
  520          continue
  530          continue
c mode numbers kz = 0, nz/2
               do 540 i = 1, ndim
               vpot(i,1,j,k1,m) = vpott(it,i,1,j,k1,m)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpot(i,l1,j,k1,m) = vpott(it,i,nz,j,k1,m)
               endif
  540          continue
            endif
  550       continue
c mode numbers kx = 0, nx/2
            n1 = joff
            if (n1.eq.0) then
c kx = 0
               do 570 l = 2, lmax
               do 560 i = 1, ndim
               vpot(i,l,1,k1,m) = vpott(it,i,2*l-2,1,k1,m)
  560          continue
  570          continue
c mode numbers kz = 0, nz/2
               do 580 i = 1, ndim
               vpot(i,1,1,k1,m) = cmplx(real(vpott(it,i,1,1,k1,m)),0.)
               if (modesz.gt.nzh) then
                  l1 = nzh + 1
                  vpot(i,l1,1,k1,m) = cmplx(real(vpott(it,i,nz,1,k1,m)),
     10.)
               endif
  580          continue
c kx  = nx/2
               if (modesx.gt.nxh) then
                  do 600 l = 2, lmax
                  l1 = nz2 - l
                  do 590 i = 1, ndim
                  vpot(i,l1,1,k1,m) = conjg(vpott(it,i,2*l-2,j1,k1,m))
  590             continue
  600             continue
c mode numbers kz = 0, nz/2
                  do 610 i = 1, ndim
                  vpot(i,1,1,k1,m) = cmplx(real(vpot(i,1,1,k1,m)),real(v
     1pott(it,i,1,j1,k1,m)))

                  if (modesz.gt.nzh) then
                     l1 = nzh + 1
                     vpot(i,l1,1,k1,m) = cmplx(real(vpot(i,l1,1,k1,m)),r
     1eal(vpott(it,i,nz,j1,k1,m)))
                  endif
  610             continue
               endif
            endif
         endif
      endif
  620 continue
  630 continue
  640 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PADDQEI32(qe,qi,nyzp,nx,nxe,nypmx,nzpmx,idds,mnblok)
c adds electron and ion densities
c assumes guard cells have already been added
c qe/qi = charge density for electrons/ions
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c nx = system length in x direction
c nxe = first dimension of field array, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
c mnblok = number of particle partitions.
      implicit none
      integer nx, nxe, nypmx, nzpmx, idds, mnblok
      real qe, qi
      integer nyzp
      dimension qe(nxe,nypmx,nzpmx,mnblok), qi(nxe,nypmx,nzpmx,mnblok)
      dimension nyzp(idds,mnblok)
      integer j, k, l, m
      do 40 m = 1, mnblok
      do 30 l = 1, nyzp(2,m)
      do 20 k = 1, nyzp(1,m)
      do 10 j = 1, nx
      qe(j,k,l,m) = qe(j,k,l,m) + qi(j,k,l,m)
   10 continue
   20 continue
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PADDQEI32X(qe,qi,nyzp,qbme,qbmi,wpmax,wpmin,nx,nxe,nypm
     1x,nzpmx,idds,mnblok)
c adds electron and ion densities, and calculates maximum and minimum
c plasma frequency.  assumes guard cells have already been added
c qe/qi = charge density for electrons/ions
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c qbme/qbmi = charge/mass ratio for electrons/ions
c wpmax/wpmin = maximum/minimum plasma frequency
c nx = system length in x direction
c nxe = first dimension of field array, must be >= nx
c nypmx = maximum size of particle partition in y, including guard cells
c nzpmx = maximum size of particle partition in z, including guard cells
c idds = dimensionality of domain decomposition
c mnblok = number of particle partitions.
      implicit none
      integer nx, nxe, nypmx, nzpmx, idds, mnblok
      real qbme, qbmi, wpmax, wpmin
      real qe, qi
      integer nyzp
      dimension qe(nxe,nypmx,nzpmx,mnblok), qi(nxe,nypmx,nzpmx,mnblok)
      dimension nyzp(idds,mnblok)
      integer j, k, l, m
      real at1
      double precision sum1
      sum1 = 0.0d0
      wpmax = qbme*qe(1,1,1,1) + qbmi*qi(1,1,1,1)
      wpmin = wpmax
      do 40 m = 1, mnblok
      do 30 l = 1, nyzp(2,m)
      do 20 k = 1, nyzp(1,m)
      do 10 j = 1, nx
      at1 = qbme*qe(j,k,l,m) + qbmi*qi(j,k,l,m)
      qe(j,k,l,m) = qe(j,k,l,m) + qi(j,k,l,m)
      wpmax = max(wpmax,at1)
      wpmin = min(wpmin,at1)
   10 continue
   20 continue
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PIMOMENT32(qi,fxyz,nyzp,pxi,pyi,pzi,dt,nx,nxe,nypmx,nzp
     1mx,idds,mnblok)
c calculate ion momentum for 3d code from integral of qi*fxyz
c assumes guard cells have already been added
      implicit none
      integer nx, nxe, nypmx, nzpmx, idds, mnblok
      real pxi, pyi, pzi, dt
      real qi, fxyz
      integer nyzp
      dimension qi(nxe,nypmx,nzpmx,mnblok)
      dimension fxyz(3,nxe,nypmx,nzpmx,mnblok)
      dimension nyzp(idds,mnblok)
      integer j, k, l, m
      double precision dt1, sum1, sum2, sum3
      sum1 = 0.0d0
      sum2 = 0.0d0
      sum3 = 0.0d0
      do 40 m = 1, mnblok
      do 30 l = 1, nyzp(2,m)
      do 20 k = 1, nyzp(1,m)
      do 10 j = 1, nx
      dt1 = dble(qi(j,k,l,m))
      sum1 = sum1 + dt1*fxyz(1,j,k,l,m)
      sum2 = sum2 + dt1*fxyz(2,j,k,l,m)
      sum3 = sum3 + dt1*fxyz(3,j,k,l,m)
   10 continue
   20 continue
   30 continue
   40 continue
      pxi = pxi + sum1*dt
      pyi = pyi + sum2*dt
      pzi = pzi + sum3*dt
      return
      end
