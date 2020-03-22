c 3d parallel PIC library for solving field equations with open (vacuum)
c boundary conditions and 2D domain decomposition
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: february 6, 2008
c-----------------------------------------------------------------------
      subroutine PFORMC32(ffg,f,fs,ft,bs,br,g,h,fpotc,mixup3,sct3,affp,a
     1r,indx1,indy1,indz1,kstrt,nxv,ny2d,nz2d,kxyp2,kyp2,kyzp2,kzp2,kyp2
     2d,kzyp2,j2blok,k2blok,jk2blok,l2blok,m2blok,ml2blok,kxyp2d,kyzp2d,
     3nz1d,nxhyz2,nxyzh2)
c this subroutine calculates the form factor array ffg needed by field
c solvers with open (vacuum) boundary conditions using hockney's method.
c the five green's functions calculated are:
c g(kx,ky,kz) = affp*inverse FFT of potr
c s(kx,ky,kz) = inverse FFT of the density of a finite-sized particle
c gx(kx,ky,kz) = affp*s(kx,ky,kz)*inverse FFT of (x/r)*Er
c gy(kx,ky,kz) = affp*s(kx,ky,kz)*inverse FFT of (y/r)*Er
c gz(kx,ky,kz) = affp*s(kx,ky,kz)*inverse FFT of (z/r)*Er
c where the fields due to the finite-sized particles are given by fpotc
c input: fpotc,mixup3,sct3,affp,ar,indx1,indy1,indz1,kstrt,nxv,ny2d,nz2d
c        kxyp2,kyp2,kyzp2,kzp2,kyp2d,j2blok,k2blok,l2blok,m2blok,kxyp2d
c        nz1d,nxhyz2,nxyzh2
c ffg(1,l,j,k) = potential green's function g
c ffg(2,l,j,k) = finite-size particle shape factor s
c ffg(3,l,j,k) = x component of electric field green's function gx
c ffg(4,l,j,k) = y component of electric field green's function gy
c ffg(5,l,j,k) = z component of electric field green's function gz
c on processor 0, ffg(i,l,k,kxp2+1,m) = ffg(i,l,k,NX+1,m)
c on other processors, ffg(i,l,kxp2+1,k,m) = ffg(i,l,1,k,m) on proc 0
c f, fs, ft = scratch arrays used by FFT
c g, h = scratch arrays used for message-passing
c fpotc = a function which calculates green's function
c mixup3/sct3 = bit-reverse and sine-cosine table used by FFT
c affp = normalization constant = nx*ny*nz/np, where np=number of
c particles
c ar = half-width of particle in r direction
c indx1/indy1/indz1 = exponent which determines FFT length in x/y/z
c direction, where 2*nx=2**indx1, 2*ny=2**indy1, 2*nz=2**indz1
c kstrt = starting data block number
c nxv = half of first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
c nz2d = third dimension of field arrays, must be >= 2*nz
c kyp2d = second dimension of f
c kxyp2/kyp2,kyzp2/kzp2 = number of data values per block in x/y/z
c kzyp2 = maximum(kyzp2,kyp2)
c j2blok/k2blok,m2blok/l2blok = number of data blocks in x/y/z
c jk2blok = maximum(j2blok,k2blok)
c ml2blok = maximum(m2blok,l2blok)
c kxyp2d = third dimension of ffg arrays, must be >= kxyp2+1
c kyzp2d = third dimension of ffg arrays, must be >= kyzp2+1
c nz1d = second dimension of ffg arrays, must be >= nz+1
c nxhyz2 = maximum of (nx,2*ny,2*nz)
c nxyzh2 = maximum of (nx,ny,nz)
      implicit none
      real ffg, f, g, h
      complex fs, ft, bs, br
      integer mixup3
      complex sct3
      real affp, ar
      integer indx1, indy1, indz1, kstrt, nxv, ny2d, nz2d, kxyp2, kyp2
      integer kyzp2, kzp2, kyp2d, kzyp2, j2blok, k2blok, jk2blok
      integer l2blok, m2blok, ml2blok, kxyp2d, kyzp2d, nz1d
      integer nxhyz2, nxyzh2
      dimension ffg(5,nz1d,kxyp2d,kyzp2d,j2blok*m2blok)
      dimension f(2*nxv,kyp2d,kzp2,k2blok*l2blok)
      dimension ft(nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension fs(ny2d,kxyp2,kzp2,j2blok*l2blok)
      dimension bs(kxyp2,kzyp2,kzp2,jk2blok*l2blok)
      dimension br(kxyp2,kzyp2,kzp2,j2blok*ml2blok)
      dimension g(5,nz1d,kyzp2,j2blok*m2blok)
      dimension h(5,nz1d,kyzp2,j2blok*m2blok)
      dimension mixup3(nxhyz2), sct3(nxyzh2)
      real fpotc
      external fpotc
c local data
      integer ntpose, nx, ny, nz, nz1, nx2, ny2, nz2, kxb2, kyb2, kzb2
      integer kyzb2, isign, i, j, k, l, mx, my, mz, m, j1, k1, l1, n2
      integer js, ks, joff, koff, loff, moff, ifun
      real an, ari, at1, at2, x, y, z, r, ttp
      real POTC3
      external POTC3
      data ntpose /1/
      nx2 = 2**indx1
      ny2 = 2**indy1
      nz2 = 2**indz1
      nx = nx2/2
      ny = ny2/2
      nz = nz2/2
      nz1 = nz + 1
      kyb2 = ny2/kyp2
      kzb2 = nz2/kzp2
      kxb2 = nx/kxyp2
      kyzb2 = ny2/kyzp2
      ari = 0.0
      if (ar.gt.0.) ari = 1.0/ar
      an = float(nx2*ny2*nz2)
c calculate potential green's function
      ifun = 1
      if (kstrt.gt.(kyb2*kzb2)) go to 60
      ks = (kstrt - 1)/kyb2
      js = kstrt - kyb2*ks - 2
      ks = ks - 1
      do 50 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      loff = kzp2*(mz + ks) - 1
      do 40 my = 1, k2blok
      koff = kyp2*(my + js) - 1
      m = my + moff
      do 30 l = 1, kzp2
      l1 = l + loff
      if (l1.gt.nz) l1 = l1 - nz2
      at1 = float(l1)**2
      do 20 k = 1, kyp2
      k1 = k + koff
      if (k1.gt.ny) k1 = k1 - ny2
      at2 = float(k1)**2 + at1
      do 10 j = 1, nx2
      j1 = j - 1
      if (j1.gt.nx) j1 = j1 - nx2
      r = sqrt(at2 + float(j1)**2)
      f(j,k,l,m) = fpotc(r,affp,ari,1)
   10 continue
   20 continue
   30 continue
   40 continue
   50 continue
      isign = -1
      call WPFFT32R(f,fs,ft,bs,br,isign,ntpose,mixup3,sct3,ttp,indx1,ind
     1y1,indz1,kstrt,nxv,ny2d,nz2d,kxyp2,kyp2,kyzp2,kzp2,kxyp2,kyp2d,kyz
     2p2,kzp2,kzyp2,j2blok,k2blok,jk2blok,l2blok,m2blok,ml2blok,nxhyz2,n
     3xyzh2)
   60 if (kstrt.gt.(kxb2*kyzb2)) go to 150
      ks = (kstrt - 1)/kxb2
      js = kstrt - kxb2*ks - 2
      ks = ks - 1
      do 140 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 130 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js)
      m = mx + moff
      do 100 k = 1, kyzp2
      k1 = k + koff
      do 80 j = 1, kxyp2
      do 70 l = 1, nz1
      ffg(ifun,l,j,k,m) = an*real(ft(l,j,k,m))
   70 continue
   80 continue
      do 90 l = 1, nz1
      ffg(ifun,l,kxyp2+1,k,m) = 0.0
   90 continue
  100 continue
c mode numbers kx = 0, nx
      if (joff.eq.0) then
c mode number ky = 0
         n2 = koff + 1
         if (n2.eq.0) then
            do 110 l = 2, nz
            l1 = nz2 + 2 - l
            ffg(ifun,l,kxyp2+1,1,m) = an*real(ft(l1,1,1,m))
  110       continue
            ffg(ifun,1,kxyp2+1,1,m) = an*aimag(ft(1,1,1,m))
            ffg(ifun,nz1,kxyp2+1,1,m) = an*aimag(ft(nz1,1,1,m))
         endif
c mode number ky = ny
         k1 = (ny/kyzp2)*kyzp2
         if (n2.eq.k1) then
            k1 = ny - k1 + 1
            do 120 l = 2, nz
            l1 = nz2 + 2 - l
            ffg(ifun,l,kxyp2+1,k1,m) = an*real(ft(l1,1,k1,m))
  120       continue
            ffg(ifun,1,kxyp2+1,k1,m) = an*aimag(ft(1,1,k1,m))
            ffg(ifun,nz1,kxyp2+1,k1,m) = an*aimag(ft(nz1,1,k1,m))
         endif
      endif
  130 continue
  140 continue
c calculate particle smoothing function
  150 ifun = ifun + 1
      if (kstrt.gt.(kyb2*kzb2)) go to 210
      ks = (kstrt - 1)/kyb2
      js = kstrt - kyb2*ks - 2
      ks = ks - 1
      do 200 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      loff = kzp2*(mz + ks) - 1
      do 190 my = 1, k2blok
      koff = kyp2*(my + js) - 1
      m = my + moff
      do 180 l = 1, kzp2
      l1 = l + loff
      if (l1.gt.nz) l1 = l1 - nz2
      at1 = float(l1)**2
      do 170 k = 1, kyp2
      k1 = k + koff
      if (k1.gt.ny) k1 = k1 - ny2
      at2 = float(k1)**2 + at1
      do 160 j = 1, nx2
      j1 = j - 1
      if (j1.gt.nx) j1 = j1 - nx2
      r = sqrt(at2 + float(j1)**2)
      f(j,k,l,m) = fpotc(r,affp,ari,2)
  160 continue
  170 continue
  180 continue
  190 continue
  200 continue
      isign = -1
      call WPFFT32R(f,fs,ft,bs,br,isign,ntpose,mixup3,sct3,ttp,indx1,ind
     1y1,indz1,kstrt,nxv,ny2d,nz2d,kxyp2,kyp2,kyzp2,kzp2,kxyp2,kyp2d,kyz
     2p2,kzp2,kzyp2,j2blok,k2blok,jk2blok,l2blok,m2blok,ml2blok,nxhyz2,n
     3xyzh2)
  210 if (kstrt.gt.(kxb2*kyzb2)) go to 300
      ks = (kstrt - 1)/kxb2
      js = kstrt - kxb2*ks - 2
      ks = ks - 1
      do 290 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 280 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js)
      m = mx + moff
      do 250 k = 1, kyzp2
      k1 = k + koff
      do 230 j = 1, kxyp2
      do 220 l = 1, nz1
      ffg(ifun,l,j,k,m) = an*real(ft(l,j,k,m))
  220 continue
  230 continue
      do 240 l = 1, nz1
      ffg(ifun,l,kxyp2+1,k,m) = 0.0
  240 continue
  250 continue
c mode numbers kx = 0, nx
      if (joff.eq.0) then
c mode number ky = 0
         n2 = koff + 1
         if (n2.eq.0) then
            do 260 l = 2, nz
            l1 = nz2 + 2 - l
            ffg(ifun,l,kxyp2+1,1,m) = an*real(ft(l1,1,1,m))
  260       continue
            ffg(ifun,1,kxyp2+1,1,m) = an*aimag(ft(1,1,1,m))
            ffg(ifun,nz1,kxyp2+1,1,m) = an*aimag(ft(nz1,1,1,m))
         endif
c mode number ky = ny
         k1 = (ny/kyzp2)*kyzp2
         if (n2.eq.k1) then
            k1 = ny - k1 + 1
            do 270 l = 2, nz
            l1 = nz2 + 2 - l
            ffg(ifun,l,kxyp2+1,k1,m) = an*real(ft(l1,1,k1,m))
  270       continue
            ffg(ifun,1,kxyp2+1,k1,m) = an*aimag(ft(1,1,k1,m))
            ffg(ifun,nz1,kxyp2+1,k1,m) = an*aimag(ft(nz1,1,k1,m))
         endif
      endif
  280 continue
  290 continue
c calculate green's function for x component of electric field
  300 ifun = ifun + 1
      if (kstrt.gt.(kyb2*kzb2)) go to 360
      ks = (kstrt - 1)/kyb2
      js = kstrt - kyb2*ks - 2
      ks = ks - 1
      do 350 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      loff = kzp2*(mz + ks) - 1
      do 340 my = 1, k2blok
      koff = kyp2*(my + js) - 1
      m = my + moff
      do 330 l = 1, kzp2
      l1 = l + loff
      if (l1.gt.nz) l1 = l1 - nz2
      at1 = float(l1)**2
      do 320 k = 1, kyp2
      k1 = k + koff
      if (k1.gt.ny) k1 = k1 - ny2
      at2 = float(k1)**2 + at1
      do 310 j = 1, nx2
      j1 = j - 1
      if (j1.gt.nx) j1 = j1 - nx2
      x = float(j1)
      r = sqrt(at2 + x*x)
      f(j,k,l,m) = fpotc(r,affp,ari,3)
      if (r.gt.0.) f(j,k,l,m) = f(j,k,l,m)*(x/r)
  310 continue
  320 continue
  330 continue
  340 continue
  350 continue
      isign = -1
      call WPFFT32R(f,fs,ft,bs,br,isign,ntpose,mixup3,sct3,ttp,indx1,ind
     1y1,indz1,kstrt,nxv,ny2d,nz2d,kxyp2,kyp2,kyzp2,kzp2,kxyp2,kyp2d,kyz
     2p2,kzp2,kzyp2,j2blok,k2blok,jk2blok,l2blok,m2blok,ml2blok,nxhyz2,n
     3xyzh2)
  360 if (kstrt.gt.(kxb2*kyzb2)) go to 460
      ks = (kstrt - 1)/kxb2
      js = kstrt - kxb2*ks - 2
      ks = ks - 1
      do 450 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 440 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js)
      m = mx + moff
      do 410 k = 1, kyzp2
      k1 = k + koff
      do 380 j = 1, kxyp2
      if ((j+joff).gt.0) then
         do 370 l = 1, nz1
         ffg(ifun,l,j,k,m) = an*aimag(ft(l,j,k,m))
  370    continue
      endif
  380 continue
c mode numbers kx = 0, nx
      if (joff.eq.0) then
         do 390 l = 1, nz1
         ffg(ifun,l,1,k,m) = an*real(ft(l,1,k,m))
  390    continue
      endif
      do 400 l = 1, nz1
      ffg(ifun,l,kxyp2+1,k,m) = 0.0
  400 continue
  410 continue
c mode numbers kx = 0, nx
      if (joff.eq.0) then
c mode number ky = 0
         n2 = koff + 1
         if (n2.eq.0) then
            do 420 l = 2, nz
            l1 = nz2 + 2 - l
            ffg(ifun,l,kxyp2+1,1,m) = an*real(ft(l1,1,1,m))
  420       continue
            ffg(ifun,1,kxyp2+1,1,m) = an*aimag(ft(1,1,1,m))
            ffg(ifun,nz1,kxyp2+1,1,m) = an*aimag(ft(nz1,1,1,m))
         endif
c mode number ky = ny
         k1 = (ny/kyzp2)*kyzp2
         if (n2.eq.k1) then
            k1 = ny - k1 + 1
            do 430 l = 2, nz
            l1 = nz2 + 2 - l
            ffg(ifun,l,kxyp2+1,k1,m) = an*real(ft(l1,1,k1,m))
  430       continue
            ffg(ifun,1,kxyp2+1,k1,m) = an*aimag(ft(1,1,k1,m))
            ffg(ifun,nz1,kxyp2+1,k1,m) = an*aimag(ft(nz1,1,k1,m))
         endif
      endif
  440 continue
  450 continue
c calculate green's function for y component of electric field
  460 ifun = ifun + 1
      if (kstrt.gt.(kyb2*kzb2)) go to 520
      ks = (kstrt - 1)/kyb2
      js = kstrt - kyb2*ks - 2
      ks = ks - 1
      do 510 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      loff = kzp2*(mz + ks) - 1
      do 500 my = 1, k2blok
      koff = kyp2*(my + js) - 1
      m = my + moff
      do 490 l = 1, kzp2
      l1 = l + loff
      if (l1.gt.nz) l1 = l1 - nz2
      at1 = float(l1)**2
      do 480 k = 1, kyp2
      k1 = k + koff
      if (k1.gt.ny) k1 = k1 - ny2
      y = float(k1)
      at2 = y*y + at1
      do 470 j = 1, nx2
      j1 = j - 1
      if (j1.gt.nx) j1 = j1 - nx2
      r = sqrt(at2 + float(j1)**2)
      f(j,k,l,m) = fpotc(r,affp,ari,3)
      if (r.gt.0.) f(j,k,l,m) = f(j,k,l,m)*(y/r)
  470 continue
  480 continue
  490 continue
  500 continue
  510 continue
      isign = -1
      call WPFFT32R(f,fs,ft,bs,br,isign,ntpose,mixup3,sct3,ttp,indx1,ind
     1y1,indz1,kstrt,nxv,ny2d,nz2d,kxyp2,kyp2,kyzp2,kzp2,kxyp2,kyp2d,kyz
     2p2,kzp2,kzyp2,j2blok,k2blok,jk2blok,l2blok,m2blok,ml2blok,nxhyz2,n
     3xyzh2)
  520 if (kstrt.gt.(kxb2*kyzb2)) go to 650
      ks = (kstrt - 1)/kxb2
      js = kstrt - kxb2*ks - 2
      ks = ks - 1
      do 640 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 630 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js)
      m = mx + moff
      do 560 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         do 540 j = 1, kxyp2
         do 530 l = 1, nz1
         ffg(ifun,l,j,k,m) = an*aimag(ft(l,j,k,m))
  530    continue
  540    continue
      endif
      do 550 l = 1, nz1
      ffg(ifun,l,kxyp2+1,k,m) = 0.0
  550 continue
  560 continue
c mode number ky = 0
      n2 = koff + 1
      if (n2.eq.0) then
         do 580 j = 1, kxyp2
         do 570 l = 1, nz1
         ffg(ifun,l,j,1,m) = an*real(ft(l,j,1,m))
  570    continue
  580    continue
      endif
c mode numbers ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 600 j = 1, kxyp2
         do 590 l = 1, nz1
         ffg(ifun,l,j,k1,m) = an*real(ft(l,j,k1,m))
  590    continue
  600    continue
      endif
c mode numbers kx = 0, nx
      if (joff.eq.0) then
c mode number ky = 0
         n2 = koff + 1
         if (n2.eq.0) then
            do 610 l = 2, nz
            l1 = nz2 + 2 - l
            ffg(ifun,l,kxyp2+1,1,m) = an*real(ft(l1,1,1,m))
  610       continue
            ffg(ifun,1,kxyp2+1,1,m) = an*aimag(ft(1,1,1,m))
            ffg(ifun,nz1,kxyp2+1,1,m) = an*aimag(ft(nz1,1,1,m))
         endif
c mode number ky = ny
         k1 = (ny/kyzp2)*kyzp2
         if (n2.eq.k1) then
            k1 = ny - k1 + 1
            do 620 l = 2, nz
            l1 = nz2 + 2 - l
            ffg(ifun,l,kxyp2+1,k1,m) = an*real(ft(l1,1,k1,m))
  620       continue
            ffg(ifun,1,kxyp2+1,k1,m) = an*aimag(ft(1,1,k1,m))
            ffg(ifun,nz1,kxyp2+1,k1,m) = an*aimag(ft(nz1,1,k1,m))
         endif
      endif
  630 continue
  640 continue
c calculate green's function for z component of electric field
  650 ifun = ifun + 1
      if (kstrt.gt.(kyb2*kzb2)) go to 710
      ks = (kstrt - 1)/kyb2
      js = kstrt - kyb2*ks - 2
      ks = ks - 1
      do 700 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      loff = kzp2*(mz + ks) - 1
      do 690 my = 1, k2blok
      koff = kyp2*(my + js) - 1
      m = my + moff
      do 680 l = 1, kzp2
      l1 = l + loff
      if (l1.gt.nz) l1 = l1 - nz2
      z = float(l1)
      at1 = z*z
      do 670 k = 1, kyp2
      k1 = k + koff
      if (k1.gt.ny) k1 = k1 - ny2
      at2 = float(k1)**2 + at1
      do 660 j = 1, nx2
      j1 = j - 1
      if (j1.gt.nx) j1 = j1 - nx2
      r = sqrt(at2 + float(j1)**2)
      f(j,k,l,m) = fpotc(r,affp,ari,3)
      if (r.gt.0.) f(j,k,l,m) = f(j,k,l,m)*(z/r)
  660 continue
  670 continue
  680 continue
  690 continue
  700 continue
      isign = -1
      call WPFFT32R(f,fs,ft,bs,br,isign,ntpose,mixup3,sct3,ttp,indx1,ind
     1y1,indz1,kstrt,nxv,ny2d,nz2d,kxyp2,kyp2,kyzp2,kzp2,kxyp2,kyp2d,kyz
     2p2,kzp2,kzyp2,j2blok,k2blok,jk2blok,l2blok,m2blok,ml2blok,nxhyz2,n
     3xyzh2)
  710 if (kstrt.gt.(kxb2*kyzb2)) go to 840
      ks = (kstrt - 1)/kxb2
      js = kstrt - kxb2*ks - 2
      ks = ks - 1
      do 830 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 820 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js)
      m = mx + moff
      do 750 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         do 730 j = 1, kxyp2
         do 720 l = 2, nz
         l1 = nz2 + 2 - l
         ffg(ifun,l,j,k,m) = an*aimag(ft(l,j,k,m))
  720    continue
         ffg(ifun,1,j,k,m) = an*real(ft(1,j,k,m))
         ffg(ifun,nz1,j,k,m) = an*real(ft(nz1,j,k,m))
  730    continue
      endif
      do 740 l = 1, nz1
      ffg(ifun,l,kxyp2+1,k,m) = 0.0
  740 continue
  750 continue
c mode number ky = 0
      n2 = koff + 1
      if (n2.eq.0) then
         do 770 j = 1, kxyp2
         do 760 l = 2, nz
         l1 = nz2 + 2 - l
         ffg(ifun,l,j,1,m) = an*aimag(ft(l,j,1,m))
  760    continue
         ffg(ifun,1,j,1,m) = an*real(ft(1,j,1,m))
         ffg(ifun,nz1,j,1,m) = an*real(ft(nz1,j,1,m))
  770    continue
      endif
c mode numbers ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 790 j = 1, kxyp2
         do 780 l = 2, nz
         l1 = nz2 + 2 - l
         ffg(ifun,l,j,k1,m) = an*aimag(ft(l,j,k1,m))
  780    continue
         ffg(ifun,1,j,k1,m) = an*real(ft(1,j,k1,m))
         ffg(ifun,nz1,j,k1,m) = an*real(ft(nz1,j,k1,m))
  790    continue
      endif
c mode numbers kx = 0, nx
      if (joff.eq.0) then
c mode number ky = 0
         n2 = koff + 1
         if (n2.eq.0) then
            do 800 l = 2, nz
            l1 = nz2 + 2 - l
            ffg(ifun,l,kxyp2+1,1,m) = an*aimag(ft(l1,1,1,m))
  800       continue
            ffg(ifun,1,kxyp2+1,1,m) = an*aimag(ft(1,1,1,m))
            ffg(ifun,nz1,kxyp2+1,1,m) = an*aimag(ft(nz1,1,1,m))
         endif
c mode number ky = ny
         k1 = (ny/kyzp2)*kyzp2
         if (n2.eq.k1) then
            k1 = ny - k1 + 1
            do 810 l = 2, nz
            l1 = nz2 + 2 - l
            ffg(ifun,l,kxyp2+1,k1,m) = an*aimag(ft(l1,1,k1,m))
  810       continue
            ffg(ifun,1,kxyp2+1,k1,m) = an*aimag(ft(1,1,k1,m))
            ffg(ifun,nz1,kxyp2+1,k1,m) = an*aimag(ft(nz1,1,k1,m))
         endif
      endif
  820 continue
  830 continue
c copy ffg(i,l,1,k,m) on node 0 to ffg(i,l,kxyp2+1,k,m) on other nodes
  840 do 910 m = 1, j2blok*m2blok
      do 870 k = 1, kyzp2
      do 860 l = 1, nz1d
      do 850 i = 1, 5
      g(i,l,k,m) = ffg(i,l,1,k,m)
      h(i,l,k,m) = ffg(i,l,kxyp2+1,k,m)
  850 continue
  860 continue
  870 continue
      call P0COPY2(g,h,5*kyzp2*nz1d,kstrt,kxb2,kyzb2,1)
      do 900 k = 1, kyzp2
      do 890 l = 1, nz1d
      do 880 i = 1, 5
      ffg(i,l,kxyp2+1,k,m) = h(i,l,k,m)
  880 continue
  890 continue
  900 continue
  910 continue
c copy ffg(i,l,k,1,m) on node 0 to ffg(i,l,k,kxp2+1,m) on other nodes
      do 920 m = 1, j2blok*m2blok
      call P0COPY2(ffg(1,1,1,1,m),ffg(1,1,1,kyzp2+1,m),5*kxyp2d*nz1d,kst
     1rt,kxb2,kyzb2,2)
  920 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFORMC32X(ffg,f,fs,ft,g,h,fpotc,mixup3,sct3,affp,ar,ind
     1x1,indy1,indz1,kstrt,nxv,ny2d,nz2d,kxyp2,kyp2,kyzp2,kzp2,kyp2d,j2b
     2lok,k2blok,l2blok,m2blok,kxyp2d,kyzp2d,nz1d,nxhyz2,nxyzh2)
c this subroutine calculates the form factor array ffg needed by field
c solvers with open (vacuum) boundary conditions using hockney's method.
c the five green's functions calculated are:
c g(kx,ky,kz) = affp*inverse FFT of potr
c s(kx,ky,kz) = inverse FFT of the density of a finite-sized particle
c gx(kx,ky,kz) = affp*s(kx,ky,kz)*inverse FFT of (x/r)*Er
c gy(kx,ky,kz) = affp*s(kx,ky,kz)*inverse FFT of (y/r)*Er
c gz(kx,ky,kz) = affp*s(kx,ky,kz)*inverse FFT of (z/r)*Er
c where the fields due to the finite-sized particles are given by fpotc
c input: fpotc,mixup3,sct3,affp,ar,indx1,indy1,indz1,kstrt,nxv,ny2d,nz2d
c        kxyp2,kyp2,kyzp2,kzp2,kyp2d,j2blok,k2blok,l2blok,m2blok,kxyp2d
c        nz1d,nxhyz2,nxyzh2
c ffg(1,l,j,k) = potential green's function g
c ffg(2,l,j,k) = finite-size particle shape factor s
c ffg(3,l,j,k) = x component of electric field green's function gx
c ffg(4,l,j,k) = y component of electric field green's function gy
c ffg(5,l,j,k) = z component of electric field green's function gz
c on processor 0, ffg(i,l,k,kxp2+1,m) = ffg(i,l,k,NX+1,m)
c on other processors, ffg(i,l,kxp2+1,k,m) = ffg(i,l,1,k,m) on proc 0
c f, fs, ft = scratch arrays used by FFT
c g, h = scratch arrays used for message-passing
c fpotc = a function which calculates green's function
c mixup3/sct3 = bit-reverse and sine-cosine table used by FFT
c affp = normalization constant = nx*ny*nz/np, where np=number of
c particles
c ar = half-width of particle in r direction
c indx1/indy1/indz1 = exponent which determines FFT length in x/y/z
c direction, where 2*nx=2**indx1, 2*ny=2**indy1, 2*nz=2**indz1
c kstrt = starting data block number
c nxv = half of first dimension of field arrays, must be >= nx
c ny2d = second dimension of field arrays, must be >= 2*ny
c nz2d = third dimension of field arrays, must be >= 2*nz
c kyp2d = second dimension of f
c kxyp2/kyp2,kyzp2/kzp2 = number of data values per block in x/y/z
c j2blok/k2blok,m2blok/l2blok = number of data blocks in x/y/z
c kxyp2d = third dimension of ffg arrays, must be >= kxyp2+1
c kyzp2d = third dimension of ffg arrays, must be >= kyzp2+1
c nz1d = second dimension of ffg arrays, must be >= nz+1
c nxhyz2 = maximum of (nx,2*ny,2*nz)
c nxyzh2 = maximum of (nx,ny,nz)
      implicit none
      real ffg, f, g, h
      complex fs, ft
      integer mixup3
      complex sct3
      real affp, ar
      integer indx1, indy1, indz1, kstrt, nxv, ny2d, nz2d, kxyp2, kyp2
      integer kyzp2, kzp2, kyp2d, j2blok, k2blok, l2blok, m2blok
      integer kxyp2d, kyzp2d, nz1d, nxhyz2, nxyzh2
      dimension ffg(5,nz1d,kxyp2d,kyzp2d,j2blok*m2blok)
      dimension f(2*nxv,kyp2d,kzp2,k2blok*l2blok)
      dimension ft(nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension fs(ny2d,kxyp2,kzp2,j2blok*l2blok)
      dimension g(5,nz1d,kyzp2,j2blok*m2blok)
      dimension h(5,nz1d,kyzp2,j2blok*m2blok)
      dimension mixup3(nxhyz2), sct3(nxyzh2)
      real fpotc
      external fpotc
c local data
      integer ntpose, nx, ny, nz, nz1, nx2, ny2, nz2, kxb2, kyb2, kzb2
      integer kyzb2, isign, i, j, k, l, mx, my, mz, m, j1, k1, l1, n2
      integer js, ks, joff, koff, loff, moff, ifun
      real an, ari, at1, at2, x, y, z, r, ttp
      real POTC3
      external POTC3
      data ntpose /1/
      nx2 = 2**indx1
      ny2 = 2**indy1
      nz2 = 2**indz1
      nx = nx2/2
      ny = ny2/2
      nz = nz2/2
      nz1 = nz + 1
      kyb2 = ny2/kyp2
      kzb2 = nz2/kzp2
      kxb2 = nx/kxyp2
      kyzb2 = ny2/kyzp2
      ari = 0.0
      if (ar.gt.0.) ari = 1.0/ar
      an = float(nx2*ny2*nz2)
c calculate potential green's function
      ifun = 1
      if (kstrt.gt.(kyb2*kzb2)) go to 60
      ks = (kstrt - 1)/kyb2
      js = kstrt - kyb2*ks - 2
      ks = ks - 1
      do 50 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      loff = kzp2*(mz + ks) - 1
      do 40 my = 1, k2blok
      koff = kyp2*(my + js) - 1
      m = my + moff
      do 30 l = 1, kzp2
      l1 = l + loff
      if (l1.gt.nz) l1 = l1 - nz2
      at1 = float(l1)**2
      do 20 k = 1, kyp2
      k1 = k + koff
      if (k1.gt.ny) k1 = k1 - ny2
      at2 = float(k1)**2 + at1
      do 10 j = 1, nx2
      j1 = j - 1
      if (j1.gt.nx) j1 = j1 - nx2
      r = sqrt(at2 + float(j1)**2)
      f(j,k,l,m) = fpotc(r,affp,ari,1)
   10 continue
   20 continue
   30 continue
   40 continue
   50 continue
      isign = -1
      call WPFFT32RX(f,fs,ft,isign,ntpose,mixup3,sct3,ttp,indx1,indy1,in
     1dz1,kstrt,nxv,ny2d,nz2d,kxyp2,kyp2,kyzp2,kzp2,kxyp2,kyp2d,kyzp2,kz
     2p2,j2blok,k2blok,l2blok,m2blok,nxhyz2,nxyzh2)
   60 if (kstrt.gt.(kxb2*kyzb2)) go to 150
      ks = (kstrt - 1)/kxb2
      js = kstrt - kxb2*ks - 2
      ks = ks - 1
      do 140 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 130 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js)
      m = mx + moff
      do 100 k = 1, kyzp2
      k1 = k + koff
      do 80 j = 1, kxyp2
      do 70 l = 1, nz1
      ffg(ifun,l,j,k,m) = an*real(ft(l,j,k,m))
   70 continue
   80 continue
      do 90 l = 1, nz1
      ffg(ifun,l,kxyp2+1,k,m) = 0.0
   90 continue
  100 continue
c mode numbers kx = 0, nx
      if (joff.eq.0) then
c mode number ky = 0
         n2 = koff + 1
         if (n2.eq.0) then
            do 110 l = 2, nz
            l1 = nz2 + 2 - l
            ffg(ifun,l,kxyp2+1,1,m) = an*real(ft(l1,1,1,m))
  110       continue
            ffg(ifun,1,kxyp2+1,1,m) = an*aimag(ft(1,1,1,m))
            ffg(ifun,nz1,kxyp2+1,1,m) = an*aimag(ft(nz1,1,1,m))
         endif
c mode number ky = ny
         k1 = (ny/kyzp2)*kyzp2
         if (n2.eq.k1) then
            k1 = ny - k1 + 1
            do 120 l = 2, nz
            l1 = nz2 + 2 - l
            ffg(ifun,l,kxyp2+1,k1,m) = an*real(ft(l1,1,k1,m))
  120       continue
            ffg(ifun,1,kxyp2+1,k1,m) = an*aimag(ft(1,1,k1,m))
            ffg(ifun,nz1,kxyp2+1,k1,m) = an*aimag(ft(nz1,1,k1,m))
         endif
      endif
  130 continue
  140 continue
c calculate particle smoothing function
  150 ifun = ifun + 1
      if (kstrt.gt.(kyb2*kzb2)) go to 210
      ks = (kstrt - 1)/kyb2
      js = kstrt - kyb2*ks - 2
      ks = ks - 1
      do 200 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      loff = kzp2*(mz + ks) - 1
      do 190 my = 1, k2blok
      koff = kyp2*(my + js) - 1
      m = my + moff
      do 180 l = 1, kzp2
      l1 = l + loff
      if (l1.gt.nz) l1 = l1 - nz2
      at1 = float(l1)**2
      do 170 k = 1, kyp2
      k1 = k + koff
      if (k1.gt.ny) k1 = k1 - ny2
      at2 = float(k1)**2 + at1
      do 160 j = 1, nx2
      j1 = j - 1
      if (j1.gt.nx) j1 = j1 - nx2
      r = sqrt(at2 + float(j1)**2)
      f(j,k,l,m) = fpotc(r,affp,ari,2)
  160 continue
  170 continue
  180 continue
  190 continue
  200 continue
      isign = -1
      call WPFFT32RX(f,fs,ft,isign,ntpose,mixup3,sct3,ttp,indx1,indy1,in
     1dz1,kstrt,nxv,ny2d,nz2d,kxyp2,kyp2,kyzp2,kzp2,kxyp2,kyp2d,kyzp2,kz
     2p2,j2blok,k2blok,l2blok,m2blok,nxhyz2,nxyzh2)
  210 if (kstrt.gt.(kxb2*kyzb2)) go to 300
      ks = (kstrt - 1)/kxb2
      js = kstrt - kxb2*ks - 2
      ks = ks - 1
      do 290 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 280 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js)
      m = mx + moff
      do 250 k = 1, kyzp2
      k1 = k + koff
      do 230 j = 1, kxyp2
      do 220 l = 1, nz1
      ffg(ifun,l,j,k,m) = an*real(ft(l,j,k,m))
  220 continue
  230 continue
      do 240 l = 1, nz1
      ffg(ifun,l,kxyp2+1,k,m) = 0.0
  240 continue
  250 continue
c mode numbers kx = 0, nx
      if (joff.eq.0) then
c mode number ky = 0
         n2 = koff + 1
         if (n2.eq.0) then
            do 260 l = 2, nz
            l1 = nz2 + 2 - l
            ffg(ifun,l,kxyp2+1,1,m) = an*real(ft(l1,1,1,m))
  260       continue
            ffg(ifun,1,kxyp2+1,1,m) = an*aimag(ft(1,1,1,m))
            ffg(ifun,nz1,kxyp2+1,1,m) = an*aimag(ft(nz1,1,1,m))
         endif
c mode number ky = ny
         k1 = (ny/kyzp2)*kyzp2
         if (n2.eq.k1) then
            k1 = ny - k1 + 1
            do 270 l = 2, nz
            l1 = nz2 + 2 - l
            ffg(ifun,l,kxyp2+1,k1,m) = an*real(ft(l1,1,k1,m))
  270       continue
            ffg(ifun,1,kxyp2+1,k1,m) = an*aimag(ft(1,1,k1,m))
            ffg(ifun,nz1,kxyp2+1,k1,m) = an*aimag(ft(nz1,1,k1,m))
         endif
      endif
  280 continue
  290 continue
c calculate green's function for x component of electric field
  300 ifun = ifun + 1
      if (kstrt.gt.(kyb2*kzb2)) go to 360
      ks = (kstrt - 1)/kyb2
      js = kstrt - kyb2*ks - 2
      ks = ks - 1
      do 350 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      loff = kzp2*(mz + ks) - 1
      do 340 my = 1, k2blok
      koff = kyp2*(my + js) - 1
      m = my + moff
      do 330 l = 1, kzp2
      l1 = l + loff
      if (l1.gt.nz) l1 = l1 - nz2
      at1 = float(l1)**2
      do 320 k = 1, kyp2
      k1 = k + koff
      if (k1.gt.ny) k1 = k1 - ny2
      at2 = float(k1)**2 + at1
      do 310 j = 1, nx2
      j1 = j - 1
      if (j1.gt.nx) j1 = j1 - nx2
      x = float(j1)
      r = sqrt(at2 + x*x)
      f(j,k,l,m) = fpotc(r,affp,ari,3)
      if (r.gt.0.) f(j,k,l,m) = f(j,k,l,m)*(x/r)
  310 continue
  320 continue
  330 continue
  340 continue
  350 continue
      isign = -1
      call WPFFT32RX(f,fs,ft,isign,ntpose,mixup3,sct3,ttp,indx1,indy1,in
     1dz1,kstrt,nxv,ny2d,nz2d,kxyp2,kyp2,kyzp2,kzp2,kxyp2,kyp2d,kyzp2,kz
     2p2,j2blok,k2blok,l2blok,m2blok,nxhyz2,nxyzh2)
  360 if (kstrt.gt.(kxb2*kyzb2)) go to 460
      ks = (kstrt - 1)/kxb2
      js = kstrt - kxb2*ks - 2
      ks = ks - 1
      do 450 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 440 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js)
      m = mx + moff
      do 410 k = 1, kyzp2
      k1 = k + koff
      do 380 j = 1, kxyp2
      if ((j+joff).gt.0) then
         do 370 l = 1, nz1
         ffg(ifun,l,j,k,m) = an*aimag(ft(l,j,k,m))
  370    continue
      endif
  380 continue
c mode numbers kx = 0, nx
      if (joff.eq.0) then
         do 390 l = 1, nz1
         ffg(ifun,l,1,k,m) = an*real(ft(l,1,k,m))
  390    continue
      endif
      do 400 l = 1, nz1
      ffg(ifun,l,kxyp2+1,k,m) = 0.0
  400 continue
  410 continue
c mode numbers kx = 0, nx
      if (joff.eq.0) then
c mode number ky = 0
         n2 = koff + 1
         if (n2.eq.0) then
            do 420 l = 2, nz
            l1 = nz2 + 2 - l
            ffg(ifun,l,kxyp2+1,1,m) = an*real(ft(l1,1,1,m))
  420       continue
            ffg(ifun,1,kxyp2+1,1,m) = an*aimag(ft(1,1,1,m))
            ffg(ifun,nz1,kxyp2+1,1,m) = an*aimag(ft(nz1,1,1,m))
         endif
c mode number ky = ny
         k1 = (ny/kyzp2)*kyzp2
         if (n2.eq.k1) then
            k1 = ny - k1 + 1
            do 430 l = 2, nz
            l1 = nz2 + 2 - l
            ffg(ifun,l,kxyp2+1,k1,m) = an*real(ft(l1,1,k1,m))
  430       continue
            ffg(ifun,1,kxyp2+1,k1,m) = an*aimag(ft(1,1,k1,m))
            ffg(ifun,nz1,kxyp2+1,k1,m) = an*aimag(ft(nz1,1,k1,m))
         endif
      endif
  440 continue
  450 continue
c calculate green's function for y component of electric field
  460 ifun = ifun + 1
      if (kstrt.gt.(kyb2*kzb2)) go to 520
      ks = (kstrt - 1)/kyb2
      js = kstrt - kyb2*ks - 2
      ks = ks - 1
      do 510 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      loff = kzp2*(mz + ks) - 1
      do 500 my = 1, k2blok
      koff = kyp2*(my + js) - 1
      m = my + moff
      do 490 l = 1, kzp2
      l1 = l + loff
      if (l1.gt.nz) l1 = l1 - nz2
      at1 = float(l1)**2
      do 480 k = 1, kyp2
      k1 = k + koff
      if (k1.gt.ny) k1 = k1 - ny2
      y = float(k1)
      at2 = y*y + at1
      do 470 j = 1, nx2
      j1 = j - 1
      if (j1.gt.nx) j1 = j1 - nx2
      r = sqrt(at2 + float(j1)**2)
      f(j,k,l,m) = fpotc(r,affp,ari,3)
      if (r.gt.0.) f(j,k,l,m) = f(j,k,l,m)*(y/r)
  470 continue
  480 continue
  490 continue
  500 continue
  510 continue
      isign = -1
      call WPFFT32RX(f,fs,ft,isign,ntpose,mixup3,sct3,ttp,indx1,indy1,in
     1dz1,kstrt,nxv,ny2d,nz2d,kxyp2,kyp2,kyzp2,kzp2,kxyp2,kyp2d,kyzp2,kz
     2p2,j2blok,k2blok,l2blok,m2blok,nxhyz2,nxyzh2)
  520 if (kstrt.gt.(kxb2*kyzb2)) go to 650
      ks = (kstrt - 1)/kxb2
      js = kstrt - kxb2*ks - 2
      ks = ks - 1
      do 640 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 630 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js)
      m = mx + moff
      do 560 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         do 540 j = 1, kxyp2
         do 530 l = 1, nz1
         ffg(ifun,l,j,k,m) = an*aimag(ft(l,j,k,m))
  530    continue
  540    continue
      endif
      do 550 l = 1, nz1
      ffg(ifun,l,kxyp2+1,k,m) = 0.0
  550 continue
  560 continue
c mode number ky = 0
      n2 = koff + 1
      if (n2.eq.0) then
         do 580 j = 1, kxyp2
         do 570 l = 1, nz1
         ffg(ifun,l,j,1,m) = an*real(ft(l,j,1,m))
  570    continue
  580    continue
      endif
c mode numbers ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 600 j = 1, kxyp2
         do 590 l = 1, nz1
         ffg(ifun,l,j,k1,m) = an*real(ft(l,j,k1,m))
  590    continue
  600    continue
      endif
c mode numbers kx = 0, nx
      if (joff.eq.0) then
c mode number ky = 0
         n2 = koff + 1
         if (n2.eq.0) then
            do 610 l = 2, nz
            l1 = nz2 + 2 - l
            ffg(ifun,l,kxyp2+1,1,m) = an*real(ft(l1,1,1,m))
  610       continue
            ffg(ifun,1,kxyp2+1,1,m) = an*aimag(ft(1,1,1,m))
            ffg(ifun,nz1,kxyp2+1,1,m) = an*aimag(ft(nz1,1,1,m))
         endif
c mode number ky = ny
         k1 = (ny/kyzp2)*kyzp2
         if (n2.eq.k1) then
            k1 = ny - k1 + 1
            do 620 l = 2, nz
            l1 = nz2 + 2 - l
            ffg(ifun,l,kxyp2+1,k1,m) = an*real(ft(l1,1,k1,m))
  620       continue
            ffg(ifun,1,kxyp2+1,k1,m) = an*aimag(ft(1,1,k1,m))
            ffg(ifun,nz1,kxyp2+1,k1,m) = an*aimag(ft(nz1,1,k1,m))
         endif
      endif
  630 continue
  640 continue
c calculate green's function for z component of electric field
  650 ifun = ifun + 1
      if (kstrt.gt.(kyb2*kzb2)) go to 710
      ks = (kstrt - 1)/kyb2
      js = kstrt - kyb2*ks - 2
      ks = ks - 1
      do 700 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      loff = kzp2*(mz + ks) - 1
      do 690 my = 1, k2blok
      koff = kyp2*(my + js) - 1
      m = my + moff
      do 680 l = 1, kzp2
      l1 = l + loff
      if (l1.gt.nz) l1 = l1 - nz2
      z = float(l1)
      at1 = z*z
      do 670 k = 1, kyp2
      k1 = k + koff
      if (k1.gt.ny) k1 = k1 - ny2
      at2 = float(k1)**2 + at1
      do 660 j = 1, nx2
      j1 = j - 1
      if (j1.gt.nx) j1 = j1 - nx2
      r = sqrt(at2 + float(j1)**2)
      f(j,k,l,m) = fpotc(r,affp,ari,3)
      if (r.gt.0.) f(j,k,l,m) = f(j,k,l,m)*(z/r)
  660 continue
  670 continue
  680 continue
  690 continue
  700 continue
      isign = -1
      call WPFFT32RX(f,fs,ft,isign,ntpose,mixup3,sct3,ttp,indx1,indy1,in
     1dz1,kstrt,nxv,ny2d,nz2d,kxyp2,kyp2,kyzp2,kzp2,kxyp2,kyp2d,kyzp2,kz
     2p2,j2blok,k2blok,l2blok,m2blok,nxhyz2,nxyzh2)
  710 if (kstrt.gt.(kxb2*kyzb2)) go to 840
      ks = (kstrt - 1)/kxb2
      js = kstrt - kxb2*ks - 2
      ks = ks - 1
      do 830 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 820 mx = 1, j2blok
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js)
      m = mx + moff
      do 750 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         do 730 j = 1, kxyp2
         do 720 l = 2, nz
         l1 = nz2 + 2 - l
         ffg(ifun,l,j,k,m) = an*aimag(ft(l,j,k,m))
  720    continue
         ffg(ifun,1,j,k,m) = an*real(ft(1,j,k,m))
         ffg(ifun,nz1,j,k,m) = an*real(ft(nz1,j,k,m))
  730    continue
      endif
      do 740 l = 1, nz1
      ffg(ifun,l,kxyp2+1,k,m) = 0.0
  740 continue
  750 continue
c mode number ky = 0
      n2 = koff + 1
      if (n2.eq.0) then
         do 770 j = 1, kxyp2
         do 760 l = 2, nz
         l1 = nz2 + 2 - l
         ffg(ifun,l,j,1,m) = an*aimag(ft(l,j,1,m))
  760    continue
         ffg(ifun,1,j,1,m) = an*real(ft(1,j,1,m))
         ffg(ifun,nz1,j,1,m) = an*real(ft(nz1,j,1,m))
  770    continue
      endif
c mode numbers ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 790 j = 1, kxyp2
         do 780 l = 2, nz
         l1 = nz2 + 2 - l
         ffg(ifun,l,j,k1,m) = an*aimag(ft(l,j,k1,m))
  780    continue
         ffg(ifun,1,j,k1,m) = an*real(ft(1,j,k1,m))
         ffg(ifun,nz1,j,k1,m) = an*real(ft(nz1,j,k1,m))
  790    continue
      endif
c mode numbers kx = 0, nx
      if (joff.eq.0) then
c mode number ky = 0
         n2 = koff + 1
         if (n2.eq.0) then
            do 800 l = 2, nz
            l1 = nz2 + 2 - l
            ffg(ifun,l,kxyp2+1,1,m) = an*aimag(ft(l1,1,1,m))
  800       continue
            ffg(ifun,1,kxyp2+1,1,m) = an*aimag(ft(1,1,1,m))
            ffg(ifun,nz1,kxyp2+1,1,m) = an*aimag(ft(nz1,1,1,m))
         endif
c mode number ky = ny
         k1 = (ny/kyzp2)*kyzp2
         if (n2.eq.k1) then
            k1 = ny - k1 + 1
            do 810 l = 2, nz
            l1 = nz2 + 2 - l
            ffg(ifun,l,kxyp2+1,k1,m) = an*aimag(ft(l1,1,k1,m))
  810       continue
            ffg(ifun,1,kxyp2+1,k1,m) = an*aimag(ft(1,1,k1,m))
            ffg(ifun,nz1,kxyp2+1,k1,m) = an*aimag(ft(nz1,1,k1,m))
         endif
      endif
  820 continue
  830 continue
c copy ffg(i,l,1,k,m) on node 0 to ffg(i,l,kxyp2+1,k,m) on other nodes
  840 do 910 m = 1, j2blok*m2blok
      do 870 k = 1, kyzp2
      do 860 l = 1, nz1d
      do 850 i = 1, 5
      g(i,l,k,m) = ffg(i,l,1,k,m)
      h(i,l,k,m) = ffg(i,l,kxyp2+1,k,m)
  850 continue
  860 continue
  870 continue
      call P0COPY2(g,h,5*kyzp2*nz1d,kstrt,kxb2,kyzb2,1)
      do 900 k = 1, kyzp2
      do 890 l = 1, nz1d
      do 880 i = 1, 5
      ffg(i,l,kxyp2+1,k,m) = h(i,l,k,m)
  880 continue
  890 continue
  900 continue
  910 continue
c copy ffg(i,l,k,1,m) on node 0 to ffg(i,l,k,kxp2+1,m) on other nodes
      do 920 m = 1, j2blok*m2blok
      call P0COPY2(ffg(1,1,1,1,m),ffg(1,1,1,kyzp2+1,m),5*kxyp2d*nz1d,kst
     1rt,kxb2,kyzb2,2)
  920 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISC32(q,fx,fy,fz,isign,ffg,we,nx,ny,nz,kstrt,nz2d,kx
     1yp2,kyzp2,j2blok,m2blok,kxyp2d,kyzp2d,nz1d)
c this subroutine solves 3d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c or for potential, or provides smoothing a smoothing function,
c with open (vacuum) boundary conditions using hockney's method,
c for distributed data with 2D spatial decomposition
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate convolution
c for isign = -1, input: q,ffg,isign,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2,
c                        j2blok,m2blok,kxyp2d,nz1d
c output: fx,fy,fz,we
c approximate flop count is:
c        57*nxc*nyc*nzc + 50*(nxc*nyc + nxc*nzc + nyc*nzc)
c for isign = 1, input: q,ffg,isign,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2,
c                       j2blok,m2blok,kxyp2d,nz1d
c output: fx,we
c approximate flop count is:
c        22*nxc*nyc*nzc + 22*(nxc*nyc + nxc*nzc + nyc*nzc)
c for isign = 2, input: q,ffg,isign,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2,
c                       j2blok,m2blok,kxyp2d,nz1d
c output: fy
c approximate flop count is:
c        4*nxc*nyc*nzc + 4*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx-1)/nvpy, nyc = (ny-1)/nvpz, nzc = nz - 1, and
c nvpy/nvpz = number of procs in y/z
c if isign < 0, force/charge is calculated using the equations:
c fx(kx,ky,kz) = gx(kx,ky,kz)*s(kx,ky,kz)*q(kx,ky,kz),
c fy(kx,ky,kz) = gy(kx,ky,kz)*s(kx,ky,kz)*q(kx,ky,kz),
c fz(kx,ky,kz) = gz(kx,ky,kz)*s(kx,ky,kz)*q(kx,ky,kz),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz,
c and j,k,l = fourier mode numbers,
c gx(kx,ky,kz) = s(kx,ky,kz)*inverse FFT of (x/r)*Er
c gy(kx,ky,kz) = s(kx,ky,kz)*inverse FFT of (y/r)*Er
c gz(kx,ky,kz) = s(kx,ky,kz)*inverse FFT of (z/r)*Er
c where Er is the electric field of a single finite-sized particle
c s(kx,ky,kz) = inverse FFT of the density of a finite-sized particle
c if isign = 1, potential is calculated using the equation:
c fx(kx,ky,kz) = g(kx,ky,kz)*q(kx,ky,kz)
c where g(kx,ky,kz) = affp*inverse FFT of potr
c where potr is the potential of a single finite-sized particle
c if isign = 2, smoothing is calculated using the equation:
c fy(kx,ky,kz) = q(kx,ky,kz)*s(kx,ky,kz)
c q(l,j,k,m) = complex charge density for fourier mode jj-1,kk-1,l-1
c fx(l,j,k,m) = x component of force/charge
c fy(l,j,k,m) = y component of force/charge
c fz(l,j,k,m) = z component of force/charge
c ffg(1,l,j,k) = potential green's function g
c ffg(2,l,j,k) = finite-size particle shape factor s
c ffg(3,l,j,k) = x component of electric field green's function gx
c ffg(4,l,j,k) = y component of electric field green's function gy
c ffg(5,l,j,k) = z component of electric field green's function gz
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp2*(mx - 1) and
c kk = k + kyzp2*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c the ffg array is calculated by the subroutine PFORMC32
c electric field energy is also calculated and returned in we
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c nz2d = first dimension of field arrays, must be >= 2*nz
c kxyp2/kyzp2 = number of complex grids in each field partition in
c x/y direction
c j2blok/m2blok = number of field partitions in x/y
c kxyp2d = third dimension of ffg arrays, must be >= kxyp2+1
c kyzp2d = third dimension of ffg arrays, must be >= kyzp2+1
c nz1d = second dimension of ffg arrays, must be >= nz+1
      implicit none
      real ffg
      complex q, fx, fy, fz
      integer isign, nx, ny, nz, kstrt, nz2d, kxyp2, kyzp2
      integer j2blok, m2blok, kxyp2d, kyzp2d, nz1d
      real we
      dimension q(nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension fx(nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension fy(nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension fz(nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension ffg(5,nz1d,kxyp2d,kyzp2d,j2blok*m2blok)
c local data
      double precision wp
      integer j, k, l, mx, my, m, nz22, kxb2, kyzb2, js, ks, kx1, ky1
      integer kxy1, joff, koff, moff, k1, l1, n1, n2, nz1
      real at1, at2, at3, at4, at5, at6, at7, at8
      complex zt1, zt2, zt3
      if (isign.eq.0) return
      nz1 = nz + 1
      nz22 = nz + nz + 2
      kxb2 = nx/kxyp2
      kyzb2 = (ny + ny)/kyzp2
      ks = (kstrt - 1)/kxb2
      js = kstrt - kxb2*ks - 2
      ks = ks - 1
      kx1 = 1
      ky1 = 1
      kxy1 = 1
      if (isign.gt.0) go to 140
c calculate force/charge and sum field energy
      wp = 0.0d0
      if (kstrt.gt.(kxb2*kyzb2)) go to 130
      do 120 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      if ((my+ks).gt.0) ky1 = kyzp2 + 1
      do 110 mx = 1, j2blok
      m = mx + moff
      if ((mx+js).gt.0) kx1 = kxyp2 + 1
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      at5 = -1.0
      do 40 k = 1, kyzp2
      k1 = k + koff
      if (k1.ge.ny) kxy1 = kxyp2 + 1
      at5 = -at5
      if ((k1.gt.0).and.(k1.ne.ny)) then
         at2 = -1.0
         do 20 j = 1, kxyp2
         at2 = -at2
         if ((j+joff).gt.0) then
            at7 = ffg(5,1,j,k,m)
            do 10 l = 2, nz
            l1 = nz22 - l
            at1 = at2*ffg(3,l,kx1,k,m)
            zt1 = cmplx(at1,ffg(3,l,j,k,m))
            fx(l,j,k,m) = zt1*q(l,j,k,m)
            fx(l1,j,k,m) = zt1*q(l1,j,k,m)
            at3 = at5*ffg(4,l,j,ky1,m)
            zt2 = cmplx(at3,ffg(4,l,j,k,m))
            fy(l,j,k,m) = zt2*q(l,j,k,m)
            fy(l1,j,k,m) = zt2*q(l1,j,k,m)
            at7 = -at7
            zt3 = cmplx(at7,ffg(5,l,j,k,m))
            fz(l,j,k,m) = zt3*q(l,j,k,m)
            fz(l1,j,k,m) = conjg(zt3)*q(l1,j,k,m)
            wp = wp + ffg(1,l,j,k,m)*(q(l,j,k,m)*conjg(q(l,j,k,m)) + q(l
     11,j,k,m)*conjg(q(l1,j,k,m)))
   10       continue
c mode numbers kz = 0, nz
            at1 = at2*ffg(3,1,kx1,k,m)
            zt1 = cmplx(at1,ffg(3,1,j,k,m))
            fx(1,j,k,m) = zt1*q(1,j,k,m)
            at3 = at5*ffg(4,1,j,ky1,m)
            zt2 = cmplx(at3,ffg(4,1,j,k,m))
            fy(1,j,k,m) = zt2*q(1,j,k,m)
            at8 = ffg(5,1,j,k,m)
            fz(1,j,k,m) = at8*q(1,j,k,m)
            wp = wp + ffg(1,1,j,k,m)*q(1,j,k,m)*conjg(q(1,j,k,m))
            at1 = at2*ffg(3,nz1,kx1,k,m)
            zt1 = cmplx(at1,ffg(3,nz1,j,k,m))
            fx(nz1,j,k,m) = zt1*q(nz1,j,k,m)
            at3 = at5*ffg(4,nz1,j,ky1,m)
            zt2 = cmplx(at3,ffg(4,nz1,j,k,m))
            fy(nz1,j,k,m) = zt2*q(nz1,j,k,m)
            at8 = ffg(5,nz1,j,k,m)
            fz(nz1,j,k,m) = at8*q(nz1,j,k,m)
            wp = wp + ffg(1,nz1,j,k,m)*q(nz1,j,k,m)*conjg(q(nz1,j,k,m))
         endif
   20    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            at7 = ffg(5,1,1,k,m)
            do 30 l = 2, nz
            l1 = nz22 - l
            at1 = ffg(3,l,1,k,m)
            fx(l,1,k,m) = at1*q(l,1,k,m)
            fx(l1,1,k,m) = at1*q(l1,1,k,m)
            at3 = at5*ffg(4,l,kxy1,ky1,m)
            zt2 = cmplx(at3,ffg(4,l,1,k,m))
            fy(l,1,k,m) = zt2*q(l,1,k,m)
            fy(l1,1,k,m) = zt2*q(l1,1,k,m)
            at7 = -at7
            zt3 = cmplx(at7,ffg(5,l,1,k,m))
            fz(l,1,k,m) = zt3*q(l,1,k,m)
            fz(l1,1,k,m) = conjg(zt3)*q(l1,1,k,m)
            wp = wp + ffg(1,l,1,k,m)*(q(l,1,k,m)*conjg(q(l,1,k,m)) + q(l
     11,1,k,m)*conjg(q(l1,1,k,m)))
   30       continue
c mode numbers kz = 0
            at2 = ffg(3,1,1,k,m)
            fx(1,1,k,m) = at2*q(1,1,k,m)
            at3 = at5*ffg(4,1,kxy1,ky1,m)
            zt2 = cmplx(at3,ffg(4,1,1,k,m))
            fy(1,1,k,m) = zt2*q(1,1,k,m)
            at8 = ffg(5,1,1,k,m)
            fz(1,1,k,m) = at8*q(1,1,k,m)
            wp = wp + ffg(1,1,1,k,m)*q(1,1,k,m)*conjg(q(1,1,k,m))
c mode numbers kz = nz
            at2 = ffg(3,nz1,1,k,m)
            fx(nz1,1,k,m) = at2*q(nz1,1,k,m)
            at3 = at5*ffg(4,nz1,kxy1,ky1,m)
            zt2 = cmplx(at3,ffg(4,nz1,1,k,m))
            fy(nz1,1,k,m) = zt2*q(nz1,1,k,m)
            at8 = ffg(5,nz1,1,k,m)
            fz(nz1,1,k,m) = at8*q(nz1,1,k,m)
            wp = wp + ffg(1,nz1,1,k,m)*q(nz1,1,k,m)*conjg(q(nz1,1,k,m))
         endif
      endif
   40 continue
c mode numbers ky = 0
      n2 = koff + 1
      if (n2.eq.0) then
         at2 = -1.0
         do 60 j = 1, kxyp2
         at2 = -at2
         if ((j+joff).gt.0) then
            at7 = ffg(5,1,j,1,m)
            do 50 l = 2, nz
            l1 = nz22 - l
            at1 = at2*ffg(3,l,kx1,1,m)
            zt1 = cmplx(at1,ffg(3,l,j,1,m))
            fx(l,j,1,m) = zt1*q(l,j,1,m)
            fx(l1,j,1,m) = zt1*q(l1,j,1,m)
            at3 = ffg(4,l,j,1,m)
            fy(l,j,1,m) = at3*q(l,j,1,m)
            fy(l1,j,1,m) = at3*q(l1,j,1,m)
            at7 = -at7
            zt3 = cmplx(at7,ffg(5,l,j,1,m))
            fz(l,j,1,m) = zt3*q(l,j,1,m)
            fz(l1,j,1,m) = conjg(zt3)*q(l1,j,1,m)
            wp = wp + ffg(1,l,j,1,m)*(q(l,j,1,m)*conjg(q(l,j,1,m)) + q(l
     11,j,1,m)*conjg(q(l1,j,1,m)))
   50       continue
c mode numbers kz = 0
            at1 = at2*ffg(3,1,kx1,1,m)
            zt1 = cmplx(at1,ffg(3,1,j,1,m))
            fx(1,j,1,m) = zt1*q(1,j,1,m)
            at4 = ffg(4,1,j,1,m)
            fy(1,j,1,m) = at4*q(1,j,1,m)
            at8 = ffg(5,1,j,1,m)
            fz(1,j,1,m) = at8*q(1,j,1,m)
            wp = wp + ffg(1,1,j,1,m)*q(1,j,1,m)*conjg(q(1,j,1,m))
c mode numbers kz = nz
            at1 = at2*ffg(3,nz1,kx1,1,m)
            zt1 = cmplx(at1,ffg(3,nz1,j,1,m))
            fx(nz1,j,1,m) = zt1*q(nz1,j,1,m)
            at4 = ffg(4,nz1,j,1,m)
            fy(nz1,j,1,m) = at4*q(nz1,j,1,m)
            at8 = ffg(5,nz1,j,1,m)
            fz(nz1,j,1,m) = at8*q(nz1,j,1,m)
            wp = wp + ffg(1,nz1,j,1,m)*q(nz1,j,1,m)*conjg(q(nz1,j,1,m))
         endif
   60    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            kx1 = kxyp2 + 1
            at7 = ffg(5,1,1,1,m)
            at8 = ffg(5,1,kx1,1,m)
            do 70 l = 2, nz
            l1 = nz22 - l
c mode numbers kx = 0
            at2 = ffg(3,l,1,1,m)
            fx(l,1,1,m) = at2*q(l,1,1,m)
            at4 = ffg(4,l,1,1,m)
            fy(l,1,1,m) = at4*q(l,1,1,m)
            at7 = -at7
            zt3 = cmplx(at7,ffg(5,l,1,1,m))
            fz(l,1,1,m) = zt3*q(l,1,1,m)
            wp = wp + ffg(1,l,1,1,m)*q(l,1,1,m)*conjg(q(l,1,1,m))
c mode numbers kx = nx
            at2 = ffg(3,l,kx1,1,m)
            fx(l1,1,1,m) = at2*q(l1,1,1,m)
            at4 = ffg(4,l,kx1,1,m)
            fy(l1,1,1,m) = at4*q(l1,1,1,m)
            at8 = -at8
            zt3 = cmplx(at8,ffg(5,l,kx1,1,m))
            fz(l1,1,1,m) = zt3*q(l1,1,1,m)
            wp = wp + ffg(1,l,kx1,1,m)*q(l1,1,1,m)*conjg(q(l1,1,1,m))
   70       continue
c mode numbers kx = 0, nx/2, kz = 0
            fx(1,1,1,m) = cmplx(ffg(3,1,1,1,m)*real(q(1,1,1,m)),ffg(3,1,
     1kx1,1,m)*aimag(q(1,1,1,m)))
            fy(1,1,1,m) = cmplx(ffg(4,1,1,1,m)*real(q(1,1,1,m)),ffg(4,1,
     1kx1,1,m)*aimag(q(1,1,1,m)))
            fz(1,1,1,m) = cmplx(ffg(5,1,1,1,m)*real(q(1,1,1,m)),ffg(5,1,
     1kx1,1,m)*aimag(q(1,1,1,m)))
            wp = wp + .5*(ffg(1,1,1,1,m)*real(q(1,1,1,m))**2 + ffg(1,1,k
     1x1,1,m)*aimag(q(1,1,1,m))**2)
c mode numbers kx = 0, nx/2, kz = nz/2
            fx(nz1,1,1,m) = cmplx(ffg(3,nz1,1,1,m)*real(q(nz1,1,1,m)),ff
     1g(3,nz1,kx1,1,m)*aimag(q(nz1,1,1,m)))
            fy(nz1,1,1,m) = cmplx(ffg(4,nz1,1,1,m)*real(q(nz1,1,1,m)),ff
     1g(4,nz1,kx1,1,m)*aimag(q(nz1,1,1,m)))
            fz(nz1,1,1,m) = cmplx(ffg(5,nz1,1,1,m)*real(q(nz1,1,1,m)),ff
     1g(5,nz1,kx1,1,m)*aimag(q(nz1,1,1,m)))
            wp = wp + .5*(ffg(1,nz1,1,1,m)*real(q(nz1,1,1,m))**2 + ffg(1
     1,nz1,kx1,1,m)*aimag(q(nz1,1,1,m))**2)
         endif
      endif
c mode numbers ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         at2 = -1.0
         do 90 j = 1, kxyp2
         at2 = -at2
         if ((j+joff).gt.0) then
            at8 = ffg(5,1,j,k1,m)
            do 80 l = 2, nz
            l1 = nz22 - l
            at1 = at2*ffg(3,l,kx1,k1,m)
            zt1 = cmplx(at1,ffg(3,l,j,k1,m))
            fx(l,j,k1,m) = zt1*q(l,j,k1,m)
            fx(l1,j,k1,m) = zt1*q(l1,j,k1,m)
            at3 = ffg(4,l,j,k1,m)
            fy(l,j,k1,m) = at3*q(l,j,k1,m)
            fy(l1,j,k1,m) = at3*q(l1,j,k1,m)
            at8 = -at8
            zt3 = cmplx(at8,ffg(5,l,j,k1,m))
            fz(l,j,k1,m) = zt3*q(l,j,k1,m)
            fz(l1,j,k1,m) = conjg(zt3)*q(l1,j,k1,m)
            wp = wp + ffg(1,l,j,k1,m)*(q(l,j,k1,m)*conjg(q(l,j,k1,m)) + 
     1q(l1,j,k1,m)*conjg(q(l1,j,k1,m)))
   80       continue
c mode numbers kz = 0
            at1 = at2*ffg(3,1,kx1,k1,m)
            zt1 = cmplx(at1,ffg(3,1,j,k1,m))
            fx(1,j,k1,m) = zt1*q(1,j,k1,m)
            at4 = ffg(4,1,j,k1,m)
            fy(1,j,k1,m) = at4*q(1,j,k1,m)
            at8 = ffg(5,1,j,k1,m)
            fz(1,j,k1,m) = at8*q(1,j,k1,m)
            wp = wp + ffg(1,1,j,k1,m)*q(1,j,k1,m)*conjg(q(1,j,k1,m))
c mode numbers kz = nz
            at1 = at2*ffg(3,nz1,kx1,k1,m)
            zt1 = cmplx(at1,ffg(3,nz1,j,k1,m))
            fx(nz1,j,k1,m) = zt1*q(nz1,j,k1,m)
            at4 = ffg(4,nz1,j,k1,m)
            fy(nz1,j,k1,m) = at4*q(nz1,j,k1,m)
            at8 = ffg(5,nz1,j,k1,m)
            fz(nz1,j,k1,m) = at8*q(nz1,j,k1,m)
            wp = wp + ffg(1,nz1,j,k1,m)*q(nz1,j,k1,m)*conjg(q(nz1,j,k1,m
     1))
         endif
   90    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            kx1 = kxyp2 + 1
            at5 = ffg(5,1,1,k1,m)
            at6 = ffg(5,1,kx1,k1,m)
            do 100 l = 2, nz
            l1 = nz22 - l
c mode numbers ky = ny/2, kx = 0
            at2 = ffg(3,l,1,k1,m)
            fx(l,1,k1,m) = at2*q(l,1,k1,m)
            at4 = ffg(4,l,1,k1,m)
            fy(l,1,k1,m) = at4*q(l,1,k1,m)
            at5 = -at5
            zt3 = cmplx(at5,ffg(5,l,1,k1,m))
            fz(l,1,k1,m) = zt3*q(l,1,k1,m)
            wp = wp + ffg(1,l,1,k1,m)*q(l,1,k1,m)*conjg(q(l,1,k1,m))
c mode numbers ky = ny/2, kx = nx
            at2 = ffg(3,l,kx1,k1,m)
            fx(l1,1,k1,m) = at2*q(l1,1,k1,m)
            at4 = ffg(4,l,kx1,k1,m)
            fy(l1,1,k1,m) = at4*q(l1,1,k1,m)
            at6 = -at6
            zt3 = cmplx(at6,ffg(5,l,kx1,k1,m))
            fz(l1,1,k1,m) = zt3*q(l1,1,k1,m)
            wp = wp + ffg(1,l,kx1,k1,m)*q(l1,1,k1,m)*conjg(q(l1,1,k1,m))
  100    continue
c mode numbers ky = ny/2, kx = 0, nx/2, kz = 0
            fx(1,1,k1,m) = cmplx(ffg(3,1,1,k1,m)*real(q(1,1,k1,m)),ffg(3
     1,1,kx1,k1,m)*aimag(q(1,1,k1,m)))
            fy(1,1,k1,m) = cmplx(ffg(4,1,1,k1,m)*real(q(1,1,k1,m)),ffg(4
     1,1,kx1,k1,m)*aimag(q(1,1,k1,m)))
            fz(1,1,k1,m) = cmplx(ffg(5,1,1,k1,m)*real(q(1,1,k1,m)),ffg(5
     1,1,kx1,k1,m)*aimag(q(1,1,k1,m)))
            wp = wp + .5*(ffg(1,1,1,k1,m)*real(q(1,1,k1,m))**2 + ffg(1,1
     1,kx1,k1,m)*aimag(q(1,1,k1,m))**2)
c mode numbers ky = ny/2, kx = 0, nx/2, kz = nz/2
            fx(nz1,1,k1,m) = cmplx(ffg(3,nz1,1,k1,m)*real(q(nz1,1,k1,m))
     1,ffg(3,nz1,kx1,k1,m)*aimag(q(nz1,1,k1,m)))
            fy(nz1,1,k1,m) = cmplx(ffg(4,nz1,1,k1,m)*real(q(nz1,1,k1,m))
     1,ffg(4,nz1,kx1,k1,m)*aimag(q(nz1,1,k1,m)))
            fz(nz1,1,k1,m) = cmplx(ffg(5,nz1,1,k1,m)*real(q(nz1,1,k1,m))
     1,ffg(5,nz1,kx1,k1,m)*aimag(q(nz1,1,k1,m)))
            wp = wp + .5*(ffg(1,nz1,1,k1,m)*real(q(nz1,1,k1,m))**2 + ffg
     1(1,nz1,kx1,k1,m)*aimag(q(nz1,1,k1,m))**2)
         endif
      endif
  110 continue
  120 continue
  130 continue
      we = 8.0*float(nx*ny*nz)*wp
      return
c calculate potential and sum field energy
  140 if (isign.gt.1) go to 280
      wp = 0.0d0
      if (kstrt.gt.(kxb2*kyzb2)) go to 270
      do 260 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 250 mx = 1, j2blok
      m = mx + moff
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      do 180 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         do 160 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 150 l = 2, nz
            l1 = nz22 - l
            at2 = ffg(1,l,j,k,m)
c           at1 = at2*ffg(2,l,j,k,m)
            at1 = at2
            fx(l,j,k,m) = at2*q(l,j,k,m)
            fx(l1,j,k,m) = at2*q(l1,j,k,m)
            wp = wp + at1*(q(l,j,k,m)*conjg(q(l,j,k,m)) + q(l1,j,k,m)*co
     1njg(q(l1,j,k,m)))
  150       continue
c mode numbers kz = 0, nz
            at2 = ffg(1,1,j,k,m)
c           at1 = at2*ffg(2,1,j,k,m)
            at1 = at2
            fx(1,j,k,m) = at2*q(1,j,k,m)
            wp = wp + at1*q(1,j,k,m)*conjg(q(1,j,k,m))
            at2 = ffg(1,nz1,j,k,m)
c           at1 = at2*ffg(2,nz1,j,k,m)
            at1 = at2
            fx(nz1,j,k,m) = at2*q(nz1,j,k,m)
            wp = wp + at1*q(nz1,j,k,m)*conjg(q(nz1,j,k,m))
         endif
  160    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 170 l = 2, nz
            l1 = nz22 - l
            at2 = ffg(1,l,1,k,m)
c           at1 = at2*ffg(2,l,1,k,m)
            at1 = at2
            fx(l,1,k,m) = at2*q(l,1,k,m)
            fx(l1,1,k,m) = at2*q(l1,1,k,m)
            wp = wp + at1*(q(l,1,k,m)*conjg(q(l,1,k,m)) + q(l1,1,k,m)*co
     1njg(q(l1,1,k,m)))
  170       continue
c mode numbers kz = 0
            at2 = ffg(1,1,1,k,m)
c           at1 = at2*ffg(2,1,1,k,m)
            at1 = at2
            fx(1,1,k,m) = at2*q(1,1,k,m)
            wp = wp + at1*q(1,1,k,m)*conjg(q(1,1,k,m))
c mode numbers kz = nz
            at2 = ffg(1,nz1,1,k,m)
c           at1 = at2*ffg(2,nz1,1,k,m)
            at1 = at2
            fx(nz1,1,k,m) = at2*q(nz1,1,k,m)
            wp = wp + at1*q(nz1,1,k,m)*conjg(q(nz1,1,k,m))
         endif
      endif
  180 continue
c mode numbers ky = 0
      n2 = koff + 1
      if (n2.eq.0) then
         do 200 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 190 l = 2, nz
            l1 = nz22 - l
            at2 = ffg(1,l,j,1,m)
c           at1 = at2*ffg(2,l,j,1,m)
            at1 = at2
            fx(l,j,1,m) = at2*q(l,j,1,m)
            fx(l1,j,1,m) = at2*q(l1,j,1,m)
            wp = wp + at1*(q(l,j,1,m)*conjg(q(l,j,1,m)) + q(l1,j,1,m)*co
     1njg(q(l1,j,1,m)))
  190       continue
c mode numbers kz = 0
            at2 = ffg(1,1,j,1,m)
c           at1 = at2*ffg(2,1,j,1,m)
            at1 = at2
            fx(1,j,1,m) = at2*q(1,j,1,m)
            wp = wp + at1*q(1,j,1,m)*conjg(q(1,j,1,m))
c mode numbers kz = nz
            at2 = ffg(1,nz1,j,1,m)
c           at1 = at2*ffg(2,nz1,j,1,m)
            at1 = at2
            fx(nz1,j,1,m) = at2*q(nz1,j,1,m)
            wp = wp + at1*q(nz1,j,1,m)*conjg(q(nz1,j,1,m))
         endif
  200    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            kx1 = kxyp2 + 1
            do 210 l = 2, nz
            l1 = nz22 - l
c mode numbers kx = 0
            at2 = ffg(1,l,1,1,m)
c           at1 = at2*ffg(2,l,1,1,m)
            at1 = at2
            fx(l,1,1,m) = at2*q(l,1,1,m)
            wp = wp + at1*q(l,1,1,m)*conjg(q(l,1,1,m))
c mode numbers kx = nx
            at2 = ffg(1,l,kx1,1,m)
c           at1 = at2*ffg(2,l,kx1,1,m)
            at1 = at2
            fx(l1,1,1,m) = at2*q(l1,1,1,m)
            wp = wp + at1*q(l1,1,1,m)*conjg(q(l1,1,1,m))
  210       continue
c mode numbers kx = 0, nx/2, kz = 0
            at2 = ffg(1,1,1,1,m)
c           at1 = at2*ffg(2,1,1,1,m)
            at1 = at2
            at4 = ffg(1,1,kx1,1,m)
c           at3 = at4*ffg(2,1,kx1,1,m)
            at3 = at4
            fx(1,1,1,m) = cmplx(at2*real(q(1,1,1,m)),at4*aimag(q(1,1,1,m
     1)))
            wp = wp + .5*(at1*real(q(1,1,1,m))**2 + at3*aimag(q(1,1,1,m)
     1)**2)
c mode numbers kx = 0, nx/2, kz = nz/2
            at2 = ffg(1,nz1,1,1,m)
c           at1 = at2*ffg(2,nz1,1,1,m)
            at1 = at2
            at4 = ffg(1,nz1,kx1,1,m)
c           at3 = at4*ffg(2,nz1,kx1,1,m)
            at3 = at4
            fx(nz1,1,1,m) = cmplx(at2*real(q(nz1,1,1,m)),at4*aimag(q(nz1
     1,1,1,m)))
            wp = wp + .5*(at1*real(q(nz1,1,1,m))**2 + at3*aimag(q(nz1,1,
     11,m))**2)
         endif
      endif
c mode numbers ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 230 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 220 l = 2, nz
            l1 = nz22 - l
            at2 = ffg(1,l,j,k1,m)
c           at1 = at2*ffg(2,l,j,k1,m)
            at1 = at2
            fx(l,j,k1,m) = at2*q(l,j,k1,m)
            fx(l1,j,k1,m) = at2*q(l1,j,k1,m)
            wp = wp + at1*(q(l,j,k1,m)*conjg(q(l,j,k1,m)) + q(l1,j,k1,m)
     1*conjg(q(l1,j,k1,m)))
  220      continue
c mode numbers kz = 0
            at2 = ffg(1,1,j,k1,m)
c           at1 = at2*ffg(2,1,j,k1,m)
            at1 = at2
            fx(1,j,k1,m) = at2*q(1,j,k1,m)
            wp = wp + at1*q(1,j,k1,m)*conjg(q(1,j,k1,m))
c mode numbers kz = nz
            at2 = ffg(1,nz1,j,k1,m)
c           at1 = at2*ffg(2,nz1,j,k1,m)
            at1 = at2
            fx(nz1,j,k1,m) = at2*q(nz1,j,k1,m)
            wp = wp + at1*q(nz1,j,k1,m)*conjg(q(nz1,j,k1,m))
         endif
  230    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            kx1 = kxyp2 + 1
            do 240 l = 2, nz
            l1 = nz22 - l
c mode numbers ky = ny/2, kx = 0
            at2 = ffg(1,l,1,k1,m)
c           at1 = at2*ffg(2,l,1,k1,m)
            at1 = at2
            fx(l,1,k1,m) = at2*q(l,1,k1,m)
            wp = wp + at1*q(l,1,k1,m)*conjg(q(l,1,k1,m))
c mode numbers ky = ny/2, kx = nx
            at2 = ffg(1,l,kx1,k1,m)
c           at1 = at2*ffg(2,l,kx1,k1,m)
            at1 = at2
            fx(l1,1,k1,m) = at2*q(l1,1,k1,m)
            wp = wp + at1*q(l1,1,k1,m)*conjg(q(l1,1,k1,m))
  240    continue
c mode numbers ky = ny/2, kx = 0, nx/2, kz = 0
            at2 = ffg(1,1,1,k1,m)
c           at1 = at2*ffg(2,1,1,k1,m)
            at1 = at2
            at4 = ffg(1,1,kx1,k1,m)
c           at3 = at4*ffg(2,1,kx1,k1,m)
            at3 = at4
            fx(1,1,k1,m) = cmplx(at2*real(q(1,1,k1,m)),at4*aimag(q(1,1,k
     11,m)))
            wp = wp + .5*(at1*real(q(1,1,k1,m))**2 + at3*aimag(q(1,1,k1,
     1m))**2)
c mode numbers ky = ny/2, kx = 0, nx/2, kz = nz/2
            at2 = ffg(1,nz1,1,k1,m)
c           at1 = at2*ffg(2,nz1,1,k1,m)
            at1 = at2
            at4 = ffg(1,nz1,kx1,k1,m)
c           at3 = at4*ffg(2,nz1,kx1,k1,m)
            at3 = at4
            fx(nz1,1,k1,m) = cmplx(at2*real(q(nz1,1,k1,m)),at4*aimag(q(n
     1z1,1,k1,m)))
            wp = wp + .5*(at1*real(q(nz1,1,k1,m))**2 + at3*aimag(q(nz1,1
     1,k1,m))**2)
         endif
      endif
  250 continue
  260 continue
  270 continue
      we = 8.0*float(nx*ny*nz)*wp
      return
c calculate smoothing
  280 if (kstrt.gt.(kxb2*kyzb2)) go to 410
      do 400 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      do 390 mx = 1, j2blok
      m = mx + moff
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      do 320 k = 1, kyzp2
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.ny)) then
         do 300 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 290 l = 2, nz
            l1 = nz22 - l
            at1 = ffg(2,l,j,k,m)
            fy(l,j,k,m) = at1*q(l,j,k,m)
            fy(l1,j,k,m) = at1*q(l1,j,k,m)
  290       continue
c mode numbers kz = 0, nz
            fy(1,j,k,m) = ffg(2,1,j,k,m)*q(1,j,k,m)
            fy(nz1,j,k,m) = ffg(2,nz1,j,k,m)*q(nz1,j,k,m)
         endif
  300    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            do 310 l = 2, nz
            l1 = nz22 - l
            at1 = ffg(2,l,1,k,m)
            fy(l,1,k,m) = at1*q(l,1,k,m)
            fy(l1,1,k,m) = at1*q(l1,1,k,m)
  310       continue
c mode numbers kz = 0
            fy(1,1,k,m) = ffg(2,1,1,k,m)*q(1,1,k,m)
c mode numbers kz = nz
            fy(nz1,1,k,m) = ffg(2,nz1,1,k,m)*q(nz1,1,k,m)
         endif
      endif
  320 continue
c mode numbers ky = 0
      n2 = koff + 1
      if (n2.eq.0) then
         do 340 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 330 l = 2, nz
            l1 = nz22 - l
            at1 = ffg(2,l,j,1,m)
            fy(l,j,1,m) = at1*q(l,j,1,m)
            fy(l1,j,1,m) = at1*q(l1,j,1,m)
  330       continue
c mode numbers kz = 0
            fy(1,j,1,m) = ffg(2,1,j,1,m)*q(1,j,1,m)
c mode numbers kz = nz
            fy(nz1,j,1,m) = ffg(2,nz1,j,1,m)*q(nz1,j,1,m)
         endif
  340    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            kx1 = kxyp2 + 1
            do 350 l = 2, nz
            l1 = nz22 - l
c mode numbers kx = 0
            fy(l,1,1,m) = ffg(2,l,1,1,m)*q(l,1,1,m)
c mode numbers kx = nx
            fy(l1,1,1,m) = ffg(2,l,kx1,1,m)*q(l1,1,1,m)
  350       continue
c mode numbers kx = 0, nx/2, kz = 0
            at1 = ffg(2,1,1,1,m)
            at3 = ffg(2,1,kx1,1,m)
            fy(1,1,1,m) = cmplx(at1*real(q(1,1,1,m)),at3*aimag(q(1,1,1,m
     1)))
c mode numbers kx = 0, nx/2, kz = nz/2
            at1 = ffg(2,nz1,1,1,m)
            at3 = ffg(2,nz1,kx1,1,m)
            fy(nz1,1,1,m) = cmplx(at1*real(q(nz1,1,1,m)),at3*aimag(q(nz1
     1,1,1,m)))
         endif
      endif
c mode numbers ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         do 370 j = 1, kxyp2
         if ((j+joff).gt.0) then
            do 360 l = 2, nz
            l1 = nz22 - l
            at1 = ffg(2,l,j,k1,m)
            fy(l,j,k1,m) = at1*q(l,j,k1,m)
            fy(l1,j,k1,m) = at1*q(l1,j,k1,m)
  360      continue
c mode numbers kz = 0
            fy(1,j,k1,m) = ffg(2,1,j,k1,m)*q(1,j,k1,m)
c mode numbers kz = nz
            fy(nz1,j,k1,m) = ffg(2,nz1,j,k1,m)*q(nz1,j,k1,m)
         endif
  370    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            kx1 = kxyp2 + 1
            do 380 l = 2, nz
            l1 = nz22 - l
c mode numbers ky = ny/2, kx = 0
            fy(l,1,k1,m) = ffg(2,l,1,k1,m)*q(l,1,k1,m)
c mode numbers ky = ny/2, kx = nx
            fy(l1,1,k1,m) = ffg(2,l,kx1,k1,m)*q(l1,1,k1,m)
  380    continue
c mode numbers ky = ny/2, kx = 0, nx/2, kz = 0
            at1 = ffg(2,1,1,k1,m)
            at3 = ffg(2,1,kx1,k1,m)
            fy(1,1,k1,m) = cmplx(at1*real(q(1,1,k1,m)),at3*aimag(q(1,1,k
     11,m)))
c mode numbers ky = ny/2, kx = 0, nx/2, kz = nz/2
            at1 = ffg(2,nz1,1,k1,m)
            at3 = ffg(2,nz1,kx1,k1,m)
            fy(nz1,1,k1,m) = cmplx(at1*real(q(nz1,1,k1,m)),at3*aimag(q(n
     1z1,1,k1,m)))
         endif
      endif
  390 continue
  400 continue
  410 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOISC332(q,fxyz,ffg,we,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2
     1,j2blok,m2blok,kxyp2d,kyzp2d,nz1d)
c this subroutine solves 3d poisson's equation in fourier space for
c force/charge (or convolution of electric field over particle shape)
c with open (vacuum) boundary conditions using hockney's method,
c for distributed data with 2D spatial decomposition
c fourier coefficients are constructed so that a real to complex fft
c will perform the appropriate convolution
c input: q,ffg,nx,ny,nz,kstrt,nz2d,kxyp2,kyzp2,j2blok,m2blok,kxyp2d,nz1d
c output: fxyz,we
c        57*nxc*nyc*nzc + 50*(nxc*nyc + nxc*nzc + nyc*nzc)
c where nxc = (nx-1)/nvpy, nyc = (ny-1)/nvpz, nzc = nz - 1, and
c nvpy/nvpz = number of procs in y/z
c force/charge is calculated using the equations:
c fx(kx,ky,kz) = gx(kx,ky,kz)*s(kx,ky,kz)*q(kx,ky,kz),
c fy(kx,ky,kz) = gy(kx,ky,kz)*s(kx,ky,kz)*q(kx,ky,kz),
c fz(kx,ky,kz) = gz(kx,ky,kz)*s(kx,ky,kz)*q(kx,ky,kz),
c where kx = 2pi*j/nx, ky = 2pi*k/ny, kz = 2pi*l/nz,
c and j,k,l = fourier mode numbers,
c gx(kx,ky,kz) = s(kx,ky,kz)*inverse FFT of (x/r)*Er
c gy(kx,ky,kz) = s(kx,ky,kz)*inverse FFT of (y/r)*Er
c gz(kx,ky,kz) = s(kx,ky,kz)*inverse FFT of (z/r)*Er
c where Er is the electric field of a single finite-sized particle
c s(kx,ky,kz) = inverse FFT of the density of a finite-sized particle
c q(l,j,k,m) = complex charge density for fourier mode jj-1,kk-1,l-1
c fxyz(1,l,j,k,m) = x component of force/charge
c fxyz(2,,j,k,m) = y component of force/charge
c fxyz(3,,j,k,m) = z component of force/charge
c ffg(1,l,j,k) = potential green's function g
c ffg(2,l,j,k) = finite-size particle shape factor s
c ffg(3,l,j,k) = x component of electric field green's function gx
c ffg(4,l,j,k) = y component of electric field green's function gy
c ffg(5,l,j,k) = z component of electric field green's function gz
c all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp2*(mx - 1) and
c kk = k + kyzp2*(my - 1), 0 <= mx <= nvpy, 0 <= my < nvpz
c the ffg array is calculated by the subroutine PFORMC32
c electric field energy is also calculated and returned in we
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c nz2d = first dimension of field arrays, must be >= 2*nz
c kxyp2/kyzp2 = number of complex grids in each field partition in
c x/y direction
c j2blok/m2blok = number of field partitions in x/y
c kxyp2d = third dimension of ffg arrays, must be >= kxyp2+1
c kyzp2d = third dimension of ffg arrays, must be >= kyzp2+1
c nz1d = second dimension of ffg arrays, must be >= nz+1
      implicit none
      real ffg
      complex q, fxyz
      integer nx, ny, nz, kstrt, nz2d, kxyp2, kyzp2, j2blok, m2blok
      integer kxyp2d, kyzp2d, nz1d
      real we
      dimension q(nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension fxyz(3,nz2d,kxyp2,kyzp2,j2blok*m2blok)
      dimension ffg(5,nz1d,kxyp2d,kyzp2d,j2blok*m2blok)
c local data
      double precision wp
      integer j, k, l, mx, my, m, nz22, kxb2, kyzb2, js, ks, kx1, ky1
      integer kxy1, joff, koff, moff, k1, l1, n1, n2, nz1
      real at1, at2, at3, at4, at5, at6, at7, at8
      complex zt1, zt2, zt3
      nz1 = nz + 1
      nz22 = nz + nz + 2
      kxb2 = nx/kxyp2
      kyzb2 = (ny + ny)/kyzp2
      ks = (kstrt - 1)/kxb2
      js = kstrt - kxb2*ks - 2
      ks = ks - 1
      kx1 = 1
      ky1 = 1
      kxy1 = 1
c calculate force/charge and sum field energy
      wp = 0.0d0
      if (kstrt.gt.(kxb2*kyzb2)) go to 130
      do 120 my = 1, m2blok
      koff = kyzp2*(my + ks) - 1
      moff = j2blok*(my - 1)
      if ((my+ks).gt.0) ky1 = kyzp2 + 1
      do 110 mx = 1, j2blok
      m = mx + moff
      if ((mx+js).gt.0) kx1 = kxyp2 + 1
c mode numbers 0 < kx < nx, 0 < ky < ny, and 0 < kz < nz
      joff = kxyp2*(mx + js) - 1
      at5 = -1.0
      do 40 k = 1, kyzp2
      k1 = k + koff
      if (k1.ge.ny) kxy1 = kxyp2 + 1
      at5 = -at5
      if ((k1.gt.0).and.(k1.ne.ny)) then
         at2 = -1.0
         do 20 j = 1, kxyp2
         at2 = -at2
         if ((j+joff).gt.0) then
            at7 = ffg(5,1,j,k,m)
            do 10 l = 2, nz
            l1 = nz22 - l
            at1 = at2*ffg(3,l,kx1,k,m)
            zt1 = cmplx(at1,ffg(3,l,j,k,m))
            fxyz(1,l,j,k,m) = zt1*q(l,j,k,m)
            fxyz(1,l1,j,k,m) = zt1*q(l1,j,k,m)
            at3 = at5*ffg(4,l,j,ky1,m)
            zt2 = cmplx(at3,ffg(4,l,j,k,m))
            fxyz(2,l,j,k,m) = zt2*q(l,j,k,m)
            fxyz(2,l1,j,k,m) = zt2*q(l1,j,k,m)
            at7 = -at7
            zt3 = cmplx(at7,ffg(5,l,j,k,m))
            fxyz(3,l,j,k,m) = zt3*q(l,j,k,m)
            fxyz(3,l1,j,k,m) = conjg(zt3)*q(l1,j,k,m)
            wp = wp + ffg(1,l,j,k,m)*(q(l,j,k,m)*conjg(q(l,j,k,m)) + q(l
     11,j,k,m)*conjg(q(l1,j,k,m)))
   10       continue
c mode numbers kz = 0, nz
            at1 = at2*ffg(3,1,kx1,k,m)
            zt1 = cmplx(at1,ffg(3,1,j,k,m))
            fxyz(1,1,j,k,m) = zt1*q(1,j,k,m)
            at3 = at5*ffg(4,1,j,ky1,m)
            zt2 = cmplx(at3,ffg(4,1,j,k,m))
            fxyz(2,1,j,k,m) = zt2*q(1,j,k,m)
            at8 = ffg(5,1,j,k,m)
            fxyz(3,1,j,k,m) = at8*q(1,j,k,m)
            wp = wp + ffg(1,1,j,k,m)*q(1,j,k,m)*conjg(q(1,j,k,m))
            at1 = at2*ffg(3,nz1,kx1,k,m)
            zt1 = cmplx(at1,ffg(3,nz1,j,k,m))
            fxyz(1,nz1,j,k,m) = zt1*q(nz1,j,k,m)
            at3 = at5*ffg(4,nz1,j,ky1,m)
            zt2 = cmplx(at3,ffg(4,nz1,j,k,m))
            fxyz(2,nz1,j,k,m) = zt2*q(nz1,j,k,m)
            at8 = ffg(5,nz1,j,k,m)
            fxyz(3,nz1,j,k,m) = at8*q(nz1,j,k,m)
            wp = wp + ffg(1,nz1,j,k,m)*q(nz1,j,k,m)*conjg(q(nz1,j,k,m))
         endif
   20    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            at7 = ffg(5,1,1,k,m)
            do 30 l = 2, nz
            l1 = nz22 - l
            at1 = ffg(3,l,1,k,m)
            fxyz(1,l,1,k,m) = at1*q(l,1,k,m)
            fxyz(1,l1,1,k,m) = at1*q(l1,1,k,m)
            at3 = at5*ffg(4,l,kxy1,ky1,m)
            zt2 = cmplx(at3,ffg(4,l,1,k,m))
            fxyz(2,l,1,k,m) = zt2*q(l,1,k,m)
            fxyz(2,l1,1,k,m) = zt2*q(l1,1,k,m)
            at7 = -at7
            zt3 = cmplx(at7,ffg(5,l,1,k,m))
            fxyz(3,l,1,k,m) = zt3*q(l,1,k,m)
            fxyz(3,l1,1,k,m) = conjg(zt3)*q(l1,1,k,m)
            wp = wp + ffg(1,l,1,k,m)*(q(l,1,k,m)*conjg(q(l,1,k,m)) + q(l
     11,1,k,m)*conjg(q(l1,1,k,m)))
   30       continue
c mode numbers kz = 0
            at2 = ffg(3,1,1,k,m)
            fxyz(1,1,1,k,m) = at2*q(1,1,k,m)
            at3 = at5*ffg(4,1,kxy1,ky1,m)
            zt2 = cmplx(at3,ffg(4,1,1,k,m))
            fxyz(2,1,1,k,m) = zt2*q(1,1,k,m)
            at8 = ffg(5,1,1,k,m)
            fxyz(3,1,1,k,m) = at8*q(1,1,k,m)
            wp = wp + ffg(1,1,1,k,m)*q(1,1,k,m)*conjg(q(1,1,k,m))
c mode numbers kz = nz
            at2 = ffg(3,nz1,1,k,m)
            fxyz(1,nz1,1,k,m) = at2*q(nz1,1,k,m)
            at3 = at5*ffg(4,nz1,kxy1,ky1,m)
            zt2 = cmplx(at3,ffg(4,nz1,1,k,m))
            fxyz(2,nz1,1,k,m) = zt2*q(nz1,1,k,m)
            at8 = ffg(5,nz1,1,k,m)
            fxyz(3,nz1,1,k,m) = at8*q(nz1,1,k,m)
            wp = wp + ffg(1,nz1,1,k,m)*q(nz1,1,k,m)*conjg(q(nz1,1,k,m))
         endif
      endif
   40 continue
c mode numbers ky = 0
      n2 = koff + 1
      if (n2.eq.0) then
         at2 = -1.0
         do 60 j = 1, kxyp2
         at2 = -at2
         if ((j+joff).gt.0) then
            at7 = ffg(5,1,j,1,m)
            do 50 l = 2, nz
            l1 = nz22 - l
            at1 = at2*ffg(3,l,kx1,1,m)
            zt1 = cmplx(at1,ffg(3,l,j,1,m))
            fxyz(1,l,j,1,m) = zt1*q(l,j,1,m)
            fxyz(1,l1,j,1,m) = zt1*q(l1,j,1,m)
            at3 = ffg(4,l,j,1,m)
            fxyz(2,l,j,1,m) = at3*q(l,j,1,m)
            fxyz(2,l1,j,1,m) = at3*q(l1,j,1,m)
            at7 = -at7
            zt3 = cmplx(at7,ffg(5,l,j,1,m))
            fxyz(3,l,j,1,m) = zt3*q(l,j,1,m)
            fxyz(3,l1,j,1,m) = conjg(zt3)*q(l1,j,1,m)
            wp = wp + ffg(1,l,j,1,m)*(q(l,j,1,m)*conjg(q(l,j,1,m)) + q(l
     11,j,1,m)*conjg(q(l1,j,1,m)))
   50       continue
c mode numbers kz = 0
            at1 = at2*ffg(3,1,kx1,1,m)
            zt1 = cmplx(at1,ffg(3,1,j,1,m))
            fxyz(1,1,j,1,m) = zt1*q(1,j,1,m)
            at4 = ffg(4,1,j,1,m)
            fxyz(2,1,j,1,m) = at4*q(1,j,1,m)
            at8 = ffg(5,1,j,1,m)
            fxyz(3,1,j,1,m) = at8*q(1,j,1,m)
            wp = wp + ffg(1,1,j,1,m)*q(1,j,1,m)*conjg(q(1,j,1,m))
c mode numbers kz = nz
            at1 = at2*ffg(3,nz1,kx1,1,m)
            zt1 = cmplx(at1,ffg(3,nz1,j,1,m))
            fxyz(1,nz1,j,1,m) = zt1*q(nz1,j,1,m)
            at4 = ffg(4,nz1,j,1,m)
            fxyz(2,nz1,j,1,m) = at4*q(nz1,j,1,m)
            at8 = ffg(5,nz1,j,1,m)
            fxyz(3,nz1,j,1,m) = at8*q(nz1,j,1,m)
            wp = wp + ffg(1,nz1,j,1,m)*q(nz1,j,1,m)*conjg(q(nz1,j,1,m))
         endif
   60    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            kx1 = kxyp2 + 1
            at7 = ffg(5,1,1,1,m)
            at8 = ffg(5,1,kx1,1,m)
            do 70 l = 2, nz
            l1 = nz22 - l
c mode numbers kx = 0
            at2 = ffg(3,l,1,1,m)
            fxyz(1,l,1,1,m) = at2*q(l,1,1,m)
            at4 = ffg(4,l,1,1,m)
            fxyz(2,l,1,1,m) = at4*q(l,1,1,m)
            at7 = -at7
            zt3 = cmplx(at7,ffg(5,l,1,1,m))
            fxyz(3,l,1,1,m) = zt3*q(l,1,1,m)
            wp = wp + ffg(1,l,1,1,m)*q(l,1,1,m)*conjg(q(l,1,1,m))
c mode numbers kx = nx
            at2 = ffg(3,l,kx1,1,m)
            fxyz(1,l1,1,1,m) = at2*q(l1,1,1,m)
            at4 = ffg(4,l,kx1,1,m)
            fxyz(2,l1,1,1,m) = at4*q(l1,1,1,m)
            at8 = -at8
            zt3 = cmplx(at8,ffg(5,l,kx1,1,m))
            fxyz(3,l1,1,1,m) = zt3*q(l1,1,1,m)
            wp = wp + ffg(1,l,kx1,1,m)*q(l1,1,1,m)*conjg(q(l1,1,1,m))
   70       continue
c mode numbers kx = 0, nx/2, kz = 0
            fxyz(1,1,1,1,m) = cmplx(ffg(3,1,1,1,m)*real(q(1,1,1,m)),ffg(
     13,1,kx1,1,m)*aimag(q(1,1,1,m)))
            fxyz(2,1,1,1,m) = cmplx(ffg(4,1,1,1,m)*real(q(1,1,1,m)),ffg(
     14,1,kx1,1,m)*aimag(q(1,1,1,m)))
            fxyz(3,1,1,1,m) = cmplx(ffg(5,1,1,1,m)*real(q(1,1,1,m)),ffg(
     15,1,kx1,1,m)*aimag(q(1,1,1,m)))
            wp = wp + .5*(ffg(1,1,1,1,m)*real(q(1,1,1,m))**2 + ffg(1,1,k
     1x1,1,m)*aimag(q(1,1,1,m))**2)
c mode numbers kx = 0, nx/2, kz = nz/2
            fxyz(1,nz1,1,1,m) = cmplx(ffg(3,nz1,1,1,m)*real(q(nz1,1,1,m)
     1),ffg(3,nz1,kx1,1,m)*aimag(q(nz1,1,1,m)))
            fxyz(2,nz1,1,1,m) = cmplx(ffg(4,nz1,1,1,m)*real(q(nz1,1,1,m)
     1),ffg(4,nz1,kx1,1,m)*aimag(q(nz1,1,1,m)))
            fxyz(3,nz1,1,1,m) = cmplx(ffg(5,nz1,1,1,m)*real(q(nz1,1,1,m)
     1),ffg(5,nz1,kx1,1,m)*aimag(q(nz1,1,1,m)))
            wp = wp + .5*(ffg(1,nz1,1,1,m)*real(q(nz1,1,1,m))**2 + ffg(1
     1,nz1,kx1,1,m)*aimag(q(nz1,1,1,m))**2)
         endif
      endif
c mode numbers ky = ny
      k1 = (ny/kyzp2)*kyzp2
      if (n2.eq.k1) then
         k1 = ny - k1 + 1
         at2 = -1.0
         do 90 j = 1, kxyp2
         at2 = -at2
         if ((j+joff).gt.0) then
            at8 = ffg(5,1,j,k1,m)
            do 80 l = 2, nz
            l1 = nz22 - l
            at1 = at2*ffg(3,l,kx1,k1,m)
            zt1 = cmplx(at1,ffg(3,l,j,k1,m))
            fxyz(1,l,j,k1,m) = zt1*q(l,j,k1,m)
            fxyz(1,l1,j,k1,m) = zt1*q(l1,j,k1,m)
            at3 = ffg(4,l,j,k1,m)
            fxyz(2,l,j,k1,m) = at3*q(l,j,k1,m)
            fxyz(2,l1,j,k1,m) = at3*q(l1,j,k1,m)
            at8 = -at8
            zt3 = cmplx(at8,ffg(5,l,j,k1,m))
            fxyz(3,l,j,k1,m) = zt3*q(l,j,k1,m)
            fxyz(3,l1,j,k1,m) = conjg(zt3)*q(l1,j,k1,m)
            wp = wp + ffg(1,l,j,k1,m)*(q(l,j,k1,m)*conjg(q(l,j,k1,m)) + 
     1q(l1,j,k1,m)*conjg(q(l1,j,k1,m)))
   80       continue
c mode numbers kz = 0
            at1 = at2*ffg(3,1,kx1,k1,m)
            zt1 = cmplx(at1,ffg(3,1,j,k1,m))
            fxyz(1,1,j,k1,m) = zt1*q(1,j,k1,m)
            at4 = ffg(4,1,j,k1,m)
            fxyz(2,1,j,k1,m) = at4*q(1,j,k1,m)
            at8 = ffg(5,1,j,k1,m)
            fxyz(3,1,j,k1,m) = at8*q(1,j,k1,m)
            wp = wp + ffg(1,1,j,k1,m)*q(1,j,k1,m)*conjg(q(1,j,k1,m))
c mode numbers kz = nz
            at1 = at2*ffg(3,nz1,kx1,k1,m)
            zt1 = cmplx(at1,ffg(3,nz1,j,k1,m))
            fxyz(1,nz1,j,k1,m) = zt1*q(nz1,j,k1,m)
            at4 = ffg(4,nz1,j,k1,m)
            fxyz(2,nz1,j,k1,m) = at4*q(nz1,j,k1,m)
            at8 = ffg(5,nz1,j,k1,m)
            fxyz(3,nz1,j,k1,m) = at8*q(nz1,j,k1,m)
            wp = wp + ffg(1,nz1,j,k1,m)*q(nz1,j,k1,m)*conjg(q(nz1,j,k1,m
     1))
         endif
   90    continue
c mode numbers kx = 0, nx
         n1 = joff + 1
         if (n1.eq.0) then
            kx1 = kxyp2 + 1
            at5 = ffg(5,1,1,k1,m)
            at6 = ffg(5,1,kx1,k1,m)
            do 100 l = 2, nz
            l1 = nz22 - l
c mode numbers ky = ny/2, kx = 0
            at2 = ffg(3,l,1,k1,m)
            fxyz(1,l,1,k1,m) = at2*q(l,1,k1,m)
            at4 = ffg(4,l,1,k1,m)
            fxyz(2,l,1,k1,m) = at4*q(l,1,k1,m)
            at5 = -at5
            zt3 = cmplx(at5,ffg(5,l,1,k1,m))
            fxyz(3,l,1,k1,m) = zt3*q(l,1,k1,m)
            wp = wp + ffg(1,l,1,k1,m)*q(l,1,k1,m)*conjg(q(l,1,k1,m))
c mode numbers ky = ny/2, kx = nx
            at2 = ffg(3,l,kx1,k1,m)
            fxyz(1,l1,1,k1,m) = at2*q(l1,1,k1,m)
            at4 = ffg(4,l,kx1,k1,m)
            fxyz(2,l1,1,k1,m) = at4*q(l1,1,k1,m)
            at6 = -at6
            zt3 = cmplx(at6,ffg(5,l,kx1,k1,m))
            fxyz(3,l1,1,k1,m) = zt3*q(l1,1,k1,m)
            wp = wp + ffg(1,l,kx1,k1,m)*q(l1,1,k1,m)*conjg(q(l1,1,k1,m))
  100    continue
c mode numbers ky = ny/2, kx = 0, nx/2, kz = 0
            fxyz(1,1,1,k1,m) = cmplx(ffg(3,1,1,k1,m)*real(q(1,1,k1,m)),f
     1fg(3,1,kx1,k1,m)*aimag(q(1,1,k1,m)))
            fxyz(2,1,1,k1,m) = cmplx(ffg(4,1,1,k1,m)*real(q(1,1,k1,m)),f
     1fg(4,1,kx1,k1,m)*aimag(q(1,1,k1,m)))
            fxyz(3,1,1,k1,m) = cmplx(ffg(5,1,1,k1,m)*real(q(1,1,k1,m)),f
     1fg(5,1,kx1,k1,m)*aimag(q(1,1,k1,m)))
            wp = wp + .5*(ffg(1,1,1,k1,m)*real(q(1,1,k1,m))**2 + ffg(1,1
     1,kx1,k1,m)*aimag(q(1,1,k1,m))**2)
c mode numbers ky = ny/2, kx = 0, nx/2, kz = nz/2
            fxyz(1,nz1,1,k1,m) = cmplx(ffg(3,nz1,1,k1,m)*real(q(nz1,1,k1
     1,m)),ffg(3,nz1,kx1,k1,m)*aimag(q(nz1,1,k1,m)))
            fxyz(2,nz1,1,k1,m) = cmplx(ffg(4,nz1,1,k1,m)*real(q(nz1,1,k1
     1,m)),ffg(4,nz1,kx1,k1,m)*aimag(q(nz1,1,k1,m)))
            fxyz(3,nz1,1,k1,m) = cmplx(ffg(5,nz1,1,k1,m)*real(q(nz1,1,k1
     1,m)),ffg(5,nz1,kx1,k1,m)*aimag(q(nz1,1,k1,m)))
            wp = wp + .5*(ffg(1,nz1,1,k1,m)*real(q(nz1,1,k1,m))**2 + ffg
     1(1,nz1,kx1,k1,m)*aimag(q(nz1,1,k1,m))**2)
         endif
      endif
  110 continue
  120 continue
  130 continue
      we = 8.0*float(nx*ny*nz)*wp
      return
      end
c-----------------------------------------------------------------------
      function POTC3(r,affp,ari,ifun)
c this function calculates the fields for finite-size gaussian particles
c in 3D:
c if ifun = 1, calculate potential function
c POTC3 = (affp/(4*pi))*erfn(r/(ar*sqrt(2.)))/r, for r > 0.
c POTC3 = (affp/(4*pi))*sqrt(2./3.14159265358979)/ar, for r = 0.
c if ifun = 2, calculate particle shape function
c POTC3 = exp(-(r/(sqrt(2.)*ar))**2)/(sqrt(2.*pi)*ar)**3, for r > 0.
c POTC3 = 1./(sqrt(2.*pi)*ar)**3, for r = 0.
c if ifun = 3, calculate radial electric field
c POTC3 = (affp/(4*pi))*(1/r)*(erf(r/(sqrt(2.)*ar))/r -
c exp(-(r/(sqrt(2.)*ar))**2)*sqrt(2./3.14159265358979)/ar, for r > 0.
c POTC3 = 0.0, for r = 0.
c where erfn is the error function
c and where the finite-size particle density is given by:
c rho(r) = exp(-(r/sqrt(2)*ar)**2)/(sqrt(2*pi)*ar)**3
c affp = 4*pi*e**2/(me*(omega0**2)*delta**3) = 1/(n0*delta**3)
c where n0*delta**3 = number density per grid
c r = radial coordinate
c affp = normalization constant
c ari = 1/ar = inverse of particle size function
c (ari = 0., means use point particle result)
c ifun = (1,2,3) = calculate (potential,shape,electric field)
      implicit none
      real r, affp, ari
      integer ifun
c local data
c pi4i = 1/4*pi, sqt2i = 1./sqrt(2.), sqt2pi = sqrt(2./pi)
      real pi4i, sqt2i, sqt2pi
      parameter(pi4i=0.5/6.28318530717959)
      parameter(sqt2i=0.707106781186548,sqt2pi=0.797884560802865)
      real POTC3, erfn
      external erfn
      real anorm, at1, ri
      anorm = affp*pi4i
c calculate potential function
      if (ifun.eq.1) then
c finite-size particles
         if (ari.gt.0.) then
            if (r.eq.0.) then
               POTC3 = anorm*sqt2pi*ari
            else
               POTC3 = anorm*erfn(r*sqt2i*ari)/r
            endif
c point particles
         else
            if (r.eq.0.) then
               POTC3 = 0.0
            else
               POTC3 = anorm/r
            endif
         endif
c calculate particle shape function
      else if (ifun.eq.2) then
         anorm = affp*(.5*sqt2pi*ari)**3
c finite-size particles
         if (ari.gt.0.) then
            if (r.eq.0.) then
               POTC3 = anorm
            else
               at1 = amin1(r*sqt2i*ari,8.0)
               POTC3 = anorm*exp(-(at1*at1))
            endif
c point particles
         else
            if (r.eq.0.) then
               POTC3 = affp
            else
               POTC3 = 0.0
            endif
         endif
c calculate radial electric field
      else if (ifun.eq.3) then
c finite-size particles
         if (ari.gt.0.) then
            if (r.eq.0.) then
               POTC3 = 0.0
            else
               ri = 1.0/r
               at1 = amin1(r*sqt2i*ari,8.0)
               POTC3 = anorm*ri*(erfn(at1)*ri - sqt2pi*ari*exp(-(at1*at1
     1)))
            endif
c point particles
         else
            if (r.eq.0.) then
               POTC3 = 0.0
            else
               POTC3 = anorm/(r*r)
            endif
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine PPOTC32(f,nyzp,noff,fpotc,ifun,affp,ar,kstrt,ix0,iy0,iz
     10,nx,nxv,nypmx,nzpmx,idds,mnblok,epsmax)
c test program for comparing fields due to finite-size particles
c calculated numerically (f) and analytically (g)
      implicit none
      real f
      real affp, ar, epsmax
      integer ifun, kstrt, ix0, iy0, iz0, nx, nxv, nypmx, nzpmx, mnblok
      integer idds
      integer nyzp, noff
      dimension f(nxv,nypmx,nzpmx,mnblok)
      dimension nyzp(idds,mnblok), noff(idds,mnblok)
      real fpotc
      external fpotc
c local data
      integer mnoff, lnoff, j, k, l, m, ifn
      real ari, at1, at2, x, y, z, r, g, eps
      ari = 0.0
      epsmax = 0.0
      if (ar.gt.0.) ari = 1.0/ar
      ifn = ifun
      if (ifun.ge.4) ifn = 3
c calculate analytic function
      do 40 m = 1, mnblok
      mnoff = noff(1,m) - 1 - iy0
      lnoff = noff(2,m) - 1 - iz0
      do 30 l = 1, nyzp(2,m) + 1
      z = float(l + lnoff)
      at1 = z*z
      do 20 k = 1, nyzp(1,m) + 1
      y = float(k + mnoff)
      at2 = y*y + at1
      do 10 j = 1, nx + 1
      x = float(j - 1 - ix0)
      r = sqrt(at2 + x*x)
      g = fpotc(r,affp,ari,ifn)
      if (ifun.eq.3) then
         if (r.gt.0.) g = g*(x/r)
      else if (ifun.eq.4) then
         if (r.gt.0.) g = g*(y/r)
      else if (ifun.eq.5) then
         if (r.gt.0.) g = g*(z/r)
      endif
c compare data
      eps = abs(f(j,k,l,m)-g)
      if (eps.gt.epsmax) then
         write (70+kstrt,*) j,k,l,f(j,k,l,m),g,eps
         epsmax = eps
      endif
   10 continue
   20 continue
   30 continue
   40 continue
      write (70+kstrt,*) 'local error = ', epsmax
      call PSUM(epsmax,eps,1,1)
      return
      end
