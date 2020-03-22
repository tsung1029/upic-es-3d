c 3d parallel PIC library for fast fourier transforms
c with 2D domain decomposition
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: april 10, 2008
c-----------------------------------------------------------------------
      subroutine PFFT32R(f,g,h,bs,br,isign,ntpose,mixup,sct,indx,indy,in
     1dz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,kzyp
     2,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd)
c this subroutine performs a three dimensional real to complex fast
c fourier transform and its inverse, using complex arithmetic,
c for data which is distributed in blocks, with 2D spatial decomposition
c for isign = 0, input: isign, indx, indy, indz, kstrt, nxhyzd, nxyzhd
c output: mixup, sct
c for isign = (-1,1), input: all, output: f, g, h, bs, br
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
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
c nxvh/nyv/nzv = first dimension of f/g/h
c kypd/kxypd = second dimension of f/g,h
c kzpd/kyzpd = third dimension of f,g/h
c kxyp/kyp,kyzp/kzp = number of data values per block in x/y/z
c kzyp = maximum(kyzp,kyp)
c jblok/kblok,mblok/lblok = number of data blocks in x/y/z
c jkblok = maximum(jblok,kblok)
c mlblok = maximum(mblok,lblok)
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c the real data is stored in a complex array of length nx/2, ny, nz
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c if ntpose = 0,
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
c if ntpose = 1,
c h(l,j,k,i) = mode jj-1,kk-1,l-1, where jj = j + kxyp*(ix - 1),
c kk = k + kyzp*(iz - 1), and i = iy + nprocy*(iz - 1)
c 1 <= jj <= nx/2, 1 <= kk <= ny, and 1 <= l <= nz, except for
c h(l,1,k,1) = mode nx/2,k-1,l-1, where
c ny/2+2 <= kk <= ny, 1 <= l <= nz, and
c h(l,1,1,1) = mode nx/2,0,l-1
c h(l,1,ny/2+1,1) = mode nx/2,ny/2,l-1, where nz/2+2 <= l <= nz, and
c imaginary part of h(1,1,1,1) = real part of mode nx/2,0,0
c imaginary part of h(1,1,ny/2+1,1) = real part of mode nx/2,ny/2,0
c imaginary part of h(nz/2+1,1,1,1) =
c    real part of mode nx/2,0,nz/2
c imaginary part of h(nz/2+1,1,ny/2+1,1) =
c    real part of mode nx/2,ny/2,nz/2
c using jpl storage convention
c written by viktor k. decyk, ucla
c parallel version
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd, kzyp
      integer jblok, kblok, jkblok, lblok, mblok, mlblok, nxhyzd, nxyzhd
      integer mixup
      complex f, g, h, bs, br, sct
      dimension f(nxvh,kypd,kzpd,kblok*lblok)
      dimension g(nyv,kxypd,kzpd,jblok*lblok)
      dimension h(nzv,kxypd,kyzpd,jblok*mblok)
      dimension bs(kxyp,kzyp,kzp,jkblok*lblok)
      dimension br(kxyp,kzyp,kzp,jblok*mlblok)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nyh, ny2, nz, nzh
      integer nz2, nxyz, nxhyz, j, k, lb, ll, jb, it, nrx, nry, nrz
      integer nxyzh, l, i, m, n, ns, ns2, km, kmr, l1, k1, k2, j1, j2
      integer js, ks, kyb, kzb, kxb, kyzb, klblok, jlblok, jmblok
      integer mx, my, mz, moff
      real dnxyz, arg, ani
      complex s, t, t1
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nz = 2**indz
      nzh = nz/2
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      if (isign) 50, 10, 520
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
   10 do 30 j = 1, nxhyz
      lb = j - 1
      ll = 0
      do 20 k = 1, ndx1yz
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
   50 kxb = nxh/kxyp
      kyb = ny/kyp
      kzb = nz/kzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      if (kstrt.gt.(kyb*kzb)) go to 210
      klblok = kblok*lblok
      nrx = nxhyz/nxh
      do 90 m = 1, klblok
c bit-reverse array elements in x
      do 80 n = 1, kzp
      do 70 j = 1, nxh
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
      nrx = nxyz/nxh
      do 150 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
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
c unscramble coefficients and normalize
      kmr = nxyz/nx
      ani = 0.5/(float(nx)*float(ny)*float(nz))
      nry = nxhyz/ny
      do 200 m = 1, klblok
      do 190 n = 1, kzp
      do 170 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 160 k = 1, kyp
      t = conjg(f(nxh2-j,k,n,m))
      s = f(j,k,n,m) + t
      t = (f(j,k,n,m) - t)*t1
      f(j,k,n,m) = ani*(s + t)
      f(nxh2-j,k,n,m) = ani*conjg(s - t)
  160 continue
  170 continue
      do 180 k = 1, kyp
      f(nxhh+1,k,n,m) = 2.*ani*conjg(f(nxhh+1,k,n,m))
      f(1,k,n,m) = 2.*ani*cmplx(real(f(1,k,n,m)) + aimag(f(1,k,n,m)),rea
     1l(f(1,k,n,m)) - aimag(f(1,k,n,m)))
  180 continue
  190 continue
  200 continue
c transpose f array to g
  210 call PTPOS3A(f,g,bs,br,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,kxypd
     1,kypd,kzpd,jblok,kblok,lblok)
      if (kstrt.gt.(kxb*kzb)) go to 360
      jlblok = jblok*lblok
      nry = nxhyz/ny
      do 250 m = 1, jlblok
c bit-reverse array elements in y
      do 240 n = 1, kzp
      do 230 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 230
      do 220 i = 1, kxyp
      t = g(k1,i,n,m)
      g(k1,i,n,m) = g(k,i,n,m)
      g(k,i,n,m) = t
  220 continue
  230 continue
  240 continue
  250 continue
c then transform in y
      nry = nxyz/ny
      do 310 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 300 m = 1, jlblok
      do 290 n = 1, kzp
      do 280 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 270 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 260 i = 1, kxyp
      t = s*g(j2,i,n,m)
      g(j2,i,n,m) = g(j1,i,n,m) - t
      g(j1,i,n,m) = g(j1,i,n,m) + t
  260 continue
  270 continue
  280 continue
  290 continue
  300 continue
  310 continue
c unscramble modes kx = 0, nx/2
      do 350 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 340 mx = 1, jblok
      if ((mx+js).gt.0) go to 340
      m = mx + moff
      do 330 n = 1, kzp
      do 320 k = 2, nyh
      s = g(ny2-k,1,n,m)
      g(ny2-k,1,n,m) = .5*cmplx(aimag(g(k,1,n,m) + s),real(g(k,1,n,m) -
     1s))
      g(k,1,n,m) = .5*cmplx(real(g(k,1,n,m) + s),aimag(g(k,1,n,m) - s))
  320 continue
  330 continue
  340 continue
  350 continue
c transpose g array to h
  360 call PTPOS3B(g,h,bs,br,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxypd
     1,kyzpd,kzpd,jblok,mblok,lblok)
      kyzb = ny/kyzp
      if (kstrt.gt.(kxb*kyzb)) go to 510
      jmblok = jblok*mblok
      nrz = nxhyz/nz
      do 400 m = 1, jmblok
c bit-reverse array elements in z
      do 390 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 390
      do 380 n = 1, kyzp
      do 370 i = 1, kxyp
      t = h(l1,i,n,m)
      h(l1,i,n,m) = h(l,i,n,m)
      h(l,i,n,m) = t
  370 continue
  380 continue
  390 continue
  400 continue
c finally transform in z
      nrz = nxyz/nz
      do 460 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 450 m = 1, jmblok
      do 440 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 430 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 420 n = 1, kyzp
      do 410 i = 1, kxyp
      t = s*h(j2,i,n,m)
      h(j2,i,n,m) = h(j1,i,n,m) - t
      h(j1,i,n,m) = h(j1,i,n,m) + t
  410 continue
  420 continue
  430 continue
  440 continue
  450 continue
  460 continue
c unscramble modes kx = 0, nx/2
      do 500 my = 1, mblok
      moff = jblok*(my - 1)
      do 490 mx = 1, jblok
      if ((mx+js).gt.0) go to 490
      m = mx + moff
      if ((my+ks).eq.0) then
         do 470 n = 2, nzh
         s = h(nz2-n,1,1,m)
         h(nz2-n,1,1,m) = .5*cmplx(aimag(h(n,1,1,m) + s),real(h(n,1,1,m)
     1 - s))
         h(n,1,1,m) = .5*cmplx(real(h(n,1,1,m) + s),aimag(h(n,1,1,m) - s
     1))
  470    continue
      endif
      kyzb = nyh/kyzp
      if ((my+ks).eq.kyzb) then
         k1 = nyh - kyzb*kyzp + 1
         do 480 n = 2, nzh
         s = h(nz2-n,1,k1,m)
         h(nz2-n,1,k1,m) = .5*cmplx(aimag(h(n,1,k1,m) + s),real(h(n,1,k1
     1,m) - s))
         h(n,1,k1,m) = .5*cmplx(real(h(n,1,k1,m) + s),aimag(h(n,1,k1,m) 
     1- s))
  480    continue
      endif
  490 continue
  500 continue
c transpose h array to f
  510 if (ntpose.eq.0) then
         call PTPOS3B(h,g,br,bs,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kx
     1ypd,kzpd,kyzpd,jblok,lblok,mblok)
         call PTPOS3A(g,f,br,bs,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,ky
     1pd,kxypd,kzpd,kblok,jblok,lblok)
      endif
      return
c forward fourier transform
c transpose f array to h
  520 if (ntpose.eq.0) then
         call PTPOS3A(f,g,bs,br,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,kx
     1ypd,kypd,kzpd,jblok,kblok,lblok)
         call PTPOS3B(g,h,bs,br,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kx
     1ypd,kyzpd,kzpd,jblok,mblok,lblok)
      endif
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      if (kstrt.gt.(kxb*kyzb)) go to 670
      jmblok = jblok*mblok
      nrz = nxhyz/nz
c scramble modes kx = 0, nx/2
      do 560 my = 1, mblok
      moff = jblok*(my - 1)
      do 550 mx = 1, jblok
      if ((mx+js).gt.0) go to 550
      m = mx + moff
      if ((my+ks).eq.0) then
         do 530 n = 2, nzh
         s = cmplx(aimag(h(nz2-n,1,1,m)),real(h(nz2-n,1,1,m)))
         h(nz2-n,1,1,m) = conjg(h(n,1,1,m) - s)
         h(n,1,1,m) = h(n,1,1,m) + s
  530    continue
      endif
      kyzb = nyh/kyzp
      if ((my+ks).eq.kyzb) then
         k1 = nyh - kyzb*kyzp + 1
         do 540 n = 2, nzh
         s = cmplx(aimag(h(nz2-n,1,k1,m)),real(h(nz2-n,1,k1,m)))
         h(nz2-n,1,k1,m) = conjg(h(n,1,k1,m) - s)
         h(n,1,k1,m) = h(n,1,k1,m) + s
  540    continue
      endif
  550 continue
  560 continue
      do 600 m = 1, jmblok
c bit-reverse array elements in z
      do 590 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 590
      do 580 n = 1, kyzp
      do 570 i = 1, kxyp
      t = h(l1,i,n,m)
      h(l1,i,n,m) = h(l,i,n,m)
      h(l,i,n,m) = t
  570 continue
  580 continue
  590 continue
  600 continue
c first transform in z
      nrz = nxyz/nz
      do 660 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 650 m = 1, jmblok
      do 640 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 630 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 620 n = 1, kyzp
      do 610 i = 1, kxyp
      t = s*h(j2,i,n,m)
      h(j2,i,n,m) = h(j1,i,n,m) - t
      h(j1,i,n,m) = h(j1,i,n,m) + t
  610 continue
  620 continue
  630 continue
  640 continue
  650 continue
  660 continue
c transpose h array to g
  670 call PTPOS3B(h,g,br,bs,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kxypd
     1,kzpd,kyzpd,jblok,lblok,mblok)
      kzb = nz/kzp
      if (kstrt.gt.(kxb*kzb)) go to 820
      jlblok = jblok*lblok
      nry = nxhyz/ny
c scramble modes kx = 0, nx/2
      do 710 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 700 mx = 1, jblok
      if ((mx+js).gt.0) go to 700
      m = mx + moff
      do 690 n = 1, kzp
      do 680 k = 2, nyh
      s = cmplx(aimag(g(ny2-k,1,n,m)),real(g(ny2-k,1,n,m)))
      g(ny2-k,1,n,m) = conjg(g(k,1,n,m) - s)
      g(k,1,n,m) = g(k,1,n,m) + s
  680 continue
  690 continue
  700 continue
  710 continue
      do 750 m = 1, jlblok
c bit-reverse array elements in y
      do 740 n = 1, kzp
      do 730 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 730
      do 720 i = 1, kxyp
      t = g(k1,i,n,m)
      g(k1,i,n,m) = g(k,i,n,m)
      g(k,i,n,m) = t
  720 continue
  730 continue
  740 continue
  750 continue
c then transform in y
      nry = nxyz/ny
      do 810 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 800 m = 1, jlblok
      do 790 n = 1, kzp
      do 780 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 770 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 760 i = 1, kxyp
      t = s*g(j2,i,n,m)
      g(j2,i,n,m) = g(j1,i,n,m) - t
      g(j1,i,n,m) = g(j1,i,n,m) + t
  760 continue
  770 continue
  780 continue
  790 continue
  800 continue
  810 continue
c transpose g array to f
  820 call PTPOS3A(g,f,br,bs,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,kypd,
     1kxypd,kzpd,kblok,jblok,lblok)
      kyb = ny/kyp
      if (kstrt.gt.(kyb*kzb)) return
      klblok = kblok*lblok
      nrx = nxhyz/nxh
c scramble coefficients
      kmr = nxyz/nx
      do 890 m = 1, klblok
      do 880 n = 1, kzp
      do 840 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 830 k = 1, kyp
      t = conjg(f(nxh2-j,k,n,m))
      s = f(j,k,n,m) + t
      t = (f(j,k,n,m) - t)*t1
      f(j,k,n,m) = s + t
      f(nxh2-j,k,n,m) = conjg(s - t)
  830 continue
  840 continue
      do 850 k = 1, kyp
      f(nxhh+1,k,n,m) = 2.*conjg(f(nxhh+1,k,n,m))
      f(1,k,n,m) = cmplx(real(f(1,k,n,m)) + aimag(f(1,k,n,m)),real(f(1,k
     1,n,m)) - aimag(f(1,k,n,m)))
  850 continue
c bit-reverse array elements in x
      do 870 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 870
      do 860 i = 1, kyp
      t = f(j1,i,n,m)
      f(j1,i,n,m) = f(j,i,n,m)
      f(j,i,n,m) = t
  860 continue
  870 continue
  880 continue
  890 continue
c finally transform in x
      nrx = nxyz/nxh
      do 950 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 940 m = 1, klblok
      do 930 n = 1, kzp
      do 920 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 910 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 900 i = 1, kyp
      t = s*f(j2,i,n,m)
      f(j2,i,n,m) = f(j1,i,n,m) - t
      f(j1,i,n,m) = f(j1,i,n,m) + t
  900 continue
  910 continue
  920 continue
  930 continue
  940 continue
  950 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT32R3(f,g,h,bs,br,isign,ntpose,mixup,sct,indx,indy,i
     1ndz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,kzy
     2p,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd)
c this subroutine performs 3 three dimensional real to complex fast
c fourier transforms, using complex arithmetic,
c for data which is distributed in blocks, with 2D spatial decomposition
c for isign = 0, input: isign, indx, indy, indz, kstrt, nxhyzd, nxyzhd
c output: mixup, sct
c for isign = (-1,1), input: all, output: f, g, bs, br
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c ntpose = (0,1) = (no,yes) input, output data are transposed
c if isign = 0, the fft tables are prepared
c if isign = -1, three inverse fourier transforms are performed
c if ntpose = 0, f is the input and output array, g, h are scratch arrays
c f(1:3,n,m,l,ld) = (1/nx*ny*nz)*(sum(f(1:3,j,k,i,id)*
c       exp(-sqrt(-1)*2pi*nn*j/nx)*exp(-sqrt(-1)*2pi*mm*kk/ny)*
c       exp(-sqrt(-1)*2pi*l*ii/nz))
c where ll = l + kzp*(ldz - 1) and ii = i + kzp*(idz - 1),
c and mm = m + kyp*(ldy - 1) and kk = k + kyp*(idy - 1),
c and ld = ldy + nprocy*(ldz - 1), id = idy + nprocy*(idz - 1),
c and nprocy = number of processors in y.
c if ntpose = 1, f is the input and h is the output, and g is scratch
c h(1:3,l,n,m,nd) = (1/nx*ny*nz)*sum(f(1:3,j,k,i,id)*
c       exp(-sqrt(-1)*2pi*nn*j/nx)*exp(-sqrt(-1)*2pi*mm*kk/ny)*
c       exp(-sqrt(-1)*2pi*ll*ii/nz))
c where nn = n + kxyp*(ndy - 1) and ii = i + kzp*(idz - 1),
c and mm = m + kyzp*(mdz - 1) and kk = k + kyp*(idy - 1),
c and nd = ndy + nprocy*(ndz - 1), and id = idy + nprocy*(idz - 1)
c if isign = 1, three forward fourier transforms are performed
c if ntpose = 0, f is the input and output array, g, h are scratch arrays
c f(1:3,n,m,l,ld) = (sum(f(1:3,j,k,i,id)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*mm*kk/ny)*exp(sqrt(-1)*2pi*ll*ii/nz))
c if ntpose = 1, h is the input and f is the output, and g is scratch
c f(1:3,j,k,i,id) = sum(h(1:3,l,n,m,nd)*exp(sqrt(-1)*2pi*nn*j/nx)*
c       exp(sqrt(-1)*2pi*mm*kk/ny)*exp(sqrt(-1)*2pi*l*ii/nz))
c bs, br = scratch arrays
c kstrt = starting data block number
c nxvh/nyv/nzv = first dimension of f/g/h
c kypd/kxypd = second dimension of f/g,h
c kzpd/kyzpd = third dimension of f,g/h
c kxyp/kyp,kyzp/kzp = number of data values per block in x/y/z
c kzyp = maximum(kyzp,kyp)
c jblok/kblok,mblok/lblok = number of data blocks in x/y/z
c jkblok = maximum(jblok,kblok)
c mlblok = maximum(mblok,lblok)
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c the real data is stored in a complex array of length nx/2, ny, nz
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c if ntpose = 0,
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
c if ntpose = 1,
c h(l,j,k,i) = mode jj-1,kk-1,l-1, where jj = j + kxyp*(ix - 1),
c kk = k + kyzp*(iz - 1), and i = iy + nprocy*(iz - 1)
c 1 <= jj <= nx/2, 1 <= kk <= ny, and 1 <= l <= nz, except for
c h(l,1,k,1) = mode nx/2,k-1,l-1, where
c ny/2+2 <= kk <= ny, 1 <= l <= nz, and
c h(l,1,1,1) = mode nx/2,0,l-1
c h(l,1,ny/2+1,1) = mode nx/2,ny/2,l-1, where nz/2+2 <= l <= nz, and
c imaginary part of h(1,1,1,1) = real part of mode nx/2,0,0
c imaginary part of h(1,1,ny/2+1,1) = real part of mode nx/2,ny/2,0
c imaginary part of h(nz/2+1,1,1,1) =
c    real part of mode nx/2,0,nz/2
c imaginary part of h(nz/2+1,1,ny/2+1,1) =
c    real part of mode nx/2,ny/2,nz/2
c using jpl storage convention
c written by viktor k. decyk, ucla
c parallel version
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd, kzyp
      integer jblok, kblok, jkblok, lblok, mblok, mlblok, nxhyzd, nxyzhd
      integer mixup
      complex f, g, h, bs, br, sct
      dimension f(3,nxvh,kypd,kzpd,kblok*lblok)
      dimension g(3,nyv,kxypd,kzpd,jblok*lblok)
      dimension h(3,nzv,kxypd,kyzpd,jblok*mblok)
      dimension bs(3,kxyp,kzyp,kzp,jkblok*lblok)
      dimension br(3,kxyp,kzyp,kzp,jblok*mlblok)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nyh, ny2, nz, nzh
      integer nz2, nxyz, nxhyz, j, k, lb, ll, jb, it, nrx, nry, nrz
      integer nxyzh, l, i, m, n, ns, ns2, km, kmr, l1, k1, k2, j1, j2
      integer js, ks, kyb, kzb, kxb, kyzb, klblok, jlblok, jmblok
      integer mx, my, mz, moff, jj
      real dnxyz, arg, ani, at1, at2
      complex s, t, t1, t2, t3
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nz = 2**indz
      nzh = nz/2
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      if (isign) 50, 10, 610
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
   10 do 30 j = 1, nxhyz
      lb = j - 1
      ll = 0
      do 20 k = 1, ndx1yz
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
   50 kxb = nxh/kxyp
      kyb = ny/kyp
      kzb = nz/kzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      if (kstrt.gt.(kyb*kzb)) go to 270
      klblok = kblok*lblok
c swap complex components
      do 90 m = 1, klblok
      do 80 n = 1, kzp
      do 70 i = 1, kyp
      do 60 j = 1, nxh
      at1 = real(f(3,j,i,n,m))
      f(3,j,i,n,m) = cmplx(real(f(2,j,i,n,m)),aimag(f(3,j,i,n,m)))
      at2 = aimag(f(2,j,i,n,m))
      f(2,j,i,n,m) = cmplx(aimag(f(1,j,i,n,m)),at1)
      f(1,j,i,n,m) = cmplx(real(f(1,j,i,n,m)),at2)
   60 continue
   70 continue
   80 continue
   90 continue
      nrx = nxhyz/nxh
      do 130 m = 1, klblok
c bit-reverse array elements in x
      do 120 n = 1, kzp
      do 110 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 110
      do 100 i = 1, kyp
      t1 = f(1,j1,i,n,m)
      t2 = f(2,j1,i,n,m)
      t3 = f(3,j1,i,n,m)
      f(1,j1,i,n,m) = f(1,j,i,n,m)
      f(2,j1,i,n,m) = f(2,j,i,n,m)
      f(3,j1,i,n,m) = f(3,j,i,n,m)
      f(1,j,i,n,m) = t1
      f(2,j,i,n,m) = t2
      f(3,j,i,n,m) = t3
  100 continue
  110 continue
  120 continue
  130 continue
c first transform in x
      nrx = nxyz/nxh
      do 190 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 180 m = 1, klblok
      do 170 n = 1, kzp
      do 160 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 150 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 140 i = 1, kyp
      t1 = s*f(1,j2,i,n,m)
      t2 = s*f(2,j2,i,n,m)
      t3 = s*f(3,j2,i,n,m)
      f(1,j2,i,n,m) = f(1,j1,i,n,m) - t1
      f(2,j2,i,n,m) = f(2,j1,i,n,m) - t2
      f(3,j2,i,n,m) = f(3,j1,i,n,m) - t3
      f(1,j1,i,n,m) = f(1,j1,i,n,m) + t1
      f(2,j1,i,n,m) = f(2,j1,i,n,m) + t2
      f(3,j1,i,n,m) = f(3,j1,i,n,m) + t3
  140 continue
  150 continue
  160 continue
  170 continue
  180 continue
  190 continue
c unscramble coefficients and normalize
      kmr = nxyz/nx
      ani = 0.5/(float(nx)*float(ny)*float(nz))
      nry = nxhyz/ny
      do 260 m = 1, klblok
      do 250 n = 1, kzp
      do 220 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 210 k = 1, kyp
      do 200 jj = 1, 3
      t = conjg(f(jj,nxh2-j,k,n,m))
      s = f(jj,j,k,n,m) + t
      t = (f(jj,j,k,n,m) - t)*t1
      f(jj,j,k,n,m) = ani*(s + t)
      f(jj,nxh2-j,k,n,m) = ani*conjg(s - t)
  200 continue
  210 continue
  220 continue
      do 240 k = 1, kyp
      do 230 jj = 1, 3
      f(jj,nxhh+1,k,n,m) = 2.*ani*conjg(f(jj,nxhh+1,k,n,m))
      f(jj,1,k,n,m) = 2.*ani*cmplx(real(f(jj,1,k,n,m)) + aimag(f(jj,1,k,
     1n,m)),real(f(jj,1,k,n,m)) - aimag(f(jj,1,k,n,m)))
  230 continue
  240 continue
  250 continue
  260 continue
c transpose f array to g
  270 call P3TPOS3A(f,g,bs,br,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,kxyp
     1d,kypd,kzpd,jblok,kblok,lblok)
      if (kstrt.gt.(kxb*kzb)) go to 430
      jlblok = jblok*lblok
      nry = nxhyz/ny
      do 310 m = 1, jlblok
c bit-reverse array elements in y
      do 300 n = 1, kzp
      do 290 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 290
      do 280 i = 1, kxyp
      t1 = g(1,k1,i,n,m)
      t2 = g(2,k1,i,n,m)
      t3 = g(3,k1,i,n,m)
      g(1,k1,i,n,m) = g(1,k,i,n,m)
      g(2,k1,i,n,m) = g(2,k,i,n,m)
      g(3,k1,i,n,m) = g(3,k,i,n,m)
      g(1,k,i,n,m) = t1
      g(2,k,i,n,m) = t2
      g(3,k,i,n,m) = t3
  280 continue
  290 continue
  300 continue
  310 continue
c then transform in y
      nry = nxyz/ny
      do 370 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 360 m = 1, jlblok
      do 350 n = 1, kzp
      do 340 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 330 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      do 320 i = 1, kxyp
      t1 = s*g(1,j2,i,n,m)
      t2 = s*g(2,j2,i,n,m)
      t3 = s*g(3,j2,i,n,m)
      g(1,j2,i,n,m) = g(1,j1,i,n,m) - t1
      g(2,j2,i,n,m) = g(2,j1,i,n,m) - t2
      g(3,j2,i,n,m) = g(3,j1,i,n,m) - t3
      g(1,j1,i,n,m) = g(1,j1,i,n,m) + t1
      g(2,j1,i,n,m) = g(2,j1,i,n,m) + t2
      g(3,j1,i,n,m) = g(3,j1,i,n,m) + t3
  320 continue
  330 continue
  340 continue
  350 continue
  360 continue
  370 continue
c unscramble modes kx = 0, nx/2
      do 420 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 410 mx = 1, jblok
      if ((mx+js).gt.0) go to 410
      m = mx + moff
      do 400 n = 1, kzp
      do 390 k = 2, nyh
      do 380 jj = 1, 3
      s = g(jj,ny2-k,1,n,m)
      g(jj,ny2-k,1,n,m) = .5*cmplx(aimag(g(jj,k,1,n,m) + s),real(g(jj,k,
     11,n,m) - s))
      g(jj,k,1,n,m) = .5*cmplx(real(g(jj,k,1,n,m) + s),aimag(g(jj,k,1,n,
     1m) - s))
  380 continue
  390 continue
  400 continue
  410 continue
  420 continue
c transpose g array to h
  430 call P3TPOS3B(g,h,bs,br,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxyp
     1d,kyzpd,kzpd,jblok,mblok,lblok)
      kyzb = ny/kyzp
      if (kstrt.gt.(kxb*kyzb)) go to 600
      jmblok = jblok*mblok
      nrz = nxhyz/nz
      do 470 m = 1, jmblok
c bit-reverse array elements in z
      do 460 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 460
      do 450 n = 1, kyzp
      do 440 i = 1, kxyp
      t1 = h(1,l1,i,n,m)
      t2 = h(2,l1,i,n,m)
      t3 = h(3,l1,i,n,m)
      h(1,l1,i,n,m) = h(1,l,i,n,m)
      h(2,l1,i,n,m) = h(2,l,i,n,m)
      h(3,l1,i,n,m) = h(3,l,i,n,m)
      h(1,l,i,n,m) = t1
      h(2,l,i,n,m) = t2
      h(3,l,i,n,m) = t3
  440 continue
  450 continue
  460 continue
  470 continue
c finally transform in z
      nrz = nxyz/nz
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
      s = sct(1+kmr*(j-1))
      do 490 n = 1, kyzp
      do 480 i = 1, kxyp
      t1 = s*h(1,j2,i,n,m)
      t2 = s*h(2,j2,i,n,m)
      t3 = s*h(3,j2,i,n,m)
      h(1,j2,i,n,m) = h(1,j1,i,n,m) - t1
      h(2,j2,i,n,m) = h(2,j1,i,n,m) - t2
      h(3,j2,i,n,m) = h(3,j1,i,n,m) - t3
      h(1,j1,i,n,m) = h(1,j1,i,n,m) + t1
      h(2,j1,i,n,m) = h(2,j1,i,n,m) + t2
      h(3,j1,i,n,m) = h(3,j1,i,n,m) + t3
  480 continue
  490 continue
  500 continue
  510 continue
  520 continue
  530 continue
c unscramble modes kx = 0, nx/2
      do 590 my = 1, mblok
      moff = jblok*(my - 1)
      do 580 mx = 1, jblok
      if ((mx+js).gt.0) go to 580
      m = mx + moff
      if ((my+ks).eq.0) then
         do 550 n = 2, nzh
         do 540 jj = 1, 3
         s = h(jj,nz2-n,1,1,m)
         h(jj,nz2-n,1,1,m) = .5*cmplx(aimag(h(jj,n,1,1,m) + s),real(h(jj
     1,n,1,1,m) - s))

         h(jj,n,1,1,m) = .5*cmplx(real(h(jj,n,1,1,m) + s),aimag(h(jj,n,1
     1,1,m) - s))
  540    continue
  550    continue
      endif
      kyzb = nyh/kyzp
      if ((my+ks).eq.kyzb) then
         k1 = nyh - kyzb*kyzp + 1
         do 570 n = 2, nzh
         do 560 jj = 1, 3
         s = h(jj,nz2-n,1,k1,m)
         h(jj,nz2-n,1,k1,m) = .5*cmplx(aimag(h(jj,n,1,k1,m) + s),real(h(
     1jj,n,1,k1,m) - s))
         h(jj,n,1,k1,m) = .5*cmplx(real(h(jj,n,1,k1,m) + s),aimag(h(jj,n
     1,1,k1,m) - s))
  560    continue
  570    continue
      endif
  580 continue
  590 continue
c transpose h array to f
  600 if (ntpose.eq.0) then
         call P3TPOS3B(h,g,br,bs,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,k
     1xypd,kzpd,kyzpd,jblok,lblok,mblok)
         call P3TPOS3A(g,f,br,bs,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,k
     1ypd,kxypd,kzpd,kblok,jblok,lblok)
      endif
      return
c forward fourier transform
c transpose f array to h
  610 if (ntpose.eq.0) then
         call P3TPOS3A(f,g,bs,br,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,k
     1xypd,kypd,kzpd,jblok,kblok,lblok)
         call P3TPOS3B(g,h,bs,br,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,k
     1xypd,kyzpd,kzpd,jblok,mblok,lblok)
      endif
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      if (kstrt.gt.(kxb*kyzb)) go to 780
      jmblok = jblok*mblok
      nrz = nxhyz/nz
c scramble modes kx = 0, nx/2
      do 670 my = 1, mblok
      moff = jblok*(my - 1)
      do 660 mx = 1, jblok
      if ((mx+js).gt.0) go to 660
      m = mx + moff
      if ((my+ks).eq.0) then
         do 630 n = 2, nzh
         do 620 jj = 1, 3
         s = cmplx(aimag(h(jj,nz2-n,1,1,m)),real(h(jj,nz2-n,1,1,m)))
         h(jj,nz2-n,1,1,m) = conjg(h(jj,n,1,1,m) - s)
         h(jj,n,1,1,m) = h(jj,n,1,1,m) + s
  620    continue
  630    continue
      endif
      kyzb = nyh/kyzp
      if ((my+ks).eq.kyzb) then
         k1 = nyh - kyzb*kyzp + 1
         do 650 n = 2, nzh
         do 640 jj = 1, 3
         s = cmplx(aimag(h(jj,nz2-n,1,k1,m)),real(h(jj,nz2-n,1,k1,m)))
         h(jj,nz2-n,1,k1,m) = conjg(h(jj,n,1,k1,m) - s)
         h(jj,n,1,k1,m) = h(jj,n,1,k1,m) + s
  640    continue
  650    continue
      endif
  660 continue
  670 continue
      do 710 m = 1, jmblok
c bit-reverse array elements in z
      do 700 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 700
      do 690 n = 1, kyzp
      do 680 i = 1, kxyp
      t1 = h(1,l1,i,n,m)
      t2 = h(2,l1,i,n,m)
      t3 = h(3,l1,i,n,m)
      h(1,l1,i,n,m) = h(1,l,i,n,m)
      h(2,l1,i,n,m) = h(2,l,i,n,m)
      h(3,l1,i,n,m) = h(3,l,i,n,m)
      h(1,l,i,n,m) = t1
      h(2,l,i,n,m) = t2
      h(3,l,i,n,m) = t3
  680 continue
  690 continue
  700 continue
  710 continue
c first transform in z
      nrz = nxyz/nz
      do 770 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 760 m = 1, jmblok
      do 750 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 740 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 730 n = 1, kyzp
      do 720 i = 1, kxyp
      t1 = s*h(1,j2,i,n,m)
      t2 = s*h(2,j2,i,n,m)
      t3 = s*h(3,j2,i,n,m)
      h(1,j2,i,n,m) = h(1,j1,i,n,m) - t1
      h(2,j2,i,n,m) = h(2,j1,i,n,m) - t2
      h(3,j2,i,n,m) = h(3,j1,i,n,m) - t3
      h(1,j1,i,n,m) = h(1,j1,i,n,m) + t1
      h(2,j1,i,n,m) = h(2,j1,i,n,m) + t2
      h(3,j1,i,n,m) = h(3,j1,i,n,m) + t3
  720 continue
  730 continue
  740 continue
  750 continue
  760 continue
  770 continue
c transpose h array to g
  780 call P3TPOS3B(h,g,br,bs,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kxyp
     1d,kzpd,kyzpd,jblok,lblok,mblok)
      kzb = nz/kzp
      if (kstrt.gt.(kxb*kzb)) go to 940
      jlblok = jblok*lblok
      nry = nxhyz/ny
c scramble modes kx = 0, nx/2
      do 830 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 820 mx = 1, jblok
      if ((mx+js).gt.0) go to 820
      m = mx + moff
      do 810 n = 1, kzp
      do 800 k = 2, nyh
      do 790 jj = 1, 3
      s = cmplx(aimag(g(jj,ny2-k,1,n,m)),real(g(jj,ny2-k,1,n,m)))
      g(jj,ny2-k,1,n,m) = conjg(g(jj,k,1,n,m) - s)
      g(jj,k,1,n,m) = g(jj,k,1,n,m) + s
  790 continue
  800 continue
  810 continue
  820 continue
  830 continue
      do 870 m = 1, jlblok
c bit-reverse array elements in y
      do 860 n = 1, kzp
      do 850 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 850
      do 840 i = 1, kxyp
      t1 = g(1,k1,i,n,m)
      t2 = g(2,k1,i,n,m)
      t3 = g(3,k1,i,n,m)
      g(1,k1,i,n,m) = g(1,k,i,n,m)
      g(2,k1,i,n,m) = g(2,k,i,n,m)
      g(3,k1,i,n,m) = g(3,k,i,n,m)
      g(1,k,i,n,m) = t1
      g(2,k,i,n,m) = t2
      g(3,k,i,n,m) = t3
  840 continue
  850 continue
  860 continue
  870 continue
c then transform in y
      nry = nxyz/ny
      do 930 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 920 m = 1, jlblok
      do 910 n = 1, kzp
      do 900 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 890 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 880 i = 1, kxyp
      t1 = s*g(1,j2,i,n,m)
      t2 = s*g(2,j2,i,n,m)
      t3 = s*g(3,j2,i,n,m)
      g(1,j2,i,n,m) = g(1,j1,i,n,m) - t1
      g(2,j2,i,n,m) = g(2,j1,i,n,m) - t2
      g(3,j2,i,n,m) = g(3,j1,i,n,m) - t3
      g(1,j1,i,n,m) = g(1,j1,i,n,m) + t1
      g(2,j1,i,n,m) = g(2,j1,i,n,m) + t2
      g(3,j1,i,n,m) = g(3,j1,i,n,m) + t3
  880 continue
  890 continue
  900 continue
  910 continue
  920 continue
  930 continue
c transpose g array to f
  940 call P3TPOS3A(g,f,br,bs,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,kypd
     1,kxypd,kzpd,kblok,jblok,lblok)
      kyb = ny/kyp
      if (kstrt.gt.(kyb*kzb)) return
      klblok = kblok*lblok
      nrx = nxhyz/nxh
c scramble coefficients
      kmr = nxyz/nx
      do 1030 m = 1, klblok
      do 1020 n = 1, kzp
      do 970 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 960 k = 1, kyp
      do 950 jj = 1, 3
      t = conjg(f(jj,nxh2-j,k,n,m))
      s = f(jj,j,k,n,m) + t
      t = (f(jj,j,k,n,m) - t)*t1
      f(jj,j,k,n,m) = s + t
      f(jj,nxh2-j,k,n,m) = conjg(s - t)
  950 continue
  960 continue
  970 continue
      do 990 k = 1, kyp
      do 980 jj = 1, 3
      f(jj,nxhh+1,k,n,m) = 2.*conjg(f(jj,nxhh+1,k,n,m))
      f(jj,1,k,n,m) = cmplx(real(f(jj,1,k,n,m)) + aimag(f(jj,1,k,n,m)),r
     1eal(f(jj,1,k,n,m)) - aimag(f(jj,1,k,n,m)))
  980 continue
  990 continue
c bit-reverse array elements in x
      do 1010 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 1010
      do 1000 i = 1, kyp
      t1 = f(1,j1,i,n,m)
      t2 = f(2,j1,i,n,m)
      t3 = f(3,j1,i,n,m)
      f(1,j1,i,n,m) = f(1,j,i,n,m)
      f(2,j1,i,n,m) = f(2,j,i,n,m)
      f(3,j1,i,n,m) = f(3,j,i,n,m)
      f(1,j,i,n,m) = t1
      f(2,j,i,n,m) = t2
      f(3,j,i,n,m) = t3
 1000 continue
 1010 continue
 1020 continue
 1030 continue
c finally transform in x
      nrx = nxyz/nxh
      do 1090 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 1080 m = 1, klblok
      do 1070 n = 1, kzp
      do 1060 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 1050 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 1040 i = 1, kyp
      t1 = s*f(1,j2,i,n,m)
      t2 = s*f(2,j2,i,n,m)
      t3 = s*f(3,j2,i,n,m)
      f(1,j2,i,n,m) = f(1,j1,i,n,m) - t1
      f(2,j2,i,n,m) = f(2,j1,i,n,m) - t2
      f(3,j2,i,n,m) = f(3,j1,i,n,m) - t3
      f(1,j1,i,n,m) = f(1,j1,i,n,m) + t1
      f(2,j1,i,n,m) = f(2,j1,i,n,m) + t2
      f(3,j1,i,n,m) = f(3,j1,i,n,m) + t3
 1040 continue
 1050 continue
 1060 continue
 1070 continue
 1080 continue
 1090 continue
c swap complex components
      do 1130 m = 1, klblok
      do 1120 n = 1, kzp 
      do 1110 i = 1, kyp
      do 1100 j = 1, nxh
      at1 = real(f(3,j,i,n,m))
      f(3,j,i,n,m) = cmplx(aimag(f(2,j,i,n,m)),aimag(f(3,j,i,n,m)))
      at2 = real(f(2,j,i,n,m))
      f(2,j,i,n,m) = cmplx(at1,aimag(f(1,j,i,n,m)))
      f(1,j,i,n,m) = cmplx(real(f(1,j,i,n,m)),at2)
 1100 continue
 1110 continue
 1120 continue
 1130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT32RX(f,g,h,isign,ntpose,mixup,sct,indx,indy,indz,ks
     1trt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,jblok,kbl
     2ok,lblok,mblok,nxhyzd,nxyzhd)
c this subroutine performs a three dimensional real to complex fast
c fourier transform and its inverse, using complex arithmetic,
c for data which is distributed in blocks, with 2D spatial decomposition
c for isign = 0, input: isign, indx, indy, indz, kstrt, nxhyzd, nxyzhd
c output: mixup, sct
c for isign = (-1,1), input: all, output: f, g, h, bs, br
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
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
c kstrt = starting data block number
c nxvh/nyv/nzv = first dimension of f/g/h
c kypd/kxypd = second dimension of f/g,h
c kzpd/kyzpd = third dimension of f,g/h
c kxyp/kyp,kyzp/kzp = number of data values per block in x/y/z
c jblok/kblok,mblok/lblok = number of data blocks in x/y/z
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c the real data is stored in a complex array of length nx/2, ny, nz
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c if ntpose = 0,
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
c if ntpose = 1,
c h(l,j,k,i) = mode jj-1,kk-1,l-1, where jj = j + kxyp*(ix - 1),
c kk = k + kyzp*(iz - 1), and i = iy + nprocy*(iz - 1)
c 1 <= jj <= nx/2, 1 <= kk <= ny, and 1 <= l <= nz, except for
c h(l,1,k,1) = mode nx/2,k-1,l-1, where
c ny/2+2 <= kk <= ny, 1 <= l <= nz, and
c h(l,1,1,1) = mode nx/2,0,l-1
c h(l,1,ny/2+1,1) = mode nx/2,ny/2,l-1, where nz/2+2 <= l <= nz, and
c imaginary part of h(1,1,1,1) = real part of mode nx/2,0,0
c imaginary part of h(1,1,ny/2+1,1) = real part of mode nx/2,ny/2,0
c imaginary part of h(nz/2+1,1,1,1) =
c    real part of mode nx/2,0,nz/2
c imaginary part of h(nz/2+1,1,ny/2+1,1) =
c    real part of mode nx/2,ny/2,nz/2
c using jpl storage convention
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd
      integer jblok, kblok, lblok, mblok, nxhyzd, nxyzhd
      integer mixup
      complex f, g, h, sct
      dimension f(nxvh,kypd,kzpd,kblok*lblok)
      dimension g(nyv,kxypd,kzpd,jblok*lblok)
      dimension h(nzv,kxypd,kyzpd,jblok*mblok)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nyh, ny2, nz, nzh
      integer nz2, nxyz, nxhyz, j, k, lb, ll, jb, it, nrx, nry, nrz
      integer nxyzh, l, i, m, n, ns, ns2, km, kmr, l1, k1, k2, j1, j2
      integer js, ks, kyb, kzb, kxb, kyzb, klblok, jlblok, jmblok
      integer mx, my, mz, moff
      real dnxyz, arg, ani
      complex s, t, t1
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nz = 2**indz
      nzh = nz/2
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      if (isign) 50, 10, 440
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
   10 do 30 j = 1, nxhyz
      lb = j - 1
      ll = 0
      do 20 k = 1, ndx1yz
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
   50 kxb = nxh/kxyp
      kyb = ny/kyp
      kzb = nz/kzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      if (kstrt.gt.(kyb*kzb)) go to 170
      klblok = kblok*lblok
      ani = 0.5/(float(nx)*float(ny)*float(nz))
      do 160 m = 1, klblok
      do 150 n = 1, kzp
      nrx = nxhyz/nxh
c bit-reverse array elements in x
      do 70 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 70
      do 60 i = 1, kyp
      t = f(j1,i,n,m)
      f(j1,i,n,m) = f(j,i,n,m)
      f(j,i,n,m) = t
   60 continue
   70 continue
c first transform in x
      nrx = nxyz/nxh
      do 110 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 100 i = 1, kyp
      do 90 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 80 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t = s*f(j2,i,n,m)
      f(j2,i,n,m) = f(j1,i,n,m) - t
      f(j1,i,n,m) = f(j1,i,n,m) + t
   80 continue
   90 continue
  100 continue
  110 continue
c unscramble coefficients and normalize
      kmr = nxyz/nx
      nry = nxhyz/ny
      do 130 k = 1, kyp
      do 120 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      t = conjg(f(nxh2-j,k,n,m))
      s = f(j,k,n,m) + t
      t = (f(j,k,n,m) - t)*t1
      f(j,k,n,m) = ani*(s + t)
      f(nxh2-j,k,n,m) = ani*conjg(s - t)
  120 continue
  130 continue
      do 140 k = 1, kyp
      f(nxhh+1,k,n,m) = 2.*ani*conjg(f(nxhh+1,k,n,m))
      f(1,k,n,m) = 2.*ani*cmplx(real(f(1,k,n,m)) + aimag(f(1,k,n,m)),rea
     1l(f(1,k,n,m)) - aimag(f(1,k,n,m)))
  140 continue
  150 continue
  160 continue
c transpose f array to g
  170 call PTPOS3AX(f,g,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,kxypd,kypd
     1,kzpd,jblok,kblok,lblok)
      if (kstrt.gt.(kxb*kzb)) go to 300
      jlblok = jblok*lblok
      do 250 m = 1, jlblok
      do 240 n = 1, kzp
      nry = nxhyz/ny
c bit-reverse array elements in y
      do 190 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 190
      do 180 i = 1, kxyp
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
      do 220 i = 1, kxyp
      do 210 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 200 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t = s*g(j2,i,n,m)
      g(j2,i,n,m) = g(j1,i,n,m) - t
      g(j1,i,n,m) = g(j1,i,n,m) + t
  200 continue
  210 continue
  220 continue
  230 continue
  240 continue
  250 continue
c unscramble modes kx = 0, nx/2
      do 290 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 280 mx = 1, jblok
      if ((mx+js).gt.0) go to 280
      m = mx + moff
      do 270 n = 1, kzp
      do 260 k = 2, nyh
      s = g(ny2-k,1,n,m)
      g(ny2-k,1,n,m) = .5*cmplx(aimag(g(k,1,n,m) + s),real(g(k,1,n,m) -
     1s))
      g(k,1,n,m) = .5*cmplx(real(g(k,1,n,m) + s),aimag(g(k,1,n,m) - s))
  260 continue
  270 continue
  280 continue
  290 continue
c transpose g array to h
  300 call PTPOS3BX(g,h,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxypd,kyzp
     1d,kzpd,jblok,mblok,lblok)
      kyzb = ny/kyzp
      if (kstrt.gt.(kxb*kyzb)) go to 430
      jmblok = jblok*mblok
      do 380 m = 1, jmblok
      do 370 n = 1, kyzp
      nrz = nxhyz/nz
c bit-reverse array elements in z
      do 320 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 320
      do 310 i = 1, kxyp
      t = h(l1,i,n,m)
      h(l1,i,n,m) = h(l,i,n,m)
      h(l,i,n,m) = t
  310 continue
  320 continue
c finally transform in z
      nrz = nxyz/nz
      do 360 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 350 i = 1, kxyp
      do 340 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 330 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = sct(1+kmr*(j-1))
      t = s*h(j2,i,n,m)
      h(j2,i,n,m) = h(j1,i,n,m) - t
      h(j1,i,n,m) = h(j1,i,n,m) + t
  330 continue
  340 continue
  350 continue
  360 continue
  370 continue
  380 continue
c unscramble modes kx = 0, nx/2
      do 420 my = 1, mblok
      moff = jblok*(my - 1)
      do 410 mx = 1, jblok
      if ((mx+js).gt.0) go to 410
      m = mx + moff
      if ((my+ks).eq.0) then
         do 390 n = 2, nzh
         s = h(nz2-n,1,1,m)
         h(nz2-n,1,1,m) = .5*cmplx(aimag(h(n,1,1,m) + s),real(h(n,1,1,m)
     1 - s))
         h(n,1,1,m) = .5*cmplx(real(h(n,1,1,m) + s),aimag(h(n,1,1,m) - s
     1))
  390    continue
      endif
      kyzb = nyh/kyzp
      if ((my+ks).eq.kyzb) then
         k1 = nyh - kyzb*kyzp + 1
         do 400 n = 2, nzh
         s = h(nz2-n,1,k1,m)
         h(nz2-n,1,k1,m) = .5*cmplx(aimag(h(n,1,k1,m) + s),real(h(n,1,k1
     1,m) - s))
         h(n,1,k1,m) = .5*cmplx(real(h(n,1,k1,m) + s),aimag(h(n,1,k1,m) 
     1- s))
  400 continue
      endif
  410 continue
  420 continue
c transpose h array to f
  430 if (ntpose.eq.0) then
         call PTPOS3BX(h,g,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kxypd,k
     1zpd,kyzpd,jblok,lblok,mblok)
         call PTPOS3AX(g,f,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,kypd,kx
     1ypd,kzpd,kblok,jblok,lblok)
      endif
      return
c forward fourier transform
c transpose f array to h
  440 if (ntpose.eq.0) then
         call PTPOS3AX(f,g,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,kxypd,k
     1ypd,kzpd,jblok,kblok,lblok)
         call PTPOS3BX(g,h,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxypd,k
     1yzpd,kzpd,jblok,mblok,lblok)
      endif
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      if (kstrt.gt.(kxb*kyzb)) go to 570
      jmblok = jblok*mblok
c scramble modes kx = 0, nx/2
      do 480 my = 1, mblok
      moff = jblok*(my - 1)
      do 470 mx = 1, jblok
      if ((mx+js).gt.0) go to 470
      m = mx + moff
      if ((my+ks).eq.0) then
         do 450 n = 2, nzh
         s = cmplx(aimag(h(nz2-n,1,1,m)),real(h(nz2-n,1,1,m)))
         h(nz2-n,1,1,m) = conjg(h(n,1,1,m) - s)
         h(n,1,1,m) = h(n,1,1,m) + s
  450    continue
      endif
      kyzb = nyh/kyzp
      if ((my+ks).eq.kyzb) then
         k1 = nyh - kyzb*kyzp + 1
         do 460 n = 2, nzh
         s = cmplx(aimag(h(nz2-n,1,k1,m)),real(h(nz2-n,1,k1,m)))
         h(nz2-n,1,k1,m) = conjg(h(n,1,k1,m) - s)
         h(n,1,k1,m) = h(n,1,k1,m) + s
  460    continue
      endif
  470 continue
  480 continue
      do 560 m = 1, jmblok
      do 550 n = 1, kyzp
      nrz = nxhyz/nz
c bit-reverse array elements in z
      do 500 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 500
      do 490 i = 1, kxyp
      t = h(l1,i,n,m)
      h(l1,i,n,m) = h(l,i,n,m)
      h(l,i,n,m) = t
  490 continue
  500 continue
c first transform in z
      nrz = nxyz/nz
      do 540 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 530 i = 1, kxyp
      do 520 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 510 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t = s*h(j2,i,n,m)
      h(j2,i,n,m) = h(j1,i,n,m) - t
      h(j1,i,n,m) = h(j1,i,n,m) + t
  510 continue
  520 continue
  530 continue
  540 continue
  550 continue
  560 continue
c transpose h array to g
  570 call PTPOS3BX(h,g,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kxypd,kzpd
     1,kyzpd,jblok,lblok,mblok)
      kzb = nz/kzp
      if (kstrt.gt.(kxb*kzb)) go to 700
      jlblok = jblok*lblok
c scramble modes kx = 0, nx/2
      do 610 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 600 mx = 1, jblok
      if ((mx+js).gt.0) go to 600
      m = mx + moff
      do 590 n = 1, kzp
      do 580 k = 2, nyh
      s = cmplx(aimag(g(ny2-k,1,n,m)),real(g(ny2-k,1,n,m)))
      g(ny2-k,1,n,m) = conjg(g(k,1,n,m) - s)
      g(k,1,n,m) = g(k,1,n,m) + s
  580 continue
  590 continue
  600 continue
  610 continue
      do 690 m = 1, jlblok
      do 680 n = 1, kzp
      nry = nxhyz/ny
c bit-reverse array elements in y
      do 630 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 630
      do 620 i = 1, kxyp
      t = g(k1,i,n,m)
      g(k1,i,n,m) = g(k,i,n,m)
      g(k,i,n,m) = t
  620 continue
  630 continue
c then transform in y
      nry = nxyz/ny
      do 670 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 660 i = 1, kxyp
      do 650 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 640 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t = s*g(j2,i,n,m)
      g(j2,i,n,m) = g(j1,i,n,m) - t
      g(j1,i,n,m) = g(j1,i,n,m) + t
  640 continue
  650 continue
  660 continue
  670 continue
  680 continue
  690 continue
c transpose g array to f
  700 call PTPOS3AX(g,f,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,kypd,kxypd
     1,kzpd,kblok,jblok,lblok)
      kyb = ny/kyp
      if (kstrt.gt.(kyb*kzb)) return
      klblok = kblok*lblok
      do 810 m = 1, klblok
      do 800 n = 1, kzp
c scramble coefficients
      kmr = nxyz/nx
      do 720 k = 1, kyp
      do 710 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      t = conjg(f(nxh2-j,k,n,m))
      s = f(j,k,n,m) + t
      t = (f(j,k,n,m) - t)*t1
      f(j,k,n,m) = s + t
      f(nxh2-j,k,n,m) = conjg(s - t)
  710 continue
  720 continue
      do 730 k = 1, kyp
      f(nxhh+1,k,n,m) = 2.*conjg(f(nxhh+1,k,n,m))
      f(1,k,n,m) = cmplx(real(f(1,k,n,m)) + aimag(f(1,k,n,m)),real(f(1,k
     1,n,m)) - aimag(f(1,k,n,m)))
  730 continue
      nrx = nxhyz/nxh
c bit-reverse array elements in x
      do 750 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 750
      do 740 i = 1, kyp
      t = f(j1,i,n,m)
      f(j1,i,n,m) = f(j,i,n,m)
      f(j,i,n,m) = t
  740 continue
  750 continue
c finally transform in x
      nrx = nxyz/nxh
      do 790 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 780 i = 1, kyp
      do 770 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 760 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      t = s*f(j2,i,n,m)
      f(j2,i,n,m) = f(j1,i,n,m) - t
      f(j1,i,n,m) = f(j1,i,n,m) + t
  760 continue
  770 continue
  780 continue
  790 continue
  800 continue
  810 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT32RX3(f,g,h,isign,ntpose,mixup,sct,indx,indy,indz,k
     1strt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,jblok,kb
     2lok,lblok,mblok,nxhyzd,nxyzhd)
c this subroutine performs 3 three dimensional real to complex fast
c fourier transforms, using complex arithmetic,
c for data which is distributed in blocks, with 2D spatial decomposition
c for isign = 0, input: isign, indx, indy, indz, kstrt, nxhyzd, nxyzhd
c output: mixup, sct
c for isign = (-1,1), input: all, output: f, g
c for isign = -1, approximate flop count: N*(5*log2(N) + 19/2)/nvp
c for isign = 1,  approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c ntpose = (0,1) = (no,yes) input, output data are transposed
c if isign = 0, the fft tables are prepared
c if isign = -1, three inverse fourier transforms are performed
c if ntpose = 0, f is the input and output array, g, h are scratch arrays
c f(1:3,n,m,l,ld) = (1/nx*ny*nz)*(sum(f(1:3,j,k,i,id)*
c       exp(-sqrt(-1)*2pi*nn*j/nx)*exp(-sqrt(-1)*2pi*mm*kk/ny)*
c       exp(-sqrt(-1)*2pi*l*ii/nz))
c where ll = l + kzp*(ldz - 1) and ii = i + kzp*(idz - 1),
c and mm = m + kyp*(ldy - 1) and kk = k + kyp*(idy - 1),
c and ld = ldy + nprocy*(ldz - 1), id = idy + nprocy*(idz - 1),
c and nprocy = number of processors in y.
c if ntpose = 1, f is the input and h is the output, and g is scratch
c h(1:3,l,n,m,nd) = (1/nx*ny*nz)*sum(f(1:3,j,k,i,id)*
c       exp(-sqrt(-1)*2pi*nn*j/nx)*exp(-sqrt(-1)*2pi*mm*kk/ny)*
c       exp(-sqrt(-1)*2pi*ll*ii/nz))
c where nn = n + kxyp*(ndy - 1) and ii = i + kzp*(idz - 1),
c and mm = m + kyzp*(mdz - 1) and kk = k + kyp*(idy - 1),
c and nd = ndy + nprocy*(ndz - 1), and id = idy + nprocy*(idz - 1)
c if isign = 1, three forward fourier transforms are performed
c if ntpose = 0, f is the input and output array, g, h are scratch arrays
c f(1:3,n,m,l,ld) = (sum(f(1:3,j,k,i,id)*exp(sqrt(-1)*2pi*n*j/nx)*
c       exp(sqrt(-1)*2pi*mm*kk/ny)*exp(sqrt(-1)*2pi*ll*ii/nz))
c if ntpose = 1, h is the input and f is the output, and g is scratch
c f(1:3,j,k,i,id) = sum(h(1:3,l,n,m,nd)*exp(sqrt(-1)*2pi*nn*j/nx)*
c       exp(sqrt(-1)*2pi*mm*kk/ny)*exp(sqrt(-1)*2pi*l*ii/nz))
c kstrt = starting data block number
c nxvh/nyv/nzv = first dimension of f/g/h
c kypd/kxypd = second dimension of f/g,h
c kzpd/kyzpd = third dimension of f,g/h
c kxyp/kyp,kyzp/kzp = number of data values per block in x/y/z
c kzyp = maximum(kyzp,kyp)
c jblok/kblok,mblok/lblok = number of data blocks in x/y/z
c jkblok = maximum(jblok,kblok)
c mlblok = maximum(mblok,lblok)
c mixup = array of bit reversed addresses
c sct = sine/cosine table
c nxhyzd = maximum of (nx/2,ny,nz)
c nxyzhd = one half of maximum of (nx,ny,nz)
c the real data is stored in a complex array of length nx/2, ny, nz
c with the odd/even x points stored in the real/imaginary parts.
c in complex notation, fourier coefficients are stored as follows:
c if ntpose = 0,
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
c if ntpose = 1,
c h(l,j,k,i) = mode jj-1,kk-1,l-1, where jj = j + kxyp*(ix - 1),
c kk = k + kyzp*(iz - 1), and i = iy + nprocy*(iz - 1)
c 1 <= jj <= nx/2, 1 <= kk <= ny, and 1 <= l <= nz, except for
c h(l,1,k,1) = mode nx/2,k-1,l-1, where
c ny/2+2 <= kk <= ny, 1 <= l <= nz, and
c h(l,1,1,1) = mode nx/2,0,l-1
c h(l,1,ny/2+1,1) = mode nx/2,ny/2,l-1, where nz/2+2 <= l <= nz, and
c imaginary part of h(1,1,1,1) = real part of mode nx/2,0,0
c imaginary part of h(1,1,ny/2+1,1) = real part of mode nx/2,ny/2,0
c imaginary part of h(nz/2+1,1,1,1) =
c    real part of mode nx/2,0,nz/2
c imaginary part of h(nz/2+1,1,ny/2+1,1) =
c    real part of mode nx/2,ny/2,nz/2
c using jpl storage convention
c written by viktor k. decyk, ucla
c parallel, RISC optimized version
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd
      integer jblok, kblok, lblok, mblok, nxhyzd, nxyzhd
      integer mixup
      complex f, g, h, sct
      dimension f(3,nxvh,kypd,kzpd,kblok*lblok)
      dimension g(3,nyv,kxypd,kzpd,jblok*lblok)
      dimension h(3,nzv,kxypd,kyzpd,jblok*mblok)
      dimension mixup(nxhyzd), sct(nxyzhd)
c local data
      integer indx1, ndx1yz, nx, nxh, nxhh, nxh2, ny, nyh, ny2, nz, nzh
      integer nz2, nxyz, nxhyz, j, k, lb, ll, jb, it, nrx, nry, nrz
      integer nxyzh, l, i, m, n, ns, ns2, km, kmr, l1, k1, k2, j1, j2
      integer js, ks, kyb, kzb, kxb, kyzb, klblok, jlblok, jmblok
      integer mx, my, mz, moff, jj
      real dnxyz, arg, ani, at1, at2
      complex s, t, t1, t2, t3
      indx1 = indx - 1
      ndx1yz = max0(indx1,indy,indz)
      nx = 2**indx
      nxh = nx/2
      nxhh = nx/4
      nxh2 = nxh + 2
      ny = 2**indy
      nyh = ny/2
      ny2 = ny + 2
      nz = 2**indz
      nzh = nz/2
      nz2 = nz + 2
      nxyz = max0(nx,ny,nz)
      nxhyz = 2**ndx1yz
      if (isign) 50, 10, 530
c prepare fft tables
c bit-reverse index table: mixup(j) = 1 + reversed bits of (j - 1)
   10 do 30 j = 1, nxhyz
      lb = j - 1
      ll = 0
      do 20 k = 1, ndx1yz
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
   50 kxb = nxh/kxyp
      kyb = ny/kyp
      kzb = nz/kzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      if (kstrt.gt.(kyb*kzb)) go to 230
      klblok = kblok*lblok
c swap complex components
      do 90 m = 1, klblok
      do 80 n = 1, kzp
      do 70 i = 1, kyp
      do 60 j = 1, nxh
      at1 = real(f(3,j,i,n,m))
      f(3,j,i,n,m) = cmplx(real(f(2,j,i,n,m)),aimag(f(3,j,i,n,m)))
      at2 = aimag(f(2,j,i,n,m))
      f(2,j,i,n,m) = cmplx(aimag(f(1,j,i,n,m)),at1)
      f(1,j,i,n,m) = cmplx(real(f(1,j,i,n,m)),at2)
   60 continue
   70 continue
   80 continue
   90 continue
      ani = 0.5/(float(nx)*float(ny)*float(nz))
      do 220 m = 1, klblok
      do 210 n = 1, kzp
      nrx = nxhyz/nxh
c bit-reverse array elements in x
      do 110 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 110
      do 100 i = 1, kyp
      t1 = f(1,j1,i,n,m)
      t2 = f(2,j1,i,n,m)
      t3 = f(3,j1,i,n,m)
      f(1,j1,i,n,m) = f(1,j,i,n,m)
      f(2,j1,i,n,m) = f(2,j,i,n,m)
      f(3,j1,i,n,m) = f(3,j,i,n,m)
      f(1,j,i,n,m) = t1
      f(2,j,i,n,m) = t2
      f(3,j,i,n,m) = t3
  100 continue
  110 continue
c first transform in x
      nrx = nxyz/nxh
      do 150 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 140 i = 1, kyp
      do 130 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 120 j = 1, ns
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
  120 continue
  130 continue
  140 continue
  150 continue
c unscramble coefficients and normalize
      kmr = nxyz/nx
      nry = nxhyz/ny
      do 180 k = 1, kyp
      do 170 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),-real(sct(1+kmr*(j-1))))
      do 160 jj = 1, 3
      t = conjg(f(jj,nxh2-j,k,n,m))
      s = f(jj,j,k,n,m) + t
      t = (f(jj,j,k,n,m) - t)*t1
      f(jj,j,k,n,m) = ani*(s + t)
      f(jj,nxh2-j,k,n,m) = ani*conjg(s - t)
  160 continue
  170 continue
  180 continue
      do 200 k = 1, kyp
      do 190 jj = 1, 3
      f(jj,nxhh+1,k,n,m) = 2.*ani*conjg(f(jj,nxhh+1,k,n,m))
      f(jj,1,k,n,m) = 2.*ani*cmplx(real(f(jj,1,k,n,m)) + aimag(f(jj,1,k,
     1n,m)),real(f(jj,1,k,n,m)) - aimag(f(jj,1,k,n,m)))
  190 continue
  200 continue
  210 continue
  220 continue
c transpose f array to g
  230 call P3TPOS3AX(f,g,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,kxypd,kyp
     1d,kzpd,jblok,kblok,lblok)
      if (kstrt.gt.(kxb*kzb)) go to 370
      jlblok = jblok*lblok
      do 310 m = 1, jlblok
      do 300 n = 1, kzp
      nry = nxhyz/ny
c bit-reverse array elements in y
      do 250 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 250
      do 240 i = 1, kxyp
      t1 = g(1,k1,i,n,m)
      t2 = g(2,k1,i,n,m)
      t3 = g(3,k1,i,n,m)
      g(1,k1,i,n,m) = g(1,k,i,n,m)
      g(2,k1,i,n,m) = g(2,k,i,n,m)
      g(3,k1,i,n,m) = g(3,k,i,n,m)
      g(1,k,i,n,m) = t1
      g(2,k,i,n,m) = t2
      g(3,k,i,n,m) = t3
  240 continue
  250 continue
c then transform in y
      nry = nxyz/ny
      do 290 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 280 i = 1, kxyp
      do 270 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 260 j = 1, ns
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
  260 continue
  270 continue
  280 continue
  290 continue
  300 continue
  310 continue
c unscramble modes kx = 0, nx/2
      do 360 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 350 mx = 1, jblok
      if ((mx+js).gt.0) go to 350
      m = mx + moff
      do 340 n = 1, kzp
      do 330 k = 2, nyh
      do 320 jj = 1, 3
      s = g(jj,ny2-k,1,n,m)
      g(jj,ny2-k,1,n,m) = .5*cmplx(aimag(g(jj,k,1,n,m) + s),real(g(jj,k,
     11,n,m) - s))
      g(jj,k,1,n,m) = .5*cmplx(real(g(jj,k,1,n,m) + s),aimag(g(jj,k,1,n,
     1m) - s))
  320 continue
  330 continue
  340 continue
  350 continue
  360 continue
c transpose g array to h
  370 call P3TPOS3BX(g,h,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxypd,kyz
     1pd,kzpd,jblok,mblok,lblok)
      kyzb = ny/kyzp
      if (kstrt.gt.(kxb*kyzb)) go to 520
      jmblok = jblok*mblok
      do 450 m = 1, jmblok
      do 440 n = 1, kyzp
      nrz = nxhyz/nz
c bit-reverse array elements in z
      do 390 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 390
      do 380 i = 1, kxyp
      t1 = h(1,l1,i,n,m)
      t2 = h(2,l1,i,n,m)
      t3 = h(3,l1,i,n,m)
      h(1,l1,i,n,m) = h(1,l,i,n,m)
      h(2,l1,i,n,m) = h(2,l,i,n,m)
      h(3,l1,i,n,m) = h(3,l,i,n,m)
      h(1,l,i,n,m) = t1
      h(2,l,i,n,m) = t2
      h(3,l,i,n,m) = t3
  380 continue
  390 continue
c finally transform in z
      nrz = nxyz/nz
      do 430 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 420 i = 1, kxyp
      do 410 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 400 j = 1, ns
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
  400 continue
  410 continue
  420 continue
  430 continue
  440 continue
  450 continue
c unscramble modes kx = 0, nx/2
      do 510 my = 1, mblok
      moff = jblok*(my - 1)
      do 500 mx = 1, jblok
      if ((mx+js).gt.0) go to 500
      m = mx + moff
      if ((my+ks).eq.0) then
         do 470 n = 2, nzh
         do 460 jj = 1, 3
         s = h(jj,nz2-n,1,1,m)
         h(jj,nz2-n,1,1,m) = .5*cmplx(aimag(h(jj,n,1,1,m) + s),real(h(jj
     1,n,1,1,m) - s))
         h(jj,n,1,1,m) = .5*cmplx(real(h(jj,n,1,1,m) + s),aimag(h(jj,n,1
     1,1,m) - s))
  460    continue
  470    continue
      endif
      kyzb = nyh/kyzp
      if ((my+ks).eq.kyzb) then
         k1 = nyh - kyzb*kyzp + 1
         do 490 n = 2, nzh
         do 480 jj = 1, 3
         s = h(jj,nz2-n,1,k1,m)
         h(jj,nz2-n,1,k1,m) = .5*cmplx(aimag(h(jj,n,1,k1,m) + s),real(h(
     1jj,n,1,k1,m) - s))
         h(jj,n,1,k1,m) = .5*cmplx(real(h(jj,n,1,k1,m) + s),aimag(h(jj,n
     1,1,k1,m) - s))
  480    continue
  490    continue
      endif
  500 continue
  510 continue
c transpose h array to f
  520 if (ntpose.eq.0) then
         call P3TPOS3BX(h,g,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kxypd,
     1kzpd,kyzpd,jblok,lblok,mblok)
         call P3TPOS3AX(g,f,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,kypd,k
     1xypd,kzpd,kblok,jblok,lblok)
      endif
      return
c forward fourier transform
c transpose f array to h
  530 if (ntpose.eq.0) then
         call P3TPOS3AX(f,g,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,kxypd,
     1kypd,kzpd,jblok,kblok,lblok)
         call P3TPOS3BX(g,h,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxypd,
     1kyzpd,kzpd,jblok,mblok,lblok)
      endif
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      if (kstrt.gt.(kxb*kyzb)) go to 680
      jmblok = jblok*mblok
c scramble modes kx = 0, nx/2
      do 590 my = 1, mblok
      moff = jblok*(my - 1)
      do 580 mx = 1, jblok
      if ((mx+js).gt.0) go to 580
      m = mx + moff
      if ((my+ks).eq.0) then
         do 550 n = 2, nzh
         do 540 jj = 1, 3
         s = cmplx(aimag(h(jj,nz2-n,1,1,m)),real(h(jj,nz2-n,1,1,m)))
         h(jj,nz2-n,1,1,m) = conjg(h(jj,n,1,1,m) - s)
         h(jj,n,1,1,m) = h(jj,n,1,1,m) + s
  540    continue
  550    continue
      endif
      kyzb = nyh/kyzp
      if ((my+ks).eq.kyzb) then
         k1 = nyh - kyzb*kyzp + 1
         do 570 n = 2, nzh
         do 560 jj = 1, 3
         s = cmplx(aimag(h(jj,nz2-n,1,k1,m)),real(h(jj,nz2-n,1,k1,m)))
         h(jj,nz2-n,1,k1,m) = conjg(h(jj,n,1,k1,m) - s)
         h(jj,n,1,k1,m) = h(jj,n,1,k1,m) + s
  560    continue
  570    continue
      endif
  580 continue
  590 continue
      do 670 m = 1, jmblok
      do 660 n = 1, kyzp
      nrz = nxhyz/nz
c bit-reverse array elements in z
      do 610 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 610
      do 600 i = 1, kxyp
      t1 = h(1,l1,i,n,m)
      t2 = h(2,l1,i,n,m)
      t3 = h(3,l1,i,n,m)
      h(1,l1,i,n,m) = h(1,l,i,n,m)
      h(2,l1,i,n,m) = h(2,l,i,n,m)
      h(3,l1,i,n,m) = h(3,l,i,n,m)
      h(1,l,i,n,m) = t1
      h(2,l,i,n,m) = t2
      h(3,l,i,n,m) = t3
  600 continue
  610 continue
c first transform in z
      nrz = nxyz/nz
      do 650 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 640 i = 1, kxyp
      do 630 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 620 j = 1, ns
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
  620 continue
  630 continue
  640 continue
  650 continue
  660 continue
  670 continue
c transpose h array to g
  680 call P3TPOS3BX(h,g,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kxypd,kzp
     1d,kyzpd,jblok,lblok,mblok)
      kzb = nz/kzp
      if (kstrt.gt.(kxb*kzb)) go to 820
      jlblok = jblok*lblok
      nry = nxhyz/ny
c scramble modes kx = 0, nx/2
      do 730 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 720 mx = 1, jblok
      if ((mx+js).gt.0) go to 720
      m = mx + moff
      do 710 n = 1, kzp
      do 700 k = 2, nyh
      do 690 jj = 1, 3
      s = cmplx(aimag(g(jj,ny2-k,1,n,m)),real(g(jj,ny2-k,1,n,m)))
      g(jj,ny2-k,1,n,m) = conjg(g(jj,k,1,n,m) - s)
      g(jj,k,1,n,m) = g(jj,k,1,n,m) + s
  690 continue
  700 continue
  710 continue
  720 continue
  730 continue
      do 810 m = 1, jlblok
      do 800 n = 1, kzp
      nry = nxhyz/ny
c bit-reverse array elements in y
      do 750 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 750
      do 740 i = 1, kxyp
      t1 = g(1,k1,i,n,m)
      t2 = g(2,k1,i,n,m)
      t3 = g(3,k1,i,n,m)
      g(1,k1,i,n,m) = g(1,k,i,n,m)
      g(2,k1,i,n,m) = g(2,k,i,n,m)
      g(3,k1,i,n,m) = g(3,k,i,n,m)
      g(1,k,i,n,m) = t1
      g(2,k,i,n,m) = t2
      g(3,k,i,n,m) = t3
  740 continue
  750 continue
c then transform in y
      nry = nxyz/ny
      do 790 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 780 i = 1, kxyp
      do 770 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 760 j = 1, ns
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
  760 continue
  770 continue
  780 continue
  790 continue
  800 continue
  810 continue
c transpose g array to f
  820 call P3TPOS3AX(g,f,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,kypd,kxyp
     1d,kzpd,kblok,jblok,lblok)
      kyb = ny/kyp
      if (kstrt.gt.(kyb*kzb)) return
      klblok = kblok*lblok
      do 950 m = 1, klblok
      do 940 n = 1, kzp
c scramble coefficients
      kmr = nxyz/nx
      do 850 k = 1, kyp
      do 840 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 830 jj = 1, 3
      t = conjg(f(jj,nxh2-j,k,n,m))
      s = f(jj,j,k,n,m) + t
      t = (f(jj,j,k,n,m) - t)*t1
      f(jj,j,k,n,m) = s + t
      f(jj,nxh2-j,k,n,m) = conjg(s - t)
  830 continue
  840 continue
  850 continue
      do 870 k = 1, kyp
      do 860 jj = 1, 3
      f(jj,nxhh+1,k,n,m) = 2.*conjg(f(jj,nxhh+1,k,n,m))
      f(jj,1,k,n,m) = cmplx(real(f(jj,1,k,n,m)) + aimag(f(jj,1,k,n,m)),r
     1eal(f(jj,1,k,n,m)) - aimag(f(jj,1,k,n,m)))
  860 continue
  870 continue
      nrx = nxhyz/nxh
c bit-reverse array elements in x
      do 890 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 890
      do 880 i = 1, kyp
      t1 = f(1,j1,i,n,m)
      t2 = f(2,j1,i,n,m)
      t3 = f(3,j1,i,n,m)
      f(1,j1,i,n,m) = f(1,j,i,n,m)
      f(2,j1,i,n,m) = f(2,j,i,n,m)
      f(3,j1,i,n,m) = f(3,j,i,n,m)
      f(1,j,i,n,m) = t1
      f(2,j,i,n,m) = t2
      f(3,j,i,n,m) = t3
  880 continue
  890 continue
c finally transform in x
      nrx = nxyz/nxh
      do 930 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 920 i = 1, kyp
      do 910 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 900 j = 1, ns
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
  900 continue
  910 continue
  920 continue
  930 continue
  940 continue
  950 continue
c swap complex components
      do 990 m = 1, klblok
      do 980 n = 1, kzp 
      do 970 i = 1, kyp
      do 960 j = 1, nxh
      at1 = real(f(3,j,i,n,m))
      f(3,j,i,n,m) = cmplx(aimag(f(2,j,i,n,m)),aimag(f(3,j,i,n,m)))
      at2 = real(f(2,j,i,n,m))
      f(2,j,i,n,m) = cmplx(at1,aimag(f(1,j,i,n,m)))
      f(1,j,i,n,m) = cmplx(real(f(1,j,i,n,m)),at2)
  960 continue
  970 continue
  980 continue
  990 continue
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
      subroutine WPFFT32RX(f,g,h,isign,ntpose,mixup,sct,ttp,indx,indy,in
     1dz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,jblo
     2k,kblok,lblok,mblok,nxhyzd,nxyzhd)
c wrapper function for real to complex fft
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd
      integer jblok, kblok, lblok, mblok, nxhyzd, nxyzhd
      integer mixup
      real ttp
      complex f, g, h, sct
      dimension f(nxvh,kypd,kzpd,kblok*lblok)
      dimension g(nyv,kxypd,kzpd,jblok*lblok)
      dimension h(nzv,kxypd,kyzpd,jblok*mblok)
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
         call PTPOS3AX(f,g,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,kxypd,k
     1ypd,kzpd,jblok,kblok,lblok)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PFFT32RXY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kxy
     1p,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,nxhyzd,nxyzhd)
c transpose g array to h
         call PWTIMERA(-1,tp,dtime)
         call PTPOS3BX(g,h,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxypd,k
     1yzpd,kzpd,jblok,mblok,lblok)
         call PWTIMERA(1,tp,dtime)
c perform z fft
         call PFFT32RXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kxy
     1p,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,nxhyzd,nxyzhd)
c transpose h array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PTPOS3BX(h,g,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kxyp
     1d,kzpd,kyzpd,jblok,lblok,mblok)
            call PTPOS3AX(g,f,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,kypd
     1,kxypd,kzpd,kblok,jblok,lblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to h
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PTPOS3AX(f,g,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,kxyp
     1d,kypd,kzpd,jblok,kblok,lblok)
            call PTPOS3BX(g,h,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxyp
     1d,kyzpd,kzpd,jblok,mblok,lblok)
            call PWTIMERA(1,tf,dtime)
         endif
c perform z fft
         call PFFT32RXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kxy
     1p,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,nxhyzd,nxyzhd)
c transpose h array to g
         call PWTIMERA(-1,tp,dtime)
         call PTPOS3BX(h,g,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kxypd,k
     1zpd,kyzpd,jblok,lblok,mblok)
         call PWTIMERA(1,tp,dtime)
c perform y fft
         call PFFT32RXY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kxy
     1p,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,nxhyzd,nxyzhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PTPOS3AX(g,f,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,kypd,kx
     1ypd,kzpd,kblok,jblok,lblok)
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
      subroutine WPFFT32RX3(f,g,h,isign,ntpose,mixup,sct,ttp,indx,indy,i
     1ndz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,jbl
     2ok,kblok,lblok,mblok,nxhyzd,nxyzhd)
c wrapper function for real to complex fft
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd
      integer jblok, kblok, lblok, mblok, nxhyzd, nxyzhd
      integer mixup
      real ttp
      complex f, g, h, sct
      dimension f(3,nxvh,kypd,kzpd,kblok*lblok)
      dimension g(3,nyv,kxypd,kzpd,jblok*lblok)
      dimension h(3,nzv,kxypd,kyzpd,jblok*mblok)
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
         call P3TPOS3AX(f,g,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,kxypd,
     1kypd,kzpd,jblok,kblok,lblok)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PFFT32R3XY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1yp,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,nxhyzd,nxyzhd)
c transpose g array to h
         call PWTIMERA(-1,tp,dtime)
         call P3TPOS3BX(g,h,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxypd,
     1kyzpd,kzpd,jblok,mblok,lblok)
         call PWTIMERA(1,tp,dtime)
c perform z fft
         call PFFT32R3XZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1yp,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,nxhyzd,nxyzhd)
c transpose h array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call P3TPOS3BX(h,g,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kxy
     1pd,kzpd,kyzpd,jblok,lblok,mblok)
            call P3TPOS3AX(g,f,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,kyp
     1d,kxypd,kzpd,kblok,jblok,lblok)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to h
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call P3TPOS3AX(f,g,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,kxy
     1pd,kypd,kzpd,jblok,kblok,lblok)
            call P3TPOS3BX(g,h,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxy
     1pd,kyzpd,kzpd,jblok,mblok,lblok)
            call PWTIMERA(-1,tf,dtime)
         endif
c perform z fft
         call PFFT32R3XZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1yp,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,nxhyzd,nxyzhd)
c transpose h array to g
         call PWTIMERA(-1,tp,dtime)
         call P3TPOS3BX(h,g,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kxypd,
     1kzpd,kyzpd,jblok,lblok,mblok)
         call PWTIMERA(1,tp,dtime)
c perform y fft
         call PFFT32R3XY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1yp,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,nxhyzd,nxyzhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call P3TPOS3AX(g,f,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,kypd,k
     1xypd,kzpd,kblok,jblok,lblok)
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
      subroutine WPFFT32RN(f,g,h,bs,br,ss,isign,ntpose,mixup,sct,ttp,ind
     1x,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,
     2kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,ndim,nxhyzd,nxyzhd
     3)
c wrapper function for real to complex fft
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd, kzyp
      integer jblok, kblok, jkblok, lblok, mblok, mlblok, ndim
      integer nxhyzd, nxyzhd
      integer mixup
      real ttp
      complex f, g, h, bs, br, ss, sct
      dimension f(ndim,nxvh,kypd,kzpd,kblok*lblok)
      dimension g(ndim,nyv,kxypd,kzpd,jblok*lblok)
      dimension h(ndim,nzv,kxypd,kyzpd,jblok*mblok)
      dimension bs(ndim,kxyp,kzyp,kzp,jkblok*lblok)
      dimension br(ndim,kxyp,kzyp,kzp,jblok*mlblok)
      dimension ss(ndim,nxvh)
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
         call PFFT32RNXX(f,ss,isign,mixup,sct,indx,indy,indz,kstrt,kypi,
     1kyp,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,ndim,nxhyzd,nxyzhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PNTPOS3A(f,g,bs,br,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,k
     1xypd,kypd,kzpd,jblok,kblok,lblok,ndim)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PFFT32RNXY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1yp,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,ndim,nxhyzd,nxyzhd)
c transpose g array to h
         call PWTIMERA(-1,tp,dtime)
         call PNTPOS3B(g,h,bs,br,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,k
     1xypd,kyzpd,kzpd,jblok,mblok,lblok,ndim)
         call PWTIMERA(1,tp,dtime)
c perform z fft
         call PFFT32RNXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1yp,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,ndim,nxhyzd,nxyzhd)
c transpose h array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PNTPOS3B(h,g,br,bs,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyz
     1p,kxypd,kzpd,kyzpd,jblok,lblok,mblok,ndim)
            call PNTPOS3A(g,f,br,bs,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kz
     1p,kypd,kxypd,kzpd,kblok,jblok,lblok,ndim)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to h
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PNTPOS3A(f,g,bs,br,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kz
     1p,kxypd,kypd,kzpd,jblok,kblok,lblok,ndim)
            call PNTPOS3B(g,h,bs,br,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kz
     1p,kxypd,kyzpd,kzpd,jblok,mblok,lblok,ndim)
            call PWTIMERA(1,tf,dtime)
         endif
c perform z fft
         call PFFT32RNXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1yp,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,ndim,nxhyzd,nxyzhd)
c transpose h array to g
         call PWTIMERA(-1,tp,dtime)
         call PNTPOS3B(h,g,br,bs,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,k
     1xypd,kzpd,kyzpd,jblok,lblok,mblok,ndim)
         call PWTIMERA(1,tp,dtime)
c perform y fft
         call PFFT32RNXY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1yp,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,ndim,nxhyzd,nxyzhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PNTPOS3A(g,f,br,bs,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,k
     1ypd,kxypd,kzpd,kblok,jblok,lblok,ndim)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PFFT32RNXX(f,ss,isign,mixup,sct,indx,indy,indz,kstrt,kypi,
     1kyp,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,ndim,nxhyzd,nxyzhd)
      endif
      ttp = ttp + tp
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine WPFFT32RXN(f,g,h,ss,isign,ntpose,mixup,sct,ttp,indx,ind
     1y,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,
     2jblok,kblok,lblok,mblok,ndim,nxhyzd,nxyzhd)
c wrapper function for real to complex vector fft
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd
      integer jblok, kblok, lblok, mblok, ndim, nxhyzd, nxyzhd
      integer mixup
      real ttp
      complex f, g, h, ss, sct
      dimension f(ndim,nxvh,kypd,kzpd,kblok*lblok)
      dimension g(ndim,nyv,kxypd,kzpd,jblok*lblok)
      dimension h(ndim,nzv,kxypd,kyzpd,jblok*mblok)
      dimension ss(ndim,nxvh)
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
         call PFFT32RNXX(f,ss,isign,mixup,sct,indx,indy,indz,kstrt,kypi,
     1kyp,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,ndim,nxhyzd,nxyzhd)
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PNTPOS3AX(f,g,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,kxypd,
     1kypd,kzpd,jblok,kblok,lblok,ndim)
         call PWTIMERA(1,ttp,dtime)
c perform y fft
         call PFFT32RNXY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1yp,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,ndim,nxhyzd,nxyzhd)
c transpose g array to h
         call PWTIMERA(-1,tp,dtime)
         call PNTPOS3BX(g,h,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxypd,
     1kyzpd,kzpd,jblok,mblok,lblok,ndim)
         call PWTIMERA(1,tp,dtime)
c perform z fft
         call PFFT32RNXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1yp,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,ndim,nxhyzd,nxyzhd)
c transpose h array to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PNTPOS3BX(h,g,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kxy
     1pd,kzpd,kyzpd,jblok,lblok,mblok,ndim)
            call PNTPOS3AX(g,f,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,kyp
     1d,kxypd,kzpd,kblok,jblok,lblok,ndim)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f array to h
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PNTPOS3AX(f,g,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,kxy
     1pd,kypd,kzpd,jblok,kblok,lblok,ndim)
            call PNTPOS3BX(g,h,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxy
     1pd,kyzpd,kzpd,jblok,mblok,lblok,ndim)
            call PWTIMERA(1,tf,dtime)
         endif
c perform z fft
         call PFFT32RNXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1yp,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,ndim,nxhyzd,nxyzhd)
c transpose h array to g
         call PWTIMERA(-1,tp,dtime)
         call PNTPOS3BX(h,g,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kxypd,
     1kzpd,kyzpd,jblok,lblok,mblok,ndim)
         call PWTIMERA(1,tp,dtime)
c perform y fft
         call PFFT32RNXY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1yp,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,ndim,nxhyzd,nxyzhd)
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PNTPOS3AX(g,f,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,kypd,k
     1xypd,kzpd,kblok,jblok,lblok,ndim)
         call PWTIMERA(1,ttp,dtime)
c perform x fft
         call PFFT32RNXX(f,ss,isign,mixup,sct,indx,indy,indz,kstrt,kypi,
     1kyp,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,ndim,nxhyzd,nxyzhd)
      endif
      ttp = ttp + tp
      if (ntpose.eq.0) ttp = ttp + tf
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
