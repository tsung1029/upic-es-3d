c-----------------------------------------------------------------------
      program ptest32
      implicit none
      integer indx, indy, indz, indnvpy, indnvpz, mshare, nx, ny, nz
      integer nxh, nyh, nzh, nvpy, nvpz, nvp, kyp, kzp, nypmx, nzpmx
      integer kxyp, kyzp, kzyp, kyzmx, kyb, kzb, kxb, jkmx, kxyb, kyzb
      integer jlmx, kzyb, kbmin, kblok, jbmin, jblok, lbmin, lblok
      integer mbmin, mblok, jkbmx, jkblok, mlbmx, mlblok, klblok
      integer iblok, nblok, inblok
      integer nmyz, nyz, nmx, nxyz, nmxh
      integer nxv, nxvh, nyv, nyvh, nzv, nxyzh, nxhyz
      parameter( indx =   3, indy =   3, indz =   3)
c     parameter( indx =   2, indy =   2, indz =   2)
c     parameter( indx =   6, indy =   3, indz =   4)
c     parameter( indx =   7, indy =   8, indz =   6)
c indnvpy/indnvpz = exponent determining number of real or virtual
c processors in y/z, indnvpy must be <= indy, indnvpz must be <= indz
c mshare = (0,1) = (no,yes) architecture is shared memory
      parameter( indnvpy =   1, indnvpz =   1, mshare =   0)
      parameter(nx=2**indx,ny=2**indy,nz=2**indz)
      parameter(nxh=nx/2,nyh=ny/2,nzh=nz/2)
      parameter(nvpy=2**indnvpy,nvpz=2**indnvpz,nvp=nvpy*nvpz)
      parameter(iblok=1+mshare*(nvpy-1),nblok=1+mshare*(nvpz-1))
      parameter(inblok=iblok*nblok)
      parameter(kyp=(ny-1)/nvpy+1,kzp=(nz-1)/nvpz+1)
      parameter(nypmx=kyp+3,nzpmx=kzp+3)
      parameter(kxyp=(nxh-1)/nvpy+1,kyzp=(ny-1)/nvpz+1)
      parameter(kyzmx=kyzp*(kyp/kyzp)+kyp*(kyzp/kyp))
      parameter(kzyp=kyzmx/(2-kyzp/kyzmx-kyp/kyzmx))
      parameter(kyb=ny/kyp,kxb=nxh/kxyp,kzb=nz/kzp,kyzb=ny/kyzp)
      parameter(jkmx=kxb*(kyb/kxb)+kyb*(kxb/kyb))
      parameter(kxyb=jkmx/(2-kxb/jkmx-kyb/jkmx))
      parameter(jlmx=kyzb*(kzb/kyzb)+kzb*(kyzb/kzb))
      parameter(kzyb=jlmx/(2-kyzb/jlmx-kzb/jlmx))
      parameter(kbmin=1+(1-mshare)*(kxyb/kxb-1))
      parameter(kblok=1+mshare*(ny/kyp-1))
      parameter(jbmin=1+(1-mshare)*(kxyb/kyb-1))
      parameter(jblok=1+mshare*(nxh/kxyp-1))
      parameter(jkbmx=jblok*(kblok/jblok)+kblok*(jblok/kblok))
      parameter(jkblok=jkbmx/(2-jblok/jkbmx-kblok/jkbmx))
      parameter(lbmin=1+(1-mshare)*(kzyb/kyzb-1))
      parameter(lblok=1+mshare*(nz/kzp-1))
      parameter(klblok=kblok*lblok)
      parameter(mbmin=1+(1-mshare)*(kzyb/kzb-1))
      parameter(mblok=1+mshare*(ny/kyzp-1))
      parameter(mlbmx=mblok*(lblok/mblok)+lblok*(mblok/lblok))
      parameter(mlblok=mlbmx/(2-mblok/mlbmx-lblok/mlbmx))
      parameter(nmyz=ny*(nz/ny)+nz*(ny/nz))
      parameter(nyz=nmyz/(2-ny/nmyz-nz/nmyz))
      parameter(nmx=nx*(nyz/nx)+nyz*(nx/nyz))
      parameter(nxyz=nmx/(2-nx/nmx-nyz/nmx))
      parameter(nmxh=nxh*(nyz/nxh)+nyz*(nxh/nyz))
      parameter(nxhyz=nmxh/(2-nxh/nmxh-nyz/nmxh))
      parameter(nxv=nx+2,nxvh=nxv/2,nyv=ny+2,nyvh=nyv/2,nzv=nz+2)
      parameter(nxyzh=nxyz/2)
c
      integer isign, ntpose, idproc, kstrt, ks, js, ls, ms, j, k, l, m
      integer mx, my, mz, ll, l1, kk, k1, j1, moff, loff, koff, joff, i
      integer mixup
      real f, g, h, gh, s, t
      real at1, time, epsmax, eps
      double precision ranorm
      complex sct, br, bs
      dimension f(nxv,nypmx,nzpmx*kbmin,kblok*lblok)
      dimension s(2*nyv,kxyp,nzpmx*jbmin*lbmin,jblok*lblok)
      dimension t(2*nzv,kxyp,kyzp*mbmin,jblok*mblok)
      dimension bs(3,kxyp,kzyp,kzp,jkblok*lblok)
      dimension br(3,kxyp,kzyp,kzp,jblok*mlblok)
      dimension g(nx,ny,nz), gh(ny,nx,nz), h(nz,nx,ny)
      dimension mixup(nxhyz), sct(nxyzh)
c
c     integer nt, nt2, modesz, modeszd
c     parameter(nt=1,nt2=2*nt,modesz=4,modeszd=2*modesz-1)
c     integer it, k2
c     complex pottg, pott
c     dimension pottg(nt,nxh,ny,2*modesz-1)
c     dimension pott(nt,modeszd,kxyp,kyzp,jblok*mblok)
c     complex cg, ct
c     dimension cg(nxh,ny,nz)
c     dimension ct(nzv,kxyp,kyzp*mbmin,jblok*mblok)
c     equivalence (g,cg), (t,ct)
c
      integer ngds, nxe, nye, nze, idps, idds, nyzp, noff
      integer ngx, ngy, ngz
      parameter(nxe=nx+4,nye=ny+3,nze=nz+3,idps=4,idds=2)
      parameter(ngds=3*((idds-1)/2+1))
      dimension nyzp(idds,inblok), noff(idds,inblok)
      real xj0, yj0, zj0, edges, qu, qug, sum1, sum2, sum3, scr, scs
      dimension edges(idps,inblok)
      dimension scs(3,nxe,nzpmx,2*ngds,kblok*lblok)
      dimension scr(3,nxe,nypmx,ngds,kblok*lblok)
      dimension qu(3,nxe,nypmx,nzpmx,kblok*lblok), qug(3,nxe,nye,nze)
c     dimension qu(nxe,nypmx,nzpmx,kblok*lblok), qug(nxe,nye,nze)

      integer indx1, indy1, indz1, nx2, ny2, nz2, nx2e, kyp2, kzp2
      integer k2blok, l2blok, kyb2, kzb2
      real q3, q3g, scb, scd, q2, q2g
      parameter(indx1=indx+1,indy1=indy+1,nx2=2*nx,ny2=2*ny,nx2e=2*nxe)
      parameter(indz1=indz+1,nz2=2*nz)
      parameter(kyp2=(ny2-1)/nvpy+1,k2blok=1+mshare*(ny2/kyp2-1))
      parameter(kzp2=(nz2-1)/nvpz+1,l2blok=1+mshare*(nz2/kzp2-1))
      parameter(kyb2=ny2/kyp2,kzb2=nz2/kzp2)
      dimension q3(3,nx2e,2*nypmx,kzp2,k2blok*l2blok)
      dimension q3g(3,nx2e,ny2,nz2)
c     dimension q3(nx2e,2*nypmx,kzp2,k2blok*l2blok)
c     dimension q3g(nx2e,ny2,nz2)
      dimension scb(3,nxe,nypmx,nzpmx,kblok*lblok)
      dimension scd(3,nxe,nzpmx,2,kblok*lblok)
      dimension q2(3,nx2e,2*nypmx,kzp,k2blok*lblok)
      dimension q2g(3,nx2e,ny2,nz)
c     dimension q2(nx2e,2*nypmx,kzp,k2blok*lblok)
c     dimension q2g(nx2e,ny2,nz)
c
      real sumg, workg
      dimension sumg(3), workg(3)
c
c ntpose = (0,1) = (no,yes) input, output data are transposed in PFFT32R
      data ntpose /1/
c initialize for parallel processing
      call ppinit(idproc,nvp)
      kstrt = idproc + 1
      ks = (kstrt - 1)/kyb
      js = kstrt - kyb*ks - 2
      ks = ks - 1
      ls = (kstrt - 1)/kxb
      ms = kstrt - kxb*ls - 2
      ls = ls - 1
c prepare fft tables
      isign = 0
      call PFFT32R(f,s,t,bs,br,isign,ntpose,mixup,sct,indx,indy,indz,kst
     1rt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxyp,nypmx,kyzp,nzpmx,kzyp,jblok
     2,kblok,jkblok,lblok,mblok,mlblok,nxhyz,nxyzh)
c     call PFFT32RX(f,s,t,isign,ntpose,mixup,sct,indx,indy,indz,kstrt,nx
c    1vh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxyp,nypmx,kyzp,nzpmx,kzyp,jblok,kblo
c    2k,lblok,mblok,nxhyz,nxyzh)
c
      call DCOMP32(edges,nyzp,noff,ny,nz,kstrt,nvpy,nvpz,idps,idds,iblok
     1,nblok)
c create test function
      qug = -9.0
      qu = -9.0
      do 50 l = 1, nz
      ll = (l - 1)/kzp
      l1 = l - kzp*ll
      do 40 k = 1, ny
      kk = (k - 1)/kyp
      k1 = k - kyp*kk
      do 30 j = 1, nx
      at1 = ranorm()
      g(j,k,l) = at1
      gh(k,j,l) = g(j,k,l)
      h(l,j,k) = g(j,k,l)
c
c     qug(j+1,k+1,l+1) = ranorm()
      qug(1,j+1,k+1,l+1) = ranorm()
      qug(2,j+1,k+1,l+1) = ranorm()
      qug(3,j+1,k+1,l+1) = ranorm()
c
      do 20 mz = 1, lblok
      moff = kblok*(mz - 1)
      do 10 my = 1, kblok
      m = my + moff
      if ((kk.eq.(my+js)).and.(ll.eq.(mz+ks))) f(j,k1,l1,m) = g(j,k,l)
c
      if ((kk.eq.(my+js)).and.(ll.eq.(mz+ks))) then
c        qu(j+1,k1+1,l1+1,m) = qug(j+1,k+1,l+1)
         qu(1,j+1,k1+1,l1+1,m) = qug(1,j+1,k+1,l+1)
         qu(2,j+1,k1+1,l1+1,m) = qug(2,j+1,k+1,l+1)
         qu(3,j+1,k1+1,l1+1,m) = qug(3,j+1,k+1,l+1)
      endif
c
   10 continue
   20 continue
   30 continue
   40 continue
   50 continue
c
c begin test
      q3 = 9.0; q3g = -9.0
      q2 = 9.0; q2g = -9.0
c-----------------------------------------------------------------------
c     call TRPSIN3C(qug(1,2,2,2),q3g,nx,ny,nz,nxe,nye,nze,nx2e,ny2,nz2)
c     call TRPSIN3D(qug(2,2,2),q3g,nx,ny,nz,nxe,nye,nze,nx2e,ny2,nz2)
c     call PTRPSIN32C(qu(1,2,2,2,1),q3,scb,scd,nx,ny,nz,kstrt,nvpy,nvpz,
c    1nxe,kyp,kzp,nypmx,nzpmx,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
c     call PTRPSIN32D(qu(2,2,2,1),q3,scb,scd,nx,ny,nz,kstrt,nvpy,nvpz,nx
c    1e,kyp,kzp,nypmx,nzpmx,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
c-----------------------------------------------------------------------
      call DBLSIN3C(qug(1,2,2,2),q2g,nx,ny,nz,nxe,nye,nze,nx2e,ny2)
c     call DBLSIN3D(qug(2,2,2),q2g,nx,ny,nz,nxe,nye,nze,nx2e,ny2)
      call PDBLSIN32C(qu(1,2,2,2,1),q2,scb,scd,nx,ny,kstrt,nvpy,nxe,kyp,
     1kzp,nypmx,kzp,kyp2,kblok,lblok,k2blok)
c     call PDBLSIN32D(qu(2,2,2,1),q2,scb,scd,nx,ny,kstrt,nvpy,nxe,kyp,kz
c    1p,nypmx,kzp,kyp2,kblok,lblok,k2blok)
c-----------------------------------------------------------------------
      qug = 99.0
      qu = -99.0
c     do k = 1, ny+1
c     do j = 1, nx+1
c     do i = 1, 3
c     qug(i,j+1,k+1,nz+2) = real(k)
c     q3g(i,j,k,nz+1) = real(k)
c     q2g(i,j,k,nz+1) = real(k)
c     enddo
c     enddo
c     enddo
      do l = 1, nz
      do j = 1, nx+1
      do i = 1, 3
c     qug(i,j+1,ny+2,l+1) = -real(l)
c     q3g(i,j,ny+1,l) = -real(l)
      q2g(i,j,ny+1,l) = -real(l)
      enddo
      enddo
      enddo
      do mz = 1, lblok
      moff = kblok*(mz - 1)
      loff = kzp*(mz + ks)
c     do mz = 1, l2blok
c     moff = k2blok*(mz - 1)
c     loff = kzp2*(mz + ks)
      do my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
c
c     do l = 1, kzp2
c     l1 = l + loff
c     if (l1.eq.(nz+1)) then
c        do k = 1, kyp2
c        k1 = k + koff
c        if (k1.le.(ny+1)) then
c           do j = 1, nx+1
c           do i = 1, 3
c           q3(i,j,k,l,m) = real(k1)
c           enddo
c           enddo
c        endif
c        enddo
c     endif
c     enddo
c
      do l = 1, kzp
c     do l = 1, kzp2
      l1 = l + loff
      if (l1.lt.(nz+1)) then
         do k = 1, kyp2
         k1 = k + koff
         if (k1.eq.(ny+1)) then
            do j = 1, nx+1
            do i = 1, 3
c           q3(i,j,k,l,m) = -real(l1)
            q2(i,j,k,l,m) = -real(l1)
            enddo
            enddo
         endif
         enddo
      endif
      enddo
c
      enddo
      enddo
c
c-----------------------------------------------------------------------
c     call HAFTRP3C(qug(1,2,2,2),q3g,nx,ny,nz,nxe,nye,nze,nx2e,ny2,nz2)
c     call HAFTRP3D(qug(2,2,2),q3g,nx,ny,nz,nxe,nye,nze,nx2e,ny2,nz2)
c     call PHAFTRP32C(qu(1,2,2,2,1),q3,scb,scd,nx,ny,nz,kstrt,nvpy,nvpz,
c    1nxe,kyp,kzp,nypmx,nzpmx,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
c     call PHAFTRP32D(qu(2,2,2,1),q3,scb,scd,nx,ny,nz,kstrt,nvpy,nvpz,nx
c    1e,kyp,kzp,nypmx,nzpmx,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
c-----------------------------------------------------------------------
      call HAFDBL3C(qug(1,2,2,2),q2g,nx,ny,nz,nxe,nye,nze,nx2e,ny2)
c     call HAFDBL3D(qug(2,2,2),q2g,nx,ny,nz,nxe,nye,nze,nx2e,ny2)
      call SET_MON(2)
      call PHAFDBL32C(qu(1,2,2,2,1),q2,scb,scd,nx,ny,kstrt,nvpy,nvpz,nxe
     1,kyp,kzp,nypmx,kzp,kyp2,kblok,lblok,k2blok)
c     call PHAFDBL32D(qu(2,2,2,1),q2,scb,scd,nx,ny,kstrt,nvpy,nvpz,nxe,k
c    1yp,kzp,nypmx,kzp,kyp2,kblok,lblok,k2blok)
      call SET_MON(1)
c-----------------------------------------------------------------------
      epsmax = 0.
c     if (kstrt.gt.(kyb*kzb)) go to 120
      do 110 mz = 1, lblok
      moff = kblok*(mz - 1)
      loff = kzp*(mz + ks)
      do 100 my = 1, kblok
      koff = kyp*(my + js)
      m = my + moff
      do 90 l = 1, kzp+1
      l1 = l + loff
      do 80 k = 1, kyp+1
      k1 = k + koff
      do 70 j = 1, nx+1
      do 60 i = 1, 3
      eps = abs(qu(i,j+1,k+1,l+1,m) - qug(i,j+1,k1+1,l1+1))
      if (eps.gt.epsmax) then
         write (71,*) i,j,k,l,':',k1,l1,qu(i,j+1,k+1,l+1,m),
     1qug(i,j+1,k1+1,l1+1),eps
         epsmax = eps
      endif
   60 continue
   70 continue
   80 continue
   90 continue
  100 continue
  110 continue
  120 continue
      write (71,*) 'local epsmax=',epsmax
      call psum(epsmax,eps,1,1)
      write (71,*) 'global epsmax=',epsmax
      if (kstrt.ge.0) then
         call ppexit
         stop
      endif
c-----------------------------------------------------------------------
      epsmax = 0.
c     if (kstrt.gt.nz) go to 100
c     if (kstrt.gt.(kyb*kzb)) go to 120
      do 200 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      loff = kzp2*(mz + ks)
      do 190 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      do 180 l = 1, kzp2
      l1 = l + loff
c     if (l1.gt.(nz+1)) go to 180
      do 170 k = 1, kyp2
      k1 = k + koff
c     if (k1.gt.(ny+1)) go to 170
      do 160 j = 1, nx2
c     if (j.gt.nx) go to 160
c     write (74,*) j,k,l,':',k1,l1,q3(j,k,l,m)
c     eps = abs(q3(j,k,l,m) - q3g(j,k1,l1))
      do 150 i = 1, 3
      eps = abs(q3(i,j,k,l,m) - q3g(i,j,k1,l1))
      if (eps.gt.epsmax) then
c        write (71,*) j,k,l,':',k1,l1,q3(j,k,l,m),q3g(j,k1,l1),eps
         write (71,*) i,j,k,l,':',k1,l1,q3(i,j,k,l,m),q3g(i,j,k1,l1),eps
         epsmax = eps
c        write (71,*) qu(j+1,k+1,l+1,m)
      endif
  150 continue
  160 continue
  170 continue
  180 continue
  190 continue
  200 continue
      write (71,*) 'local epsmax=',epsmax
      call psum(epsmax,eps,1,1)
      write (71,*) 'global epsmax=',epsmax
      if (kstrt.ge.0) then
         call ppexit
         stop
      endif
c-----------------------------------------------------------------------
      epsmax = 0.
      do 300 mz = 1, lblok
      moff = k2blok*(mz - 1)
      loff = kzp*(mz + ks)
      do 290 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      do 280 l = 1, kzp
      l1 = l + loff
      do 270 k = 1, kyp2
      k1 = k + koff
c     if (k1.gt.(ny+1)) go to 270
      do 260 j = 1, nx2
c     if (j.gt.nx) go to 260
c     write (74,*) j,k,l,':',k1,l1,q2(j,k,l,m)
c     eps = abs(q2(j,k,l,m) - q2g(j,k1,l1))
      do 250 i = 1, 3
c     write (74,*) i,j,k,l,':',i,k1,l1,q2(i,j,k,l,m)
      eps = abs(q2(i,j,k,l,m) - q2g(i,j,k1,l1))
      if (eps.gt.epsmax) then
c        write (71,*) j,k,l,':',k1,l1,q2(j,k,l,m),q2g(j,k1,l1),eps
         write (71,*) i,j,k,l,':',k1,l1,q2(i,j,k,l,m),q2g(i,j,k1,l1),eps
         epsmax = eps
      endif
  250 continue
  260 continue
  270 continue
  280 continue
  290 continue
  300 continue
      write (71,*) 'local epsmax=',epsmax
      call psum(epsmax,eps,1,1)
      write (71,*) 'global epsmax=',epsmax
      if (kstrt.ge.0) then
         call ppexit
         stop
      endif
c-----------------------------------------------------------------------
c     call PTPOS3A(f,s,bs,br,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,kxyp,
c    1nypmx,nzpmx,jblok,kblok,lblok)
c     call PTPOS3B(s,t,bs,br,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxyp,
c    1kyzp,nzpmx,jblok,mblok,lblok)
c     dimension pottg(nt2,nxh,ny,2*modesz-1)
c     dimension pott(nt,modeszd,kxyp,kyzp,jblok*mblok)
c     dimension t(2*nzv,kxyp,kyzp*mbmin,jblok*mblok)
c     it = 1
c     pottg = cmplx(999.,999.)
c     pott = cmplx(999.,999.)
c-----------------------------------------------------------------------
c     call GTMODES3(g,pottg,nx,ny,nz,it,nxh,nyh,modesz,nx,ny,nz,nt2,nxh,
c    1ny,modeszd)
c     call PGTMODES33(t,pott,nx,ny,nz,it,kstrt,modesz,nzv,kxyp,kyzp,jblo
c    1k,mblok,nt,modeszd)
c     call PTMODES3(g,pottg,nx,ny,nz,it,nxh,nyh,modesz,nx,ny,nz,nt2,nxh,
c    1ny,modeszd)
c     call PPTMODES33(t,pott,nx,ny,nz,it,kstrt,modesz,nzv,kxyp,kyzp,jblo
c    1k,mblok,nt,modeszd)
c
c     epsmax = 0.
c     if (kstrt.gt.(kyb*kzb)) go to 120
c     do 110 mz = 1, lblok
c     moff = kblok*(mz - 1)
c     loff = kzp*(mz + ks)
c     do 100 my = 1, kblok
c     koff = kyp*(my + js)
c     m = my + moff
c     do 90 l = 1, kzp
c     l1 = l + loff
c     do 80 k = 1, kyp
c     k1 = k + koff
c     do 70 j = 1, nx
c     eps = abs(f(j,k,l,m) - g(j,k1,l1))
c     if (eps.gt.epsmax) then
c        write (61,*) 'f',j,k,l,m,f(j,k,l,m)
c        write (61,*) 'g',j,k1,l1,g(j,k1,l1)
c        write (61,*) 'eps=',eps
c        epsmax = eps
c     endif
c  70 continue
c  80 continue
c  90 continue
c 100 continue
c 110 continue
c 120 continue
c     do 230 my = 1, mblok
c     moff = jblok*(my - 1)
c     koff = kyzp*(my + ls)
c     do 220 mx = 1, jblok
c     joff = kxyp*(mx + ms)
c     m = mx + moff
c     do 210 k = 1, kyzp
c     k1 = k + koff
c     do 200 j = 1, kxyp
c     j1 = j + joff
c     do 190 l = 1, nz
c     eps = abs(ct(l,j,k,m) - cg(j1,k1,l))
c     if (eps.gt.epsmax) then
c        write (71,*) l,j,j1,k,k1,ct(l,j,k,m),cg(j1,k1,l),eps
c        epsmax = eps
c     endif
c 190 continue
c 200 continue
c 210 continue
c 220 continue
c 230 continue
c 240 continue
c     write (71,*) 'local epsmax=',epsmax
c     call psum (epsmax,eps,1,1)
c     write (71,*) 'global epsmax=',epsmax
      call ppexit
      stop
      end
c-----------------------------------------------------------------------
      subroutine ppinit(idproc,nvp)
c this subroutine initializes parallel processing
c input: nvp, output: idproc
c idproc = processor id
c nvp = number of real or virtual processors requested
      implicit none
      integer idproc, nvp
c get definition of MPI constants
      include 'mpif.h'
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=8)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mreal = default datatype for reals
c mint = default datatype for integers
c mcplx = default datatype for complex type
      common /pparms/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer ierror, ndprec
      save /pparms/
c ndprec = (0,1) = (no,yes) use (normal,autodouble) precision
      data ndprec /1/
c this segment is used for shared memory computers
c     nproc = nvp
c     idproc = 0
c this segment is used for mpi computers
      if (MPI_STATUS_SIZE.gt.lstat) then
         write (2,*) ' status size too small, actual/required = ', lstat
     1, MPI_STATUS_SIZE
         stop
      endif
c initialize the MPI execution environment
      call MPI_INIT(ierror)
      if (ierror.ne.0) stop
      lgrp = MPI_COMM_WORLD
c determine the rank of the calling process in the communicator
      call MPI_COMM_RANK(lgrp,idproc,ierror)
c determine the size of the group associated with a communicator
      call MPI_COMM_SIZE(lgrp,nproc,ierror)
c set default datatypes
         mint = MPI_INTEGER
c single precision
      if (ndprec.eq.0) then
         mreal = MPI_REAL
         mcplx = MPI_COMPLEX
c double precision
      else
         mreal = MPI_DOUBLE_PRECISION
         mcplx = MPI_DOUBLE_COMPLEX
      endif
c requested number of processors not obtained
      if (nproc.ne.nvp) then
         write (2,*) ' processor number error: nvp, nproc=', nvp, nproc
         call ppexit
         stop
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine ppexit
c this subroutine terminates parallel processing
      implicit none
c common block for parallel processing
      integer nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c lgrp = current communicator
      common /pparms/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
      integer ierror
c synchronize processes
      call MPI_BARRIER(lgrp,ierror)
c terminate MPI execution environment
      call MPI_FINALIZE(ierror)
      return
      end
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
      ani = 1./float(2*nx*ny*nz)
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
  480 continue
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
c for isign =/ 0 , input: all, output: f, g, h, bs, br
c approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c ntpose = (0,1) = (no,yes) input, output data are transposed
c if isign = 0, the fft tables are prepared
c if isign =/ 0, 3 forward fourier transforms are performed
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
      real dnxyz, arg, at1, at2
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
      if (isign.ne.0) go to 50
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
c forward fourier transform
c transpose f array to h
   50 if (ntpose.eq.0) then
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
      if (kstrt.gt.(kxb*kyzb)) go to 220
      jmblok = jblok*mblok
      nrz = nxhyz/nz
c scramble modes kx = 0, nx/2
      do 110 my = 1, mblok
      moff = jblok*(my - 1)
      do 100 mx = 1, jblok
      if ((mx+js).gt.0) go to 100
      m = mx + moff
      if ((my+ks).eq.0) then
         do 70 n = 2, nzh
         do 60 jj = 1, 3
         s = cmplx(aimag(h(jj,nz2-n,1,1,m)),real(h(jj,nz2-n,1,1,m)))
         h(jj,nz2-n,1,1,m) = conjg(h(jj,n,1,1,m) - s)
         h(jj,n,1,1,m) = h(jj,n,1,1,m) + s
   60    continue
   70    continue
      endif
      kyzb = nyh/kyzp
      if ((my+ks).eq.kyzb) then
         k1 = nyh - kyzb*kyzp + 1
         do 90 n = 2, nzh
         do 80 jj = 1, 3
         s = cmplx(aimag(h(jj,nz2-n,1,k1,m)),real(h(jj,nz2-n,1,k1,m)))
         h(jj,nz2-n,1,k1,m) = conjg(h(jj,n,1,k1,m) - s)
         h(jj,n,1,k1,m) = h(jj,n,1,k1,m) + s
   80    continue
   90    continue
      endif
  100 continue
  110 continue
      do 150 m = 1, jmblok
c bit-reverse array elements in z
      do 140 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 140
      do 130 n = 1, kyzp
      do 120 i = 1, kxyp
      t1 = h(1,l1,i,n,m)
      t2 = h(2,l1,i,n,m)
      t3 = h(3,l1,i,n,m)
      h(1,l1,i,n,m) = h(1,l,i,n,m)
      h(2,l1,i,n,m) = h(2,l,i,n,m)
      h(3,l1,i,n,m) = h(3,l,i,n,m)
      h(1,l,i,n,m) = t1
      h(2,l,i,n,m) = t2
      h(3,l,i,n,m) = t3
  120 continue
  130 continue
  140 continue
  150 continue
c first transform in z
      nrz = nxyz/nz
      do 210 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 200 m = 1, jmblok
      do 190 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 180 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 170 n = 1, kyzp
      do 160 i = 1, kxyp
      t1 = s*h(1,j2,i,n,m)
      t2 = s*h(2,j2,i,n,m)
      t3 = s*h(3,j2,i,n,m)
      h(1,j2,i,n,m) = h(1,j1,i,n,m) - t1
      h(2,j2,i,n,m) = h(2,j1,i,n,m) - t2
      h(3,j2,i,n,m) = h(3,j1,i,n,m) - t3
      h(1,j1,i,n,m) = h(1,j1,i,n,m) + t1
      h(2,j1,i,n,m) = h(2,j1,i,n,m) + t2
      h(3,j1,i,n,m) = h(3,j1,i,n,m) + t3
  160 continue
  170 continue
  180 continue
  190 continue
  200 continue
  210 continue
c transpose h array to g
  220 call P3TPOS3B(h,g,br,bs,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kxyp
     1d,kzpd,kyzpd,jblok,lblok,mblok)
      kzb = nz/kzp
      if (kstrt.gt.(kxb*kzb)) go to 380
      jlblok = jblok*lblok
      nry = nxhyz/ny
c scramble modes kx = 0, nx/2
      do 270 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 260 mx = 1, jblok
      if ((mx+js).gt.0) go to 260
      m = mx + moff
      do 250 n = 1, kzp
      do 240 k = 2, nyh
      do 230 jj = 1, 3
      s = cmplx(aimag(g(jj,ny2-k,1,n,m)),real(g(jj,ny2-k,1,n,m)))
      g(jj,ny2-k,1,n,m) = conjg(g(jj,k,1,n,m) - s)
      g(jj,k,1,n,m) = g(jj,k,1,n,m) + s
  230 continue
  240 continue
  250 continue
  260 continue
  270 continue
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
      s = conjg(sct(1+kmr*(j-1)))
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
c transpose g array to f
  380 call P3TPOS3A(g,f,br,bs,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,kypd
     1,kxypd,kzpd,kblok,jblok,lblok)
      kyb = ny/kyp
      if (kstrt.gt.(kyb*kzb)) return
      klblok = kblok*lblok
      nrx = nxhyz/nxh
c scramble coefficients
      kmr = nxyz/nx
      do 470 m = 1, klblok
      do 460 n = 1, kzp
      do 410 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 400 k = 1, kyp
      do 390 jj = 1, 3
      t = conjg(f(jj,nxh2-j,k,n,m))
      s = f(jj,j,k,n,m) + t
      t = (f(jj,j,k,n,m) - t)*t1
      f(jj,j,k,n,m) = s + t
      f(jj,nxh2-j,k,n,m) = conjg(s - t)
  390 continue
  400 continue
  410 continue
      do 430 k = 1, kyp
      do 420 jj = 1, 3
      f(jj,nxhh+1,k,n,m) = 2.*conjg(f(jj,nxhh+1,k,n,m))
      f(jj,1,k,n,m) = cmplx(real(f(jj,1,k,n,m)) + aimag(f(jj,1,k,n,m)),r
     1eal(f(jj,1,k,n,m)) - aimag(f(jj,1,k,n,m)))
  420 continue
  430 continue
c bit-reverse array elements in x
      do 450 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 450
      do 440 i = 1, kyp
      t1 = f(1,j1,i,n,m)
      t2 = f(2,j1,i,n,m)
      t3 = f(3,j1,i,n,m)
      f(1,j1,i,n,m) = f(1,j,i,n,m)
      f(2,j1,i,n,m) = f(2,j,i,n,m)
      f(3,j1,i,n,m) = f(3,j,i,n,m)
      f(1,j,i,n,m) = t1
      f(2,j,i,n,m) = t2
      f(3,j,i,n,m) = t3
  440 continue
  450 continue
  460 continue
  470 continue
c finally transform in x
      nrx = nxyz/nxh
      do 530 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 520 m = 1, klblok
      do 510 n = 1, kzp
      do 500 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 490 j = 1, ns
      j1 = j + k1
      j2 = j + k2
      s = conjg(sct(1+kmr*(j-1)))
      do 480 i = 1, kyp
      t1 = s*f(1,j2,i,n,m)
      t2 = s*f(2,j2,i,n,m)
      t3 = s*f(3,j2,i,n,m)
      f(1,j2,i,n,m) = f(1,j1,i,n,m) - t1
      f(2,j2,i,n,m) = f(2,j1,i,n,m) - t2
      f(3,j2,i,n,m) = f(3,j1,i,n,m) - t3
      f(1,j1,i,n,m) = f(1,j1,i,n,m) + t1
      f(2,j1,i,n,m) = f(2,j1,i,n,m) + t2
      f(3,j1,i,n,m) = f(3,j1,i,n,m) + t3
  480 continue
  490 continue
  500 continue
  510 continue
  520 continue
  530 continue
c swap complex components
      do 570 m = 1, klblok
      do 560 n = 1, kzp 
      do 550 i = 1, kyp
      do 540 j = 1, nxh
      at1 = real(f(3,j,i,n,m))
      f(3,j,i,n,m) = cmplx(aimag(f(2,j,i,n,m)),aimag(f(3,j,i,n,m)))
      at2 = real(f(2,j,i,n,m))
      f(2,j,i,n,m) = cmplx(at1,aimag(f(1,j,i,n,m)))
      f(1,j,i,n,m) = cmplx(real(f(1,j,i,n,m)),at2)
  540 continue
  550 continue
  560 continue
  570 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PFFT32RX(f,g,h,isign,ntpose,mixup,sct,indx,indy,indz,ks
     1trt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,kzyp,jblo
     2k,kblok,lblok,mblok,nxhyzd,nxyzhd)
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
c kzyp = maximum(kyzp,kyp)
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
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd, kzyp
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
      ani = 1./float(2*nx*ny*nz)
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
     1strt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,kzyp,jbl
     2ok,kblok,lblok,mblok,nxhyzd,nxyzhd)
c this subroutine performs 3 three dimensional real to complex fast
c fourier transforms, using complex arithmetic,
c for data which is distributed in blocks, with 2D spatial decomposition
c for isign = 0, input: isign, indx, indy, indz, kstrt, nxhyzd, nxyzhd
c output: mixup, sct
c for isign =/ 0 , input: all, output: f, g, h, bs, br
c approximate flop count: N*(5*log2(N) + 15/2)/nvp
c where N = (nx/2)*ny*nz, and nvp = number of procs
c indx/indy/indz = exponent which determines length in x/y/z direction,
c where nx=2**indx, ny=2**indy, nz=2**indz
c ntpose = (0,1) = (no,yes) input, output data are transposed
c if isign = 0, the fft tables are prepared
c if isign =/ 0, 3 forward fourier transforms are performed
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
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd, kzyp
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
      real dnxyz, arg, at1, at2
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
      if (isign.ne.0) go to 50
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
c forward fourier transform
c transpose f array to h
   50 if (ntpose.eq.0) then
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
      if (kstrt.gt.(kxb*kyzb)) go to 200
      jmblok = jblok*mblok
c scramble modes kx = 0, nx/2
      do 110 my = 1, mblok
      moff = jblok*(my - 1)
      do 100 mx = 1, jblok
      if ((mx+js).gt.0) go to 100
      m = mx + moff
      if ((my+ks).eq.0) then
         do 70 n = 2, nzh
         do 60 jj = 1, 3
         s = cmplx(aimag(h(jj,nz2-n,1,1,m)),real(h(jj,nz2-n,1,1,m)))
         h(jj,nz2-n,1,1,m) = conjg(h(jj,n,1,1,m) - s)
         h(jj,n,1,1,m) = h(jj,n,1,1,m) + s
   60    continue
   70    continue
      endif
      kyzb = nyh/kyzp
      if ((my+ks).eq.kyzb) then
         k1 = nyh - kyzb*kyzp + 1
         do 90 n = 2, nzh
         do 80 jj = 1, 3
         s = cmplx(aimag(h(jj,nz2-n,1,k1,m)),real(h(jj,nz2-n,1,k1,m)))
         h(jj,nz2-n,1,k1,m) = conjg(h(jj,n,1,k1,m) - s)
         h(jj,n,1,k1,m) = h(jj,n,1,k1,m) + s
   80    continue
   90    continue
      endif
  100 continue
  110 continue
      do 190 m = 1, jmblok
      do 180 n = 1, kyzp
      nrz = nxhyz/nz
c bit-reverse array elements in z
      do 130 l = 1, nz
      l1 = (mixup(l) - 1)/nrz + 1
      if (l.ge.l1) go to 130
      do 120 i = 1, kxyp
      t1 = h(1,l1,i,n,m)
      t2 = h(2,l1,i,n,m)
      t3 = h(3,l1,i,n,m)
      h(1,l1,i,n,m) = h(1,l,i,n,m)
      h(2,l1,i,n,m) = h(2,l,i,n,m)
      h(3,l1,i,n,m) = h(3,l,i,n,m)
      h(1,l,i,n,m) = t1
      h(2,l,i,n,m) = t2
      h(3,l,i,n,m) = t3
  120 continue
  130 continue
c first transform in z
      nrz = nxyz/nz
      do 170 l = 1, indz
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nzh/ns
      kmr = km*nrz
      do 160 i = 1, kxyp
      do 150 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 140 j = 1, ns
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
  140 continue
  150 continue
  160 continue
  170 continue
  180 continue
  190 continue
c transpose h array to g
  200 call P3TPOS3BX(h,g,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kxypd,kzp
     1d,kyzpd,jblok,lblok,mblok)
      kzb = nz/kzp
      if (kstrt.gt.(kxb*kzb)) go to 340
      jlblok = jblok*lblok
      nry = nxhyz/ny
c scramble modes kx = 0, nx/2
      do 250 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 240 mx = 1, jblok
      if ((mx+js).gt.0) go to 240
      m = mx + moff
      do 230 n = 1, kzp
      do 220 k = 2, nyh
      do 210 jj = 1, 3
      s = cmplx(aimag(g(jj,ny2-k,1,n,m)),real(g(jj,ny2-k,1,n,m)))
      g(jj,ny2-k,1,n,m) = conjg(g(jj,k,1,n,m) - s)
      g(jj,k,1,n,m) = g(jj,k,1,n,m) + s
  210 continue
  220 continue
  230 continue
  240 continue
  250 continue
      do 330 m = 1, jlblok
      do 320 n = 1, kzp
      nry = nxhyz/ny
c bit-reverse array elements in y
      do 270 k = 1, ny
      k1 = (mixup(k) - 1)/nry + 1
      if (k.ge.k1) go to 270
      do 260 i = 1, kxyp
      t1 = g(1,k1,i,n,m)
      t2 = g(2,k1,i,n,m)
      t3 = g(3,k1,i,n,m)
      g(1,k1,i,n,m) = g(1,k,i,n,m)
      g(2,k1,i,n,m) = g(2,k,i,n,m)
      g(3,k1,i,n,m) = g(3,k,i,n,m)
      g(1,k,i,n,m) = t1
      g(2,k,i,n,m) = t2
      g(3,k,i,n,m) = t3
  260 continue
  270 continue
c then transform in y
      nry = nxyz/ny
      do 310 l = 1, indy
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nyh/ns
      kmr = km*nry
      do 300 i = 1, kxyp
      do 290 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 280 j = 1, ns
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
  280 continue
  290 continue
  300 continue
  310 continue
  320 continue
  330 continue
c transpose g array to f
  340 call P3TPOS3AX(g,f,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,kypd,kxyp
     1d,kzpd,kblok,jblok,lblok)
      kyb = ny/kyp
      if (kstrt.gt.(kyb*kzb)) return
      klblok = kblok*lblok
      do 470 m = 1, klblok
      do 460 n = 1, kzp
c scramble coefficients
      kmr = nxyz/nx
      do 370 k = 1, kyp
      do 360 j = 2, nxhh
      t1 = cmplx(aimag(sct(1+kmr*(j-1))),real(sct(1+kmr*(j-1))))
      do 350 jj = 1, 3
      t = conjg(f(jj,nxh2-j,k,n,m))
      s = f(jj,j,k,n,m) + t
      t = (f(jj,j,k,n,m) - t)*t1
      f(jj,j,k,n,m) = s + t
      f(jj,nxh2-j,k,n,m) = conjg(s - t)
  350 continue
  360 continue
  370 continue
      do 390 k = 1, kyp
      do 380 jj = 1, 3
      f(jj,nxhh+1,k,n,m) = 2.*conjg(f(jj,nxhh+1,k,n,m))
      f(jj,1,k,n,m) = cmplx(real(f(jj,1,k,n,m)) + aimag(f(jj,1,k,n,m)),r
     1eal(f(jj,1,k,n,m)) - aimag(f(jj,1,k,n,m)))
  380 continue
  390 continue
      nrx = nxhyz/nxh
c bit-reverse array elements in x
      do 410 j = 1, nxh
      j1 = (mixup(j) - 1)/nrx + 1
      if (j.ge.j1) go to 410
      do 400 i = 1, kyp
      t1 = f(1,j1,i,n,m)
      t2 = f(2,j1,i,n,m)
      t3 = f(3,j1,i,n,m)
      f(1,j1,i,n,m) = f(1,j,i,n,m)
      f(2,j1,i,n,m) = f(2,j,i,n,m)
      f(3,j1,i,n,m) = f(3,j,i,n,m)
      f(1,j,i,n,m) = t1
      f(2,j,i,n,m) = t2
      f(3,j,i,n,m) = t3
  400 continue
  410 continue
c finally transform in x
      nrx = nxyz/nxh
      do 450 l = 1, indx1
      ns = 2**(l - 1)
      ns2 = ns + ns
      km = nxhh/ns
      kmr = km*nrx
      do 440 i = 1, kyp
      do 430 k = 1, km
      k1 = ns2*(k - 1)
      k2 = k1 + ns
      do 420 j = 1, ns
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
  420 continue
  430 continue
  440 continue
  450 continue
  460 continue
  470 continue
c swap complex components
      do 510 m = 1, klblok
      do 500 n = 1, kzp 
      do 490 i = 1, kyp
      do 480 j = 1, nxh
      at1 = real(f(3,j,i,n,m))
      f(3,j,i,n,m) = cmplx(aimag(f(2,j,i,n,m)),aimag(f(3,j,i,n,m)))
      at2 = real(f(2,j,i,n,m))
      f(2,j,i,n,m) = cmplx(at1,aimag(f(1,j,i,n,m)))
      f(1,j,i,n,m) = cmplx(real(f(1,j,i,n,m)),at2)
  480 continue
  490 continue
  500 continue
  510 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PTPOS3A(f,g,s,t,nx,ny,nz,kstrt,nxv,nyv,kxyp,kyp,kzp,kxy
     1pd,kypd,kzpd,jblok,kblok,lblok)
c this subroutine performs a transpose of a matrix f, distributed in y
c and z to a matrix g, distributed in x and z, that is,
c g(k+kyp*(m-1),j,l,mx,mz) = f(j+kxyp*(l-1),k,l,my,mz), where
c 1 <= j <= kxyp, 1 <= k <= kyp, 1 <= l <= kzp, and
c 1 <= mx <= nx/kxyp, 1 <= my <= ny/kyp, 1 <= mz <= nz/kzp
c and where indices mx, my, and mz can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = complex input array
c g = complex output array
c s, t = complex scratch arrays
c nx/ny/nz = number of points in x/y/z
c kstrt = starting data block number
c nxv/nyv = first dimension of f/g
c kypd/kxypd = second dimension of f/g
c kzpd = third dimension of f and g
c kxyp/kyp/kzp = number of data values per block in x/y/z
c jblok/kblok/lblok = number of data blocks in x/y/z
      implicit none
      integer nx, ny, nz, kstrt, nxv, nyv, kxyp, kyp, kzp
      integer kxypd, kypd, kzpd, jblok, kblok, lblok
      complex f, g, s, t
      dimension f(nxv,kypd,kzpd,kblok*lblok)
      dimension g(nyv,kxypd,kzpd,jblok*lblok)
      dimension s(kxyp,kyp,kzp,kblok*lblok), t(kxyp,kyp,kzp,jblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=8)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /pparms/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, ls, kxb, kyb, kzb
      integer jkblok, kxym, mtr, ntr, mntr
      integer m, mx, mz, l, i, moff, ioff, joff, koff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyb = ny/kyp
      kzb = nz/kzp
      ks = (kstrt - 1)/kyb
      js = kstrt - kyb*ks - 2
      ks = ks - 1
      ls = (kstrt - 1)/kxb - 1
c this segment is used for shared memory computers
c     do 60 mz = 1, lblok
c     moff = jblok*(mz - 1)
c     ioff = kblok*(mz - 1)
c     do 50 mx = 1, jblok
c     joff = kxyp*(mx + js)
c     m = mx + moff
c     do 40 i = 1, kyb
c     koff = kyp*(i - 1)
c     do 30 l = 1, kzp
c     do 20 k = 1, kyp
c     do 10 j = 1, kxyp
c     g(k+koff,j,l,m) = f(j+joff,k,l,i+ioff)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      do 100 mz = 1, lblok
      moff = jkblok*(mz - 1)
      do 90 mx = 1, jkblok
      m = mx + moff
      do 80 i = 1, kxym
      ir0 = iand(kxym-1,ieor(mx+js,i-1))
      is0 = ir0
      do 70 ii = 1, mntr
c post receive
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         koff = kyp*ir
         ir = ir + kyb*(mz + ls) + 1
         call MPI_IRECV(t(1,1,1,m),kxyp*kyp*kzp,mcplx,ir-1,ir+kxym+1,lgr
     1p,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kyb*kzb)).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         joff = kxyp*is
         is = is + kxb*(mz + ks) + 1
         do 30 l = 1, kzp
         do 20 k = 1, kyp
         do 10 j = 1, kxyp
         s(j,k,l,m) = f(j+joff,k,l,m)
   10    continue
   20    continue
   30    continue
         call MPI_SEND(s(1,1,1,m),kxyp*kyp*kzp,mcplx,is-1,m+kstrt+kxym,l
     1grp,ierr)
      endif
c receive data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
         do 60 l = 1, kzp
         do 50 k = 1, kyp
         do 40 j = 1, kxyp
         g(k+koff,j,l,m) = t(j,k,l,m)
   40    continue
   50    continue
   60    continue
      endif
   70 continue
   80 continue
   90 continue
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PTPOS3B(g,h,s,t,nx,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kx
     1ypd,kyzpd,kzpd,jblok,mblok,lblok)
c this subroutine performs a transpose of a matrix g, distributed in x
c and z to a matrix h, distributed in x and y, that is,
c h(j+kzp*(l-1),k,l,my,mz) = g(k+kyzp*(m-1),j,l,mx,mz), where
c 1 <= j <= kxyp, 1 <= k <= kyzp, 1 <= l <= kzp, and
c 1 <= mx <= nx/kxyp, 1 <= my <= ny/kyzp, 1 <= mz <= nz/kzp
c and where indices mx, my, and mz can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c g = complex input array
c h = complex output array
c s, t = complex scratch arrays
c nx/ny/nz = number of points in x/y/z
c kstrt = starting data block number
c nyv/nzv = first dimension of g/h
c kxypd = second dimension of f and g
c kzpd/kyzpd = third dimension of g/h
c kxyp/kyzp/kzp = number of data values per block in x/y/z
c jblok/mblok/lblok = number of data blocks in x/y/z
      implicit none
      integer nx, ny, nz, kstrt, nyv, nzv, kxyp, kyzp, kzp
      integer kxypd, kyzpd, kzpd, jblok, mblok, lblok
      complex g, h, s, t
      dimension g(nyv,kxypd,kzpd,jblok*lblok)
      dimension h(nzv,kxypd,kyzpd,jblok*mblok)
      dimension s(kyzp,kxyp,kzp,jblok*lblok)
      dimension t(kyzp,kxyp,kzp,jblok*mblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=8)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /pparms/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, kxb, kyzb, kzb
      integer mlblok, kyzm, mtr, ntr, mntr
      integer m, mx, my, l, i, moff, ioff, koff, loff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyzb = ny/kyzp
      kzb = nz/kzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
c this segment is used for shared memory computers
c     do 60 my = 1, mblok
c     moff = jblok*(my - 1)
c     koff = kyzp*(my + ks)
c     do 50 mx = 1, jblok
c     m = mx + moff
c     do 40 i = 1, kzb
c     loff = kzp*(i - 1)
c     ioff = jblok*(i - 1)
c     do 30 l = 1, kzp
c     do 20 j = 1, kxyp
c     do 10 k = 1, kyzp
c     h(l+loff,j,k,m) = g(k+koff,j,l,mx+ioff)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c this segment is used for mpi computers
      mlblok = max0(mblok,lblok)
      kyzm = min0(kyzb,kzb)
      mtr = kzb/kyzm
      ntr = kyzb/kyzm
      mntr = max0(mtr,ntr)
      do 100 my = 1, mlblok
      moff = jblok*(my - 1)
      do 90 mx = 1, jblok
      m = mx + moff
      do 80 i = 1, kyzm
      ir0 = iand(kyzm-1,ieor(my+ks,i-1))
      is0 = ir0
      do 70 ii = 1, mntr
c post receive
      if ((kstrt.le.(kxb*kyzb)).and.(ii.le.mtr)) then
         ir = ir0 + kyzm*(ii - 1)
         loff = kzp*ir
         ir = mx + js + kxb*ir + 1
         call MPI_IRECV(t(1,1,1,m),kxyp*kyzp*kzp,mcplx,ir-1,ir+kyzm+1,lg
     1rp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.ntr)) then
         is = is0 + kyzm*(ii - 1)
         koff = kyzp*is
         is = mx + js + kxb*is + 1
         do 30 l = 1, kzp
         do 20 j = 1, kxyp
         do 10 k = 1, kyzp
         s(k,j,l,m) = g(k+koff,j,l,m)
   10    continue
   20    continue
   30    continue
         call MPI_SEND(s(1,1,1,m),kxyp*kyzp*kzp,mcplx,is-1,m+kstrt+kyzm,
     1lgrp,ierr)
      endif
c receive data
      if ((kstrt.le.(kxb*kyzb)).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
         do 60 l = 1, kzp
         do 50 j = 1, kxyp
         do 40 k = 1, kyzp
         h(l+loff,j,k,m) = t(k,j,l,m)
   40    continue
   50    continue
   60    continue
      endif
   70 continue
   80 continue
   90 continue
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine P3TPOS3A(f,g,s,t,nx,ny,nz,kstrt,nxv,nyv,kxyp,kyp,kzp,kx
     1ypd,kypd,kzpd,jblok,kblok,lblok)
c this subroutine performs a transpose of a matrix f, distributed in y
c and z to a matrix g, distributed in x and z, that is,
c g(k+kyp*(m-1),j,l,mx,mz) = f(j+kxyp*(l-1),k,l,my,mz), where
c 1 <= j <= kxyp, 1 <= k <= kyp, 1 <= l <= kzp, and
c 1 <= mx <= nx/kxyp, 1 <= my <= ny/kyp, 1 <= mz <= nz/kzp
c and where indices mx, my, and mz can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c f = complex input array
c g = complex output array
c s, t = complex scratch arrays
c nx/ny/nz = number of points in x/y/z
c kstrt = starting data block number
c nxv/nyv = first dimension of f/g
c kypd/kxypd = second dimension of f/g
c kzpd = third dimension of f and g
c kxyp/kyp/kzp = number of data values per block in x/y/z
c jblok/kblok/lblok = number of data blocks in x/y/z
      implicit none
      integer nx, ny, nz, kstrt, nxv, nyv, kxyp, kyp, kzp
      integer kxypd, kypd, kzpd, jblok, kblok, lblok
      complex f, g, s, t
      dimension f(3,nxv,kypd,kzpd,kblok*lblok)
      dimension g(3,nyv,kxypd,kzpd,jblok*lblok)
      dimension s(3,kxyp,kyp,kzp,kblok*lblok)
      dimension t(3,kxyp,kyp,kzp,jblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=8)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /pparms/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, ls, kxb, kyb, kzb
      integer jkblok, kxym, mtr, ntr, mntr
      integer m, mx, mz, l, i, moff, ioff, joff, koff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyb = ny/kyp
      kzb = nz/kzp
      ks = (kstrt - 1)/kyb
      js = kstrt - kyb*ks - 2
      ks = ks - 1
      ls = (kstrt - 1)/kxb - 1
c this segment is used for shared memory computers
c     do 60 mz = 1, lblok
c     moff = jblok*(mz - 1)
c     ioff = kblok*(mz - 1)
c     do 50 mx = 1, jblok
c     joff = kxyp*(mx + js)
c     m = mx + moff
c     do 40 i = 1, kyb
c     koff = kyp*(i - 1)
c     do 30 l = 1, kzp
c     do 20 k = 1, kyp
c     do 10 j = 1, kxyp
c     g(1,k+koff,j,l,m) = f(1,j+joff,k,l,i+ioff)
c     g(2,k+koff,j,l,m) = f(2,j+joff,k,l,i+ioff)
c     g(3,k+koff,j,l,m) = f(3,j+joff,k,l,i+ioff)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
      do 100 mz = 1, lblok
      moff = jkblok*(mz - 1)
      do 90 mx = 1, jkblok
      m = mx + moff
      do 80 i = 1, kxym
      ir0 = iand(kxym-1,ieor(mx+js,i-1))
      is0 = ir0
      do 70 ii = 1, mntr
c post receive
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.mtr)) then
         ir = ir0 + kxym*(ii - 1)
         koff = kyp*ir
         ir = ir + kyb*(mz + ls) + 1
         call MPI_IRECV(t(1,1,1,1,m),3*kxyp*kyp*kzp,mcplx,ir-1,ir+kxym+1
     1,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kyb*kzb)).and.(ii.le.ntr)) then
         is = is0 + kxym*(ii - 1)
         joff = kxyp*is
         is = is + kxb*(mz + ks) + 1
         do 30 l = 1, kzp
         do 20 k = 1, kyp
         do 10 j = 1, kxyp
         s(1,j,k,l,m) = f(1,j+joff,k,l,m)
         s(2,j,k,l,m) = f(2,j+joff,k,l,m)
         s(3,j,k,l,m) = f(3,j+joff,k,l,m)
   10    continue
   20    continue
   30    continue
         call MPI_SEND(s(1,1,1,1,m),3*kxyp*kyp*kzp,mcplx,is-1,m+kstrt+kx
     1ym,lgrp,ierr)
      endif
c receive data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
         do 60 l = 1, kzp
         do 50 k = 1, kyp
         do 40 j = 1, kxyp
         g(1,k+koff,j,l,m) = t(1,j,k,l,m)
         g(2,k+koff,j,l,m) = t(2,j,k,l,m)
         g(3,k+koff,j,l,m) = t(3,j,k,l,m)
   40    continue
   50    continue
   60    continue
      endif
   70 continue
   80 continue
   90 continue
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine P3TPOS3B(g,h,s,t,nx,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,k
     1xypd,kyzpd,kzpd,jblok,mblok,lblok)
c this subroutine performs a transpose of a matrix g, distributed in x
c and z to a matrix h, distributed in x and y, that is,
c h(j+kzp*(l-1),k,l,my,mz) = g(k+kyzp*(m-1),j,l,mx,mz), where
c 1 <= j <= kxyp, 1 <= k <= kyzp, 1 <= l <= kzp, and
c 1 <= mx <= nx/kxyp, 1 <= my <= ny/kyzp, 1 <= mz <= nz/kzp
c and where indices mx, my, and mz can be distributed across processors.
c this subroutine sends and receives one message at a time, either
c synchronously or asynchronously. it uses a minimum of system resources
c g = complex input array
c h = complex output array
c s, t = complex scratch arrays
c nx/ny/nz = number of points in x/y/z
c kstrt = starting data block number
c nyv/nzv = first dimension of g/h
c kxypd = second dimension of f and g
c kzpd/kyzpd = third dimension of g/h
c kxyp/kyzp/kzp = number of data values per block in x/y/z
c jblok/mblok/lblok = number of data blocks in x/y/z
      implicit none
      integer nx, ny, nz, kstrt, nyv, nzv, kxyp, kyzp, kzp
      integer kxypd, kyzpd, kzpd, jblok, mblok, lblok
      complex g, h, s, t
      dimension g(3,nyv,kxypd,kzpd,jblok*lblok)
      dimension h(3,nzv,kxypd,kyzpd,jblok*mblok)
      dimension s(3,kyzp,kxyp,kzp,jblok*lblok)
      dimension t(3,kyzp,kxyp,kzp,jblok*mblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=8)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /pparms/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, kxb, kyzb, kzb
      integer mlblok, kyzm, mtr, ntr, mntr
      integer m, mx, my, l, i, moff, ioff, koff, loff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyzb = ny/kyzp
      kzb = nz/kzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
c this segment is used for shared memory computers
c     do 60 my = 1, mblok
c     moff = jblok*(my - 1)
c     koff = kyzp*(my + ks)
c     do 50 mx = 1, jblok
c     m = mx + moff
c     do 40 i = 1, kzb
c     loff = kzp*(i - 1)
c     ioff = jblok*(i - 1)
c     do 30 l = 1, kzp
c     do 20 j = 1, kxyp
c     do 10 k = 1, kyzp
c     h(1,l+loff,j,k,m) = g(1,k+koff,j,l,mx+ioff)
c     h(2,l+loff,j,k,m) = g(2,k+koff,j,l,mx+ioff)
c     h(3,l+loff,j,k,m) = g(3,k+koff,j,l,mx+ioff)
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c this segment is used for mpi computers
      mlblok = max0(mblok,lblok)
      kyzm = min0(kyzb,kzb)
      mtr = kzb/kyzm
      ntr = kyzb/kyzm
      mntr = max0(mtr,ntr)
      do 100 my = 1, mlblok
      moff = jblok*(my - 1)
      do 90 mx = 1, jblok
      m = mx + moff
      do 80 i = 1, kyzm
      ir0 = iand(kyzm-1,ieor(my+ks,i-1))
      is0 = ir0
      do 70 ii = 1, mntr
c post receive
      if ((kstrt.le.(kxb*kyzb)).and.(ii.le.mtr)) then
         ir = ir0 + kyzm*(ii - 1)
         loff = kzp*ir
         ir = mx + js + kxb*ir + 1
         call MPI_IRECV(t(1,1,1,1,m),3*kxyp*kyzp*kzp,mcplx,ir-1,ir+kyzm+
     11,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.ntr)) then
         is = is0 + kyzm*(ii - 1)
         koff = kyzp*is
         is = mx + js + kxb*is + 1
         do 30 l = 1, kzp
         do 20 j = 1, kxyp
         do 10 k = 1, kyzp
         s(1,k,j,l,m) = g(1,k+koff,j,l,m)
         s(2,k,j,l,m) = g(2,k+koff,j,l,m)
         s(3,k,j,l,m) = g(3,k+koff,j,l,m)
   10    continue
   20    continue
   30    continue
         call MPI_SEND(s(1,1,1,1,m),3*kxyp*kyzp*kzp,mcplx,is-1,m+kstrt+k
     1yzm,lgrp,ierr)
      endif
c receive data
      if ((kstrt.le.(kxb*kyzb)).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
         do 60 l = 1, kzp
         do 50 j = 1, kxyp
         do 40 k = 1, kyzp
         h(1,l+loff,j,k,m) = t(1,k,j,l,m)
         h(2,l+loff,j,k,m) = t(2,k,j,l,m)
         h(3,l+loff,j,k,m) = t(3,k,j,l,m)
   40    continue
   50    continue
   60    continue
      endif
   70 continue
   80 continue
   90 continue
  100 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PTPOS3AX(f,g,nx,ny,nz,kstrt,nxv,nyv,kxyp,kyp,kzp,kxypd,
     1kypd,kzpd,jblok,kblok,lblok)
c this subroutine performs a transpose of a matrix f, distributed in y
c and z to a matrix g, distributed in x and z, that is,
c g(k+kyp*(m-1),j,l,mx,mz) = f(j+kxyp*(l-1),k,l,my,mz), where
c 1 <= j <= kxyp, 1 <= k <= kyp, 1 <= l <= kzp, and
c 1 <= mx <= nx/kxyp, 1 <= my <= ny/kyp, 1 <= mz <= nz/kzp
c and where indices mx, my, and mz can be distributed across processors.
c this subroutine sends and receives multiple asynchronous messages
c f = complex input array
c g = complex output array
c nx/ny/nz = number of points in x/y/z
c kstrt = starting data block number
c nxv/nyv = first dimension of f/g
c kypd/kxypd = second dimension of f/g
c kzpd = third dimension of f and g
c kxyp/kyp/kzp = number of data values per block in x/y/z
c jblok/kblok/lblok = number of data blocks in x/y/z
      implicit none
      integer nx, ny, nz, kstrt, nxv, nyv, kxyp, kyp, kzp
      integer kxypd, kypd, kzpd, jblok, kblok, lblok
      complex f, g
      dimension f(nxv*kypd*kzpd*kblok*lblok)
      dimension g(nyv*kxypd*kzpd*jblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=8)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /pparms/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, ls, kxb, kyb, kzb
      integer jkblok, kxym, mtr, ntr, mntr
      integer m, mx, mz, l, i, moff, ioff, joff, koff, loff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyb = ny/kyp
      kzb = nz/kzp
      ks = (kstrt - 1)/kyb
      js = kstrt - kyb*ks - 2
      ks = ks - 1
      ls = (kstrt - 1)/kxb - 1
c this segment is used for shared memory computers
c     do 60 mz = 1, lblok
c     moff = jblok*(mz - 1)
c     ioff = kblok*(mz - 1) - 1
c     do 50 mx = 1, jblok
c     joff = kxyp*(mx + js)
c     m = mx + moff
c     do 40 i = 1, kyb
c     koff = kyp*(i - 1)
c     do 30 l = 1, kzp
c     do 20 k = 1, kyp
c     do 10 j = 1, kxyp
c     g(k+koff+nyv*(j-1+kxypd*(l-1+kzpd*(m-1)))) = f(j+joff+nxv*(k-1+kyp
c    1d*(l-1+kzpd*(i+ioff))))
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
c transpose local data
      do 70 mz = 1, lblok
      moff = kblok*(mz - 1)
      do 60 mx = 1, kblok
      m = mx + moff
      ioff = kxb*(m - 1)
      loff = kzpd*(m - 1) - 1
      do 50 i = 1, kxym
      is0 = iand(kxym-1,ieor(mx+js,i-1))
      do 40 ii = 1, ntr
      if (kstrt.le.(kyb*kzb)) then
         is = is0 + kxym*(ii - 1)
         joff = kxyp*is
         is = kzp*(is + ioff) - 1
         do 30 l = 1, kzp
         ir = kyp*(l + is) - 1
         koff = kypd*(l + loff) - 1
         do 20 k = 1, kyp
         do 10 j = 1, kxyp
         g(j+kxyp*(k+ir)) = f(j+joff+nxv*(k+koff))
   10    continue
   20    continue
   30    continue
      endif
   40 continue
   50 continue
   60 continue
   70 continue
c exchange data
      do 110 mz = 1, lblok
      moff = jkblok*(mz - 1)
      do 100 mx = 1, jkblok
      m = mx + moff
      do 90 i = 1, kxym
      ir0 = iand(kxym-1,ieor(mx+js,i-1))
      is0 = ir0
      do 80 ii = 1, mntr
c post receive
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.mtr)) then
         is = ir0 + kxym*(ii - 1)
         ir = is + kyb*(mz + ls) + 1
         call MPI_IRECV(f(1+kxyp*kyp*kzp*is),kxyp*kyp*kzp,mcplx,ir-1,ir+
     1kxym+1,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kyb*kzb)).and.(ii.le.ntr)) then
         ir = is0 + kxym*(ii - 1)
         is = ir + kxb*(mz + ks) + 1
         call MPI_SEND(g(1+kxyp*kyp*kzp*ir),kxyp*kyp*kzp,mcplx,is-1,m+ks
     1trt+kxym,lgrp,ierr)
      endif
c receive data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
   80 continue
   90 continue
  100 continue
  110 continue
c transpose local data
      do 180 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 170 mx = 1, jblok
      m = mx + moff
      ioff = kyb*(m - 1) 
      loff = kzpd*(m - 1) - 1
      do 160 i = 1, kxym
      ir0 = iand(kxym-1,ieor(mx+js,i-1))
      do 150 ii = 1, mtr
      if (kstrt.le.(kxb*kzb)) then
         ir = ir0 + kxym*(ii - 1)
         koff = kyp*ir
         ir = kzp*(ir + ioff) - 1
         do 140 l = 1, kzp
         is = kyp*(l + ir) - 1
         joff = kxypd*(l + loff) - 1
         do 130 k = 1, kyp
         do 120 j = 1, kxyp
         g(k+koff+nyv*(j+joff)) = f(j+kxyp*(k+is))
  120    continue
  130    continue
  140    continue
      endif
  150 continue
  160 continue
  170 continue
  180 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PTPOS3BX(g,h,nx,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxypd
     1,kyzpd,kzpd,jblok,mblok,lblok)
c this subroutine performs a transpose of a matrix g, distributed in x
c and z to a matrix h, distributed in x and y, that is,
c h(j+kzp*(l-1),k,l,my,mz) = g(k+kyzp*(m-1),j,l,mx,mz), where
c 1 <= j <= kxyp, 1 <= k <= kyzp, 1 <= l <= kzp, and
c 1 <= mx <= nx/kxyp, 1 <= my <= ny/kyzp, 1 <= mz <= nz/kzp
c and where indices mx, my, and mz can be distributed across processors.
c this subroutine sends and receives multiple asynchronous messages
c g = complex input array
c h = complex output array
c nx/ny/nz = number of points in x/y/z
c kstrt = starting data block number
c nyv/nzv = first dimension of g/h
c kxypd = second dimension of f and g
c kzpd/kyzpd = third dimension of g/h
c kxyp/kyzp/kzp = number of data values per block in x/y/z
c jblok/mblok/lblok = number of data blocks in x/y/z
      implicit none
      integer nx, ny, nz, kstrt, nyv, nzv, kxyp, kyzp, kzp
      integer kxypd, kyzpd, kzpd, jblok, mblok, lblok
      complex g, h
      dimension g(nyv*kxypd*kzpd*jblok*lblok)
      dimension h(nzv*kxypd*kyzpd*jblok*mblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=8)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /pparms/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, kxb, kyzb, kzb
      integer mlblok, kyzm, mtr, ntr, mntr, nzyv
      integer m, mx, my, l, i, moff, ioff, joff, koff, loff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyzb = ny/kyzp
      kzb = nz/kzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
c this segment is used for shared memory computers
c     do 60 my = 1, mblok
c     moff = jblok*(my - 1)
c     koff = kyzp*(my + ks)
c     do 50 mx = 1, jblok
c     m = mx + moff
c     do 40 i = 1, kzb
c     loff = kzp*(i - 1)
c     ioff = jblok*(i - 1) - 1
c     do 30 l = 1, kzp
c     do 20 j = 1, kxyp
c     do 10 k = 1, kyzp
c     h(l+loff+nzv*(j-1+kxypd*(k-1+kyzpd*(m-1)))) = g(k+koff+nyv*(j-1+kx
c    1ypd*(l-1+kzpd*(mx+ioff))))
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c this segment is used for mpi computers
      mlblok = max0(mblok,lblok)
      kyzm = min0(kyzb,kzb)
      mtr = kzb/kyzm
      ntr = kyzb/kyzm
      mntr = max0(mtr,ntr)
c transpose local data
      do 70 my = 1, lblok
      moff = jblok*(my - 1)
      do 60 mx = 1, jblok
      m = mx + moff
      ioff = kyzb*(m - 1)
      loff = kzpd*(m - 1) - 1
      do 50 i = 1, kyzm
      is0 = iand(kyzm-1,ieor(my+ks,i-1))
      do 40 ii = 1, ntr
      if (kstrt.le.(kxb*kzb)) then
         is = is0 + kyzm*(ii - 1)
         koff = kyzp*is
         is = kzp*(is + ioff) - 1
         do 30 l = 1, kzp
         ir = kxyp*(l + is) - 1
         joff = kxypd*(l + loff) - 1
         do 20 j = 1, kxyp
         do 10 k = 1, kyzp
         h(k+kyzp*(j+ir)) = g(k+koff+nyv*(j+joff))
   10    continue
   20    continue
   30    continue
      endif
   40 continue
   50 continue
   60 continue
   70 continue
c exchange data
      do 110 my = 1, mlblok
      moff = jblok*(my - 1)
      do 100 mx = 1, jblok
      m = mx + moff
      do 90 i = 1, kyzm
      ir0 = iand(kyzm-1,ieor(my+ks,i-1))
      is0 = ir0
      do 80 ii = 1, mntr
c post receive
      if ((kstrt.le.(kxb*kyzb)).and.(ii.le.mtr)) then
         is = ir0 + kyzm*(ii - 1)
         ir = mx + js + kxb*is + 1
         call MPI_IRECV(g(1+kxyp*kyzp*kzp*is),kxyp*kyzp*kzp,mcplx,ir-1,i
     1r+kyzm+1,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.ntr)) then
         ir = is0 + kyzm*(ii - 1)
         is = mx + js + kxb*ir + 1
         call MPI_SEND(h(1+kxyp*kyzp*kzp*ir),kxyp*kyzp*kzp,mcplx,is-1,m+
     1kstrt+kyzm,lgrp,ierr)
      endif
c receive data
      if ((kstrt.le.(kxb*kyzb)).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
   80 continue
   90 continue
  100 continue
  110 continue
c transpose local data
      nzyv = nzv*kxypd
      do 180 my = 1, mblok
      moff = jblok*(my - 1)
      do 170 mx = 1, jblok
      m = mx + moff
      ioff = kzb*(m - 1)
      koff = kyzpd*(m - 1) - 1
      do 160 i = 1, kyzm
      ir0 = iand(kyzm-1,ieor(my+ks,i-1))
      do 150 ii = 1, mtr
      if (kstrt.le.(kxb*kyzb)) then
         ir = ir0 + kyzm*(ii - 1)
         joff = kzp*ir
         ir = kzp*(ir + ioff) - 1
         do 140 l = 1, kzp
         is = kxyp*(l + ir) - 1
         do 130 j = 1, kxyp
         loff = nzv*(j - 1) + joff
         do 120 k = 1, kyzp
         h(l+loff+nzyv*(k+koff)) = g(k+kyzp*(j+is))
  120    continue
  130    continue
  140    continue
      endif
  150 continue
  160 continue
  170 continue
  180 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine P3TPOS3AX(f,g,nx,ny,nz,kstrt,nxv,nyv,kxyp,kyp,kzp,kxypd
     1,kypd,kzpd,jblok,kblok,lblok)
c this subroutine performs a transpose of a matrix f, distributed in y
c and z to a matrix g, distributed in x and z, that is,
c g(1:3,k+kyp*(m-1),j,l,mx,mz) = f(1:3,j+kxyp*(l-1),k,l,my,mz), where
c 1 <= j <= kxyp, 1 <= k <= kyp, 1 <= l <= kzp, and
c 1 <= mx <= nx/kxyp, 1 <= my <= ny/kyp, 1 <= mz <= nz/kzp
c and where indices mx, my, and mz can be distributed across processors.
c this subroutine sends and receives multiple asynchronous messages
c f = complex input array
c g = complex output array
c nx/ny/nz = number of points in x/y/z
c kstrt = starting data block number
c nxv/nyv = first dimension of f/g
c kypd/kxypd = second dimension of f/g
c kzpd = third dimension of f and g
c kxyp/kyp/kzp = number of data values per block in x/y/z
c jblok/kblok/lblok = number of data blocks in x/y/z
c optimized version
      implicit none
      integer nx, ny, nz, kstrt, nxv, nyv, kxyp, kyp, kzp
      integer kxypd, kypd, kzpd, jblok, kblok, lblok
      complex f, g
      dimension f(3*nxv*kypd*kzpd*kblok*lblok)
      dimension g(3*nyv*kxypd*kzpd*jblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=8)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /pparms/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, ls, kxb, kyb, kzb
      integer jkblok, kxym, mtr, ntr, mntr
      integer m, mx, mz, l, i, moff, ioff, joff, koff, loff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyb = ny/kyp
      kzb = nz/kzp
      ks = (kstrt - 1)/kyb
      js = kstrt - kyb*ks - 2
      ks = ks - 1
      ls = (kstrt - 1)/kxb - 1
c this segment is used for shared memory computers
c     do 60 mz = 1, lblok
c     moff = jblok*(mz - 1)
c     ioff = kblok*(mz - 1) - 1
c     do 50 mx = 1, jblok
c     joff = kxyp*(mx + js) - 1
c     m = mx + moff
c     do 40 i = 1, kyb
c     koff = kyp*(i - 1) - 1
c     do 30 l = 1, kzp
c     do 20 k = 1, kyp
c     do 10 j = 1, kxyp
c     g(1+3*(k+koff+nyv*(j-1+kxypd*(l-1+kzpd*(m-1))))) = f(1+3*(j+joff+n
c    1xv*(k-1+kypd*(l-1+kzpd*(i+ioff)))))
c     g(2+3*(k+koff+nyv*(j-1+kxypd*(l-1+kzpd*(m-1))))) = f(2+3*(j+joff+n
c    1xv*(k-1+kypd*(l-1+kzpd*(i+ioff)))))
c     g(3+3*(k+koff+nyv*(j-1+kxypd*(l-1+kzpd*(m-1))))) = f(3+3*(j+joff+n
c    1xv*(k-1+kypd*(l-1+kzpd*(i+ioff)))))
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c this segment is used for mpi computers
      jkblok = max0(jblok,kblok)
      kxym = min0(kxb,kyb)
      mtr = kyb/kxym
      ntr = kxb/kxym
      mntr = max0(mtr,ntr)
c transpose local data
      do 70 mz = 1, lblok
      moff = kblok*(mz - 1)
      do 60 mx = 1, kblok
      m = mx + moff
      ioff = kxb*(m - 1)
      loff = kzpd*(m - 1) - 1
      do 50 i = 1, kxym
      is0 = iand(kxym-1,ieor(mx+js,i-1))
      do 40 ii = 1, ntr
      if (kstrt.le.(kyb*kzb)) then
         is = is0 + kxym*(ii - 1)
         joff = 3*kxyp*is
         is = kzp*(is + ioff) - 1
         do 30 l = 1, kzp
         ir = kyp*(l + is) - 1
         koff = kypd*(l + loff) - 1
         do 20 k = 1, kyp
         do 10 j = 1, 3*kxyp
         g(j+3*kxyp*(k+ir)) = f(j+joff+3*nxv*(k+koff))
   10    continue
   20    continue
   30    continue
      endif
   40 continue
   50 continue
   60 continue
   70 continue
c exchange data
      do 110 mz = 1, lblok
      moff = jkblok*(mz - 1)
      do 100 mx = 1, jkblok
      m = mx + moff
      do 90 i = 1, kxym
      ir0 = iand(kxym-1,ieor(mx+js,i-1))
      is0 = ir0
      do 80 ii = 1, mntr
c post receive
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.mtr)) then
         is = ir0 + kxym*(ii - 1)
         ir = is + kyb*(mz + ls) + 1
         call MPI_IRECV(f(1+3*kxyp*kyp*kzp*is),3*kxyp*kyp*kzp,mcplx,ir-1
     1,ir+kxym+1,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kyb*kzb)).and.(ii.le.ntr)) then
         ir = is0 + kxym*(ii - 1)
         is = ir + kxb*(mz + ks) + 1
         call MPI_SEND(g(1+3*kxyp*kyp*kzp*ir),3*kxyp*kyp*kzp,mcplx,is-1,
     1m+kstrt+kxym,lgrp,ierr)
      endif
c receive data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
   80 continue
   90 continue
  100 continue
  110 continue
c transpose local data
      do 180 mz = 1, lblok
      moff = jblok*(mz - 1)
      do 170 mx = 1, jblok
      m = mx + moff
      ioff = kyb*(m - 1) 
      loff = kzpd*(m - 1) - 1
      do 160 i = 1, kxym
      ir0 = iand(kxym-1,ieor(mx+js,i-1))
      do 150 ii = 1, mtr
      if (kstrt.le.(kxb*kzb)) then
         ir = ir0 + kxym*(ii - 1)
         koff = kyp*ir
         ir = kzp*(ir + ioff) - 1
         do 140 l = 1, kzp
         is = kyp*(l + ir) - 1
         joff = kxypd*(l + loff) - 1
         do 130 k = 1, kyp
         do 120 j = 1, kxyp
         g(3*(k+koff+nyv*(j+joff))-2) = f(3*(j+kxyp*(k+is))-2)
         g(3*(k+koff+nyv*(j+joff))-1) = f(3*(j+kxyp*(k+is))-1)
         g(3*(k+koff+nyv*(j+joff))) = f(3*(j+kxyp*(k+is)))
  120    continue
  130    continue
  140    continue
      endif
  150 continue
  160 continue
  170 continue
  180 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine P3TPOS3BX(g,h,nx,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxyp
     1d,kyzpd,kzpd,jblok,mblok,lblok)
c this subroutine performs a transpose of a matrix g, distributed in x
c and z to a matrix h, distributed in x and y, that is,
c h(1:3,j+kzp*(l-1),k,l,my,mz) = g(1:3,k+kyzp*(m-1),j,l,mx,mz), where
c 1 <= j <= kxyp, 1 <= k <= kyzp, 1 <= l <= kzp, and
c 1 <= mx <= nx/kxyp, 1 <= my <= ny/kyzp, 1 <= mz <= nz/kzp
c and where indices mx, my, and mz can be distributed across processors.
c this subroutine sends and receives multiple asynchronous messages
c g = complex input array
c h = complex output array
c nx/ny/nz = number of points in x/y/z
c kstrt = starting data block number
c nyv/nzv = first dimension of g/h
c kxypd = second dimension of f and g
c kzpd/kyzpd = third dimension of g/h
c kxyp/kyzp/kzp = number of data values per block in x/y/z
c jblok/mblok/lblok = number of data blocks in x/y/z
c optimized version
      implicit none
      integer nx, ny, nz, kstrt, nyv, nzv, kxyp, kyzp, kzp
      integer kxypd, kyzpd, kzpd, jblok, mblok, lblok
      complex g, h
      dimension g(3*nyv*kxypd*kzpd*jblok*lblok)
      dimension h(3*nzv*kxypd*kyzpd*jblok*mblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=8)
c lgrp = current communicator
c mcplx = default datatype for complex
      common /pparms/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer js, ks, kxb, kyzb, kzb
      integer mlblok, kyzm, mtr, ntr, mntr, nzyv
      integer m, mx, my, l, i, moff, ioff, joff, koff, loff, k, j
      integer ir0, is0, ii, ir, is, ierr, msid, istatus
      dimension istatus(lstat)
      kxb = nx/kxyp
      kyzb = ny/kyzp
      kzb = nz/kzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
c this segment is used for shared memory computers
c     do 60 my = 1, mblok
c     moff = jblok*(my - 1)
c     koff = kyzp*(my + ks) - 1
c     do 50 mx = 1, jblok
c     m = mx + moff
c     do 40 i = 1, kzb
c     loff = kzp*(i - 1) - 1
c     ioff = jblok*(i - 1) - 1
c     do 30 l = 1, kzp
c     do 20 j = 1, kxyp
c     do 10 k = 1, kyzp
c     h(1+3*(l+loff+nzv*(j-1+kxypd*(k-1+kyzpd*(m-1))))) = g(1+3*(k+koff+
c    1nyv*(j-1+kxypd*(l-1+kzpd*(mx+ioff)))))
c     h(2+3*(l+loff+nzv*(j-1+kxypd*(k-1+kyzpd*(m-1))))) = g(2+3*(k+koff+
c    1nyv*(j-1+kxypd*(l-1+kzpd*(mx+ioff)))))
c     h(3+3*(l+loff+nzv*(j-1+kxypd*(k-1+kyzpd*(m-1))))) = g(3+3*(k+koff+
c    1nyv*(j-1+kxypd*(l-1+kzpd*(mx+ioff)))))
c  10 continue
c  20 continue
c  30 continue
c  40 continue
c  50 continue
c  60 continue
c this segment is used for mpi computers
      mlblok = max0(mblok,lblok)
      kyzm = min0(kyzb,kzb)
      mtr = kzb/kyzm
      ntr = kyzb/kyzm
      mntr = max0(mtr,ntr)
c transpose local data
      do 70 my = 1, lblok
      moff = jblok*(my - 1)
      do 60 mx = 1, jblok
      m = mx + moff
      ioff = kyzb*(m - 1)
      loff = kzpd*(m - 1) - 1
      do 50 i = 1, kyzm
      is0 = iand(kyzm-1,ieor(my+ks,i-1))
      do 40 ii = 1, ntr
      if (kstrt.le.(kxb*kzb)) then
         is = is0 + kyzm*(ii - 1)
         koff = 3*kyzp*is
         is = kzp*(is + ioff) - 1
         do 30 l = 1, kzp
         ir = kxyp*(l + is) - 1
         joff = kxypd*(l + loff) - 1
         do 20 j = 1, kxyp
         do 10 k = 1, 3*kyzp
         h(k+3*kyzp*(j+ir)) = g(k+koff+3*nyv*(j+joff))
   10    continue
   20    continue
   30    continue
      endif
   40 continue
   50 continue
   60 continue
   70 continue
c exchange data
      do 110 my = 1, mlblok
      moff = jblok*(my - 1)
      do 100 mx = 1, jblok
      m = mx + moff
      do 90 i = 1, kyzm
      ir0 = iand(kyzm-1,ieor(my+ks,i-1))
      is0 = ir0
      do 80 ii = 1, mntr
c post receive
      if ((kstrt.le.(kxb*kyzb)).and.(ii.le.mtr)) then
         is = ir0 + kyzm*(ii - 1)
         ir = mx + js + kxb*is + 1
         call MPI_IRECV(g(1+3*kxyp*kyzp*kzp*is),3*kxyp*kyzp*kzp,mcplx,ir
     1-1,ir+kyzm+1,lgrp,msid,ierr)
      endif
c send data
      if ((kstrt.le.(kxb*kzb)).and.(ii.le.ntr)) then
         ir = is0 + kyzm*(ii - 1)
         is = mx + js + kxb*ir + 1
         call MPI_SEND(h(1+3*kxyp*kyzp*kzp*ir),3*kxyp*kyzp*kzp,mcplx,is-
     11,m+kstrt+kyzm,lgrp,ierr)
      endif
c receive data
      if ((kstrt.le.(kxb*kyzb)).and.(ii.le.mtr)) then
         call MPI_WAIT(msid,istatus,ierr)
      endif
   80 continue
   90 continue
  100 continue
  110 continue
c transpose local data
      nzyv = nzv*kxypd
      do 180 my = 1, mblok
      moff = jblok*(my - 1)
      do 170 mx = 1, jblok
      m = mx + moff
      ioff = kzb*(m - 1)
      koff = kyzpd*(m - 1) - 1
      do 160 i = 1, kyzm
      ir0 = iand(kyzm-1,ieor(my+ks,i-1))
      do 150 ii = 1, mtr
      if (kstrt.le.(kxb*kyzb)) then
         ir = ir0 + kyzm*(ii - 1)
         joff = kzp*ir
         ir = kzp*(ir + ioff) - 1
         do 140 l = 1, kzp
         is = kxyp*(l + ir) - 1
         do 130 j = 1, kxyp
         loff = nzv*(j - 1) + joff
         do 120 k = 1, kyzp
         h(3*(l+loff+nzyv*(k+koff))-2) = g(3*(k+kyzp*(j+is))-2)
         h(3*(l+loff+nzyv*(k+koff))-1) = g(3*(k+kyzp*(j+is))-1)
         h(3*(l+loff+nzyv*(k+koff))) = g(3*(k+kyzp*(j+is)))
  120    continue
  130    continue
  140    continue
      endif
  150 continue
  160 continue
  170 continue
  180 continue
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
      subroutine timera(icntrl,chr,time)
c this subroutine performs timing
c input: icntrl, chr
c icntrl = (-1,0,1) = (initialize,ignore,read) clock
c clock should be initialized before it is read!
c chr = character variable for labeling timings
c time = elapsed time in seconds
c written for mpi
      implicit none
      integer icntrl
      character*8 chr
      real time
c get definition of MPI constants
      include 'mpif.h'
c common block for parallel processing
      integer nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c lgrp = current communicator
c mreal = default datatype for reals
      common /pparms/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer idproc, ierr
      real nclock, mclock
      double precision jclock
      save jclock
   91 format (1x,a8,1x,'max/min real time = ',e14.7,1x,e14.7,1x,'sec')
      data jclock /0.0d0/
      if (icntrl.eq.0) return
      if (icntrl.eq.1) go to 10
c initialize clock
      call MPI_BARRIER(lgrp,ierr)
      jclock = MPI_WTIME()
      return
c read clock and write time difference from last clock initialization
   10 nclock = real(MPI_WTIME() - jclock)
      call MPI_ALLREDUCE(nclock,time,1,mreal,MPI_MIN,lgrp,ierr)
      mclock = time
      call MPI_ALLREDUCE(nclock,time,1,mreal,MPI_MAX,lgrp,ierr)
      call MPI_COMM_RANK(lgrp,idproc,ierr)
      if (idproc.eq.0) write (6,91) chr, time, mclock
      return
      end
c-----------------------------------------------------------------------
      subroutine psum (f,g,nxp,nblok)
c this subroutine performs a parallel sum of a vector, that is:
c f(j,k) = sum over k of f(j,k)
c assumes the number of processors nproc is a power of two.
c the algorithm performs partial sums in binary pairs, as follows:
c first, adjacent processors exchange vectors and sum them.  next,
c processors separated by 2 exchange the new vectors and sum them, then
c those separated by 4, up to processors separated by nproc/2.  at the
c end, all processors contain the same summation.
c f = input and output data
c g = scratch array
c nxp = number of data values in vector
c nblok = number of data blocks
c written by viktor k. decyk, ucla
      implicit none
      real f, g
      integer nxp, nblok
      dimension f(nxp,nblok), g(nxp,nblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=8)
c nproc = number of real or virtual processors obtained
c lgrp = current communicator
c mreal = default datatype for reals
      common /pparms/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus
      integer idproc, ierr, kstrt, ks, l, kxs, k, kb, lb, msid, j
      dimension istatus(lstat)
c find processor id
c this line is used for shared memory computers
c     idproc = 0
c this line is used for mpi computers
      call MPI_COMM_RANK(lgrp,idproc,ierr)
      kstrt = idproc + 1
      if (kstrt.gt.nproc) return
      ks = kstrt - 2
      l = 1
      kxs = 1
c main iteration loop
   10 if (kxs.ge.nproc) go to 60
c shift data
      do 30 k = 1, nblok
      kb = k + ks
      lb = kb/kxs
      kb = kb + 1
      lb = lb - 2*(lb/2)
c this loop is used for shared memory computers
c     do 20 j = 1, nxp
c     if (lb.eq.0) then
c        g(j,k) = f(j,kb+kxs)
c     else
c        g(j,k) = f(j,kb-kxs)
c     endif
c  20 continue
c this segment is used for mpi computers
      if (lb.eq.0) then
         call MPI_IRECV(g,nxp,mreal,kb+kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(f,nxp,mreal,kb+kxs-1,l+nxp,lgrp,ierr)
      else
         call MPI_IRECV(g,nxp,mreal,kb-kxs-1,l+nxp,lgrp,msid,ierr)
         call MPI_SEND(f,nxp,mreal,kb-kxs-1,l+nxp,lgrp,ierr)
      endif
      call MPI_WAIT(msid,istatus,ierr)
   30 continue
c perform sum
      do 50 k = 1, nblok
      do 40 j = 1, nxp
      f(j,k) = f(j,k) + g(j,k)
   40 continue
   50 continue
      l = l + 1
      kxs = kxs + kxs
      go to 10
   60 return
      end
c-----------------------------------------------------------------------
      subroutine DCOMP32(edges,nyzp,noff,ny,nz,kstrt,nvpy,nvpz,idps,idds
     1,mblok,nblok)
c this subroutine determines spatial boundaries for particle
c decomposition, calculates number of grid points in each spatial
c region, and the offset of these grid points from the global address
c for 2D spatial decomposition
c edges(1,m) = lower boundary in y of particle partition m
c edges(2,m) = upper boundary in y of particle partition m
c edges(3,m) = back boundary in z of particle partition m
c edges(4,m) = front boundary in z of particle partition m
c nyzp(1,m) = number of primary gridpoints in y in particle partition m
c nyzp(2,m) = number of primary gridpoints in z in particle partition m
c noff(1,m) = lowermost global gridpoint in y in particle partition m
c noff(2,m) = backmost global gridpoint in z in particle partition m
c where m = my + mblok*(mz - 1)
c ny/nz = system length in y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c idps = number of particle partition boundaries
c idds = dimensionality of domain decomposition
c mblok/nblok = number of particle partitions in y/z
      implicit none
      real edges
      integer nyzp, noff, ny, nz, kstrt, nvpy, nvpz, idps, idds
      integer mblok, nblok
      dimension edges(idps,mblok*nblok)
      dimension nyzp(idds,mblok*nblok), noff(idds,mblok*nblok)
c local data
      integer js, ks, jb, kb, kr, m, moff, my, mz
      real at1, at2
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      at1 = float(ny)/float(nvpy)
      at2 = float(nz)/float(nvpz)
      do 20 mz = 1, nblok
      moff = mblok*(mz - 1)
      kb = mz + ks
      do 10 my = 1, mblok
      m = my + moff
      jb = my + js
      edges(1,m) = at1*float(jb)
      noff(1,m) = edges(1,m) + .5
      edges(2,m) = at1*float(jb + 1)
      kr = edges(2,m) + .5
      nyzp(1,m) = kr - noff(1,m)
      edges(3,m) = at2*float(kb)
      noff(2,m) = edges(3,m) + .5
      edges(4,m) = at2*float(kb + 1)
      kr = edges(4,m) + .5
      nyzp(2,m) = kr - noff(2,m)
   10 continue
   20 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine GTMODES3(pot,pott,nx,ny,nz,it,modesx,modesy,modesz,nxe,
     1nye,nze,nt2,modesxd,modesyd,modeszd)
c this subroutine extracts lowest order modes from complex array pot
c and stores them into a location in a time history array pott
c nx/ny/nz = system length in x/y/z direction
c it = current time
c modesx/modesy/modesz = number of modes to store in x/y/z direction,
c where modesx <= nx/2, modesy <= ny/2, modesz <= nz/2
c nxe = first dimension of input array pot, nxe >= nx
c nye = second dimension of input array pot, nye >= ny
c nze = third dimension of input array pot, nze >= nz
c nt2 = first dimension of output array pott, nt2 >= 2*it
c modesxd = second dimension of output array pott, modesxd >= modesx
c modesyd = third dimension of output array pott, modesyd  = 2*modesy-1
c modeszd = fourth dimension of output array pott, modeszd  = 2*modesz-1
      dimension pot(nxe,nye,nze), pott(nt2,modesxd,modesyd,modeszd)
      i2 = it + it
      i1 = i2 - 1
      if (i2.gt.nt2) return
      if ((modesx.gt.(nx/2)).or.(modesy.gt.(ny/2)).or.(modesz.gt.(nz/2))
     1) return
      ny2 = ny + 2
      nz2 = nz + 2
      do 40 l = 2, modesz
      l1 = nz2 - l
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      do 20 k = 2, modesy
      k1 = ny2 - k
      do 10 j = 2, modesx
      pott(i1,j,2*k-2,2*l-2) = pot(2*j-1,k,l)
      pott(i2,j,2*k-2,2*l-2) = pot(2*j,k,l)
      pott(i1,j,2*k-1,2*l-2) = pot(2*j-1,k1,l)
      pott(i2,j,2*k-1,2*l-2) = pot(2*j,k1,l)
      pott(i1,j,2*k-2,2*l-1) = pot(2*j-1,k,l1)
      pott(i2,j,2*k-2,2*l-1) = pot(2*j,k,l1)
      pott(i1,j,2*k-1,2*l-1) = pot(2*j-1,k1,l1)
      pott(i2,j,2*k-1,2*l-1) = pot(2*j,k1,l1)
   10 continue
c mode numbers kx = 0, nx/2
      pott(i1,1,2*k-2,2*l-2) = pot(1,k,l)
      pott(i2,1,2*k-2,2*l-2) = pot(2,k,l)
      pott(i1,1,2*k-1,2*l-2) = pot(1,k,l1)
      pott(i2,1,2*k-1,2*l-2) = -pot(2,k,l1)
      pott(i1,1,2*k-2,2*l-1) = pot(1,k,l1)
      pott(i2,1,2*k-2,2*l-1) = pot(2,k,l1)
      pott(i1,1,2*k-1,2*l-1) = pot(1,k,l)
      pott(i2,1,2*k-1,2*l-1) = -pot(2,k,l)
   20 continue
c mode numbers ky = 0, ny/2
      do 30 j = 2, modesx
      pott(i1,j,1,2*l-2) = pot(2*j-1,1,l)
      pott(i2,j,1,2*l-2) = pot(2*j,1,l)
      pott(i1,j,1,2*l-1) = pot(2*j-1,1,l1)
      pott(i2,j,1,2*l-1) = pot(2*j,1,l1)
   30 continue
c mode numbers kx = 0, nx/2
      pott(i1,1,1,2*l-2) = pot(1,1,l)
      pott(i2,1,1,2*l-2) = pot(2,1,l)
      pott(i1,1,1,2*l-1) = pot(1,1,l)
      pott(i2,1,1,2*l-1) = -pot(2,1,l)
   40 continue
c mode numbers kz = 0, nz/2
      do 60 k = 2, modesy
      k1 = ny2 - k
      do 50 j = 2, modesx
      pott(i1,j,2*k-2,1) = pot(2*j-1,k,1)
      pott(i2,j,2*k-2,1) = pot(2*j,k,1)
      pott(i1,j,2*k-1,1) = pot(2*j-1,k1,1)
      pott(i2,j,2*k-1,1) = pot(2*j,k1,1)
   50 continue
c mode numbers kx = 0, nx/2
      pott(i1,1,2*k-2,1) = pot(1,k,1)
      pott(i2,1,2*k-2,1) = pot(2,k,1)
      pott(i1,1,2*k-1,1) = pot(1,k,1)
      pott(i2,1,2*k-1,1) = -pot(2,k,1)
   60 continue
c mode numbers ky = 0, ny/2
      do 70 j = 2, modesx
      pott(i1,j,1,1) = pot(2*j-1,1,1)
      pott(i2,j,1,1) = pot(2*j,1,1)
   70 continue
      pott(i1,1,1,1) = pot(1,1,1)
      pott(i2,1,1,1) = 0.
      return
      end
c-----------------------------------------------------------------------
      subroutine PGTMODES33(pot,pott,nx,ny,nz,it,kstrt,modesz,nzv,kxyp,k
     1yzp,jblok,mblok,nt,modeszd)
c this subroutine extracts lowest order modes from complex array pot
c and stores them into a location in a time history array pott
c nx/ny/nz = system length in x/y/z direction
c it = current time
c modesx/modesy/modesz = number of modes to store in x/y/z direction,
c where modesx <= nx/2, modesy <= ny/2, modesz <= nz/2
c nzv = first dimension of field arrays, must be >= nz
c kstrt = starting data block number
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c jblok/mblok = number of field partitions in x/y
c nt = first dimension of output array pott, nt >= it
c modeszd = second dimension of output array pott, modeszd  = 2*modesz-1
      complex pot, pott, zero
      dimension pot(nzv,kxyp,kyzp,jblok*mblok)
      dimension pott(nt,modeszd,kxyp,kyzp,jblok*mblok)
      if (it.gt.nt) return
      if (modesz.gt.(nz/2)) return
      nxh = nx/2
      nyh = ny/2
      nz2 = nz + 2
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      zero = cmplx(0.,0.)
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
         do 20 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 10 l = 2, modesz
            l1 = nz2 - l
            pott(it,2*l-2,j,k,m) = pot(l,j,k,m)
            pott(it,2*l-1,j,k,m) = pot(l1,j,k,m)
   10       continue
c mode numbers kz = 0, nz/2
            pott(it,1,j,k,m) = pot(1,j,k,m)
         endif
   20    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
c kx = 0
            if (k1.gt.0) then
               do 30 l = 2, modesz
               l1 = nz2 - l
               pott(it,2*l-2,1,k,m) = pot(l,1,k,m)
               pott(it,2*l-1,1,k,m) = pot(l1,1,k,m)
   30          continue
c mode numbers kz = 0, nz/2
               pott(it,1,1,k,m) = pot(1,1,k,m)
c kx = nx/2
            else
               do 40 l = 2, modesz
               l1 = nz2 - l
               pott(it,2*l-2,1,k,m) = zero
               pott(it,2*l-1,1,k,m) = zero
   40          continue
c mode numbers kz = 0, nz/2
               pott(it,1,1,k,m) = zero
            endif
         endif
      endif
   50 continue
c mode numbers ky = 0, ny/2
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 70 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 60 l = 2, modesz
            l1 = nz2 - l
            pott(it,2*l-2,j,1,m) = pot(l,j,1,m)
            pott(it,2*l-1,j,1,m) = pot(l1,j,1,m)
   60       continue
c mode numbers kz = 0, nz/2
            pott(it,1,j,1,m) = pot(1,j,1,m)
         endif
   70    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 80 l = 2, modesz
            pott(it,2*l-2,1,1,m) = pot(l,1,1,m)
            pott(it,2*l-1,1,1,m) = conjg(pot(l,1,1,m))
   80       continue
c mode numbers kz = 0, nz/2
            pott(it,1,1,1,m) = cmplx(real(pot(1,1,1,m)),0.)
         endif
      endif
c throw away ky = ny/2
      k1 = (nyh/kyzp)*kyzp
      if (n2.eq.k1) then
         k1 = nyh - k1 + 1
         do 100 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 90 l = 2, modesz
            pott(it,2*l-2,j,k1,m) = zero
            pott(it,2*l-1,j,k1,m) = zero
   90       continue
c mode numbers kz = 0, nz/2
            pott(it,1,j,k1,m) = zero
         endif
  100    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
         do 110 l = 2, modesz
            pott(it,2*l-2,1,k1,m) = zero
            pott(it,2*l-1,1,k1,m) = zero
  110       continue
c mode numbers kz = 0, nz/2
            pott(it,1,1,k1,m) = zero
         endif
      endif
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PTMODES3(pot,pott,nx,ny,nz,it,modesx,modesy,modesz,nxe,
     1nye,nze,nt2,modesxd,modesyd,modeszd)
c this subroutine extracts lowest order modes from a location in a time
c history array pott and stores them into complex array pot
c nx/ny/nz = system length in x/y/z direction
c it = current time
c modesx/modesy/modesz = number of modes to store in x/y/z direction,
c where modesx <= nx/2, modesy <= ny/2, modesz <= nz/2
c nxe = first dimension of input array pot, nxe >= nx
c nye = second dimension of input array pot, nye >= ny
c nze = third dimension of input array pot, nze >= nz
c nt2 = first dimension of output array pott, nt2 >= 2*it
c modesxd = second dimension of output array pott, modesxd >= modesx
c modesyd = third dimension of output array pott, modesyd  = 2*modesy-1
c modeszd = fourth dimension of output array pott, modeszd  = 2*modesz-1
      dimension pot(nxe,nye,nze), pott(nt2,modesxd,modesyd,modeszd)
      i2 = it + it
      i1 = i2 - 1
      if (i2.gt.nt2) return 
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2
      if ((modesx.gt.nxh).or.(modesy.gt.nyh).or.(modesz.gt.nzh)) return
      ny2 = ny + 2
      nz2 = nz + 2
      do 80 l = 2, modesz
      l1 = nz2 - l
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      do 30 k = 2, modesy
      k1 = ny2 - k
      do 10 j = 2, modesx
      pot(2*j-1,k,l) = pott(i1,j,2*k-2,2*l-2)
      pot(2*j,k,l) = pott(i2,j,2*k-2,2*l-2)
      pot(2*j-1,k1,l) = pott(i1,j,2*k-1,2*l-2)
      pot(2*j,k1,l) = pott(i2,j,2*k-1,2*l-2)
      pot(2*j-1,k,l1) = pott(i1,j,2*k-2,2*l-1)
      pot(2*j,k,l1) = pott(i2,j,2*k-2,2*l-1)
      pot(2*j-1,k1,l1) = pott(i1,j,2*k-1,2*l-1)
      pot(2*j,k1,l1) = pott(i2,j,2*k-1,2*l-1)
   10 continue
      do 20 j = modesx+1, nxh
      pot(2*j-1,k,l) = 0.
      pot(2*j,k,l) = 0.
      pot(2*j-1,k1,l) = 0.
      pot(2*j,k1,l) = 0.
      pot(2*j-1,k,l1) = 0.
      pot(2*j,k,l1) = 0.
      pot(2*j-1,k1,l1) = 0.
      pot(2*j,k1,l1) = 0.
   20 continue
c mode numbers kx = 0, nx/2
      pot(1,k,l) = pott(i1,1,2*k-2,2*l-2)
      pot(2,k,l) = pott(i2,1,2*k-2,2*l-2)
      pot(1,k1,l) = 0.
      pot(2,k1,l) = 0.
      pot(1,k,l1) = pott(i1,1,2*k-2,2*l-1)
      pot(2,k,l1) = pott(i2,1,2*k-2,2*l-1)
      pot(1,k1,l1) = 0.
      pot(2,k1,l1) = 0.
   30 continue
      do 50 k = modesy+1, nyh
      k1 = ny2 - k
      do 40 j = 1, nxh
      pot(2*j-1,k,l) = 0.
      pot(2*j,k,l) = 0.
      pot(2*j-1,k1,l) = 0.
      pot(2*j,k1,l) = 0.
      pot(2*j-1,k,l1) = 0.
      pot(2*j,k,l1) = 0.
      pot(2*j-1,k1,l1) = 0.
      pot(2*j,k1,l1) = 0.
   40 continue
   50 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 60 j = 2, modesx
      pot(2*j-1,1,l) = pott(i1,j,1,2*l-2)
      pot(2*j,1,l) = pott(i2,j,1,2*l-2)
      pot(2*j-1,k1,l) = 0.
      pot(2*j,k1,l) = 0.
      pot(2*j-1,1,l1) = pott(i1,j,1,2*l-1)
      pot(2*j,1,l1) = pott(i2,j,1,2*l-1)
      pot(2*j-1,k1,l1) = 0.
      pot(2*j,k1,l1) = 0.
   60 continue
c mode numbers kx = 0, nx/2
      pot(1,1,l) = pott(i1,1,1,2*l-2)
      pot(2,1,l) = pott(i2,1,1,2*l-2)
      pot(1,k1,l) = 0.
      pot(2,k1,l) = 0.
      pot(1,1,l1) = 0.
      pot(2,1,l1) = 0.
      pot(1,k1,l1) = 0.
      pot(2,k1,l1) = 0.
      do 70 j = modesx+1, nxh
      pot(2*j-1,1,l) = 0.
      pot(2*j,1,l) = 0.
      pot(2*j-1,k1,l) = 0.
      pot(2*j,k1,l) = 0.
      pot(2*j-1,1,l1) = 0.
      pot(2*j,1,l1) = 0.
      pot(2*j-1,k1,l1) = 0.
      pot(2*j,k1,l1) = 0.
   70 continue
   80 continue
      do 120 l = modesz+1, nzh
      l1 = nz2 - l
      do 100 k = 2, nyh
      k1 = ny2 - k
      do 90 j = 1, nxh
      pot(2*j-1,k,l) = 0.
      pot(2*j,k,l) = 0.
      pot(2*j-1,k1,l) = 0.
      pot(2*j,k1,l) = 0.
      pot(2*j-1,k,l1) = 0.
      pot(2*j,k,l1) = 0.
      pot(2*j-1,k1,l1) = 0.
      pot(2*j,k1,l1) = 0.
   90 continue
  100 continue
      k1 = nyh + 1
      do 110 j = 1, nxh
      pot(2*j-1,1,l) = 0.
      pot(2*j,1,l) = 0.
      pot(2*j-1,k1,l) = 0.
      pot(2*j,k1,l) = 0.
      pot(2*j-1,1,l1) = 0.
      pot(2*j,1,l1) = 0.
      pot(2*j-1,k1,l1) = 0.
      pot(2*j,k1,l1) = 0.
  110 continue
  120 continue
c mode numbers kz = 0, nz/2
      l1 = nzh + 1
      do 150 k = 2, nyh
      k1 = ny2 - k
      do 130 j = 2, modesx
      pot(2*j-1,k,1) = pott(i1,j,2*k-2,1)
      pot(2*j,k,1) = pott(i2,j,2*k-2,1)
      pot(2*j-1,k1,1) = pott(i1,j,2*k-1,1)
      pot(2*j,k1,1) = pott(i2,j,2*k-1,1)
      pot(2*j-1,k,l1) = 0.
      pot(2*j,k,l1) = 0.
      pot(2*j-1,k1,l1) = 0.
      pot(2*j,k1,l1) = 0.
  130 continue
      do 140 j = modesx+1, nxh
      pot(2*j-1,k,1) = 0.
      pot(2*j,k,1) = 0.
      pot(2*j-1,k1,1) = 0.
      pot(2*j,k1,1) = 0.
      pot(2*j-1,k,l1) = 0.
      pot(2*j,k,l1) = 0.
      pot(2*j-1,k1,l1) = 0.
      pot(2*j,k1,l1) = 0.
  140 continue
c mode numbers kx = 0, nx/2
      pot(1,k,1) = pott(i1,1,2*k-2,1)
      pot(2,k,1) = pott(i2,1,2*k-2,1)
      pot(1,k1,1) = 0.
      pot(2,k1,1) = 0.
      pot(1,k,l1) = 0.
      pot(2,k,l1) = 0.
      pot(1,k1,l1) = 0.
      pot(2,k1,l1) = 0.
  150 continue
      do 170 k = modesy+1, nyh
      k1 = ny2 - k
      do 160 j = 1, nxh
      pot(2*j-1,k,1) = 0.
      pot(2*j,k,1) = 0.
      pot(2*j-1,k1,1) = 0.
      pot(2*j,k1,1) = 0.
      pot(2*j-1,k,l1) = 0.
      pot(2*j,k,l1) = 0.
      pot(2*j-1,k1,l1) = 0.
      pot(2*j,k1,l1) = 0.
  160 continue
  170 continue
c mode numbers ky = 0, ny/2
      k1 = nyh + 1
      do 180 j = 2, modesx
      pot(2*j-1,1,1) = pott(i1,j,1,1)
      pot(2*j,1,1) = pott(i2,j,1,1)
      pot(2*j-1,k1,1) = 0.
      pot(2*j,k1,1) = 0.
      pot(2*j-1,1,l1) = 0.
      pot(2*j,1,l1) = 0.
      pot(2*j-1,k1,l1) = 0.
      pot(2*j,k1,l1) = 0.
  180 continue
      do 190 j = modesx+1, nxh
      pot(2*j-1,1,1) = 0.
      pot(2*j,1,1) = 0.
      pot(2*j-1,k1,1) = 0.
      pot(2*j,k1,1) = 0.
      pot(2*j-1,1,l1) = 0.
      pot(2*j,1,l1) = 0.
      pot(2*j-1,k1,l1) = 0.
      pot(2*j,k1,l1) = 0.
  190 continue
      pot(1,1,1) = pott(i1,1,1,1)
      pot(2,1,1) = 0.
      pot(1,k1,1) = 0.
      pot(2,k1,1) = 0.
      pot(1,1,l1) = 0.
      pot(2,1,l1) = 0.
      pot(1,k1,l1) = 0.
      pot(2,k1,l1) = 0.
      return
      end
c-----------------------------------------------------------------------
      subroutine PPTMODES33(pot,pott,nx,ny,nz,it,kstrt,modesz,nzv,kxyp,k
     1yzp,jblok,mblok,nt,modeszd)
c this subroutine extracts lowest order modes from complex array pot
c and stores them into a location in a time history array pott
c nx/ny/nz = system length in x/y/z direction
c it = current time
c modesx/modesy/modesz = number of modes to store in x/y/z direction,
c where modesx <= nx/2, modesy <= ny/2, modesz <= nz/2
c nzv = first dimension of field arrays, must be >= nz
c kstrt = starting data block number
c kxyp/kyzp = number of complex grids in each field partition in
c x/y direction
c jblok/mblok = number of field partitions in x/y
c nt = first dimension of output array pott, nt >= it
c modeszd = second dimension of output array pott, modeszd  = 2*modesz-1
      complex pot, pott, zero
      dimension pot(nzv,kxyp,kyzp,jblok*mblok)
      dimension pott(nt,modeszd,kxyp,kyzp,jblok*mblok)
      if (it.gt.nt) return
      nxh = nx/2
      nyh = ny/2
      nzh = nz/2
      if (modesz.gt.nzh) return
      nz2 = nz + 2
      kxb = nxh/kxyp
      kyzb = ny/kyzp
      ks = (kstrt - 1)/kxb
      js = kstrt - kxb*ks - 2
      ks = ks - 1
      zero = cmplx(0.,0.)
      if (kstrt.gt.(kxb*kyzb)) return
      do 170 my = 1, mblok
      koff = kyzp*(my + ks) - 1
      moff = jblok*(my - 1)
      do 160 mx = 1, jblok
c mode numbers 0 < kx < nx/2, 0 < ky < ny/2, and 0 < kz < nz/2
      joff = kxyp*(mx + js) - 1
      m = mx + moff
      do 70 k = 1, kyzp
      k1 = k + koff
      if ((k1.gt.0).and.(k1.ne.nyh)) then
         if (k1.gt.nyh) k1 = k1 - ny
         do 30 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 10 l = 2, modesz
            l1 = nz2 - l
            pot(l,j,k,m) = pott(it,2*l-2,j,k,m)
            pot(l1,j,k,m) = pott(it,2*l-1,j,k,m)
   10       continue
            do 20 l = modesz+1, nzh
            l1 = nz2 - l
            pot(l,j,k,m) = zero
            pot(l1,j,k,m) = zero
   20       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1 
            pot(1,j,k,m) = pott(it,1,j,k,m)
            pot(l1,j,k,m) = zero
         endif
   30    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
c keep kx = 0
            if (k1.gt.0) then
               do 40 l = 2, modesz
               l1 = nz2 - l
               pot(l,1,k,m) = pott(it,2*l-2,1,k,m)
               pot(l1,1,k,m) = pott(it,2*l-1,1,k,m)
   40          continue
               do 50 l = modesz+1, nzh
               l1 = nz2 - l
               pot(l,1,k,m) = zero
               pot(l1,1,k,m) = zero
   50          continue
c mode numbers kz = 0, nz/2
               l1 = nzh + 1
               pot(1,1,k,m) = pott(it,1,1,k,m)
               pot(l1,1,k,m) = zero
c throw away kx = nx/2
            else
               do 60 l = 1, nz
               pot(l,1,k,m) = zero
   60          continue
            endif
         endif
      endif
   70 continue
c mode numbers ky = 0, ny/2
      n2 = koff + 1
c keep ky = 0
      if (n2.eq.0) then
         do 100 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 80 l = 2, modesz
            l1 = nz2 - l
            pot(l,j,1,m) = pott(it,2*l-2,j,1,m)
            pot(l1,j,1,m) = pott(it,2*l-1,j,1,m)
   80       continue
            do 90 l = modesz+1, nzh
            l1 = nz2 - l
            pot(l,j,1,m) = zero
            pot(l1,j,1,m) = zero
   90       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            pot(1,j,1,m) = pott(it,1,j,1,m)
            pot(l1,j,1,m) = zero
         endif
  100    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 110 l = 2, modesz
            l1 = nz2 - l
            pot(l,1,1,m) = pott(it,2*l-2,1,1,m)
            pot(l1,1,1,m) = zero
  110       continue
            do 120 l = modesz+1, nzh
            l1 = nz2 - l
            pot(l,1,1,m) = zero
            pot(l1,1,1,m) = zero
  120       continue
c mode numbers kz = 0, nz/2
            l1 = nzh + 1
            pot(1,1,1,m) = cmplx(real(pott(it,1,1,1,m)),0.)
            pot(l1,1,1,m) = zero
         endif
      endif
c throw away ky = ny/2
      k1 = (nyh/kyzp)*kyzp
      if (n2.eq.k1) then
         k1 = nyh - k1 + 1
         do 140 j = 1, kxyp
         if ((j+joff).gt.0) then
            do 130 l = 1, nz
            pot(l,j,k1,m) = zero
  130       continue
         endif
  140    continue
c mode numbers kx = 0, nx/2
         n1 = joff + 1
         if (n1.eq.0) then
            do 150 l = 1, nz
            pot(l,1,k1,m) = zero
  150       continue
         endif
      endif
  160 continue
  170 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PTRPSIN32C(cu,cu3,scb,scd,nx,ny,nz,kstrt,nvpy,nvpz,nxv,
     1kyp,kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
c this subroutine creates a tripled vector array cu3 from a vector array
c cu, so that various 3d sine/cosine transforms can be performed with a
c 3d real to complex fft.  the x component is an odd function in y,
c y component is an odd function in x, and the z component is an odd
c function in both x and y.  Asummes vector cu vanishes at end points
c linear interpolation for distributed data with 2D domain decomposition
c cu3 array may be modified
c scb/scd = scratch arrays
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = second dimension of input array cu, must be >= nx
c kyp/kzp = number of data values per block in cu in y/z
c kypd = third dimension of input array cu, must be >= kyp
c kzpd = fourth dimension of input array cu, must be >= kzp
c kyp2/kzp2 = number of data values per block in cu3 in y/z
c kblok/lblok = number of data blocks in y/z
c k2blok/l2blok = number of data blocks in y/z for tripled data
      implicit none
      real cu, cu3, scb, scd
      integer nx, ny, nz, kstrt, nvpy, nvpz, nxv, kyp, kzp, kypd, kzpd
      integer kyp2, kzp2, kblok, lblok, k2blok, l2blok
      dimension cu(3,nxv,kypd,kzpd,kblok*lblok)
      dimension cu3(3,2*nxv,2*kypd,kzp2,k2blok*l2blok)
      dimension scb(3,nxv,kypd,kzpd,kblok*lblok)
      dimension scd(3,nxv,kzpd,2,kblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      logical lt1
      integer istatus, lsid, msid, nsid, ierr
      integer i, j, k, l, m, my, mz, nxs, nys, nzs, ks, js, ny2, nz2, kr
      integer kyb, kyb2, kzb, kzb2, nxvy, noff, moff, loff, koff, joff
      integer jm, km, jl, kl, ky, kz, kk, ll, mm, lm, l1, l2, k0, k1, k2
      integer ii
      real at1
      dimension istatus(lstat)
      nxs = nx - 1
      nys = ny - 1
      nzs = nz - 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      kyb = ny/kyp
      ny2 = ny + ny
      kyb2 = ny2/kyp2
      kzb = nz/kzp
      nz2 = nz + nz
      kzb2 = nz2/kzp2
      nxvy = nxv*kypd
      noff = kypd*kzpd + kyb*kzb
c copy to triple array in y direction
      do 200 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 190 my = 1, k2blok
      m = my + moff
      loff = kyp2*(my + js)
      jm = loff/kyp + 1
      loff = kzp2*(mz + ks)
      km = loff/kzp + 1
      loff = kyp*(my + js)
      jl = loff/kyp2 + 1
      loff = kzp*(mz + ks)
      kl = loff - kzp2*(loff/kzp2)
      kz = nvpy*(mz + ks)
c special case for one processor
      if ((kyb2*kzb2).eq.1) then
         do 60 l = 1, nzs
         do 30 k = 1, nys
         do 10 j = 1, nxs
         cu3(1,j+1,k+1,l+1,m) = cu(1,j+1,k+1,l+1,m)
         cu3(2,j+1,k+1,l+1,m) = cu(2,j+1,k+1,l+1,m)
         cu3(3,j+1,k+1,l+1,m) = cu(3,j+1,k+1,l+1,m)
         cu3(1,nx+j+1,k+1,l+1,m) = cu(1,nx-j+1,k+1,l+1,m)
         cu3(2,nx+j+1,k+1,l+1,m) = -cu(2,nx-j+1,k+1,l+1,m)
         cu3(3,nx+j+1,k+1,l+1,m) = -cu(3,nx-j+1,k+1,l+1,m)
         cu3(1,j+1,ny+k+1,l+1,m) = -cu(1,j+1,ny-k+1,l+1,m)
         cu3(2,j+1,ny+k+1,l+1,m) = cu(2,j+1,ny-k+1,l+1,m)
         cu3(3,j+1,ny+k+1,l+1,m) = -cu(3,j+1,ny-k+1,l+1,m)
         cu3(1,nx+j+1,ny+k+1,l+1,m) = -cu(1,nx-j+1,ny-k+1,l+1,m)
         cu3(2,nx+j+1,ny+k+1,l+1,m) = -cu(2,nx-j+1,ny-k+1,l+1,m)
         cu3(3,nx+j+1,ny+k+1,l+1,m) = cu(3,nx-j+1,ny-k+1,l+1,m)
         cu3(1,j+1,k+1,nz+l+1,m) = -cu(1,j+1,k+1,nz-l+1,m)
         cu3(2,j+1,k+1,nz+l+1,m) = -cu(2,j+1,k+1,nz-l+1,m)
         cu3(3,j+1,k+1,nz+l+1,m) = cu(3,j+1,k+1,nz-l+1,m)
         cu3(1,nx+j+1,k+1,nz+l+1,m) = -cu(1,nx-j+1,k+1,nz-l+1,m)
         cu3(2,nx+j+1,k+1,nz+l+1,m) = cu(2,nx-j+1,k+1,nz-l+1,m)
         cu3(3,nx+j+1,k+1,nz+l+1,m) = -cu(3,nx-j+1,k+1,nz-l+1,m)
         cu3(1,j+1,ny+k+1,nz+l+1,m) = cu(1,j+1,ny-k+1,nz-l+1,m)
         cu3(2,j+1,ny+k+1,nz+l+1,m) = -cu(2,j+1,ny-k+1,nz-l+1,m)
         cu3(3,j+1,ny+k+1,nz+l+1,m) = -cu(3,j+1,ny-k+1,nz-l+1,m)
         cu3(1,nx+j+1,ny+k+1,nz+l+1,m) = cu(1,nx-j+1,ny-k+1,nz-l+1,m)
         cu3(2,nx+j+1,ny+k+1,nz+l+1,m) = cu(2,nx-j+1,ny-k+1,nz-l+1,m)
         cu3(3,nx+j+1,ny+k+1,nz+l+1,m) = cu(3,nx-j+1,ny-k+1,nz-l+1,m)
   10    continue
         do 20 i = 1, 3
         cu3(i,1,k+1,l+1,m) = 0.
         cu3(i,nx+1,k+1,l+1,m) = 0.
         cu3(i,1,k+ny+1,l+1,m) = 0.
         cu3(i,nx+1,k+ny+1,l+1,m) = 0.
         cu3(i,1,k+1,nz+l+1,m) = 0.
         cu3(i,nx+1,k+1,nz+l+1,m) = 0.
         cu3(i,1,k+ny+1,nz+l+1,m) = 0.
         cu3(i,nx+1,k+ny+1,nz+l+1,m) = 0.
   20    continue
   30    continue
         do 50 j = 1, nx
         do 40 i = 1, 3
         cu3(i,j,1,l+1,m) = 0.
         cu3(i,j+nx,1,l+1,m) = 0.
         cu3(i,j,ny+1,l+1,m) = 0.
         cu3(i,j+nx,ny+1,l+1,m) = 0.
         cu3(i,j,1,nz+l+1,m) = 0.
         cu3(i,j+nx,1,nz+l+1,m) = 0.
         cu3(i,j,ny+1,nz+l+1,m) = 0.
         cu3(i,j+nx,ny+1,nz+l+1,m) = 0.
   40    continue
   50    continue
   60    continue
         nxs = nx + nx
         nys = ny + ny
         do 90 k = 1, nys
         do 80 j = 1, nxs
         do 70 i = 1, 3
         cu3(i,j,k,1,m) = 0.
         cu3(i,j,k,nz+1,m) = 0.
   70    continue
   80    continue
   90    continue
         return
      endif
c copy main data in y direction
      do 180 ii = 1, 2
c for odd rows in z
      if (ii.eq.1) then
         mm = jm + 2*kz
         lm = jl + kz/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in z
      else if (ii.eq.2) then
         if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 180
         mm = mm + nvpy
         lm = lm - nvpy/2
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kyb).and.(km.le.kzb)) then
c        do 130 l = 1, kzp
c        do 120 k = 1, kyp
c        do 110 j = 1, nx
c        do 100 i = 1, 3
c        cu3(i,j,k,l+loff,m) = cu(i,j,k,l,mm)
c 100    continue
c 110    continue
c 120    continue
c 130    continue
c        if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
c           do 170 l = 1, kzp
c           do 160 k = 1, kyp
c           do 150 j = 1, nx
c           do 140 i = 1, 3
c           cu3(i,j,k+kyp,l+loff,m) = cu(i,j,k,l,mm+1)
c 140       continue
c 150       continue
c 160       continue
c 170       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kyb).and.(km.le.kzb)) then
         call MPI_IRECV(cu3(1,1,1,loff+1,m),3*kzp*nxvy,mreal,mm-1,noff+i
     1i,lgrp,msid,ierr)
         if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
            call MPI_IRECV(scb(1,1,1,1,m),3*kzp*nxvy,mreal,mm,noff+ii,lg
     1rp,nsid,ierr)
         endif
      endif
      if ((jl.le.((kyb2-1)/2+1)).and.lt1) then
         call MPI_SEND(cu(1,1,1,1,m),3*kzp*nxvy,mreal,lm-1,noff+ii,lgrp,
     1ierr)
      endif
c wait for data and unpack it
      if ((jm.le.kyb).and.(km.le.kzb)) then
         call MPI_WAIT(msid,istatus,ierr)
         do 130 l = 1, kzp
         l1 = kzp - l + 1
         l2 = (l1 - 1)/4 + 1
         koff = kypd*(l1 - 4*(l2 - 1) - 1)
         do 120 k = 1, kypd
         k1 = kypd - k + 1
         k0 = k1 - 1 + koff
         k2 = k0/2 + 1
         joff = nxv*(k0 - 2*(k2 - 1))
         do 110 j = 1, nxv
         do 100 i = 1, 3
         cu3(i,j,k1,l1+loff,m) = cu3(i,j+joff,k2,l2+loff,m)
  100    continue
  110    continue
  120    continue
  130    continue
         if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 170 l = 1, kzp
            do 160 k = 1, kypd
            do 150 j = 1, nxv
            do 140 i = 1, 3
            cu3(i,j,k+kyp,l+loff,m) = scb(i,j,k,l,m)
  140       continue
  150       continue
  160       continue
  170       continue
         endif
      endif
  180 continue
  190 continue
  200 continue
c copying in y direction not needed
      if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 380
c copy reflected data in y direction
      do 370 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 360 my = 1, k2blok
      m = my + moff
      loff = kyp2*(my + js)
      jm = (ny2 - loff - 1)/kyp + 1
      loff = kzp2*(mz + ks)
      km = loff/kzp + 1
      loff = kyp*(my + js)
      jl = (ny2 - loff - 1)/kyp2 + 1
      ky = loff + kyp2*jl - ny2
      loff = kzp*(mz + ks)
      kl = loff - kzp2*(loff/kzp2)
      kz = nvpy*(mz + ks)
      do 350 ii = 3, 4
c for odd rows in z
      if (ii.eq.3) then
         mm = jm + 2*kz
         lm = jl + kz/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in z
      else if (ii.eq.4) then
         if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 350
         mm = mm + nvpy
         lm = lm - nvpy/2
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kyb).and.(km.le.kzb)) then
c        if ((jm+1).le.kyb) then
c           do 230 l = 1, kzp
c           do 220 j = 1, nx
c           do 210 i = 1, 3
c           cu3(i,j,1,l+loff,m) = cu(i,j,1,l,mm+1)
c 210       continue
c 220       continue
c 230       continue
c        endif
c        if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
c           do 270 l = 1, kzp
c           do 260 k = 1, kyp
c           do 250 j = 1, nx
c           do 240 i = 1, 3
c           cu3(i,j,k+kyp,l+loff,m) = cu(i,j,k,l,mm)
c 240       continue
c 250       continue
c 260       continue
c 270       continue
c        endif
c        if (kyp.gt.1) then
c           do 310 l = 1, kzp
c           do 300 k = 2, kyp
c           do 290 j = 1, nx
c           do 280 i = 1, 3
c           cu3(i,j,k,l+loff,m) = cu(i,j,k,l,mm-1)
c 280       continue
c 290       continue
c 300       continue
c 310       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if ((jm+1).le.kyb) then
            call MPI_IRECV(scd(1,1,1,2,m),3*kzp*nxv,mreal,mm,noff+ii,lgr
     1p,lsid,ierr)
         endif
         if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
            call MPI_IRECV(scb(1,1,1,1,m),3*kzp*nxvy,mreal,mm-1,noff+ii,
     1lgrp,msid,ierr)
         endif
         if (kyp.gt.1) then
            call MPI_IRECV(cu3(1,1,1,loff+1,m),3*kzp*nxvy,mreal,mm-2,nof
     1f+ii,lgrp,nsid,ierr)
         endif
      endif
      if ((jl.gt.((kyb2-1)/2+1)).and.(jl.le.kyb2).and.lt1) then
         if (ky.eq.0) then
            if ((jl+1).le.kyb2) then
               do 230 l = 1, kzp
               do 220 j = 1, nxv
               do 210 i = 1, 3
               scd(i,j,l,1,m) = cu(i,j,1,l,m)
  210          continue
  220          continue
  230          continue
               call MPI_SEND(scd(1,1,1,1,m),3*kzp*nxv,mreal,lm,noff+ii,l
     1grp,ierr)
            endif
            if (kyp.gt.1) then
               call MPI_SEND(cu(1,1,1,1,m),3*kzp*nxvy,mreal,lm-1,noff+ii
     1,lgrp,ierr)
            endif
         else
            call MPI_SEND(cu(1,1,1,1,m),3*kzp*nxvy,mreal,lm-1,noff+ii,lg
     1rp,ierr)
         endif
      endif
c wait for data and unpack it
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if (kyp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 270 l = 1, kzp
            l1 = kzp - l + 1
            l2 = (l1 - 1)/4 + 1
            koff = kypd*(l1 - 4*(l2 - 1) - 1)
            do 260 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1 + koff
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 250 j = 1, nxv
            do 240 i = 1, 3
            cu3(i,j,k1,l1+loff,m) = cu3(i,j+joff,k2,l2+loff,m)
  240       continue
  250       continue
  260       continue
  270       continue
         endif
         if ((jm+1).le.kyb) then
            call MPI_WAIT(lsid,istatus,ierr)
            do 300 l = 1, kzp
            do 290 j = 1, nxv
            do 280 i = 1, 3
            cu3(i,j,1,l+loff,m) = scd(i,j,l,2,m)
  280       continue
  290       continue
  300       continue
         endif
         if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 340 l = 1, kzp
            do 330 k = 1, kypd
            do 320 j = 1, nxv
            do 310 i = 1, 3
            cu3(i,j,k+kyp,l+loff,m) = scb(i,j,k,l,m)
  310       continue
  320       continue
  330       continue
  340       continue
         endif
      endif
  350 continue
  360 continue
  370 continue
c copying in z direction not needed
  380 if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 860
c copy reflected data in z direction
      do 560 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 550 my = 1, k2blok
      m = my + moff
      loff = kzp2*(mz + ks)
      jm = (nz2 - loff - 1)/kzp + 1
      loff = kyp2*(my + js)
      km = loff/kyp + 1
      loff = kzp*(mz + ks)
      jl = (nz2 - loff - 1)/kzp2 + 1
      kz = loff + kzp2*jl - nz2
      loff = kyp*(my + js)
      kl = loff - kyp2*(loff/kyp2)
      kr = my + js
      do 500 ii = 5, 6
c for odd rows in y
      if (ii.eq.5) then
         mm = nvpy*jm + 2*kr
         lm = nvpy*jl + kr/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in y
      else if (ii.eq.6) then
         if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 500
         mm = mm + 1
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kzb).and.(km.le.kyb)) then
c        if (kzp.gt.1) then
c           do 420 l = 2, kzp
c           do 410 k = 1, kypd
c           do 400 j = 1, nx
c           do 390 i = 1, 3
c           cu3(i,j,k,l+loff,m) = cu(i,j,k,l,mm-2*nvpy+1)
c 390       continue
c 400       continue
c 410       continue
c 420       continue
c        endif
c        if ((jm+1).le.kzb) then
c           do 450 k = 1, kypd
c           do 440 j = 1, nx
c           do 430 i = 1, 3
c           cu3(i,j,k,loff+1,m) = cu(i,j,k,1,mm+1)
c 430       continue
c 440       continue
c 450       continue
c        endif
c        if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
c           do 490 l = 1, kzp
c           do 480 k = 1, kypd
c           do 470 j = 1, nx
c           do 460 i = 1, 3
c           cu3(i,j,k+kyp,l+loff,m) = cu(i,j,k,l,mm-nvpy+1)
c 460       continue
c 470       continue
c 480       continue
c 490       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kzb).and.(km.le.kyb)) then
         if ((jm+1).le.kzb) then
            call MPI_IRECV(cu3(1,1,1,loff+1,m),3*nxvy,mreal,mm,noff+ii,l
     1grp,lsid,ierr)
         endif
         if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
            call MPI_IRECV(scb(1,1,1,1,m),3*kzp*nxvy,mreal,mm-nvpy,noff+
     1ii,lgrp,msid,ierr)
         endif
         if (kzp.gt.1) then
            call MPI_IRECV(cu3(1,1,1,loff+2,m),3*(kzp-1)*nxvy,mreal,mm-2
     1*nvpy,noff+ii,lgrp,nsid,ierr)
         endif
      endif
      if ((jl.gt.((kzb2-1)/2+1)).and.(jl.le.kzb2).and.lt1) then
         if (kz.eq.0) then
            if ((jl+1).le.kzb2) then
               call MPI_SEND(cu(1,1,1,1,m),3*nxvy,mreal,lm,noff+ii,lgrp,
     1ierr)
            endif
            if (kzp.gt.1) then
               call MPI_SEND(cu(1,1,1,2,m),3*(kzp-1)*nxvy,mreal,lm-nvpy,
     1noff+ii,lgrp,ierr)
            endif
         else
            call MPI_SEND(cu(1,1,1,1,m),3*kzp*nxvy,mreal,lm-nvpy,noff+ii
     1,lgrp,ierr)
         endif
      endif
c wait for data and unpack it
      if ((jm.le.kzb).and.(km.le.kyb)) then
         if (kzp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 420 l = 2, kzp
            l1 = kzp - l + 2
            l2 = (l1 + 2)/4 + 1
            koff = kypd*(l1 - 4*(l2 - 1) + 2)
            do 410 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1 + koff
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 400 j = 1, nxv
            do 390 i = 1, 3
            cu3(i,j,k1,l1+loff,m) = cu3(i,j+joff,k2,l2+loff,m)
  390       continue
  400       continue
  410       continue
  420       continue
         endif
         if ((jm+1).le.kzb) then
            call MPI_WAIT(lsid,istatus,ierr)
            do 450 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 440 j = 1, nxv
            do 430 i = 1, 3
            cu3(i,j,k1,loff+1,m) = cu3(i,j+joff,k2,loff+1,m)
  430       continue
  440       continue
  450       continue
         endif
         if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 490 l = 1, kzp
            do 480 k = 1, kypd
            do 470 j = 1, nxv
            do 460 i = 1, 3
            cu3(i,j,k+kyp,l+loff,m) = scb(i,j,k,l,m)
  460       continue
  470       continue
  480       continue
  490       continue
         endif
      endif
  500 continue
c switch internal data
      if ((jm.le.kzb).and.(km.le.kyb)) then
         do 540 l = 1, kzp
         do 530 k = 1, kyp
         do 520 j = 1, nxv
         do 510 i = 1, 3
         at1 = cu3(i,j,k+kyp,l,m)
         cu3(i,j,k+kyp,l,m) = cu3(i,j,k,l+kzp,m)
         cu3(i,j,k,l+kzp,m) = at1
  510    continue
  520    continue
  530    continue
  540    continue
      endif
  550 continue
  560 continue
c copying in y direction not needed
      if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 860
c copy reflected data in y and z direction
      do 740 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 730 my = 1, k2blok
      m = my + moff
      loff = kzp2*(mz + ks)
      jm = (nz2 - loff - 1)/kzp + 1
      loff = kyp2*(my + js)
      km = (ny2 - loff - 1)/kyp + 1
      loff = kzp*(mz + ks)
      jl = (nz2 - loff - 1)/kzp2 + 1
      kz = loff + kzp2*jl - nz2
      loff = kyp*(my + js)
      kl = loff - kyp2*(loff/kyp2)
      kr = nvpy - (my + js) - 1
      do 680 ii = 7, 8
c for odd rows in y
      if (ii.eq.7) then
         mm = nvpy*jm + 2*kr
         lm = nvpy*jl + (nvpy + kr)/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in y
      else if (ii.eq.8) then
         if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 680
         mm = mm + 1
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kzb).and.(km.le.kyb)) then
c        if (kzp.gt.1) then
c           do 600 l = 2, kzp
c           do 590 k = 1, kypd
c           do 580 j = 1, nx
c           do 570 i = 1, 3
c           cu3(i,j,k,l+loff,m) = cu(i,j,k,l,mm-2*nvpy+1)
c 570       continue
c 580       continue
c 590       continue
c 600       continue
c        endif
c        if ((jm+1).le.kzb) then
c           do 630 k = 1, kypd
c           do 620 j = 1, nx
c           do 610 i = 1, 3
c           cu3(i,j,k,loff+1,m) = cu(i,j,k,1,mm+1)
c 610       continue
c 620       continue
c 630       continue
c        endif
c        if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
c           do 670 l = 1, kzp
c           do 660 k = 1, kypd
c           do 650 j = 1, nx
c           do 640 i = 1, 3
c           cu3(i,j,k+kyp,l+loff,m) = cu(i,j,k,l,mm-nvpy+1)
c 640       continue
c 650       continue
c 660       continue
c 670       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kzb).and.(km.le.kyb)) then
         if ((jm+1).le.kzb) then
            call MPI_IRECV(cu3(1,1,1,loff+1,m),3*nxvy,mreal,mm,noff+ii,l
     1grp,lsid,ierr)
         endif
         if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
            call MPI_IRECV(scb(1,1,1,1,m),3*kzp*nxvy,mreal,mm-nvpy,noff+
     1ii,lgrp,msid,ierr)
         endif
         if (kzp.gt.1) then
            call MPI_IRECV(cu3(1,1,1,loff+2,m),3*(kzp-1)*nxvy,mreal,mm-2
     1*nvpy,noff+ii,lgrp,nsid,ierr)
         endif
      endif
      if ((jl.gt.((kzb2-1)/2+1)).and.(jl.le.kzb2).and.lt1) then
         if (kz.eq.0) then
            if ((jl+1).le.kzb2) then
               call MPI_SEND(cu(1,1,1,1,m),3*nxvy,mreal,lm,noff+ii,lgrp,
     1ierr)
            endif
            if (kzp.gt.1) then
               call MPI_SEND(cu(1,1,1,2,m),3*(kzp-1)*nxvy,mreal,lm-nvpy,
     1noff+ii,lgrp,ierr)
            endif
         else
            call MPI_SEND(cu(1,1,1,1,m),3*kzp*nxvy,mreal,lm-nvpy,noff+ii
     1,lgrp,ierr)
         endif
      endif
c wait for data and unpack it
      if ((jm.le.kzb).and.(km.le.kyb)) then
         if (kzp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 600 l = 2, kzp
            l1 = kzp - l + 2
            l2 = (l1 + 2)/4 + 1
            koff = kypd*(l1 - 4*(l2 - 1) + 2)
            do 590 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1 + koff
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 580 j = 1, nxv
            do 570 i = 1, 3
            cu3(i,j,k1,l1+loff,m) = cu3(i,j+joff,k2,l2+loff,m)
  570       continue
  580       continue
  590       continue
  600       continue
         endif
         if ((jm+1).le.kzb) then
            call MPI_WAIT(lsid,istatus,ierr)
            do 630 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 620 j = 1, nxv
            do 610 i = 1, 3
            cu3(i,j,k1,loff+1,m) = cu3(i,j+joff,k2,loff+1,m)
  610       continue
  620       continue
  630       continue
         endif
         if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 670 l = 1, kzp
            do 660 k = 1, kypd
            do 650 j = 1, nxv
            do 640 i = 1, 3
            cu3(i,j,k+kyp,l+loff,m) = scb(i,j,k,l,m)
  640       continue
  650       continue
  660       continue
  670       continue
         endif
      endif
  680 continue
c switch internal data
      if ((jm.le.kzb).and.(km.le.kyb)) then
         do 720 l = 1, kzp
         do 710 k = 1, kyp
         do 700 j = 1, nxv
         do 690 i = 1, 3
         at1 = cu3(i,j,k+kyp,l,m)
         cu3(i,j,k+kyp,l,m) = cu3(i,j,k,l+kzp,m)
         cu3(i,j,k,l+kzp,m) = at1
  690    continue
  700    continue
  710    continue
  720    continue
      endif
  730 continue
  740 continue
c finish copy reflected data in y and z direction
      do 850 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 840 my = 1, k2blok
      m = my + moff
      loff = kyp2*(my + js)
      jm = (ny2 - loff - 1)/kyp + 1
      loff = kzp2*(mz + ks)
      km = (nz2 - loff - 1)/kzp + 1
      loff = kyp*(my + js)
      jl = (ny2 - loff - 1)/kyp2 + 1
      ky = loff + kyp2*jl - ny2
      loff = kzp*(mz + ks)
      kl = loff - kzp2*(loff/kzp2)
      kz = nvpy*(nvpz - (mz + ks) - 1)
      do 810 ii = 9, 10
c for odd rows in z
      if (ii.eq.9) then
         mm = jm + 2*kz
         lm = jl + (nvpy*(nvpz - 1) + kz)/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in z
      else if (ii.eq.10) then
         mm = mm + nvpy
         lm = lm + nvpy/2
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kyb).and.(km.le.kzb)) then
c        if ((jm+1).le.kyb) then
c           do 770 l = 1, kzp
c           do 760 j = 1, nx
c           do 750 i = 1, 3
c           cu3(i,j,1,l+loff,m) = cu(i,j,1,l,mm+1)
c 750       continue
c 760       continue
c 770       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if ((jm+1).le.kyb) then
            call MPI_IRECV(scd(1,1,1,2,m),3*kzp*nxv,mreal,mm,noff+ii,lgr
     1p,lsid,ierr)
         endif
      endif
      if ((jl.gt.((kyb2-1)/2+1)).and.(jl.le.kyb2).and.lt1) then
         if (ky.eq.0) then
            if ((jl+1).le.kyb2) then
               do 770 l = 1, kzp
               do 760 j = 1, nxv
               do 750 i = 1, 3
               scd(i,j,l,1,m) = cu(i,j,1,l,m)
  750          continue
  760          continue
  770          continue
               call MPI_SEND(scd(1,1,1,1,m),3*kzp*nxv,mreal,lm,noff+ii,l
     1grp,ierr)
            endif
         endif
      endif
c wait for data and unpack it
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if ((jm+1).le.kyb) then
            call MPI_WAIT(lsid,istatus,ierr)
            do 800 l = 1, kzp
            do 790 j = 1, nxv
            do 780 i = 1, 3
            cu3(i,j,1,l+loff,m) = scd(i,j,l,2,m)
  780       continue
  790       continue
  800       continue
         endif
      endif
  810 continue
c last point
      lm = jl + nvpy*(nvpz + mz + ks)/2
      mm = jm + nvpy*(2*(mz + ks) - nvpz)
      lt1 = (jl+1).le.kyb2
      loff = kzp*(mz + ks)
      jl = (nz2 - loff - 1)/kzp2 + 1
      kz = loff + kzp2*jl - nz2
      loff = kyp*(my + js)
      kl = loff - kyp2*(loff/kyp2)
c this segment is used for shared memory computers
c     if ((jm.le.kyb).and.(km.le.kzb)) then
c        if (((jm+1).le.kyb).and.((km+1).le.kzb)) then
c           do 830 j = 1, nx
c           do 820 i = 1, 3
c           cu3(i,j,1,1,m) = cu(i,j,1,1,mm+1)
c 820       continue
c 830       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if (((jm+1).le.kyb).and.((km+1).le.kzb)) then
            call MPI_IRECV(cu3(1,1,1,1,m),3*nxv,mreal,mm,noff+11,lgrp,ls
     1id,ierr)
         endif
      endif
      if ((jl.gt.((kzb2-1)/2+1)).and.(jl.le.kzb2).and.(kl.eq.0)) then
         if (kz.eq.0) then
            if (((jl+1).le.kzb2).and.lt1) then
               call MPI_SEND(cu(1,1,1,1,m),3*nxv,mreal,lm,noff+11,lgrp,i
     1err)
            endif
         endif
      endif
c wait for data
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if (((jm+1).le.kyb).and.((km+1).le.kzb)) then
            call MPI_WAIT(lsid,istatus,ierr)
         endif
      endif
  840 continue
  850 continue
c create odd array
  860 do 1130 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      loff = kzp2*(mz + ks)
      do 1120 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      do 1110 l = 1, kzp2
      ll = l + loff
      if ((ll.eq.1).or.(ll.eq.(nz+1))) then
         do 890 k = 1, kyp2
         do 880 j = 1, nx
         do 870 i = 1, 3
         cu3(i,j,k,l,m) = 0.
         cu3(i,j+nx,k,l,m) = 0.
  870    continue
  880    continue
  890    continue
      else if (ll.le.nz) then
         do 950 k = 1, kyp2
         kk = k + koff
         if ((kk.eq.1).or.(kk.eq.(ny+1))) then
            do 910 j = 1, nx
            do 900 i = 1, 3
            cu3(i,j,k,l,m) = 0.
            cu3(i,j+nx,k,l,m) = 0.
  900       continue
  910       continue
         else if (kk.le.ny) then
            do 920 j = 1, nxs
            cu3(1,nx+j+1,k,l,m) = cu3(1,nx-j+1,k,l,m)
            cu3(2,nx+j+1,k,l,m) = -cu3(2,nx-j+1,k,l,m)
            cu3(3,nx+j+1,k,l,m) = -cu3(3,nx-j+1,k,l,m)
  920       continue
            cu3(1,1,k,l,m) = 0.
            cu3(2,1,k,l,m) = 0.
            cu3(3,1,k,l,m) = 0.
            cu3(1,nx+1,k,l,m) = 0.
            cu3(2,nx+1,k,l,m) = 0.
            cu3(3,nx+1,k,l,m) = 0.
         else if (kk.gt.(ny+1)) then
            if (k.eq.1) then
               do 930 j = 1, nxs
               cu3(1,nx+j+1,k,l,m) = -cu3(1,nx-j+1,k,l,m)
               cu3(2,nx+j+1,k,l,m) = -cu3(2,nx-j+1,k,l,m)
               cu3(3,nx+j+1,k,l,m) = cu3(3,nx-j+1,k,l,m)
  930          continue
            else
               do 940 j = 1, nxs
               cu3(1,nx+j+1,k,l,m) = -cu3(1,nx-j+1,kyp2-k+2,l,m)
               cu3(2,nx+j+1,k,l,m) = -cu3(2,nx-j+1,kyp2-k+2,l,m)
               cu3(3,nx+j+1,k,l,m) = cu3(3,nx-j+1,kyp2-k+2,l,m)
  940          continue
            endif
            cu3(1,1,k,l,m) = 0.
            cu3(2,1,k,l,m) = 0.
            cu3(3,1,k,l,m) = 0.
            cu3(1,nx+1,k,l,m) = 0.
            cu3(2,nx+1,k,l,m) = 0.
            cu3(3,nx+1,k,l,m) = 0.
         endif
  950    continue
      else if (ll.gt.(nz+1)) then
         if (l.eq.1) then
            do 1010 k = 1, kyp2
            kk = k + koff
            if (kk.le.ny) then
               do 960 j = 1, nxs
               cu3(1,nx+j+1,k,l,m) = -cu3(1,nx-j+1,k,l,m)
               cu3(2,nx+j+1,k,l,m) = cu3(2,nx-j+1,k,l,m)
               cu3(3,nx+j+1,k,l,m) = -cu3(3,nx-j+1,k,l,m)
  960          continue
            else if (kk.gt.(ny+1)) then
               if (k.eq.1) then
                  do 980 j = 1, nxs
                  do 970 i = 1, 3
                  cu3(i,nx+j+1,k,l,m) = cu3(i,nx-j+1,k,l,m)
  970             continue
  980             continue
               else
                  do 1000 j = 1, nxs
                  do 990 i = 1, 3
                  cu3(i,nx+j+1,k,l,m) = cu3(i,nx-j+1,kyp2-k+2,l,m)
  990             continue
 1000             continue
               endif
            endif
 1010       continue
         else
            do 1070 k = 1, kyp2
            kk = k + koff
            if (kk.le.ny) then
               do 1020 j = 1, nxs
               cu3(1,nx+j+1,k,l,m) = -cu3(1,nx-j+1,k,kzp2-l+2,m)
               cu3(2,nx+j+1,k,l,m) = cu3(2,nx-j+1,k,kzp2-l+2,m)
               cu3(3,nx+j+1,k,l,m) = -cu3(3,nx-j+1,k,kzp2-l+2,m)
 1020          continue
            else if (kk.gt.(ny+1)) then
               if (k.eq.1) then
                  do 1040 j = 1, nxs
                  do 1030 i = 1, 3
                  cu3(i,nx+j+1,k,l,m) = cu3(i,nx-j+1,k,kzp2-l+2,m)
 1030             continue
 1040             continue
               else
                  do 1060 j = 1, nxs
                  do 1050 i = 1, 3
                  cu3(i,nx+j+1,k,l,m) = cu3(i,nx-j+1,kyp2-k+2,kzp2-l+2,m
     1)
 1050             continue
 1060             continue
               endif
            endif
 1070       continue
         endif
         do 1100 k = 1, kyp2
         kk = k + koff
         if ((kk.eq.1).or.(kk.eq.(ny+1))) then
            do 1090 j = 1, nx
            do 1080 i = 1, 3
            cu3(i,j,k,l,m) = 0.
            cu3(i,j+nx,k,l,m) = 0.
 1080       continue
 1090       continue
         else
            cu3(1,1,k,l,m) = 0.
            cu3(2,1,k,l,m) = 0.
            cu3(3,1,k,l,m) = 0.
            cu3(1,nx+1,k,l,m) = 0.
            cu3(2,nx+1,k,l,m) = 0.
            cu3(3,nx+1,k,l,m) = 0.
         endif
 1100    continue
      endif
 1110 continue
 1120 continue
 1130 continue
c finish odd array
      do 1210 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      loff = kzp2*(mz + ks)
      do 1200 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      do 1190 l = 1, kzp2
      ll = l + loff
      if (ll.le.nz) then
         do 1160 k = 1, kyp2
         kk = k + koff
         if (kk.gt.(ny+1)) then
            if (k.eq.1) then
               do 1140 j = 1, nxs
               cu3(1,nx-j+1,k,l,m) = cu3(1,nx+j+1,k,l,m)
               cu3(2,nx-j+1,k,l,m) = -cu3(2,nx+j+1,k,l,m)
               cu3(3,nx-j+1,k,l,m) = -cu3(3,nx+j+1,k,l,m)
 1140          continue
               cu3(1,nx+1,k,l,m) = cu3(1,nx+1,k,l,m)
               cu3(2,nx+1,k,l,m) = -cu3(2,nx+1,k,l,m)
               cu3(3,nx+1,k,l,m) = -cu3(3,nx+1,k,l,m)
            else
               do 1150 j = 1, nxs
               cu3(1,nx-j+1,k,l,m) = cu3(1,nx+j+1,k,l,m)
               cu3(2,nx-j+1,k,l,m) = -cu3(2,nx+j+1,k,l,m)
               cu3(3,nx-j+1,k,l,m) = -cu3(3,nx+j+1,k,l,m)
 1150          continue
               cu3(1,nx+1,k,l,m) = cu3(1,nx+1,k,l,m)
               cu3(2,nx+1,k,l,m) = -cu3(2,nx+1,k,l,m)
               cu3(3,nx+1,k,l,m) = -cu3(3,nx+1,k,l,m)
            endif
         endif
 1160    continue
      else if (ll.gt.(nz+1)) then
         do 1180 k = 1, kyp2
         kk = k + koff
         if ((kk.le.ny).or.(kk.gt.(ny+1))) then
            do 1170 j = 1, nxs
            cu3(1,nx-j+1,k,l,m) = cu3(1,nx+j+1,k,l,m)
            cu3(2,nx-j+1,k,l,m) = -cu3(2,nx+j+1,k,l,m)
            cu3(3,nx-j+1,k,l,m) = -cu3(3,nx+j+1,k,l,m)
 1170       continue
            cu3(1,nx+1,k,l,m) = cu3(1,nx+1,k,l,m)
            cu3(2,nx+1,k,l,m) = -cu3(2,nx+1,k,l,m)
            cu3(3,nx+1,k,l,m) = -cu3(3,nx+1,k,l,m)
         endif
 1180    continue
      endif
 1190 continue
 1200 continue
 1210 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PTRPSIN32D(q,q3,scb,scd,nx,ny,nz,kstrt,nvpy,nvpz,nxv,ky
     1p,kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
c this subroutine creates an odd array q3 from an array q, so that
c a 3d sine transform can be performed with a 3d real to complex fft.
c linear interpolation for distributed data with 2D domain decomposition
c q3 array may be modified
c scb/scd = scratch arrays
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = first dimension of input array q, must be >= nx
c kyp/kzp = number of data values per block in q in y/z
c kypd = second dimension of input array q, must be >= kyp
c kzpd = third dimension of input array q, must be >= kzp
c kyp2/kzp2 = number of data values per block in q3 in y/z
c kblok/lblok = number of data blocks in y/z
c k2blok/l2blok = number of data blocks in y/z for tripled data
      implicit none
      real q, q3, scb, scd
      integer nx, ny, nz, kstrt, nvpy, nvpz, nxv, kyp, kzp, kypd, kzpd
      integer kyp2, kzp2, kblok, lblok, k2blok, l2blok
      dimension q(nxv,kypd,kzpd,kblok*lblok)
      dimension q3(2*nxv,2*kypd,kzp2,k2blok*l2blok)
      dimension scb(nxv,kypd,kzpd,kblok*lblok)
      dimension scd(nxv,kzpd,2,kblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      logical lt1
      integer istatus, lsid, msid, nsid, ierr
      integer i, j, k, l, m, my, mz, nxs, nys, nzs, ks, js, ny2, nz2, kr
      integer kyb, kyb2, kzb, kzb2, nxvy, noff, moff, loff, koff, joff
      integer jm, km, jl, kl, ky, kz, kk, ll, mm, lm, l1, l2, k0, k1, k2
      real at1
      dimension istatus(lstat)
      nxs = nx - 1
      nys = ny - 1
      nzs = nz - 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      kyb = ny/kyp
      ny2 = ny + ny
      kyb2 = ny2/kyp2
      kzb = nz/kzp
      nz2 = nz + nz
      kzb2 = nz2/kzp2
      nxvy = nxv*kypd
      noff = kypd*kzpd + kyb*kzb
c copy to triple array in y direction
      do 150 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 140 my = 1, k2blok
      m = my + moff
      loff = kyp2*(my + js)
      jm = loff/kyp + 1
      loff = kzp2*(mz + ks)
      km = loff/kzp + 1
      loff = kyp*(my + js)
      jl = loff/kyp2 + 1
      loff = kzp*(mz + ks)
      kl = loff - kzp2*(loff/kzp2)
      kz = nvpy*(mz + ks)
c special case for one processor
      if ((kyb2*kzb2).eq.1) then
         do 40 l = 1, nzs
         do 20 k = 1, nys
         do 10 j = 1, nxs
         q3(j+1,k+1,l+1,m) = q(j+1,k+1,l+1,m)
         q3(nx+j+1,k+1,l+1,m) = -q(nx-j+1,k+1,l+1,m)
         q3(j+1,ny+k+1,l+1,m) = -q(j+1,ny-k+1,l+1,m)
         q3(nx+j+1,ny+k+1,l+1,m) = q(nx-j+1,ny-k+1,l+1,m)
         q3(j+1,k+1,nz+l+1,m) = -q(j+1,k+1,nz-l+1,m)
         q3(nx+j+1,k+1,nz+l+1,m) = q(nx-j+1,k+1,nz-l+1,m)
         q3(j+1,ny+k+1,nz+l+1,m) = q(j+1,ny-k+1,nz-l+1,m)
         q3(nx+j+1,ny+k+1,nz+l+1,m) = -q(nx-j+1,ny-k+1,nz-l+1,m)
   10    continue
         q3(1,k+1,l+1,m) = 0.
         q3(nx+1,k+1,l+1,m) = 0.
         q3(1,k+ny+1,l+1,m) = 0.
         q3(nx+1,k+ny+1,l+1,m) = 0.
         q3(1,k+1,nz+l+1,m) = 0.
         q3(nx+1,k+1,nz+l+1,m) = 0.
         q3(1,k+ny+1,nz+l+1,m) = 0.
         q3(nx+1,k+ny+1,nz+l+1,m) = 0.
   20    continue
         do 30 j = 1, nx
         q3(j,1,l+1,m) = 0.
         q3(j+nx,1,l+1,m) = 0.
         q3(j,ny+1,l+1,m) = 0.
         q3(j+nx,ny+1,l+1,m) = 0.
         q3(j,1,nz+l+1,m) = 0.
         q3(j+nx,1,nz+l+1,m) = 0.
         q3(j,ny+1,nz+l+1,m) = 0.
         q3(j+nx,ny+1,nz+l+1,m) = 0.
   30    continue
   40    continue
         nxs = nx + nx
         nys = ny + ny
         do 60 k = 1, nys
         do 50 j = 1, nxs
         q3(j,k,1,m) = 0.
         q3(j,k,nz+1,m) = 0.
   50    continue
   60    continue
         return
      endif
c copy main data in y direction
      do 130 i = 1, 2
c for odd rows in z
      if (i.eq.1) then
         mm = jm + 2*kz
         lm = jl + kz/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in z
      else if (i.eq.2) then
         if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 130
         mm = mm + nvpy
         lm = lm - nvpy/2
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kyb).and.(km.le.kzb)) then
c        do 90 l = 1, kzp
c        do 80 k = 1, kyp
c        do 70 j = 1, nx
c        q3(j,k,l+loff,m) = q(j,k,l,mm)
c  70    continue
c  80    continue
c  90    continue
c        if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
c           do 120 l = 1, kzp
c           do 110 k = 1, kyp
c           do 100 j = 1, nx
c           q3(j,k+kyp,l+loff,m) = q(j,k,l,mm+1)
c 100       continue
c 110       continue
c 120       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kyb).and.(km.le.kzb)) then
         call MPI_IRECV(q3(1,1,loff+1,m),kzp*nxvy,mreal,mm-1,noff+i,lgrp
     1,msid,ierr)
         if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
            call MPI_IRECV(scb(1,1,1,m),kzp*nxvy,mreal,mm,noff+i,lgrp,ns
     1id,ierr)
         endif
      endif
      if ((jl.le.((kyb2-1)/2+1)).and.lt1) then
         call MPI_SEND(q(1,1,1,m),kzp*nxvy,mreal,lm-1,noff+i,lgrp,ierr)
      endif
c wait for data and unpack it
      if ((jm.le.kyb).and.(km.le.kzb)) then
         call MPI_WAIT(msid,istatus,ierr)
         do 90 l = 1, kzp
         l1 = kzp - l + 1
         l2 = (l1 - 1)/4 + 1
         koff = kypd*(l1 - 4*(l2 - 1) - 1)
         do 80 k = 1, kypd
         k1 = kypd - k + 1
         k0 = k1 - 1 + koff
         k2 = k0/2 + 1
         joff = nxv*(k0 - 2*(k2 - 1))
         do 70 j = 1, nxv
         q3(j,k1,l1+loff,m) = q3(j+joff,k2,l2+loff,m)
   70    continue
   80    continue
   90    continue
         if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 120 l = 1, kzp
            do 110 k = 1, kypd
            do 100 j = 1, nxv
            q3(j,k+kyp,l+loff,m) = scb(j,k,l,m)
  100       continue
  110       continue
  120       continue
         endif
      endif
  130 continue
  140 continue
  150 continue
c copying in y direction not needed
      if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 290
c copy reflected data in y direction
      do 280 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 270 my = 1, k2blok
      m = my + moff
      loff = kyp2*(my + js)
      jm = (ny2 - loff - 1)/kyp + 1
      loff = kzp2*(mz + ks)
      km = loff/kzp + 1
      loff = kyp*(my + js)
      jl = (ny2 - loff - 1)/kyp2 + 1
      ky = loff + kyp2*jl - ny2
      loff = kzp*(mz + ks)
      kl = loff - kzp2*(loff/kzp2)
      kz = nvpy*(mz + ks)
      do 260 i = 3, 4
c for odd rows in z
      if (i.eq.3) then
         mm = jm + 2*kz
         lm = jl + kz/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in z
      else if (i.eq.4) then
         if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 260
         mm = mm + nvpy
         lm = lm - nvpy/2
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kyb).and.(km.le.kzb)) then
c        if ((jm+1).le.kyb) then
c           do 170 l = 1, kzp
c           do 160 j = 1, nx
c           q3(j,1,l+loff,m) = q(j,1,l,mm+1)
c 160       continue
c 170       continue
c        endif
c        if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
c           do 200 l = 1, kzp
c           do 190 k = 1, kyp
c           do 180 j = 1, nx
c           q3(j,k+kyp,l+loff,m) = q(j,k,l,mm)
c 180       continue
c 190       continue
c 200       continue
c        endif
c        if (kyp.gt.1) then
c           do 230 l = 1, kzp
c           do 220 k = 2, kyp
c           do 210 j = 1, nx
c           q3(j,k,l+loff,m) = q(j,k,l,mm-1)
c 210       continue
c 220       continue
c 230       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if ((jm+1).le.kyb) then
            call MPI_IRECV(scd(1,1,2,m),kzp*nxv,mreal,mm,noff+i,lgrp,lsi
     1d,ierr)
         endif
         if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
            call MPI_IRECV(scb(1,1,1,m),kzp*nxvy,mreal,mm-1,noff+i,lgrp,
     1msid,ierr)
         endif
         if (kyp.gt.1) then
            call MPI_IRECV(q3(1,1,loff+1,m),kzp*nxvy,mreal,mm-2,noff+i,l
     1grp,nsid,ierr)
         endif
      endif
      if ((jl.gt.((kyb2-1)/2+1)).and.(jl.le.kyb2).and.lt1) then
         if (ky.eq.0) then
            if ((jl+1).le.kyb2) then
               do 170 l = 1, kzp
               do 160 j = 1, nxv
               scd(j,l,1,m) = q(j,1,l,m)
  160          continue
  170          continue
               call MPI_SEND(scd(1,1,1,m),kzp*nxv,mreal,lm,noff+i,lgrp,i
     1err)
            endif
            if (kyp.gt.1) then
               call MPI_SEND(q(1,1,1,m),kzp*nxvy,mreal,lm-1,noff+i,lgrp,
     1ierr)
            endif
         else
            call MPI_SEND(q(1,1,1,m),kzp*nxvy,mreal,lm-1,noff+i,lgrp,ier
     1r)
         endif
      endif
c wait for data and unpack it
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if (kyp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 200 l = 1, kzp
            l1 = kzp - l + 1
            l2 = (l1 - 1)/4 + 1
            koff = kypd*(l1 - 4*(l2 - 1) - 1)
            do 190 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1 + koff
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 180 j = 1, nxv
            q3(j,k1,l1+loff,m) = q3(j+joff,k2,l2+loff,m)
  180       continue
  190       continue
  200       continue
         endif
         if ((jm+1).le.kyb) then
            call MPI_WAIT(lsid,istatus,ierr)
            do 220 l = 1, kzp
            do 210 j = 1, nxv
            q3(j,1,l+loff,m) = scd(j,l,2,m)
  210       continue
  220       continue
         endif
         if ((kyp.lt.kyp2).and.(kyb2.gt.1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 250 l = 1, kzp
            do 240 k = 1, kypd
            do 230 j = 1, nxv
            q3(j,k+kyp,l+loff,m) = scb(j,k,l,m)
  230       continue
  240       continue
  250       continue
         endif
      endif
  260 continue
  270 continue
  280 continue
c copying in z direction not needed
  290 if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 660
c copy reflected data in z direction
      do 430 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 420 my = 1, k2blok
      m = my + moff
      loff = kzp2*(mz + ks)
      jm = (nz2 - loff - 1)/kzp + 1
      loff = kyp2*(my + js)
      km = loff/kyp + 1
      loff = kzp*(mz + ks)
      jl = (nz2 - loff - 1)/kzp2 + 1
      kz = loff + kzp2*jl - nz2
      loff = kyp*(my + js)
      kl = loff - kyp2*(loff/kyp2)
      kr = my + js
      do 380 i = 5, 6
c for odd rows in y
      if (i.eq.5) then
         mm = nvpy*jm + 2*kr
         lm = nvpy*jl + kr/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in y
      else if (i.eq.6) then
         if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 380
         mm = mm + 1
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kzb).and.(km.le.kyb)) then
c        if (kzp.gt.1) then
c           do 320 l = 2, kzp
c           do 310 k = 1, kypd
c           do 300 j = 1, nx
c           q3(j,k,l+loff,m) = q(j,k,l,mm-2*nvpy+1)
c 300       continue
c 310       continue
c 320       continue
c        endif
c        if ((jm+1).le.kzb) then
c           do 340 k = 1, kypd
c           do 330 j = 1, nx
c           q3(j,k,loff+1,m) = q(j,k,1,mm+1)
c 330       continue
c 340       continue
c        endif
c        if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
c           do 370 l = 1, kzp
c           do 360 k = 1, kypd
c           do 350 j = 1, nx
c           q3(j,k+kyp,l+loff,m) = q(j,k,l,mm-nvpy+1)
c 350       continue
c 360       continue
c 370       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kzb).and.(km.le.kyb)) then
         if ((jm+1).le.kzb) then
            call MPI_IRECV(q3(1,1,loff+1,m),nxvy,mreal,mm,noff+i,lgrp,ls
     1id,ierr)
         endif
         if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
            call MPI_IRECV(scb(1,1,1,m),kzp*nxvy,mreal,mm-nvpy,noff+i,lg
     1rp,msid,ierr)
         endif
         if (kzp.gt.1) then
            call MPI_IRECV(q3(1,1,loff+2,m),(kzp-1)*nxvy,mreal,mm-2*nvpy
     1,noff+i,lgrp,nsid,ierr)
         endif
      endif
      if ((jl.gt.((kzb2-1)/2+1)).and.(jl.le.kzb2).and.lt1) then
         if (kz.eq.0) then
            if ((jl+1).le.kzb2) then
               call MPI_SEND(q(1,1,1,m),nxvy,mreal,lm,noff+i,lgrp,ierr)
            endif
            if (kzp.gt.1) then
               call MPI_SEND(q(1,1,2,m),(kzp-1)*nxvy,mreal,lm-nvpy,noff+
     1i,lgrp,ierr)
            endif
         else
            call MPI_SEND(q(1,1,1,m),kzp*nxvy,mreal,lm-nvpy,noff+i,lgrp,
     1ierr)
         endif
      endif
c wait for data and unpack it
      if ((jm.le.kzb).and.(km.le.kyb)) then
         if (kzp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 320 l = 2, kzp
            l1 = kzp - l + 2
            l2 = (l1 + 2)/4 + 1
            koff = kypd*(l1 - 4*(l2 - 1) + 2)
            do 310 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1 + koff
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 300 j = 1, nxv
            q3(j,k1,l1+loff,m) = q3(j+joff,k2,l2+loff,m)
  300       continue
  310       continue
  320       continue
         endif
         if ((jm+1).le.kzb) then
            call MPI_WAIT(lsid,istatus,ierr)
            do 340 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 330 j = 1, nxv
            q3(j,k1,loff+1,m) = q3(j+joff,k2,loff+1,m)
  330       continue
  340       continue
         endif
         if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 370 l = 1, kzp
            do 360 k = 1, kypd
            do 350 j = 1, nxv
            q3(j,k+kyp,l+loff,m) = scb(j,k,l,m)
  350       continue
  360       continue
  370       continue
         endif
      endif
  380 continue
c switch internal data
      if ((jm.le.kzb).and.(km.le.kyb)) then
         do 410 l = 1, kzp
         do 400 k = 1, kyp
         do 390 j = 1, nxv
         at1 = q3(j,k+kyp,l,m)
         q3(j,k+kyp,l,m) = q3(j,k,l+kzp,m)
         q3(j,k,l+kzp,m) = at1
  390    continue
  400    continue
  410    continue
      endif
  420 continue
  430 continue
c copying in y direction not needed
      if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 660
c copy reflected data in y and z direction
      do 570 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 560 my = 1, k2blok
      m = my + moff
      loff = kzp2*(mz + ks)
      jm = (nz2 - loff - 1)/kzp + 1
      loff = kyp2*(my + js)
      km = (ny2 - loff - 1)/kyp + 1
      loff = kzp*(mz + ks)
      jl = (nz2 - loff - 1)/kzp2 + 1
      kz = loff + kzp2*jl - nz2
      loff = kyp*(my + js)
      kl = loff - kyp2*(loff/kyp2)
      kr = nvpy - (my + js) - 1
      do 520 i = 7, 8
c for odd rows in y
      if (i.eq.7) then
         mm = nvpy*jm + 2*kr
         lm = nvpy*jl + (nvpy + kr)/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in y
      else if (i.eq.8) then
         if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 520
         mm = mm + 1
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kzb).and.(km.le.kyb)) then
c        if (kzp.gt.1) then
c           do 460 l = 2, kzp
c           do 450 k = 1, kypd
c           do 440 j = 1, nx
c           q3(j,k,l+loff,m) = q(j,k,l,mm-2*nvpy+1)
c 440       continue
c 450       continue
c 460       continue
c        endif
c        if ((jm+1).le.kzb) then
c           do 480 k = 1, kypd
c           do 470 j = 1, nx
c           q3(j,k,loff+1,m) = q(j,k,1,mm+1)
c 470       continue
c 480       continue
c        endif
c        if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
c           do 510 l = 1, kzp
c           do 500 k = 1, kypd
c           do 490 j = 1, nx
c           q3(j,k+kyp,l+loff,m) = q(j,k,l,mm-nvpy+1)
c 490       continue
c 500       continue
c 510       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kzb).and.(km.le.kyb)) then
         if ((jm+1).le.kzb) then
            call MPI_IRECV(q3(1,1,loff+1,m),nxvy,mreal,mm,noff+i,lgrp,ls
     1id,ierr)
         endif
         if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
            call MPI_IRECV(scb(1,1,1,m),kzp*nxvy,mreal,mm-nvpy,noff+i,lg
     1rp,msid,ierr)
         endif
         if (kzp.gt.1) then
            call MPI_IRECV(q3(1,1,loff+2,m),(kzp-1)*nxvy,mreal,mm-2*nvpy
     1,noff+i,lgrp,nsid,ierr)
         endif
      endif
      if ((jl.gt.((kzb2-1)/2+1)).and.(jl.le.kzb2).and.lt1) then
         if (kz.eq.0) then
            if ((jl+1).le.kzb2) then
               call MPI_SEND(q(1,1,1,m),nxvy,mreal,lm,noff+i,lgrp,ierr)
            endif
            if (kzp.gt.1) then
               call MPI_SEND(q(1,1,2,m),(kzp-1)*nxvy,mreal,lm-nvpy,noff+
     1i,lgrp,ierr)
            endif
         else
            call MPI_SEND(q(1,1,1,m),kzp*nxvy,mreal,lm-nvpy,noff+i,lgrp,
     1ierr)
         endif
      endif
c wait for data and unpack it
      if ((jm.le.kzb).and.(km.le.kyb)) then
         if (kzp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 460 l = 2, kzp
            l1 = kzp - l + 2
            l2 = (l1 + 2)/4 + 1
            koff = kypd*(l1 - 4*(l2 - 1) + 2)
            do 450 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1 + koff
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 440 j = 1, nxv
            q3(j,k1,l1+loff,m) = q3(j+joff,k2,l2+loff,m)
  440       continue
  450       continue
  460       continue
         endif
         if ((jm+1).le.kzb) then
            call MPI_WAIT(lsid,istatus,ierr)
            do 480 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 470 j = 1, nxv
            q3(j,k1,loff+1,m) = q3(j+joff,k2,loff+1,m)
  470       continue
  480       continue
         endif
         if ((kzp.lt.kzp2).and.(kzb2.gt.1)) then
            call MPI_WAIT(msid,istatus,ierr)
            do 510 l = 1, kzp
            do 500 k = 1, kypd
            do 490 j = 1, nxv
            q3(j,k+kyp,l+loff,m) = scb(j,k,l,m)
  490       continue
  500       continue
  510       continue
         endif
      endif
  520 continue
c switch internal data
      if ((jm.le.kzb).and.(km.le.kyb)) then
         do 550 l = 1, kzp
         do 540 k = 1, kyp
         do 530 j = 1, nxv
         at1 = q3(j,k+kyp,l,m)
         q3(j,k+kyp,l,m) = q3(j,k,l+kzp,m)
         q3(j,k,l+kzp,m) = at1
  530    continue
  540    continue
  550    continue
      endif
  560 continue
  570 continue
c finish copy reflected data in y and z direction
      do 650 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 640 my = 1, k2blok
      m = my + moff
      loff = kyp2*(my + js)
      jm = (ny2 - loff - 1)/kyp + 1
      loff = kzp2*(mz + ks)
      km = (nz2 - loff - 1)/kzp + 1
      loff = kyp*(my + js)
      jl = (ny2 - loff - 1)/kyp2 + 1
      ky = loff + kyp2*jl - ny2
      loff = kzp*(mz + ks)
      kl = loff - kzp2*(loff/kzp2)
      kz = nvpy*(nvpz - (mz + ks) - 1)
      do 620 i = 9, 10
c for odd rows in z
      if (i.eq.9) then
         mm = jm + 2*kz
         lm = jl + (nvpy*(nvpz - 1) + kz)/2
         loff = 0
         lt1 = kl.eq.0
c for even rows in z
      else if (i.eq.10) then
         mm = mm + nvpy
         lm = lm + nvpy/2
         loff = kzp
         lt1 = kl.ne.0
      endif
c this segment is used for shared memory computers
c     if ((jm.le.kyb).and.(km.le.kzb)) then
c        if ((jm+1).le.kyb) then
c           do 590 l = 1, kzp
c           do 580 j = 1, nx
c           q3(j,1,l+loff,m) = q(j,1,l,mm+1)
c 580       continue
c 590       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if ((jm+1).le.kyb) then
            call MPI_IRECV(scd(1,1,2,m),kzp*nxv,mreal,mm,noff+i,lgrp,lsi
     1d,ierr)
         endif
      endif
      if ((jl.gt.((kyb2-1)/2+1)).and.(jl.le.kyb2).and.lt1) then
         if (ky.eq.0) then
            if ((jl+1).le.kyb2) then
               do 590 l = 1, kzp
               do 580 j = 1, nxv
               scd(j,l,1,m) = q(j,1,l,m)
  580          continue
  590          continue
               call MPI_SEND(scd(1,1,1,m),kzp*nxv,mreal,lm,noff+i,lgrp,i
     1err)
            endif
         endif
      endif
c wait for data and unpack it
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if ((jm+1).le.kyb) then
            call MPI_WAIT(lsid,istatus,ierr)
            do 610 l = 1, kzp
            do 600 j = 1, nxv
            q3(j,1,l+loff,m) = scd(j,l,2,m)
  600       continue
  610       continue
         endif
      endif
  620 continue
c last point
      lm = jl + nvpy*(nvpz + mz + ks)/2
      mm = jm + nvpy*(2*(mz + ks) - nvpz)
      lt1 = (jl+1).le.kyb2
      loff = kzp*(mz + ks)
      jl = (nz2 - loff - 1)/kzp2 + 1
      kz = loff + kzp2*jl - nz2
      loff = kyp*(my + js)
      kl = loff - kyp2*(loff/kyp2)
c this segment is used for shared memory computers
c     if ((jm.le.kyb).and.(km.le.kzb)) then
c        if (((jm+1).le.kyb).and.((km+1).le.kzb)) then
c           do 630 j = 1, nx
c           q3(j,1,1,m) = q(j,1,1,mm+1)
c 630       continue
c        endif
c     endif
c this segment is used for mpi computers
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if (((jm+1).le.kyb).and.((km+1).le.kzb)) then
            call MPI_IRECV(q3(1,1,1,m),nxv,mreal,mm,noff+11,lgrp,lsid,ie
     1rr)
         endif
      endif
      if ((jl.gt.((kzb2-1)/2+1)).and.(jl.le.kzb2).and.(kl.eq.0)) then
         if (kz.eq.0) then
            if (((jl+1).le.kzb2).and.lt1) then
               call MPI_SEND(q(1,1,1,m),nxv,mreal,lm,noff+11,lgrp,ierr)
            endif
         endif
      endif
c wait for data
      if ((jm.le.kyb).and.(km.le.kzb)) then
         if (((jm+1).le.kyb).and.((km+1).le.kzb)) then
            call MPI_WAIT(lsid,istatus,ierr)
         endif
      endif
  640 continue
  650 continue
c create odd array
  660 do 860 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      loff = kzp2*(mz + ks)
      do 850 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      do 840 l = 1, kzp2
      ll = l + loff
      if ((ll.eq.1).or.(ll.eq.(nz+1))) then
         do 680 k = 1, kyp2
         do 670 j = 1, nx
         q3(j,k,l,m) = 0.
         q3(j+nx,k,l,m) = 0.
  670    continue
  680    continue
      else if (ll.le.nz) then
         do 730 k = 1, kyp2
         kk = k + koff
         if ((kk.eq.1).or.(kk.eq.(ny+1))) then
            do 690 j = 1, nx
            q3(j,k,l,m) = 0.
            q3(j+nx,k,l,m) = 0.
  690       continue
         else if (kk.le.ny) then
            do 700 j = 1, nxs
            q3(nx+j+1,k,l,m) = -q3(nx-j+1,k,l,m)
  700       continue
            q3(1,k,l,m) = 0.
            q3(nx+1,k,l,m) = 0.
         else if (kk.gt.(ny+1)) then
            if (k.eq.1) then
               do 710 j = 1, nxs
               q3(nx+j+1,k,l,m) = q3(nx-j+1,k,l,m)
  710          continue
            else
               do 720 j = 1, nxs
               q3(nx+j+1,k,l,m) = q3(nx-j+1,kyp2-k+2,l,m)
  720          continue
            endif
            q3(1,k,l,m) = 0.
            q3(nx+1,k,l,m) = 0.
         endif
  730    continue
      else if (ll.gt.(nz+1)) then
         if (l.eq.1) then
            do 770 k = 1, kyp2
            kk = k + koff
            if (kk.le.ny) then
               do 740 j = 1, nxs
               q3(nx+j+1,k,l,m) = q3(nx-j+1,k,l,m)
  740          continue
            else if (kk.gt.(ny+1)) then
               if (k.eq.1) then
                  do 750 j = 1, nxs
                  q3(nx+j+1,k,l,m) = -q3(nx-j+1,k,l,m)
  750             continue
               else
                  do 760 j = 1, nxs
                  q3(nx+j+1,k,l,m) = -q3(nx-j+1,kyp2-k+2,l,m)
  760             continue
               endif
            endif
  770       continue
         else
            do 810 k = 1, kyp2
            kk = k + koff
            if (kk.le.ny) then
               do 780 j = 1, nxs
               q3(nx+j+1,k,l,m) = q3(nx-j+1,k,kzp2-l+2,m)
  780          continue
            else if (kk.gt.(ny+1)) then
               if (k.eq.1) then
                  do 790 j = 1, nxs
                  q3(nx+j+1,k,l,m) = -q3(nx-j+1,k,kzp2-l+2,m)
  790             continue
               else
                  do 800 j = 1, nxs
                  q3(nx+j+1,k,l,m) = -q3(nx-j+1,kyp2-k+2,kzp2-l+2,m)
  800             continue
               endif
            endif
  810       continue
         endif
         do 830 k = 1, kyp2
         kk = k + koff
         if ((kk.eq.1).or.(kk.eq.(ny+1))) then
            do 820 j = 1, nx
            q3(j,k,l,m) = 0.
            q3(j+nx,k,l,m) = 0.
  820       continue
         else
            q3(1,k,l,m) = 0.
            q3(nx+1,k,l,m) = 0.
         endif
  830    continue
      endif
  840 continue
  850 continue
  860 continue
c finish odd array
      do 940 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      loff = kzp2*(mz + ks)
      do 930 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      do 920 l = 1, kzp2
      ll = l + loff
      if (ll.le.nz) then
         do 890 k = 1, kyp2
         kk = k + koff
         if (kk.gt.(ny+1)) then
            if (k.eq.1) then
               do 870 j = 1, nxs
               q3(nx-j+1,k,l,m) = -q3(nx+j+1,k,l,m)
  870          continue
               q3(nx+1,k,l,m) = -q3(nx+1,k,l,m)
            else
               do 880 j = 1, nxs
               q3(nx-j+1,k,l,m) = -q3(nx+j+1,k,l,m)
  880          continue
               q3(nx+1,k,l,m) = -q3(nx+1,k,l,m)
            endif
         endif
  890    continue
      else if (ll.gt.(nz+1)) then
         do 910 k = 1, kyp2
         kk = k + koff
         if ((kk.le.ny).or.(kk.gt.(ny+1))) then
            do 900 j = 1, nxs
            q3(nx-j+1,k,l,m) = -q3(nx+j+1,k,l,m)
  900       continue
            q3(nx+1,k,l,m) = -q3(nx+1,k,l,m)
         endif
  910    continue
      endif
  920 continue
  930 continue
  940 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine TRPSIN3C(cu,cu3,nx,ny,nz,nxv,nyv,nzv,nx2v,ny2,nz2)
c this subroutine creates a tripled vector array cu3 from a vector array
c cu, so that various 3d sine/cosine transforms can be performed with a
c 3d real to complex fft.  the x component is an odd function in y,
c y component is an odd function in x, and the z component is an odd
c function in both x and y.  Asummes vector cu vanishes at end points
c linear interpolation
c nx/ny/nz = system length in x/y/z direction
c nxv = second dimension of input array cu, must be >= nx
c nyv = third dimension of input array cu, must be >= ny
c nzv = fourth dimension of input array cu, must be >= nz
c nx2v = second dimension of output array cu3, must be >= 2*nx
c ny2 = third dimension of output array cu3, must be >= 2*ny
c nz2 = fourth dimension of output array cu3, must be >= 2*nz
      implicit none
      real cu, cu3
      integer nx, ny, nz, nxv, nyv, nzv, nx2v, ny2, nz2
      dimension cu(3,nxv,nyv,nzv), cu3(3,nx2v,ny2,nz2)
c local data
      integer i, j, k, l, nxs, nys, nzs
c copy to triple array
      nxs = nx - 1
      nys = ny - 1
      nzs = nz - 1
      do 60 l = 1, nzs
      do 30 k = 1, nys
      do 10 j = 1, nxs
      cu3(1,j+1,k+1,l+1) = cu(1,j+1,k+1,l+1)
      cu3(2,j+1,k+1,l+1) = cu(2,j+1,k+1,l+1)
      cu3(3,j+1,k+1,l+1) = cu(3,j+1,k+1,l+1)
      cu3(1,nx+j+1,k+1,l+1) = cu(1,nx-j+1,k+1,l+1)
      cu3(2,nx+j+1,k+1,l+1) = -cu(2,nx-j+1,k+1,l+1)
      cu3(3,nx+j+1,k+1,l+1) = -cu(3,nx-j+1,k+1,l+1)
      cu3(1,j+1,ny+k+1,l+1) = -cu(1,j+1,ny-k+1,l+1)
      cu3(2,j+1,ny+k+1,l+1) = cu(2,j+1,ny-k+1,l+1)
      cu3(3,j+1,ny+k+1,l+1) = -cu(3,j+1,ny-k+1,l+1)
      cu3(1,nx+j+1,ny+k+1,l+1) = -cu(1,nx-j+1,ny-k+1,l+1)
      cu3(2,nx+j+1,ny+k+1,l+1) = -cu(2,nx-j+1,ny-k+1,l+1)
      cu3(3,nx+j+1,ny+k+1,l+1) = cu(3,nx-j+1,ny-k+1,l+1)
      cu3(1,j+1,k+1,nz+l+1) = -cu(1,j+1,k+1,nz-l+1)
      cu3(2,j+1,k+1,nz+l+1) = -cu(2,j+1,k+1,nz-l+1)
      cu3(3,j+1,k+1,nz+l+1) = cu(3,j+1,k+1,nz-l+1)
      cu3(1,nx+j+1,k+1,nz+l+1) = -cu(1,nx-j+1,k+1,nz-l+1)
      cu3(2,nx+j+1,k+1,nz+l+1) = cu(2,nx-j+1,k+1,nz-l+1)
      cu3(3,nx+j+1,k+1,nz+l+1) = -cu(3,nx-j+1,k+1,nz-l+1)
      cu3(1,j+1,ny+k+1,nz+l+1) = cu(1,j+1,ny-k+1,nz-l+1)
      cu3(2,j+1,ny+k+1,nz+l+1) = -cu(2,j+1,ny-k+1,nz-l+1)
      cu3(3,j+1,ny+k+1,nz+l+1) = -cu(3,j+1,ny-k+1,nz-l+1)
      cu3(1,nx+j+1,ny+k+1,nz+l+1) = cu(1,nx-j+1,ny-k+1,nz-l+1)
      cu3(2,nx+j+1,ny+k+1,nz+l+1) = cu(2,nx-j+1,ny-k+1,nz-l+1)
      cu3(3,nx+j+1,ny+k+1,nz+l+1) = cu(3,nx-j+1,ny-k+1,nz-l+1)
   10 continue
      do 20 i = 1, 3
      cu3(i,1,k+1,l+1) = 0.
      cu3(i,nx+1,k+1,l+1) = 0.
      cu3(i,1,k+ny+1,l+1) = 0.
      cu3(i,nx+1,k+ny+1,l+1) = 0.
      cu3(i,1,k+1,nz+l+1) = 0.
      cu3(i,nx+1,k+1,nz+l+1) = 0.
      cu3(i,1,k+ny+1,nz+l+1) = 0.
      cu3(i,nx+1,k+ny+1,nz+l+1) = 0.
   20 continue
   30 continue
      do 50 j = 1, nx
      do 40 i = 1, 3
      cu3(i,j,1,l+1) = 0.
      cu3(i,j+nx,1,l+1) = 0.
      cu3(i,j,ny+1,l+1) = 0.
      cu3(i,j+nx,ny+1,l+1) = 0.
      cu3(i,j,1,nz+l+1) = 0.
      cu3(i,j+nx,1,nz+l+1) = 0.
      cu3(i,j,ny+1,nz+l+1) = 0.
      cu3(i,j+nx,ny+1,nz+l+1) = 0.
   40 continue
   50 continue
   60 continue
      nxs = nx + nx
      nys = ny + ny
      do 90 k = 1, nys
      do 80 j = 1, nxs
      do 70 i = 1, 3
      cu3(i,j,k,1) = 0.
      cu3(i,j,k,nz+1) = 0.
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine TRPSIN3D(q,q3,nx,ny,nz,nxv,nyv,nzv,nx2v,ny2,nz2)
c this subroutine creates an odd array q3 from an array q, so that
c a 3d sine transform can be performed with a 3d real to complex fft.
c linear interpolation
c nx/ny/nz = system length in x/y/z direction
c nxv = first dimension of input array q, must be >= nx
c nyv = second dimension of input array q, must be >= ny
c nzv = third dimension of input array q, must be >= nz
c nx2v = first dimension of output array q3, must be >= 2*nx
c ny2 = second dimension of output array q3, must be >= 2*ny
c nz2 = third dimension of output array q3, must be >= 2*nz
      implicit none
      real q, q3
      integer nx, ny, nz, nxv, nyv, nzv, nx2v, ny2, nz2
      dimension q(nxv,nyv,nzv), q3(nx2v,ny2,nz2)
c local data
      integer j, k, l, nxs, nys, nzs
c copy to triple array
      nxs = nx - 1
      nys = ny - 1
      nzs = nz - 1
      do 40 l = 1, nzs
      do 20 k = 1, nys
      do 10 j = 1, nxs
      q3(j+1,k+1,l+1) = q(j+1,k+1,l+1)
      q3(nx+j+1,k+1,l+1) = -q(nx-j+1,k+1,l+1)
      q3(j+1,ny+k+1,l+1) = -q(j+1,ny-k+1,l+1)
      q3(nx+j+1,ny+k+1,l+1) = q(nx-j+1,ny-k+1,l+1)
      q3(j+1,k+1,nz+l+1) = -q(j+1,k+1,nz-l+1)
      q3(nx+j+1,k+1,nz+l+1) = q(nx-j+1,k+1,nz-l+1)
      q3(j+1,ny+k+1,nz+l+1) = q(j+1,ny-k+1,nz-l+1)
      q3(nx+j+1,ny+k+1,nz+l+1) = -q(nx-j+1,ny-k+1,nz-l+1)
   10 continue
      q3(1,k+1,l+1) = 0.
      q3(nx+1,k+1,l+1) = 0.
      q3(1,k+ny+1,l+1) = 0.
      q3(nx+1,k+ny+1,l+1) = 0.
      q3(1,k+1,nz+l+1) = 0.
      q3(nx+1,k+1,nz+l+1) = 0.
      q3(1,k+ny+1,nz+l+1) = 0.
      q3(nx+1,k+ny+1,nz+l+1) = 0.
   20 continue
      do 30 j = 1, nx
      q3(j,1,l+1) = 0.
      q3(j+nx,1,l+1) = 0.
      q3(j,ny+1,l+1) = 0.
      q3(j+nx,ny+1,l+1) = 0.
      q3(j,1,nz+l+1) = 0.
      q3(j+nx,1,nz+l+1) = 0.
      q3(j,ny+1,nz+l+1) = 0.
      q3(j+nx,ny+1,nz+l+1) = 0.
   30 continue
   40 continue
      nxs = nx + nx
      nys = ny + ny
      do 60 k = 1, nys
      do 50 j = 1, nxs
      q3(j,k,1) = 0.
      q3(j,k,nz+1) = 0.
   50 continue
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDBLSIN32C(cu,cu2,scb,scd,nx,ny,kstrt,nvpy,nxv,kyp,kzp,
     1kypd,kzpd,kyp2,kblok,lblok,k2blok)
c this subroutine creates a doubled vector array cu2 from a vector array
c cu, so that various 2d sine/cosine transforms can be performed with a
c 3d real to complex fft.  the x component is an odd function in y,
c and y component is an odd function in x.
c Asummes vector cu vanishes at end points
c linear interpolation for distributed data
c cu2 array may be modified
c scb/scd = scratch arrays
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nvpy = number of real or virtual processors in y
c nxv = second dimension of input array cu, must be >= nx
c kyp/kzp = number of data values per block in cu in y/z
c kypd = third dimension of input array cu, must be >= kyp
c kzpd = fourth dimension of input array cu, must be >= kzp
c kyp2 = number of data values per block in cu3 in y
c kblok/lblok = number of data blocks in y/z
c k2blok = number of data blocks in y for doubled data
      implicit none
      real cu, cu2, scb, scd
      integer nx, ny, kstrt, nvpy, nxv, kyp, kzp, kypd, kzpd, kyp2
      integer kblok, lblok, k2blok
      dimension cu(3,nxv,kypd,kzpd,kblok*lblok)
      dimension cu2(3,2*nxv,2*kypd,kzp,k2blok*lblok)
      dimension scb(3,nxv,kypd,kzpd,kblok*lblok)
      dimension scd(3,nxv,kzpd,2,kblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, lsid, msid, nsid, ierr
      integer i, j, k, l, m, my, mz, nxs, nys, ny2, kyb, kyb2, ks, js
      integer nxvy, noff, moff, koff, joff, kz, ll, lm, mm, ml, kk
      integer l1, l2, k0, k1, k2
      dimension istatus(lstat)
      nxs = nx - 1
      nys = ny - 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      kyb = ny/kyp
      ny2 = ny + ny
      kyb2 = ny2/kyp2
      nxvy = nxv*kypd
      noff = kypd*kzpd + kyb
c copy to double array in x direction
      do 140 mz = 1, lblok
      moff = k2blok*(mz - 1)
      do 130 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      ll = koff/kyp + 1
      koff = kyp*(my + js)
      lm = koff/kyp2 + 1
      kz = nvpy*(mz + ks)
      mm = ll + kz
      ml = lm + kz
c special case for one processor in y direction
      if (kyb2.eq.1) then
         do 40 l = 1, kzp
         do 20 k = 1, nys
         do 10 j = 1, nxs
         cu2(1,j+1,k+1,l,m) = cu(1,j+1,k+1,l,m)
         cu2(2,j+1,k+1,l,m) = cu(2,j+1,k+1,l,m)
         cu2(3,j+1,k+1,l,m) = cu(3,j+1,k+1,l,m)
         cu2(1,nx+j+1,k+1,l,m) = cu(1,nx-j+1,k+1,l,m)
         cu2(2,nx+j+1,k+1,l,m) = -cu(2,nx-j+1,k+1,l,m)
         cu2(3,nx+j+1,k+1,l,m) = -cu(3,nx-j+1,k+1,l,m)
         cu2(1,j+1,ny+k+1,l,m) = -cu(1,j+1,ny-k+1,l,m)
         cu2(2,j+1,ny+k+1,l,m) = cu(2,j+1,ny-k+1,l,m)
         cu2(3,j+1,ny+k+1,l,m) = -cu(3,j+1,ny-k+1,l,m)
         cu2(1,nx+j+1,ny+k+1,l,m) = -cu(1,nx-j+1,ny-k+1,l,m)
         cu2(2,nx+j+1,ny+k+1,l,m) = -cu(2,nx-j+1,ny-k+1,l,m)
         cu2(3,nx+j+1,ny+k+1,l,m) = cu(3,nx-j+1,ny-k+1,l,m)
   10    continue
         cu2(1,1,k+1,l,m) = 0.
         cu2(2,1,k+1,l,m) = 0.
         cu2(3,1,k+1,l,m) = 0.
         cu2(1,nx+1,k+1,l,m) = 0.
         cu2(2,nx+1,k+1,l,m) = 0.
         cu2(3,nx+1,k+1,l,m) = 0.
         cu2(1,1,k+ny+1,l,m) = 0.
         cu2(2,1,k+ny+1,l,m) = 0.
         cu2(3,1,k+ny+1,l,m) = 0.
         cu2(1,nx+1,k+ny+1,l,m) = 0.
         cu2(2,nx+1,k+ny+1,l,m) = 0.
         cu2(3,nx+1,k+ny+1,l,m) = 0.
   20    continue
         do 30 j = 1, nx
         cu2(1,j,1,l,m) = 0.
         cu2(2,j,1,l,m) = 0.
         cu2(3,j,1,l,m) = 0.
         cu2(1,j+nx,1,l,m) = 0.
         cu2(2,j+nx,1,l,m) = 0.
         cu2(3,j+nx,1,l,m) = 0.
         cu2(1,j,ny+1,l,m) = 0.
         cu2(2,j,ny+1,l,m) = 0.
         cu2(3,j,ny+1,l,m) = 0.
         cu2(1,j+nx,ny+1,l,m) = 0.
         cu2(2,j+nx,ny+1,l,m) = 0.
         cu2(3,j+nx,ny+1,l,m) = 0.
   30    continue
   40    continue
         return
      endif
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        do 80 l = 1, kzp
c        do 70 k = 1, kyp
c        do 60 j = 1, nx
c        do 50 i = 1, 3
c        cu2(i,j,k,l,m) = cu(i,j,k,l,mm)
c  50    continue
c  60    continue
c  70    continue
c  80   continue
c        if (kyp.lt.kyp2) then
c           do 120 l = 1, kzp
c           do 110 k = 1, kyp
c           do 100 j = 1, nx
c           do 90 i = 1, 3
c           cu2(i,j,k+kyp,l,m) = cu(i,j,k,l,mm+1)
c  90       continue
c 100       continue
c 110       continue
c 120       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         call MPI_IRECV(cu2(1,1,1,1,m),3*kzp*nxvy,mreal,mm-1,noff+1,lgrp
     1,msid,ierr)
         if (kyp.lt.kyp2) then
            call MPI_IRECV(scb(1,1,1,1,m),3*kzp*nxvy,mreal,mm,noff+1,lgr
     1p,nsid,ierr)
         endif
      endif
      if (lm.le.(kyb2/2)) then
         call MPI_SEND(cu(1,1,1,1,m),3*kzp*nxvy,mreal,ml-1,noff+1,lgrp,i
     1err)
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         call MPI_WAIT(msid,istatus,ierr)
         do 80 l = 1, kzp
         l1 = kzp - l + 1
         l2 = (l1 - 1)/4 + 1
         koff = kypd*(l1 - 4*(l2 - 1) - 1)
         do 70 k = 1, kypd
         k1 = kypd - k + 1
         k0 = k1 - 1 + koff
         k2 = k0/2 + 1
         joff = nxv*(k0 - 2*(k2 - 1))
         do 60 j = 1, nxv
         do 50 i = 1, 3
         cu2(i,j,k1,l1,m) = cu2(i,j+joff,k2,l2,m)
   50    continue
   60    continue
   70    continue
   80    continue
         if (kyp.lt.kyp2) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 120 l = 1, kzp
            do 110 k = 1, kypd
            do 100 j = 1, nxv
            do 90 i = 1, 3
            cu2(i,j,k+kyp,l,m) = scb(i,j,k,l,m)
   90       continue
  100       continue
  110       continue
  120       continue
         endif
      endif
  130 continue
  140 continue
c copy to double array in y direction
      do 360 mz = 1, lblok
      moff = k2blok*(mz - 1)
      do 350 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      ll = (ny2 - koff - 1)/kyp + 1
      koff = kyp*(my + js)
      lm = (ny2 - koff - 1)/kyp2 + 1
      koff = koff + kyp2*lm - ny2
      kz = nvpy*(mz + ks)
      mm = ll + kz
      ml = lm + kz
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        if ((ll+1).le.kyb) then
c           do 170 l = 1, kzp
c           do 160 j = 1, nx
c           do 150 i = 1, 3
c           cu2(i,j,1,l,m) = cu(i,j,1,l,mm+1)
c 150       continue
c 160       continue
c 170       continue
c        endif
c        if (kyp.lt.kyp2) then
c           do 210 l = 1, kzp
c           do 200 k = 1, kyp
c           do 190 j = 1, nx
c           do 180 i = 1, 3
c           cu2(i,j,k+kyp,l,m) = cu(i,j,k,l,mm)
c 180       continue
c 190       continue
c 200       continue
c 210       continue
c        endif
c        if (kyp.gt.1) then
c           do 250 l = 1, kzp
c           do 240 k = 2, kyp
c           do 230 j = 1, nx
c           do 220 i = 1, 3
c           cu2(i,j,k,l,m) = cu(i,j,k,l,mm-1)
c 220       continue
c 230       continue
c 240       continue
c 250       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         if ((ll+1).le.kyb) then
            call MPI_IRECV(scd(1,1,1,2,m),3*kzp*nxv,mreal,mm,noff+2,lgrp
     1,lsid,ierr)
         endif
         if (kyp.lt.kyp2) then
            call MPI_IRECV(scb(1,1,1,1,m),3*kzp*nxvy,mreal,mm-1,noff+2,l
     1grp,msid,ierr)
         endif
         if (kyp.gt.1) then
            call MPI_IRECV(cu2(1,1,1,1,m),3*kzp*nxvy,mreal,mm-2,noff+2,l
     1grp,nsid,ierr)
         endif
      endif
      if ((lm.gt.(kyb2/2)).and.(lm.le.kyb2)) then
         if (koff.eq.0) then
            if ((lm+1).le.kyb2) then
               do 170 l = 1, kzp
               do 160 j = 1, nxv
               do 150 i = 1, 3
               scd(i,j,l,1,m) = cu(i,j,1,l,m)
  150          continue
  160          continue
  170          continue
               call MPI_SEND(scd(1,1,1,1,m),3*kzp*nxv,mreal,ml,noff+2,lg
     1rp,ierr)
            endif
            if (kyp.gt.1) then
               call MPI_SEND(cu(1,1,1,1,m),3*kzp*nxvy,mreal,ml-1,noff+2,
     1lgrp,ierr)
            endif
         else
            call MPI_SEND(cu(1,1,1,1,m),3*kzp*nxvy,mreal,ml-1,noff+2,lgr
     1p,ierr)
         endif
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         if (kyp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 210 l = 1, kzp
            l1 = kzp - l + 1
            l2 = (l1 - 1)/4 + 1
            koff = kypd*(l1 - 4*(l2 - 1) - 1)
            do 200 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1 + koff
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 190 j = 1, nxv
            do 180 i = 1, 3
            cu2(i,j,k1,l1,m) = cu2(i,j+joff,k2,l2,m)
  180       continue
  190       continue
  200       continue
  210       continue
         endif
         if ((ll+1).le.kyb) then
            call MPI_WAIT(lsid,istatus,ierr)
            do 240 l = 1, kzp
            do 230 j = 1, nxv
            do 220 i = 1, 3
            cu2(i,j,1,l,m) = scd(i,j,l,2,m)
  220       continue
  230       continue
  240       continue
         endif
         if (kyp.lt.kyp2) then
            call MPI_WAIT(msid,istatus,ierr)
            do 340 l = 1, kzp
            do 330 k = 1, kypd
            do 320 j = 1, nxv
            do 310 i = 1, 3
            cu2(i,j,k+kyp,l,m) = scb(i,j,k,l,m)
  310       continue
  320       continue
  330       continue
  340       continue
         endif
      endif
  350 continue
  360 continue
c create odd array
      do 470 mz = 1, lblok
      moff = k2blok*(mz - 1)
      do 460 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      do 450 l = 1, kzp
      do 440 k = 1, kyp2
      kk = k + koff
      if ((kk.eq.1).or.(kk.eq.(ny+1))) then
         do 380 j = 1, nx
         do 370 i = 1, 3
         cu2(i,j,k,l,m) = 0.
         cu2(i,j+nx,k,l,m) = 0.
  370    continue
  380    continue
      else if (kk.le.ny) then
         do 390 j = 1, nxs
         cu2(1,nx+j+1,k,l,m) = cu2(1,nx-j+1,k,l,m)
         cu2(2,nx+j+1,k,l,m) = -cu2(2,nx-j+1,k,l,m)
         cu2(3,nx+j+1,k,l,m) = -cu2(3,nx-j+1,k,l,m)
  390    continue
         do 400 i = 1, 3
         cu2(i,1,k,l,m) = 0.
         cu2(i,nx+1,k,l,m) = 0.
  400    continue
      else if (kk.gt.(ny+1)) then
         if (k.eq.1) then
            do 410 j = 1, nxs
            cu2(1,nx+j+1,k,l,m) = -cu2(1,nx-j+1,k,l,m)
            cu2(2,nx+j+1,k,l,m) = -cu2(2,nx-j+1,k,l,m)
            cu2(3,nx+j+1,k,l,m) = cu2(3,nx-j+1,k,l,m)
  410       continue
         else
            do 420 j = 1, nxs
            cu2(1,nx+j+1,kyp2-k+2,l,m) = -cu2(1,nx-j+1,k,l,m)
            cu2(2,nx+j+1,kyp2-k+2,l,m) = -cu2(2,nx-j+1,k,l,m)
            cu2(3,nx+j+1,kyp2-k+2,l,m) = cu2(3,nx-j+1,k,l,m)
  420       continue
         endif
         do 430 i = 1, 3
         cu2(i,1,k,l,m) = 0.
         cu2(i,nx+1,k,l,m) = 0.
  430    continue
      endif
  440 continue
  450 continue
  460 continue
  470 continue
c finish odd array
      do 520 mz = 1, lblok
      moff = k2blok*(mz - 1)
      do 510 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      do 500 l = 1, kzp
      do 490 k = 1, kyp2
      kk = k + koff
      if (kk.gt.(ny+1)) then
         do 480 j = 1, nxs
         cu2(1,nx-j+1,k,l,m) = cu2(1,nx+j+1,k,l,m)
         cu2(2,nx-j+1,k,l,m) = -cu2(2,nx+j+1,k,l,m)
         cu2(3,nx-j+1,k,l,m) = -cu2(3,nx+j+1,k,l,m)
  480    continue
         cu2(1,nx+1,k,l,m) = cu2(1,nx+1,k,l,m)
         cu2(2,nx+1,k,l,m) = -cu2(2,nx+1,k,l,m)
         cu2(3,nx+1,k,l,m) = -cu2(3,nx+1,k,l,m)
      endif
  490 continue
  500 continue
  510 continue
  520 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PDBLSIN32D(q,q2,scb,scd,nx,ny,kstrt,nvpy,nxv,kyp,kzp,ky
     1pd,kzpd,kyp2,kblok,lblok,k2blok)
c this subroutine creates an odd array q2 from an array q, so that
c a 2d sine transform can be performed with a 3d real to complex fft.
c linear interpolation for distributed data with 2D domain decomposition
c Asummes q vanishes at end points
c linear interpolation for distributed data
c q2 array may be modified
c scb/scd = scratch arrays
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nvpy = number of real or virtual processors in y
c nxv = first dimension of input array q, must be >= nx
c kyp/kzp = number of data values per block in cu in y/z
c kypd = second dimension of input array q, must be >= kyp
c kzpd = third dimension of input array q, must be >= kzp
c kyp2 = number of data values per block in q2 in y
c kblok/lblok = number of data blocks in y/z
c k2blok = number of data blocks in y for doubled data
      implicit none
      real q, q2, scb, scd
      integer nx, ny, kstrt, nvpy, nxv, kyp, kzp, kypd, kzpd, kyp2
      integer kblok, lblok, k2blok
      dimension q(nxv,kypd,kzpd,kblok*lblok)
      dimension q2(2*nxv,2*kypd,kzp,k2blok*lblok)
      dimension scb(nxv,kypd,kzpd,kblok*lblok)
      dimension scd(nxv,kzpd,2,kblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, lsid, msid, nsid, ierr
      integer j, k, l, m, my, mz, nxs, nys, ny2, kyb, kyb2, ks, js
      integer nxvy, noff, moff, koff, joff, kz, ll, lm, mm, ml, kk
      integer l1, l2, k0, k1, k2
      dimension istatus(lstat)
      nxs = nx - 1
      nys = ny - 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      kyb = ny/kyp
      ny2 = ny + ny
      kyb2 = ny2/kyp2
      nxvy = nxv*kypd
      noff = kypd*kzpd + kyb
c copy to double array in x direction
      do 120 mz = 1, lblok
      moff = k2blok*(mz - 1)
      do 110 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      ll = koff/kyp + 1
      koff = kyp*(my + js)
      lm = koff/kyp2 + 1
      kz = nvpy*(mz + ks)
      mm = ll + kz
      ml = lm + kz
c special case for one processor in y direction
      if (kyb2.eq.1) then
         do 40 l = 1, kzp
         do 20 k = 1, nys
         do 10 j = 1, nxs
         q2(j+1,k+1,l,m) = q(j+1,k+1,l,m)
         q2(nx+j+1,k+1,l,m) = -q(nx-j+1,k+1,l,m)
         q2(j+1,ny+k+1,l,m) = -q(j+1,ny-k+1,l,m)
         q2(nx+j+1,ny+k+1,l,m) = q(nx-j+1,ny-k+1,l,m)
   10    continue
         q2(1,k+1,l,m) = 0.
         q2(nx+1,k+1,l,m) = 0.
         q2(1,k+ny+1,l,m) = 0.
         q2(nx+1,k+ny+1,l,m) = 0.
   20    continue
         do 30 j = 1, nx
         q2(j,1,l,m) = 0.
         q2(j+nx,1,l,m) = 0.
         q2(j,ny+1,l,m) = 0.
         q2(j+nx,ny+1,l,m) = 0.
   30    continue
   40    continue
         return
      endif
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        do 70 l = 1, kzp
c        do 60 k = 1, kyp
c        do 50 j = 1, nx
c        q2(j,k,l,m) = q(j,k,l,mm)
c  50    continue
c  60    continue
c  70   continue
c        if (kyp.lt.kyp2) then
c           do 100 l = 1, kzp
c           do 90 k = 1, kyp
c           do 80 j = 1, nx
c           q2(j,k+kyp,l,m) = q(j,k,l,mm+1)
c  80       continue
c  90       continue
c 100       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         call MPI_IRECV(q2(1,1,1,m),kzp*nxvy,mreal,mm-1,noff+1,lgrp,msid
     1,ierr)
         if (kyp.lt.kyp2) then
            call MPI_IRECV(scb(1,1,1,m),kzp*nxvy,mreal,mm,noff+1,lgrp,ns
     1id,ierr)
         endif
      endif
      if (lm.le.(kyb2/2)) then
         call MPI_SEND(q(1,1,1,m),kzp*nxvy,mreal,ml-1,noff+1,lgrp,ierr)
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         call MPI_WAIT(msid,istatus,ierr)
         do 70 l = 1, kzp
         l1 = kzp - l + 1
         l2 = (l1 - 1)/4 + 1
         koff = kypd*(l1 - 4*(l2 - 1) - 1)
         do 60 k = 1, kypd
         k1 = kypd - k + 1
         k0 = k1 - 1 + koff
         k2 = k0/2 + 1
         joff = nxv*(k0 - 2*(k2 - 1))
         do 50 j = 1, nxv
         q2(j,k1,l1,m) = q2(j+joff,k2,l2,m)
   50    continue
   60    continue
   70    continue
         if (kyp.lt.kyp2) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 100 l = 1, kzp
            do 90 k = 1, kypd
            do 80 j = 1, nxv
            q2(j,k+kyp,l,m) = scb(j,k,l,m)
   80       continue
   90       continue
  100       continue
         endif
      endif
  110 continue
  120 continue
c copy to double array in y direction
      do 240 mz = 1, lblok
      moff = k2blok*(mz - 1)
      do 230 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      ll = (ny2 - koff - 1)/kyp + 1
      koff = kyp*(my + js)
      lm = (ny2 - koff - 1)/kyp2 + 1
      koff = koff + kyp2*lm - ny2
      kz = nvpy*(mz + ks)
      mm = ll + kz
      ml = lm + kz
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        if ((ll+1).le.kyb) then
c           do 140 l = 1, kzp
c           do 130 j = 1, nx
c           q2(j,1,l,m) = q(j,1,l,mm+1)
c 130       continue
c 140       continue
c        endif
c        if (kyp.lt.kyp2) then
c           do 170 l = 1, kzp
c           do 160 k = 1, kyp
c           do 150 j = 1, nx
c           q2(j,k+kyp,l,m) = q(j,k,l,mm)
c 150       continue
c 160       continue
c 170       continue
c        endif
c        if (kyp.gt.1) then
c           do 200 l = 1, kzp
c           do 190 k = 2, kyp
c           do 180 j = 1, nx
c           q2(j,k,l,m) = q(j,k,l,mm-1)
c 180       continue
c 190       continue
c 200       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         if ((ll+1).le.kyb) then
            call MPI_IRECV(scd(1,1,2,m),kzp*nxv,mreal,mm,noff+2,lgrp,lsi
     1d,ierr)
         endif
         if (kyp.lt.kyp2) then
            call MPI_IRECV(scb(1,1,1,m),kzp*nxvy,mreal,mm-1,noff+2,lgrp,
     1msid,ierr)
         endif
         if (kyp.gt.1) then
            call MPI_IRECV(q2(1,1,1,m),kzp*nxvy,mreal,mm-2,noff+2,lgrp,n
     1sid,ierr)
         endif
      endif
      if ((lm.gt.(kyb2/2)).and.(lm.le.kyb2)) then
         if (koff.eq.0) then
            if ((lm+1).le.kyb2) then
               do 140 l = 1, kzp
               do 130 j = 1, nxv
               scd(j,l,1,m) = q(j,1,l,m)
  130          continue
  140          continue
               call MPI_SEND(scd(1,1,1,m),kzp*nxv,mreal,ml,noff+2,lgrp,i
     1err)
            endif
            if (kyp.gt.1) then
               call MPI_SEND(q(1,1,1,m),kzp*nxvy,mreal,ml-1,noff+2,lgrp,
     1ierr)
            endif
         else
            call MPI_SEND(q(1,1,1,m),kzp*nxvy,mreal,ml-1,noff+2,lgrp,ier
     1r)
         endif
      endif
c wait for data and unpack it
      if (ll.le.kyb) then
         if (kyp.gt.1) then
            call MPI_WAIT(nsid,istatus,ierr)
            do 170 l = 1, kzp
            l1 = kzp - l + 1
            l2 = (l1 - 1)/4 + 1
            koff = kypd*(l1 - 4*(l2 - 1) - 1)
            do 160 k = 1, kypd
            k1 = kypd - k + 1
            k0 = k1 - 1 + koff
            k2 = k0/2 + 1
            joff = nxv*(k0 - 2*(k2 - 1))
            do 150 j = 1, nxv
            q2(j,k1,l1,m) = q2(j+joff,k2,l2,m)
  150       continue
  160       continue
  170       continue
         endif
         if ((ll+1).le.kyb) then
            call MPI_WAIT(lsid,istatus,ierr)
            do 190 l = 1, kzp
            do 180 j = 1, nxv
            q2(j,1,l,m) = scd(j,l,2,m)
  180       continue
  190       continue
         endif
         if (kyp.lt.kyp2) then
            call MPI_WAIT(msid,istatus,ierr)
            do 220 l = 1, kzp
            do 210 k = 1, kypd
            do 200 j = 1, nxv
            q2(j,k+kyp,l,m) = scb(j,k,l,m)
  200       continue
  210       continue
  220       continue
         endif
      endif
  230 continue
  240 continue
c create odd array
      do 320 mz = 1, lblok
      moff = k2blok*(mz - 1)
      do 310 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      do 300 l = 1, kzp
      do 290 k = 1, kyp2
      kk = k + koff
      if ((kk.eq.1).or.(kk.eq.(ny+1))) then
         do 250 j = 1, nx
         q2(j,k,l,m) = 0.
         q2(j+nx,k,l,m) = 0.
  250    continue
      else if (kk.le.ny) then
         do 260 j = 1, nxs
         q2(nx+j+1,k,l,m) = -q2(nx-j+1,k,l,m)
  260    continue
         q2(1,k,l,m) = 0.
         q2(nx+1,k,l,m) = 0.
      else if (kk.gt.(ny+1)) then
         if (k.eq.1) then
            do 270 j = 1, nxs
            q2(nx+j+1,k,l,m) = q2(nx-j+1,k,l,m)
  270       continue
         else
            do 280 j = 1, nxs
            q2(nx+j+1,kyp2-k+2,l,m) = q2(nx-j+1,k,l,m)
  280       continue
         endif
         q2(1,k,l,m) = 0.
         q2(nx+1,k,l,m) = 0.
      endif
  290 continue
  300 continue
  310 continue
  320 continue
c finish odd array
      do 370 mz = 1, lblok
      moff = k2blok*(mz - 1)
      do 360 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      do 350 l = 1, kzp
      do 340 k = 1, kyp2
      kk = k + koff
      if (kk.gt.(ny+1)) then
         do 330 j = 1, nxs
         q2(nx-j+1,k,l,m) = -q2(nx+j+1,k,l,m)
  330    continue
         q2(nx+1,k,l,m) = -q2(nx+1,k,l,m)
      endif
  340 continue
  350 continue
  360 continue
  370 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DBLSIN3C(cu,cu2,nx,ny,nz,nxv,nyv,nzv,nx2v,ny2)
c this subroutine creates a doubled vector array cu2 from a vector array
c cu, so that various 2d sine/cosine transforms can be performed with a
c 3d real to complex fft.  the x component is an odd function in y,
c y component is an odd function in x, and the z component is an odd
c function in both x and y.  Asummes vector cu vanishes at end points
c linear interpolation
c nx/ny/nz = system length in x/y/z direction
c nxv = second dimension of input array cu, must be >= nx
c nyv = third dimension of input array cu, must be >= ny
c nzv = fourth dimension of input array cu, must be >= nz
c nx2v = second dimension of output array cu2, must be >= 2*nx
c ny2 = third dimension of output array cu2, must be >= 2*ny
      implicit none
      real cu, cu2
      integer nx, ny, nz, nxv, nyv, nzv, nx2v, ny2
      dimension cu(3,nxv,nyv,nzv), cu2(3,nx2v,ny2,nzv)
c local data
      integer j, k, l, nxs, nys
c copy to double array
      nxs = nx - 1
      nys = ny - 1
      do 40 l = 1, nz
      do 20 k = 1, nys
      do 10 j = 1, nxs
      cu2(1,j+1,k+1,l) = cu(1,j+1,k+1,l)
      cu2(2,j+1,k+1,l) = cu(2,j+1,k+1,l)
      cu2(3,j+1,k+1,l) = cu(3,j+1,k+1,l)
      cu2(1,nx+j+1,k+1,l) = cu(1,nx-j+1,k+1,l)
      cu2(2,nx+j+1,k+1,l) = -cu(2,nx-j+1,k+1,l)
      cu2(3,nx+j+1,k+1,l) = -cu(3,nx-j+1,k+1,l)
      cu2(1,j+1,ny+k+1,l) = -cu(1,j+1,ny-k+1,l)
      cu2(2,j+1,ny+k+1,l) = cu(2,j+1,ny-k+1,l)
      cu2(3,j+1,ny+k+1,l) = -cu(3,j+1,ny-k+1,l)
      cu2(1,nx+j+1,ny+k+1,l) = -cu(1,nx-j+1,ny-k+1,l)
      cu2(2,nx+j+1,ny+k+1,l) = -cu(2,nx-j+1,ny-k+1,l)
      cu2(3,nx+j+1,ny+k+1,l) = cu(3,nx-j+1,ny-k+1,l)
   10 continue
      cu2(1,1,k+1,l) = 0.
      cu2(2,1,k+1,l) = 0.
      cu2(3,1,k+1,l) = 0.
      cu2(1,nx+1,k+1,l) = 0.
      cu2(2,nx+1,k+1,l) = 0.
      cu2(3,nx+1,k+1,l) = 0.
      cu2(1,1,k+ny+1,l) = 0.
      cu2(2,1,k+ny+1,l) = 0.
      cu2(3,1,k+ny+1,l) = 0.
      cu2(1,nx+1,k+ny+1,l) = 0.
      cu2(2,nx+1,k+ny+1,l) = 0.
      cu2(3,nx+1,k+ny+1,l) = 0.
   20 continue
      do 30 j = 1, nx
      cu2(1,j,1,l) = 0.
      cu2(2,j,1,l) = 0.
      cu2(3,j,1,l) = 0.
      cu2(1,j+nx,1,l) = 0.
      cu2(2,j+nx,1,l) = 0.
      cu2(3,j+nx,1,l) = 0.
      cu2(1,j,ny+1,l) = 0.
      cu2(2,j,ny+1,l) = 0.
      cu2(3,j,ny+1,l) = 0.
      cu2(1,j+nx,ny+1,l) = 0.
      cu2(2,j+nx,ny+1,l) = 0.
      cu2(3,j+nx,ny+1,l) = 0.
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine DBLSIN3D(q,q2,nx,ny,nz,nxv,nyv,nzv,nx2v,ny2)
c this subroutine creates an odd array q2 from an array q, so that
c a 2d sine transform can be performed with a 3d real to complex fft.
c linear interpolation
c nx/ny/nz = system length in x/y/z direction
c nxv = first dimension of input array q, must be >= nx
c nyv = second dimension of input array q, must be >= ny
c nzv = third dimension of input array q, must be >= nz
c nx2v = first dimension of output array q2, must be >= 2*nx
c ny2 = second dimension of output array q2, must be >= 2*ny
      implicit none
      real q, q2
      integer nx, ny, nz, nxv, nyv, nzv, nx2v, ny2
      dimension q(nxv,nyv,nzv), q2(nx2v,ny2,nzv)
c local data
      integer j, k, l, nxs, nys
c copy to double array
      nxs = nx - 1
      nys = ny - 1
      do 40 l = 1, nz
      do 20 k = 1, nys
      do 10 j = 1, nxs
      q2(j+1,k+1,l) = q(j+1,k+1,l)
      q2(nx+j+1,k+1,l) = -q(nx-j+1,k+1,l)
      q2(j+1,ny+k+1,l) = -q(j+1,ny-k+1,l)
      q2(nx+j+1,ny+k+1,l) = q(nx-j+1,ny-k+1,l)
   10 continue
      q2(1,k+1,l) = 0.
      q2(nx+1,k+1,l) = 0.
      q2(1,k+ny+1,l) = 0.
      q2(nx+1,k+ny+1,l) = 0.
   20 continue
      do 30 j = 1, nx
      q2(j,1,l) = 0.
      q2(j+nx,1,l) = 0.
      q2(j,ny+1,l) = 0.
      q2(j+nx,ny+1,l) = 0.
   30 continue
   40 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PHAFTRP32C(fxyz, fxyz3,scb,scd,nx,ny,nz,kstrt,nvpy,nvpz
     1,nxv,kyp,kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
c this subroutine copies data from a triple array to regular array
c with guard cells for vector field and linear interpolation
c for distributed data with 2D domain decomposition
c fxyz array may be modified
c scb/scd = scratch arrays
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = second dimension of input array cu, must be >= nx
c kyp/kzp = number of data values per block in cu in y/z
c kypd = third dimension of input array cu, must be >= kyp
c kzpd = fourth dimension of input array cu, must be >= kzp
c kyp2/kzp2 = number of data values per block in cu3 in y/z
c kblok/lblok = number of data blocks in y/z
c k2blok/l2blok = number of data blocks in y/z for tripled data
      implicit none
      real fxyz, fxyz3, scb, scd
      integer nx, ny, nz, kstrt, nvpy, nvpz, nxv, kyp, kzp, kypd, kzpd
      integer kyp2, kzp2, kblok, lblok, k2blok, l2blok
      dimension fxyz(3,nxv,kypd,kzpd,kblok*lblok)
      dimension fxyz3(3,2*nxv,2*kypd,kzp2,k2blok*l2blok)
      dimension scb(3,nxv,kypd,kzpd,kblok*lblok)
      dimension scd(3,nxv,kzpd,2,kblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, lsid, msid, nsid, ierr
      integer i, j, k, l, m, my, mz, nx1, ny1, nz1, kzp1, ks, js, kr, ls
      integer kyb, kyb2, kzb, kzb2, nxvy, noff, moff, loff, koff, jm, km
      integer jl, kl, kz, kk, mm, lm, l2, k0, k1, k2, is, ii
      dimension istatus(lstat)
      nx1 = nx + 1
      ny1 = ny + 1
      nz1 = nz + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      kyb = ny/kyp
      kyb2 = (ny + ny)/kyp2
      kzb = nz/kzp
      kzb2 = (nz + nz)/kzp2
      kzp1 = kzp + 1
      nxvy = nxv*kypd
      noff = kypd*kzpd + kyb*kzb
c copy to triple array in y direction
      do 370 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 360 my = 1, k2blok
      m = my + moff
      loff = kyp2*(my + js)
      jl = loff/kyp + 1
      loff = kzp2*(mz + ks)
      kl = loff/kzp + 1
      loff = kyp*(my + js)
      jm = loff/kyp2 + 1
      koff = loff - kyp2*(jm - 1)
      loff = kzp*(mz + ks)
      kr = loff/kzp2
      km = loff - kzp2*kr
      mm = jm + nvpy*kr
      kz = nvpy*(mz + ks)
      lm = jl + 2*kz
c special case for one processor
      if ((kyb2*kzb2).eq.1) then
         do 40 l = 1, nz1
         do 30 k = 1, ny1
         do 20 j = 1, nx1
         do 10 i = 1, 3
         fxyz(i,j,k,l,m) = fxyz3(i,j,k,l,m)
   10    continue
   20    continue
   30    continue
   40    continue
         return
      endif
      if (km.eq.0) then
         ls = kzp1
         is = 1
      else
         ls = kzp
         is = 2
      endif
c this segment is used for shared memory computers
c     if (jm.le.kyb) then
c        if ((koff.eq.0).and.(kyp.lt.kyp2)) then
c           do 80 l = 1, kzp1
c           do 70 k = 1, kyp1
c           do 60 j = 1, nx1
c           do 50 i = 1, 3
c           fxyz(i,j,k,l,m) = fxyz3(i,j,k,l,mm)
c  50       continue
c  60       continue
c  70       continue
c  80       continue
c        else
c           do 120 l = 1, kzp1
c           do 110 k = 1, kyp
c           do 100 j = 1, nx1
c           do 90 i = 1, 3
c           fxyz(i,j,k,l,m) = fxyz3(i,j,k+koff,l,mm)
c  90       continue
c 100       continue
c 110       continue
c 120       continue
c           do 150 l = 1, kzp1
c           do 140 j = 1, nx1
c           do 130 i = 1, 3
c           fxyz(i,j,kyp+1,l,m) = fxyz3(i,j,1,l,mm+1)
c 130       continue
c 140       continue
c 150       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jm.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_IRECV(fxyz(1,1,1,1,m),3*ls*nxvy,mreal,mm-1,noff+is,
     1lgrp,msid,ierr)
         else
            call MPI_IRECV(fxyz(1,1,1,1,m),3*ls*nxvy,mreal,mm-1,noff+is,
     1lgrp,msid,ierr)
            call MPI_IRECV(scd(1,1,1,2,m),3*ls*nxv,mreal,mm,noff+is,lgrp
     1,nsid,ierr)
         endif
      endif
      do 230 ii = 1, 2
c for odd rows in z
      if (ii.eq.1) then
         loff = 0
         ls = kzp1
         lm = jl + 2*kz
c for even rows in z
      else if (ii.eq.2) then
         if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 230
         loff = kzp
         ls = kzp
         lm = lm + nvpy
      endif
c pack data and send it
      if ((jl.le.kyb).and.(kl.le.kzb)) then
         if (kyp.lt.kyp2) then
            do 80 l = 1, ls
            do 70 k = 1, kypd
            do 60 j = 1, nxv
            do 50 i = 1, 3
            scb(i,j,k,l,m) = fxyz3(i,j,k,l+loff,m)
   50       continue
   60       continue
   70       continue
   80       continue
            call MPI_SEND(scb(1,1,1,1,m),3*ls*nxvy,mreal,lm-1,noff+ii,lg
     1rp,ierr)
            if (kyb2.gt.1) then
               do 120 l = 1, ls
               do 110 k = 1, kypd
               do 100 j = 1, nxv
               do 90 i = 1, 3
               scb(i,j,k,l,m) = fxyz3(i,j,k+kyp,l+loff,m)
   90          continue
  100          continue
  110          continue
  120          continue
               call MPI_SEND(scb(1,1,1,1,m),3*ls*nxvy,mreal,lm,noff+ii,l
     1grp,ierr)
            endif
         else
            do 160 l = 1, ls
            do 150 k = 1, kypd
            do 140 j = 1, nxv
            do 130 i = 1, 3
            scb(i,j,k,l,m) = fxyz3(i,j,k,l+loff,m)
  130       continue
  140       continue
  150       continue
  160       continue
            call MPI_SEND(scb(1,1,1,1,m),3*ls*nxvy,mreal,lm-1,noff+ii,lg
     1rp,ierr)
         endif
         if (jl.gt.1) then
            do 190 l = 1, ls
            do 180 j = 1, nxv
            do 170 i = 1, 3
            scd(i,j,l,1,m) = fxyz3(i,j,1,l+loff,m)
  170       continue
  180       continue
  190       continue
            call MPI_SEND(scd(1,1,1,1,m),3*ls*nxv,mreal,lm-2,noff+ii,lgr
     1p,ierr)
         endif
      else if ((jl.eq.(kyb+1)).and.(kl.le.kzb)) then
         do 220 l = 1, ls
         do 210 j = 1, nxv
         do 200 i = 1, 3
         scd(i,j,l,1,m) = fxyz3(i,j,1,l+loff,m)
  200    continue
  210    continue
  220    continue
         call MPI_SEND(scd(1,1,1,1,m),3*ls*nxv,mreal,lm-2,noff+ii,lgrp,i
     1err)
      endif
  230 continue
c wait for data
      if (km.eq.0) then
         ls = kzp1
      else
         ls = kzp
      endif
      if (jm.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_WAIT(msid,istatus,ierr)
         else
            call MPI_WAIT(msid,istatus,ierr)
            call MPI_WAIT(nsid,istatus,ierr)
            do 260 l = 1, ls
            do 250 j = 1, nxv
            do 240 i = 1, 3
            fxyz(i,j,kyp+1,l,m) = scd(i,j,l,2,m)
  240       continue
  250       continue
  260       continue
         endif
      endif
c last point
      if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 360
      mm = mm + nvpy - 1
      if (koff.eq.0) then
         is = 3
      else
         is = 4
      endif
c this segment is used for shared memory computers
c     if ((kr.lt.kzb).and.(km.ne.0)) then
c        do 290 k = 1, kypd
c        do 280 j = 1, nx1
c        do 270 i = 1, 3
c        fxyz(i,j,k,kzp+1,m) = fxyz3(i,j,k,1,mm+1)
c 270    continue
c 280    continue
c 290    continue
c     endif
c this segment is used for mpi computers
      if ((kr.lt.kzb).and.(km.ne.0)) then
         call MPI_IRECV(fxyz(1,1,1,kzp+1,m),3*nxvy,mreal,mm,noff+is,lgrp
     1,nsid,ierr)
      endif
      do 330 ii = 3, 4
c for odd rows in z
      if (ii.eq.3) then
         loff = 0
         lm = jl + 2*kz - nvpy + 1
c for even rows in z
      else if (ii.eq.4) then
         if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 330
         loff = kyp
         lm = lm + 1
      endif
c pack data and send it
      if (jl.le.kyb) then
         if ((kl.le.kzb).and.(kl.gt.1)) then
            do 290 k = 1, kypd
            do 280 j = 1, nxv
            do 270 i = 1, 3
            scb(i,j,k,1,m) = fxyz3(i,j,k+loff,1,m)
  270       continue
  280       continue
  290       continue
            call MPI_SEND(scb(1,1,1,1,m),3*nxvy,mreal,lm-2,noff+ii,lgrp,
     1ierr)
         else if (kl.eq.(kzb+1)) then
            do 320 k = 1, kypd
            do 310 j = 1, nxv
            do 300 i = 1, 3
            scb(i,j,k,1,m) = fxyz3(i,j,k+loff,1,m)
  300       continue
  310       continue
  320       continue
            call MPI_SEND(scb(1,1,1,1,m),3*nxvy,mreal,lm-2,noff+ii,lgrp,
     1ierr)
         endif
      endif
  330 continue
c wait for data
      if ((kr.lt.kzb).and.(km.ne.0)) then
         call MPI_WAIT(nsid,istatus,ierr)
      endif
c last point
      if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 360
      mm = mm + 1
      lm = jl + 2*kz - nvpy
c this segment is used for shared memory computers
c     if ((kr.lt.kzb).and.(km.ne.0).and.(koff.ne.0)) then
c        do 350 j = 1, nx1
c        do 340 i = 1, 3
c        fxyz(i,j,kyp+1,kzp+1,m) = fxyz3(i,j,1,1,mm+1)
c 340    continue
c 350    continue
c     endif
c this segment is used for mpi computers
      if ((kr.lt.kzb).and.(km.ne.0).and.(koff.ne.0)) then
         call MPI_IRECV(fxyz(1,1,kyp+1,kzp+1,m),3*nxv,mreal,mm,noff+5,lg
     1rp,nsid,ierr)
      endif
c send data
      if ((jl.le.(kyb+1)).and.(jl.gt.1)) then
         if ((kl.le.(kzb+1)).and.(kl.gt.1)) then
            call MPI_SEND(fxyz3(1,1,1,1,m),3*nxv,mreal,lm-2,noff+5,lgrp,
     1ierr)
         endif
      endif
c wait for data
      if ((kr.lt.kzb).and.(km.ne.0).and.(koff.ne.0)) then
         call MPI_WAIT(nsid,istatus,ierr)
      endif
  360 continue
  370 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PHAFTRP32D(q,q3,scb,scd,nx,ny,nz,kstrt,nvpy,nvpz,nxv,ky
     1p,kzp,kypd,kzpd,kyp2,kzp2,kblok,lblok,k2blok,l2blok)
c this subroutine copies data from a triple array to regular array
c with guard cells for scalar field and linear interpolation
c for distributed data with 2D domain decomposition
c q array may be modified
c scb/scd = scratch arrays
c nx/ny/nz = system length in x/y/z direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = second dimension of input array cu, must be >= nx
c kyp/kzp = number of data values per block in cu in y/z
c kypd = third dimension of input array cu, must be >= kyp
c kzpd = fourth dimension of input array cu, must be >= kzp
c kyp2/kzp2 = number of data values per block in cu3 in y/z
c kblok/lblok = number of data blocks in y/z
c k2blok/l2blok = number of data blocks in y/z for tripled data
      implicit none
      real q, q3, scb, scd
      integer nx, ny, nz, kstrt, nvpy, nvpz, nxv, kyp, kzp, kypd, kzpd
      integer kyp2, kzp2, kblok, lblok, k2blok, l2blok
      dimension q(nxv,kypd,kzpd,kblok*lblok)
      dimension q3(2*nxv,2*kypd,kzp2,k2blok*l2blok)
      dimension scb(nxv,kypd,kzpd,kblok*lblok)
      dimension scd(nxv,kzpd,2,kblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, lsid, msid, nsid, ierr
      integer j, k, l, m, my, mz, nx1, ny1, nz1, kzp1, ks, js, kr, ls
      integer kyb, kyb2, kzb, kzb2, nxvy, noff, moff, loff, koff, jm, km
      integer jl, kl, kz, kk, mm, lm, l2, k0, k1, k2, is, ii
      dimension istatus(lstat)
      nx1 = nx + 1
      ny1 = ny + 1
      nz1 = nz + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      kyb = ny/kyp
      kyb2 = (ny + ny)/kyp2
      kzb = nz/kzp
      kzb2 = (nz + nz)/kzp2
      kzp1 = kzp + 1
      nxvy = nxv*kypd
      noff = kypd*kzpd + kyb*kzb
c copy to triple array in y direction
      do 270 mz = 1, l2blok
      moff = k2blok*(mz - 1)
      do 260 my = 1, k2blok
      m = my + moff
      loff = kyp2*(my + js)
      jl = loff/kyp + 1
      loff = kzp2*(mz + ks)
      kl = loff/kzp + 1
      loff = kyp*(my + js)
      jm = loff/kyp2 + 1
      koff = loff - kyp2*(jm - 1)
      loff = kzp*(mz + ks)
      kr = loff/kzp2
      km = loff - kzp2*kr
      mm = jm + nvpy*kr
      kz = nvpy*(mz + ks)
      lm = jl + 2*kz
c special case for one processor
      if ((kyb2*kzb2).eq.1) then
         do 30 l = 1, nz1
         do 20 k = 1, ny1
         do 10 j = 1, nx1
         q(j,k,l,m) = q3(j,k,l,m)
   10    continue
   20    continue
   30    continue
         return
      endif
      if (km.eq.0) then
         ls = kzp1
         is = 1
      else
         ls = kzp
         is = 2
      endif
c this segment is used for shared memory computers
c     if (jm.le.kyb) then
c        if ((koff.eq.0).and.(kyp.lt.kyp2)) then
c           do 60 l = 1, kzp1
c           do 50 k = 1, kyp1
c           do 40 j = 1, nx1
c           q(j,k,l,m) = q3(j,k,l,mm)
c  40       continue
c  50       continue
c  60       continue
c        else
c           do 90 l = 1, kzp1
c           do 80 k = 1, kyp
c           do 70 j = 1, nx1
c           q(j,k,l,m) = q3(j,k+koff,l,mm)
c  70       continue
c  80       continue
c  90       continue
c           do 110 l = 1, kzp1
c           do 100 j = 1, nx1
c           q(j,kyp+1,l,m) = q3(j,1,l,mm+1)
c 100       continue
c 110       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (jm.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_IRECV(q(1,1,1,m),ls*nxvy,mreal,mm-1,noff+is,lgrp,ms
     1id,ierr)
         else
            call MPI_IRECV(q(1,1,1,m),ls*nxvy,mreal,mm-1,noff+is,lgrp,ms
     1id,ierr)
            call MPI_IRECV(scd(1,1,2,m),ls*nxv,mreal,mm,noff+is,lgrp,nsi
     1d,ierr)
         endif
      endif
      do 170 ii = 1, 2
c for odd rows in z
      if (ii.eq.1) then
         loff = 0
         ls = kzp1
         lm = jl + 2*kz
c for even rows in z
      else if (ii.eq.2) then
         if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 170
         loff = kzp
         ls = kzp
         lm = lm + nvpy
      endif
c pack data and send it
      if ((jl.le.kyb).and.(kl.le.kzb)) then
         if (kyp.lt.kyp2) then
            do 60 l = 1, ls
            do 50 k = 1, kypd
            do 40 j = 1, nxv
            scb(j,k,l,m) = q3(j,k,l+loff,m)
   40       continue
   50       continue
   60       continue
            call MPI_SEND(scb(1,1,1,m),ls*nxvy,mreal,lm-1,noff+ii,lgrp,i
     1err)
            if (kyb2.gt.1) then
               do 90 l = 1, ls
               do 80 k = 1, kypd
               do 70 j = 1, nxv
               scb(j,k,l,m) = q3(j,k+kyp,l+loff,m)
   70          continue
   80          continue
   90          continue
               call MPI_SEND(scb(1,1,1,m),ls*nxvy,mreal,lm,noff+ii,lgrp,
     1ierr)
            endif
         else
            do 120 l = 1, ls
            do 110 k = 1, kypd
            do 100 j = 1, nxv
            scb(j,k,l,m) = q3(j,k,l+loff,m)
  100       continue
  110       continue
  120       continue
            call MPI_SEND(scb(1,1,1,m),ls*nxvy,mreal,lm-1,noff+ii,lgrp,i
     1err)
         endif
         if (jl.gt.1) then
            do 140 l = 1, ls
            do 130 j = 1, nxv
            scd(j,l,1,m) = q3(j,1,l+loff,m)
  130       continue
  140       continue
            call MPI_SEND(scd(1,1,1,m),ls*nxv,mreal,lm-2,noff+ii,lgrp,ie
     1rr)
         endif
      else if ((jl.eq.(kyb+1)).and.(kl.le.kzb)) then
         do 160 l = 1, ls
         do 150 j = 1, nxv
         scd(j,l,1,m) = q3(j,1,l+loff,m)
  150    continue
  160    continue
         call MPI_SEND(scd(1,1,1,m),ls*nxv,mreal,lm-2,noff+ii,lgrp,ierr)
      endif
  170 continue
c wait for data
      if (km.eq.0) then
         ls = kzp1
      else
         ls = kzp
      endif
      if (jm.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_WAIT(msid,istatus,ierr)
         else
            call MPI_WAIT(msid,istatus,ierr)
            call MPI_WAIT(nsid,istatus,ierr)
            do 190 l = 1, ls
            do 180 j = 1, nxv
            q(j,kyp+1,l,m) = scd(j,l,2,m)
  180       continue
  190       continue
         endif
      endif
c last point
      if ((kzp.eq.kzp2).or.(kzb2.eq.1)) go to 260
      mm = mm + nvpy - 1
      if (koff.eq.0) then
         is = 3
      else
         is = 4
      endif
c this segment is used for shared memory computers
c     if ((kr.lt.kzb).and.(km.ne.0)) then
c        do 210 k = 1, kypd
c        do 200 j = 1, nx1
c        q(j,k,kzp+1,m) = q3(j,k,1,mm+1)
c 200    continue
c 210    continue
c     endif
c this segment is used for mpi computers
      if ((kr.lt.kzb).and.(km.ne.0)) then
         call MPI_IRECV(q(1,1,kzp+1,m),nxvy,mreal,mm,noff+is,lgrp,nsid,i
     1err)
      endif
      do 240 ii = 3, 4
c for odd rows in z
      if (ii.eq.3) then
         loff = 0
         lm = jl + 2*kz - nvpy + 1
c for even rows in z
      else if (ii.eq.4) then
         if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 240
         loff = kyp
         lm = lm + 1
      endif
c pack data and send it
      if (jl.le.kyb) then
         if ((kl.le.kzb).and.(kl.gt.1)) then
            do 210 k = 1, kypd
            do 200 j = 1, nxv
            scb(j,k,1,m) = q3(j,k+loff,1,m)
  200       continue
  210       continue
            call MPI_SEND(scb(1,1,1,m),nxvy,mreal,lm-2,noff+ii,lgrp,ierr
     1)
         else if (kl.eq.(kzb+1)) then
            do 230 k = 1, kypd
            do 220 j = 1, nxv
            scb(j,k,1,m) = q3(j,k+loff,1,m)
  220       continue
  230       continue
            call MPI_SEND(scb(1,1,1,m),nxvy,mreal,lm-2,noff+ii,lgrp,ierr
     1)
         endif
      endif
  240 continue
c wait for data
      if ((kr.lt.kzb).and.(km.ne.0)) then
         call MPI_WAIT(nsid,istatus,ierr)
      endif
c last point
      if ((kyp.eq.kyp2).or.(kyb2.eq.1)) go to 260
      mm = mm + 1
      lm = jl + 2*kz - nvpy
c this segment is used for shared memory computers
c     if ((kr.lt.kzb).and.(km.ne.0).and.(koff.ne.0)) then
c        do 250 j = 1, nx1
c        q(j,kyp+1,kzp+1,m) = q3(j,1,1,mm+1)
c 250    continue
c     endif
c this segment is used for mpi computers
      if ((kr.lt.kzb).and.(km.ne.0).and.(koff.ne.0)) then
         call MPI_IRECV(q(1,kyp+1,kzp+1,m),nxv,mreal,mm,noff+5,lgrp,nsid
     1,ierr)
      endif
c send data
      if ((jl.le.(kyb+1)).and.(jl.gt.1)) then
         if ((kl.le.(kzb+1)).and.(kl.gt.1)) then
            call MPI_SEND(q3(1,1,1,m),nxv,mreal,lm-2,noff+5,lgrp,ierr)
         endif
      endif
c wait for data
      if ((kr.lt.kzb).and.(km.ne.0).and.(koff.ne.0)) then
         call MPI_WAIT(nsid,istatus,ierr)
      endif
  260 continue
  270 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine HAFTRP3C(fxyz,fxyz3,nx,ny,nz,nxe,nye,nze,nx2v,ny2,nz2)
c this subroutine copies data from a triple array to regular array
c with guard cells for vector field and linear interpolation
c nx/ny/nz = system length in x/y/z direction
c nxe = first dimension of output array fxy, must be >= nx+1
c nye = second dimension of ouput array fxy, must be >= ny+1
c nze = third dimension of ouput array fxy, must be >= nz+1
c nx2v = first dimension of input array fxy3, must be >= 2*nx
c ny2 = second dimension of input array fxy3, must be >= 2*ny
c nz2 = third dimension of input array fxy3, must be >= 2*nz
      implicit none
      real fxyz, fxyz3
      integer nx, ny, nz, nxe, nye, nze, nx2v, ny2, nz2
      dimension fxyz(3,nxe,nye,nze), fxyz3(3,nx2v,ny2,nz2)
c local data
      integer j, k, l, nx1, ny1, nz1
      nx1 = nx + 1
      ny1 = ny + 1
      nz1 = nz + 1
      do 30 l = 1, nz1
      do 20 k = 1, ny1
      do 10 j = 1, nx1
      fxyz(1,j,k,l) = fxyz3(1,j,k,l)
      fxyz(2,j,k,l) = fxyz3(2,j,k,l)
      fxyz(3,j,k,l) = fxyz3(3,j,k,l)
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine HAFTRP3D(q,q3,nx,ny,nz,nxe,nye,nze,nx2v,ny2,nz2)
c this subroutine copies data from a triple array to regular array
c with guard cells for scalar field and linear interpolation
c nx/ny/nz = system length in x/y/z direction
c nxe = first dimension of output array q, must be >= nx+1
c nye = second dimension of ouput array q, must be >= ny+1
c nze = third dimension of ouput array q, must be >= nz+1
c nx2v = first dimension of input array q3, must be >= 2*nx
c ny2 = second dimension of input array q3, must be >= 2*ny
c nz2 = third dimension of input array q3, must be >= 2*nz
      implicit none
      real q, q3
      integer nx, ny, nz, nxe, nye, nze, nx2v, ny2, nz2
      dimension q(nxe,nye,nze), q3(nx2v,ny2,nz2)
c local data
      integer j, k, l, nx1, ny1, nz1
      nx1 = nx + 1
      ny1 = ny + 1
      nz1 = nz + 1
      do 30 l = 1, nz1
      do 20 k = 1, ny1
      do 10 j = 1, nx1
      q(j,k,l) = q3(j,k,l)
   10 continue
   20 continue
   30 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PHAFDBL32C(fxyz,fxyz2,scb,scd,nx,ny,kstrt,nvpy,nvpz,nxv
     1,kyp,kzp,kypd,kzpd,kyp2,kblok,lblok,k2blok)
c this subroutine copies data from a double array to regular array
c with guard cells for vector field and linear interpolation
c for distributed data with 2D domain decomposition
c fxyz array may be modified
c scb/scd = scratch arrays
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = second dimension of input array fxyz, must be >= nx
c kyp/kzp = number of data values per block in fxyz in y/z
c kypd = third dimension of input array fxyz, must be >= kyp
c kzpd = fourth dimension of input array fxyz, must be >= kzp
c kyp2 = number of data values per block in fxyz2 in y
c kblok/lblok = number of data blocks in y/z
c k2blok = number of data blocks in y for doubled data
      implicit none
      real fxyz, fxyz2, scb, scd
      integer nx, ny, kstrt, nvpy, nvpz, nxv, kyp, kzp, kypd, kzpd, kyp2
      integer kblok, lblok, k2blok
      dimension fxyz(3,nxv,kypd,kzpd,kblok*lblok)
      dimension fxyz2(3,2*nxv,2*kypd,kzp,k2blok*lblok)
      dimension scb(3,nxv,kypd,kzpd,kblok*lblok)
      dimension scd(3,nxv,kzpd,2,kblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, lsid, msid, nsid, ierr
      integer i, j, k, l, m, my, mz, nx1, ny1, kyb, kyb2, ks, js, kr, kl
      integer nxvy, noff, moff, koff, ky, kz, ll, lm, mm, ml
      dimension istatus(lstat)
      nx1 = nx + 1
      ny1 = ny + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      kyb = ny/kyp
      kyb2 = (ny + ny)/kyp2
      nxvy = nxv*kypd
      noff = kypd*kzpd + kyb
c copy from double array in x direction
      do 270 mz = 1, lblok
      moff = k2blok*(mz - 1)
      do 260 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      lm = koff/kyp + 1
      koff = kyp*(my + js)
      ll = koff/kyp2 + 1
      koff = koff - kyp2*(ll - 1)
      kz = nvpy*(mz + ks)
      mm = ll + kz
      ml = lm + kz
c special case for one processor in y direction
      if (kyb2.eq.1) then
         do 40 l = 1, kzp
         do 30 k = 1, ny1
         do 20 j = 1, nx1
         do 10 i = 1, 3
         fxyz(i,j,k,l,m) = fxyz2(i,j,k,l,m)
   10    continue
   20    continue
   30    continue
   40    continue
         go to 260
      endif
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        if ((koff.eq.0).and.(kyp.lt.kyp2)) then
c           do 80 l = 1, kzp
c           do 70 k = 1, kyp1
c           do 60 j = 1, nx1
c           do 50 i = 1, 3
c           fxyz(i,j,k,l,m) = fxyz2(i,j,k,l,mm)
c  50       continue
c  60       continue
c  70       continue
c  80       continue
c        else
c           do 120 l = 1, kzp
c           do 110 k = 1, kyp
c           do 100 j = 1, nx1
c           do 90 i = 1, 3
c           fxyz(i,j,k,l,m) = fxyz2(i,j,k+kyp,l,mm)
c  90       continue
c 100       continue
c 110       continue
c 120       continue
c           do 150 l = 1, kzp
c           do 140 j = 1, nx1
c           do 130 i = 1, 3
c           fxyz(i,j,kyp+1,l,m) = fxyz2(i,j,1,l,mm+1)
c 130       continue
c 140       continue
c 150       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_IRECV(fxyz(1,1,1,1,m),3*kzp*nxvy,mreal,mm-1,noff+1,
     1lgrp,msid,ierr)
         else
            call MPI_IRECV(fxyz(1,1,1,1,m),3*kzp*nxvy,mreal,mm-1,noff+1,
     1lgrp,msid,ierr)
            call MPI_IRECV(scd(1,1,1,2,m),3*kzp*nxv,mreal,mm,noff+1,lgrp
     1,nsid,ierr)
         endif
      endif
c pack data and send it
      if (lm.le.kyb) then
         if (kyp.lt.kyp2) then
            do 80 l = 1, kzp
            do 70 k = 1, kypd
            do 60 j = 1, nxv
            do 50 i = 1, 3
            scb(i,j,k,l,m) = fxyz2(i,j,k,l,m)
   50       continue
   60       continue
   70       continue
   80       continue
            call MPI_SEND(scb(1,1,1,1,m),3*kzp*nxvy,mreal,ml-1,noff+1,lg
     1rp,ierr)
            do 120 l = 1, kzp
            do 110 k = 1, kypd
            do 100 j = 1, nxv
            do 90 i = 1, 3
            scb(i,j,k,l,m) = fxyz2(i,j,k+kyp,l,m)
   90       continue
  100       continue
  110       continue
  120       continue
            call MPI_SEND(scb(1,1,1,1,m),3*kzp*nxvy,mreal,ml,noff+1,lgrp
     1,ierr)
         else
            do 160 l = 1, kzp
            do 150 k = 1, kypd
            do 140 j = 1, nxv
            do 130 i = 1, 3
            scb(i,j,k,l,m) = fxyz2(i,j,k,l,m)
  130       continue
  140       continue
  150       continue
  160       continue
            call MPI_SEND(scb(1,1,1,1,m),3*kzp*nxvy,mreal,ml-1,noff+1,lg
     1rp,ierr)
         endif
         if (lm.gt.1) then
            do 190 l = 1, kzp
            do 180 j = 1, nxv
            do 170 i = 1, 3
            scd(i,j,l,1,m) = fxyz2(i,j,1,l,m)
  170       continue
  180       continue
  190       continue
            call MPI_SEND(scd(1,1,1,1,m),3*kzp*nxv,mreal,ml-2,noff+1,lgr
     1p,ierr)
         endif
      else if (lm.eq.(kyb+1)) then
         do 220 l = 1, kzp
         do 210 j = 1, nxv
         do 200 i = 1, 3
         scd(i,j,l,1,m) = fxyz2(i,j,1,l,m)
  200    continue
  210    continue
  220    continue
         call MPI_SEND(scd(1,1,1,1,m),3*kzp*nxv,mreal,ml-2,noff+1,lgrp,i
     1err)
      endif
c wait for data
      if (ll.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_WAIT(msid,istatus,ierr)
         else
            call MPI_WAIT(msid,istatus,ierr)
            call MPI_WAIT(nsid,istatus,ierr)
            do 250 l = 1, kzp
            do 240 j = 1, nxv
            do 230 i = 1, 3
            fxyz(i,j,kyp+1,l,m) = scd(i,j,l,2,m)
  230       continue
  240       continue
  250       continue
         endif
      endif
  260 continue
  270 continue
c copy to guard cells in z
      do 320 mz = 1, lblok
      moff = kblok*(mz - 1)
      do 310 my = 1, kblok
      m = my + moff
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      kl = kz - 1
      if (kl.lt.0) kl = kl + nvpz
      kr = ky + nvpy*kr
      kl = ky + nvpy*kl
c this segment is used for shared memory computers
c     do 300 k = 1, nypmx
c     do 290 j = 1, nxv
c     do 280 i = 1, 3
c     fxyz(i,j,k,kzp+1,m) = fxyz(i,j,k,1,kr)
c 280 continue
c 290 continue
c 300 continue
c this segment is used for mpi computers
      call MPI_IRECV(fxyz(1,1,1,kzp+1,m),3*nxvy,mreal,kr-1,noff+2,lgrp,m
     1sid,ierr)
      call MPI_SEND(fxyz(1,1,1,1,m),3*nxvy,mreal,kl-1,noff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
  310 continue
  320 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine PHAFDBL32D(q,q2,scb,scd,nx,ny,kstrt,nvpy,nvpz,nxv,kyp,k
     1zp,kypd,kzpd,kyp2,kblok,lblok,k2blok)
c this subroutine copies data from a double array to regular array
c with guard cells for scalar field and linear interpolation
c for distributed data with 2D domain decomposition
c q array may be modified
c scb/scd = scratch arrays
c nx/ny = system length in x/y direction
c kstrt = starting data block number
c nvpy/nvpz = number of real or virtual processors in y/z
c nxv = second dimension of input array q, must be >= nx
c kyp/kzp = number of data values per block in q in y/z
c kypd = third dimension of input array q, must be >= kyp
c kzpd = fourth dimension of input array q, must be >= kzp
c kyp2 = number of data values per block in q in y
c kblok/lblok = number of data blocks in y/z
c k2blok = number of data blocks in y for doubled data
      implicit none
      real q, q2, scb, scd
      integer nx, ny, kstrt, nvpy, nvpz, nxv, kyp, kzp, kypd, kzpd, kyp2
      integer kblok, lblok, k2blok
      dimension q(nxv,kypd,kzpd,kblok*lblok)
      dimension q2(2*nxv,2*kypd,kzp,k2blok*lblok)
      dimension scb(nxv,kypd,kzpd,kblok*lblok)
      dimension scd(nxv,kzpd,2,kblok*lblok)
c common block for parallel processing
      integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
c lstat = length of status array
      parameter(lstat=10)
c lgrp = current communicator
c mreal = default datatype for reals
      common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
c local data
      integer istatus, lsid, msid, nsid, ierr
      integer j, k, l, m, my, mz, nx1, ny1, kyb, kyb2, ks, js, kr, kl
      integer nxvy, noff, moff, koff, ky, kz, ll, lm, mm, ml
      dimension istatus(lstat)
      nx1 = nx + 1
      ny1 = ny + 1
      ks = (kstrt - 1)/nvpy
      js = kstrt - nvpy*ks - 2
      ks = ks - 1
      kyb = ny/kyp
      kyb2 = (ny + ny)/kyp2
      nxvy = nxv*kypd
      noff = kypd*kzpd + kyb
c copy from double array in x direction
      do 200 mz = 1, lblok
      moff = k2blok*(mz - 1)
      do 190 my = 1, k2blok
      m = my + moff
      koff = kyp2*(my + js)
      lm = koff/kyp + 1
      koff = kyp*(my + js)
      ll = koff/kyp2 + 1
      koff = koff - kyp2*(ll - 1)
      kz = nvpy*(mz + ks)
      mm = ll + kz
      ml = lm + kz
c special case for one processor in y direction
      if (kyb2.eq.1) then
         do 30 l = 1, kzp
         do 20 k = 1, ny1
         do 10 j = 1, nx1
         q(j,k,l,m) = q2(j,k,l,m)
   10    continue
   20    continue
   30    continue
         go to 190
      endif
c this segment is used for shared memory computers
c     if (ll.le.kyb) then
c        if ((koff.eq.0).and.(kyp.lt.kyp2)) then
c           do 60 l = 1, kzp
c           do 50 k = 1, kyp1
c           do 40 j = 1, nx1
c           q(j,k,l,m) = q2(j,k,l,mm)
c  40       continue
c  50       continue
c  60       continue
c        else
c           do 90 l = 1, kzp
c           do 80 k = 1, kyp
c           do 70 j = 1, nx1
c           q(j,k,l,m) = q2(j,k+kyp,l,mm)
c  70       continue
c  80       continue
c  90       continue
c           do 110 l = 1, kzp
c           do 100 j = 1, nx1
c           q(j,kyp+1,l,m) = q2(j,1,l,mm+1)
c 100       continue
c 110       continue
c        endif
c     endif
c this segment is used for mpi computers
      if (ll.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_IRECV(q(1,1,1,m),kzp*nxvy,mreal,mm-1,noff+1,lgrp,ms
     1id,ierr)
         else
            call MPI_IRECV(q(1,1,1,m),kzp*nxvy,mreal,mm-1,noff+1,lgrp,ms
     1id,ierr)
            call MPI_IRECV(scd(1,1,2,m),kzp*nxv,mreal,mm,noff+1,lgrp,nsi
     1d,ierr)
         endif
      endif
c pack data and send it
      if (lm.le.kyb) then
         if (kyp.lt.kyp2) then
            do 60 l = 1, kzp
            do 50 k = 1, kypd
            do 40 j = 1, nxv
            scb(j,k,l,m) = q2(j,k,l,m)
   40       continue
   50       continue
   60       continue
            call MPI_SEND(scb(1,1,1,m),kzp*nxvy,mreal,ml-1,noff+1,lgrp,i
     1err)
            do 90 l = 1, kzp
            do 80 k = 1, kypd
            do 70 j = 1, nxv
            scb(j,k,l,m) = q2(j,k+kyp,l,m)
   70       continue
   80       continue
   90       continue
            call MPI_SEND(scb(1,1,1,m),kzp*nxvy,mreal,ml,noff+1,lgrp,ier
     1r)
         else
            do 120 l = 1, kzp
            do 110 k = 1, kypd
            do 100 j = 1, nxv
            scb(j,k,l,m) = q2(j,k,l,m)
  100       continue
  110       continue
  120       continue
            call MPI_SEND(scb(1,1,1,m),kzp*nxvy,mreal,ml-1,noff+1,lgrp,i
     1err)
         endif
         if (lm.gt.1) then
            do 140 l = 1, kzp
            do 130 j = 1, nxv
            scd(j,l,1,m) = q2(j,1,l,m)
  130       continue
  140       continue
            call MPI_SEND(scd(1,1,1,m),kzp*nxv,mreal,ml-2,noff+1,lgrp,ie
     1rr)
         endif
      else if (lm.eq.(kyb+1)) then
         do 160 l = 1, kzp
         do 150 j = 1, nxv
         scd(j,l,1,m) = q2(j,1,l,m)
  150    continue
  160    continue
         call MPI_SEND(scd(1,1,1,m),kzp*nxv,mreal,ml-2,noff+1,lgrp,ierr)
      endif
c wait for data
      if (ll.le.kyb) then
         if ((koff.eq.0).and.(kyp.lt.kyp2)) then
            call MPI_WAIT(msid,istatus,ierr)
         else
            call MPI_WAIT(msid,istatus,ierr)
            call MPI_WAIT(nsid,istatus,ierr)
            do 180 l = 1, kzp
            do 170 j = 1, nxv
            q(j,kyp+1,l,m) = scd(j,l,2,m)
  170       continue
  180       continue
         endif
      endif
  190 continue
  200 continue
c copy to guard cells in z
      do 240 mz = 1, lblok
      moff = kblok*(mz - 1)
      do 230 my = 1, kblok
      m = my + moff
      ky = my + js + 1
      kz = mz + ks
      kr = kz + 1
      if (kr.ge.nvpz) kr = kr - nvpz
      kl = kz - 1
      if (kl.lt.0) kl = kl + nvpz
      kr = ky + nvpy*kr
      kl = ky + nvpy*kl
c this segment is used for shared memory computers
c     do 220 k = 1, nypmx
c     do 210 j = 1, nxv
c     q(j,k,kzp+1,m) = q(j,k,1,kr)
c 210 continue
c 220 continue
c this segment is used for mpi computers
      call MPI_IRECV(q(1,1,kzp+1,m),nxvy,mreal,kr-1,noff+2,lgrp,msid,ier
     1r)
      call MPI_SEND(q(1,1,1,m),nxvy,mreal,kl-1,noff+2,lgrp,ierr)
      call MPI_WAIT(msid,istatus,ierr)
  230 continue
  240 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine HAFDBL3C(fxyz,fxyz2,nx,ny,nz,nxe,nye,nze,nx2v,ny2,nzv)
c this subroutine copies data from a double array to regular array
c with guard cells for vector field and linear interpolation
c nx/ny/nz = system length in x/y/z direction
c nxe = first dimension of output array fxyz, must be >= nx+1
c nye = second dimension of ouput array fxyz, must be >= ny+1
c nze = third dimension of ouput array fxyz, must be >= nz+1
c nx2v = first dimension of input array fxyz2, must be >= 2*nx
c ny2 = second dimension of input array fxyz2, must be >= 2*ny
      implicit none
      real fxyz, fxyz2
      integer nx, ny, nz, nxe, nye, nze, nx2v, ny2, nzv
      dimension fxyz(3,nxe,nye,nze), fxyz2(3,nx2v,ny2,nzv)
c local data
      integer j, k, l, nx1, ny1, nz1
      nx1 = nx + 1
      ny1 = ny + 1
      nz1 = nz + 1
      do 30 l = 1, nz
      do 20 k = 1, ny1
      do 10 j = 1, nx1
      fxyz(1,j,k,l) = fxyz2(1,j,k,l)
      fxyz(2,j,k,l) = fxyz2(2,j,k,l)
      fxyz(3,j,k,l) = fxyz2(3,j,k,l)
   10 continue
   20 continue
   30 continue
c copy periodic edge
      do 50 k = 1, ny1
      do 40 j = 1, nx1
      fxyz(1,j,k,nz1) = fxyz(1,j,k,1)
      fxyz(2,j,k,nz1) = fxyz(2,j,k,1)
      fxyz(3,j,k,nz1) = fxyz(3,j,k,1)
   40 continue
   50 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine HAFDBL3D(q,q2,nx,ny,nz,nxe,nye,nze,nx2v,ny2,nzv)
c this subroutine copies data from a double array to regular array
c with guard cells for scalar field and linear interpolation
c nx/ny/nz = system length in x/y/z direction
c nxe = first dimension of output array q, must be >= nx+1
c nye = second dimension of ouput array q, must be >= ny+1
c nze = third dimension of ouput array q, must be >= nz+1
c nx2v = first dimension of input array q2, must be >= 2*nx
c ny2 = second dimension of input array q2, must be >= 2*ny
c nzv = third dimension of input array q2, must be >= nz
      implicit none
      real q, q2
      integer nx, ny, nz, nxe, nye, nze, nx2v, ny2, nzv
      dimension q(nxe,nye,nze), q2(nx2v,ny2,nzv)
c local data
      integer j, k, l, nx1, ny1, nz1
      nx1 = nx + 1
      ny1 = ny + 1
      nz1 = nz + 1
      do 30 l = 1, nz
      do 20 k = 1, ny1
      do 10 j = 1, nx1
      q(j,k,l) = q2(j,k,l)
   10 continue
   20 continue
   30 continue
c copy periodic edge
      do 50 k = 1, ny1
      do 40 j = 1, nx1
      q(j,k,nz1) = q(j,k,1)
   40 continue
   50 continue
      return
      end
