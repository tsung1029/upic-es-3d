c 3d parallel PIC multi-tasking library fast fourier transforms
c with 2D domain decomposition
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: june 16, 2008
c-----------------------------------------------------------------------
      subroutine MPFFT32R(f,g,h,bs,br,isign,ntpose,mixup,sct,ttp,indx,in
     1dy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd
     2,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd,kyzip,if
     3task,nmt,ierr)
c multi-tasking real to complex fft
c nyzip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd, kzyp
      integer jblok, kblok, jkblok, lblok, mblok, mlblok, nxhyzd, nxyzhd
      integer mixup, kyzip, iftask, nmt, ierr
      real ttp
      complex f, g, h, bs, br, sct
      dimension f(nxvh,kypd,kzpd,kblok*lblok)
      dimension g(nyv,kxypd,kzpd,jblok*lblok)
      dimension h(nzv,kxypd,kyzpd,jblok*mblok)
      dimension bs(kxyp,kzyp,kzp,jkblok*lblok)
      dimension br(kxyp,kzyp,kzp,jblok*mlblok)
      dimension mixup(nxhyzd), sct(nxyzhd)
      dimension kyzip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, nz, nxh, kypi, kxypi, kypp, kxypp
      integer kypl, kxypl
      integer i
      real tp, tf
      double precision dtime
      external PFFT32RXX, PFFT32RXY, PFFT32RXZ
      data nargs /19/
      nx = 2**indx
      ny = 2**indy
      nz = 2**indz
      nxh = nx/2
      kypp = kyp/(nmt + 1)
      kxypp = kxyp/(nmt + 1)
      kypi = kypp*nmt
      kxypi = kxypp*nmt
      kypl = kyp - kypi
      kxypl = kxyp - kxypi
      kypi = kypi + 1
      kxypi = kxypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x fft tasks
         do 10 i = 1, nmt
         kyzip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RXX,nargs,f,isign,mixup,sct,i
     1ndx,indy,indz,kstrt,kyzip(i),kypp,nxvh,kyp,kzp,kypd,kzpd,kblok,lbl
     2ok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x fft
         call PFFT32RXX(f,isign,mixup,sct,indx,indy,indz,kstrt,kypi,kypl
     1,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PTPOS3A(f,g,bs,br,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,kx
     1ypd,kypd,kzpd,jblok,kblok,lblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 30 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RXY,nargs,g,isign,mixup,sct,i
     1ndx,indy,indz,kstrt,kyzip(i),kxypp,nyv,kxyp,kzp,kxypd,kzpd,jblok,l
     2blok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call PFFT32RXY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kxy
     1pl,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to h
         call PWTIMERA(-1,tp,dtime)
         call PTPOS3B(g,h,bs,br,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kx
     1ypd,kyzpd,kzpd,jblok,mblok,lblok)
         call PWTIMERA(1,tp,dtime)
c start z fft tasks
         do 50 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RXZ,nargs,h,isign,mixup,sct,i
     1ndx,indy,indz,kstrt,kyzip(i),kxypp,nzv,kxyp,kyzp,kxypd,kyzpd,jblok
     2,mblok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish z fft
         call PFFT32RXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kxy
     1pl,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
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
c start z fft tasks
         do 70 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RXZ,nargs,h,isign,mixup,sct,i
     1ndx,indy,indz,kstrt,kyzip(i),kxypp,nzv,kxyp,kyzp,kxypd,kyzpd,jblok
     2,mblok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish z fft
         call PFFT32RXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kxy
     1pl,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
c transpose h array to g
         call PWTIMERA(-1,tp,dtime)
         call PTPOS3B(h,g,br,bs,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kx
     1ypd,kzpd,kyzpd,jblok,lblok,mblok)
         call PWTIMERA(1,tp,dtime)
c start y fft tasks
         do 90 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RXY,nargs,g,isign,mixup,sct,i
     1ndx,indy,indz,kstrt,kyzip(i),kxypp,nyv,kxyp,kzp,kxypd,kzpd,jblok,l
     2blok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   90    continue
c finish y fft
         call PFFT32RXY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kxy
     1pl,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 100 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
  100    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PTPOS3A(g,f,br,bs,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,ky
     1pd,kxypd,kzpd,kblok,jblok,lblok)
         call PWTIMERA(1,ttp,dtime)
c start x fft tasks
         do 110 i = 1, nmt
         kyzip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RXX,nargs,f,isign,mixup,sct,i
     1ndx,indy,indz,kstrt,kyzip(i),kypp,nxvh,kyp,kzp,kypd,kzpd,kblok,lbl
     2ok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  110    continue
c finish x fft
         call PFFT32RXX(f,isign,mixup,sct,indx,indy,indz,kstrt,kypi,kypl
     1,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 120 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  120    continue
      endif
      ttp = ttp + tp
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFFT32RX(f,g,h,isign,ntpose,mixup,sct,ttp,indx,indy,in
     1dz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,jblo
     2k,kblok,lblok,mblok,nxhyzd,nxyzhd,kyzip,iftask,nmt,ierr)
c multi-tasking real to complex fft
c nyzip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd
      integer jblok, kblok, lblok, mblok, nxhyzd, nxyzhd
      integer mixup, kyzip, iftask, nmt, ierr
      real ttp
      complex f, g, h, sct
      dimension f(nxvh,kypd,kzpd,kblok*lblok)
      dimension g(nyv,kxypd,kzpd,jblok*lblok)
      dimension h(nzv,kxypd,kyzpd,jblok*mblok)
      dimension mixup(nxhyzd), sct(nxyzhd)
      dimension kyzip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, nz, nxh, kypi, kxypi, kypp, kxypp
      integer kypl, kxypl
      integer i
      real tp, tf
      double precision dtime
      external PFFT32RXX, PFFT32RXY, PFFT32RXZ
      data nargs /19/
      nx = 2**indx
      ny = 2**indy
      nz = 2**indz
      nxh = nx/2
      kypp = kyp/(nmt + 1)
      kxypp = kxyp/(nmt + 1)
      kypi = kypp*nmt
      kxypi = kxypp*nmt
      kypl = kyp - kypi
      kxypl = kxyp - kxypi
      kypi = kypi + 1
      kxypi = kxypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x fft tasks
         do 10 i = 1, nmt
         kyzip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RXX,nargs,f,isign,mixup,sct,i
     1ndx,indy,indz,kstrt,kyzip(i),kypp,nxvh,kyp,kzp,kypd,kzpd,kblok,lbl
     2ok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x fft
         call PFFT32RXX(f,isign,mixup,sct,indx,indy,indz,kstrt,kypi,kypl
     1,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PTPOS3AX(f,g,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,kxypd,k
     1ypd,kzpd,jblok,kblok,lblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 30 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RXY,nargs,g,isign,mixup,sct,i
     1ndx,indy,indz,kstrt,kyzip(i),kxypp,nyv,kxyp,kzp,kxypd,kzpd,jblok,l
     2blok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call PFFT32RXY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kxy
     1pl,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to h
         call PWTIMERA(-1,tp,dtime)
         call PTPOS3BX(g,h,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxypd,k
     1yzpd,kzpd,jblok,mblok,lblok)
         call PWTIMERA(1,tp,dtime)
c start z fft tasks
         do 50 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RXZ,nargs,h,isign,mixup,sct,i
     1ndx,indy,indz,kstrt,kyzip(i),kxypp,nzv,kxyp,kyzp,kxypd,kyzpd,jblok
     2,mblok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish z fft
         call PFFT32RXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kxy
     1pl,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
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
c start z fft tasks
         do 70 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RXZ,nargs,h,isign,mixup,sct,i
     1ndx,indy,indz,kstrt,kyzip(i),kxypp,nzv,kxyp,kyzp,kxypd,kyzpd,jblok
     2,mblok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish z fft
         call PFFT32RXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kxy
     1pl,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
c transpose h array to g
         call PWTIMERA(-1,tp,dtime)
         call PTPOS3BX(h,g,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kxypd,k
     1zpd,kyzpd,jblok,lblok,mblok)
         call PWTIMERA(1,tp,dtime)
c start y fft tasks
         do 90 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RXY,nargs,g,isign,mixup,sct,i
     1ndx,indy,indz,kstrt,kyzip(i),kxypp,nyv,kxyp,kzp,kxypd,kzpd,jblok,l
     2blok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   90    continue
c finish y fft
         call PFFT32RXY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kxy
     1pl,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 100 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
  100    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PTPOS3AX(g,f,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,kypd,kx
     1ypd,kzpd,kblok,jblok,lblok)
         call PWTIMERA(1,ttp,dtime)
c start x fft tasks
         do 110 i = 1, nmt
         kyzip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RXX,nargs,f,isign,mixup,sct,i
     1ndx,indy,indz,kstrt,kyzip(i),kypp,nxvh,kyp,kzp,kypd,kzpd,kblok,lbl
     2ok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  110    continue
c finish x fft
         call PFFT32RXX(f,isign,mixup,sct,indx,indy,indz,kstrt,kypi,kypl
     1,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 120 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  120    continue
      endif
      ttp = ttp + tp
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFFT32R3(f,g,h,bs,br,isign,ntpose,mixup,sct,ttp,indx,i
     1ndy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzp
     2d,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,nxhyzd,nxyzhd,kyzip,i
     3ftask,nmt,ierr)
c multi-tasking real to complex fft
c nyzip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd, kzyp
      integer jblok, kblok, jkblok, lblok, mblok, mlblok, nxhyzd, nxyzhd
      integer mixup, kyzip, iftask, nmt, ierr
      real ttp
      complex f, g, h, bs, br, sct
      dimension f(3,nxvh,kypd,kzpd,kblok*lblok)
      dimension g(3,nyv,kxypd,kzpd,jblok*lblok)
      dimension h(3,nzv,kxypd,kyzpd,jblok*mblok)
      dimension bs(3,kxyp,kzyp,kzp,jkblok*lblok)
      dimension br(3,kxyp,kzyp,kzp,jblok*mlblok)
      dimension mixup(nxhyzd), sct(nxyzhd)
      dimension kyzip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, nz, nxh, kypi, kxypi, kypp, kxypp
      integer kypl, kxypl
      integer i
      real tp, tf
      double precision dtime
      external PFFT32R3XX, PFFT32R3XY, PFFT32R3XZ
      data nargs /19/
      nx = 2**indx
      ny = 2**indy
      nz = 2**indz
      nxh = nx/2
      kypp = kyp/(nmt + 1)
      kxypp = kxyp/(nmt + 1)
      kypi = kypp*nmt
      kxypi = kxypp*nmt
      kypl = kyp - kypi
      kxypl = kxyp - kxypi
      kypi = kypi + 1
      kxypi = kxypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x fft tasks
         do 10 i = 1, nmt
         kyzip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32R3XX,nargs,f,isign,mixup,sct,
     1indx,indy,indz,kstrt,kyzip(i),kypp,nxvh,kyp,kzp,kypd,kzpd,kblok,lb
     2lok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x fft
         call PFFT32R3XX(f,isign,mixup,sct,indx,indy,indz,kstrt,kypi,kyp
     1l,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call P3TPOS3A(f,g,bs,br,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,k
     1xypd,kypd,kzpd,jblok,kblok,lblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 30 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32R3XY,nargs,g,isign,mixup,sct,
     1indx,indy,indz,kstrt,kyzip(i),kxypp,nyv,kxyp,kzp,kxypd,kzpd,jblok,
     2lblok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call PFFT32R3XY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1ypl,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to h
         call PWTIMERA(-1,tp,dtime)
         call P3TPOS3B(g,h,bs,br,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,k
     1xypd,kyzpd,kzpd,jblok,mblok,lblok)
         call PWTIMERA(1,tp,dtime)
c start z fft tasks
         do 50 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32R3XZ,nargs,h,isign,mixup,sct,
     1indx,indy,indz,kstrt,kyzip(i),kxypp,nzv,kxyp,kyzp,kxypd,kyzpd,jblo
     2k,mblok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish z fft
         call PFFT32R3XZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1ypl,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
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
c start z fft tasks
         do 70 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32R3XZ,nargs,h,isign,mixup,sct,
     1indx,indy,indz,kstrt,kyzip(i),kxypp,nzv,kxyp,kyzp,kxypd,kyzpd,jblo
     2k,mblok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish z fft
         call PFFT32R3XZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1ypl,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
c transpose h array to g
         call PWTIMERA(-1,tp,dtime)
         call P3TPOS3B(h,g,br,bs,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,k
     1xypd,kzpd,kyzpd,jblok,lblok,mblok)
         call PWTIMERA(1,tp,dtime)
c start y fft tasks
         do 90 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32R3XY,nargs,g,isign,mixup,sct,
     1indx,indy,indz,kstrt,kyzip(i),kxypp,nyv,kxyp,kzp,kxypd,kzpd,jblok,
     2lblok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   90    continue
c finish y fft
         call PFFT32R3XY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1ypl,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 100 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
  100    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call P3TPOS3A(g,f,br,bs,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,k
     1ypd,kxypd,kzpd,kblok,jblok,lblok)
         call PWTIMERA(1,ttp,dtime)
c start x fft tasks
         do 110 i = 1, nmt
         kyzip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32R3XX,nargs,f,isign,mixup,sct,
     1indx,indy,indz,kstrt,kyzip(i),kypp,nxvh,kyp,kzp,kypd,kzpd,kblok,lb
     2lok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  110    continue
c finish x fft
         call PFFT32R3XX(f,isign,mixup,sct,indx,indy,indz,kstrt,kypi,kyp
     1l,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 120 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  120    continue
      endif
      ttp = ttp + tp
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFFT32RX3(f,g,h,isign,ntpose,mixup,sct,ttp,indx,indy,i
     1ndz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,jbl
     2ok,kblok,lblok,mblok,nxhyzd,nxyzhd,kyzip,iftask,nmt,ierr)
c multi-tasking real to complex fft
c nyzip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd
      integer jblok, kblok, lblok, mblok, nxhyzd, nxyzhd
      integer mixup, kyzip, iftask, nmt, ierr
      real ttp
      complex f, g, h, sct
      dimension f(3,nxvh,kypd,kzpd,kblok*lblok)
      dimension g(3,nyv,kxypd,kzpd,jblok*lblok)
      dimension h(3,nzv,kxypd,kyzpd,jblok*mblok)
      dimension mixup(nxhyzd), sct(nxyzhd)
      dimension kyzip(nmt), iftask(nmt)
c local data
      integer nargs, nx, ny, nz, nxh, kypi, kxypi, kypp, kxypp
      integer kypl, kxypl
      integer i
      real tp, tf
      double precision dtime
      external PFFT32R3XX, PFFT32R3XY, PFFT32R3XZ
      data nargs /19/
      nx = 2**indx
      ny = 2**indy
      nz = 2**indz
      nxh = nx/2
      kypp = kyp/(nmt + 1)
      kxypp = kxyp/(nmt + 1)
      kypi = kypp*nmt
      kxypi = kxypp*nmt
      kypl = kyp - kypi
      kxypl = kxyp - kxypi
      kypi = kypi + 1
      kxypi = kxypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x fft tasks
         do 10 i = 1, nmt
         kyzip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32R3XX,nargs,f,isign,mixup,sct,
     1indx,indy,indz,kstrt,kyzip(i),kypp,nxvh,kyp,kzp,kypd,kzpd,kblok,lb
     2lok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x fft
         call PFFT32R3XX(f,isign,mixup,sct,indx,indy,indz,kstrt,kypi,kyp
     1l,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call P3TPOS3AX(f,g,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,kxypd,
     1kypd,kzpd,jblok,kblok,lblok)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 30 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32R3XY,nargs,g,isign,mixup,sct,
     1indx,indy,indz,kstrt,kyzip(i),kxypp,nyv,kxyp,kzp,kxypd,kzpd,jblok,
     2lblok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call PFFT32R3XY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1ypl,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to h
         call PWTIMERA(-1,tp,dtime)
         call P3TPOS3BX(g,h,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxypd,
     1kyzpd,kzpd,jblok,mblok,lblok)
         call PWTIMERA(1,tp,dtime)
c start z fft tasks
         do 50 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32R3XZ,nargs,h,isign,mixup,sct,
     1indx,indy,indz,kstrt,kyzip(i),kxypp,nzv,kxyp,kyzp,kxypd,kyzpd,jblo
     2k,mblok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish z fft
         call PFFT32R3XZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1ypl,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
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
            call PWTIMERA(1,tf,dtime)
         endif
c start z fft tasks
         do 70 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32R3XZ,nargs,h,isign,mixup,sct,
     1indx,indy,indz,kstrt,kyzip(i),kxypp,nzv,kxyp,kyzp,kxypd,kyzpd,jblo
     2k,mblok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish z fft
         call PFFT32R3XZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1ypl,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
c transpose h array to g
         call PWTIMERA(-1,tp,dtime)
         call P3TPOS3BX(h,g,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kxypd,
     1kzpd,kyzpd,jblok,lblok,mblok)
         call PWTIMERA(1,tp,dtime)
c start y fft tasks
         do 90 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32R3XY,nargs,g,isign,mixup,sct,
     1indx,indy,indz,kstrt,kyzip(i),kxypp,nyv,kxyp,kzp,kxypd,kzpd,jblok,
     2lblok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   90    continue
c finish y fft
         call PFFT32R3XY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1ypl,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 100 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
  100    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call P3TPOS3AX(g,f,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,kypd,k
     1xypd,kzpd,kblok,jblok,lblok)
         call PWTIMERA(1,ttp,dtime)
c start x fft tasks
         do 110 i = 1, nmt
         kyzip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32R3XX,nargs,f,isign,mixup,sct,
     1indx,indy,indz,kstrt,kyzip(i),kypp,nxvh,kyp,kzp,kypd,kzpd,kblok,lb
     2lok,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  110    continue
c finish x fft
         call PFFT32R3XX(f,isign,mixup,sct,indx,indy,indz,kstrt,kypi,kyp
     1l,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 120 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  120    continue
      endif
      ttp = ttp + tp
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFFT32RN(f,g,h,bs,br,ss,isign,ntpose,mixup,sct,ttp,ind
     1x,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,
     2kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,ndim,nxhyzd,nxyzhd
     3,kyzip,iftask,nmt,ierr)
c multi-tasking real to complex fft
c nyzip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd, kzyp
      integer jblok, kblok, jkblok, lblok, mblok, mlblok, ndim
      integer nxhyzd, nxyzhd
      integer mixup, kyzip, iftask, nmt, ierr
      real ttp
      complex f, g, h, bs, br, ss, sct
      dimension f(ndim,nxvh,kypd,kzpd,kblok*lblok)
      dimension g(ndim,nyv,kxypd,kzpd,jblok*lblok)
      dimension h(ndim,nzv,kxypd,kyzpd,jblok*mblok)
      dimension bs(ndim,kxyp,kzyp,kzp,jkblok*lblok)
      dimension br(ndim,kxyp,kzyp,kzp,jblok*mlblok)
      dimension ss(ndim,nxvh,nmt+1)
      dimension mixup(nxhyzd), sct(nxyzhd)
      dimension kyzip(nmt), iftask(nmt)
c local data
      integer i, nargs, margs, nx, ny, nz, nxh, nmtt, kypi, kxypi, kypp
      integer kxypp, kypl, kxypl
      real tp, tf
      double precision dtime
      external PFFT32RNXX, PFFT32RNXY, PFFT32RNXZ
      data nargs, margs /20,21/
      nx = 2**indx
      ny = 2**indy
      nz = 2**indz
      nxh = nx/2
      nmtt = nmt + 1
      kypp = kyp/nmtt
      kxypp = kxyp/nmtt
      kypi = kypp*nmt
      kxypi = kxypp*nmt
      kypl = kyp - kypi
      kxypl = kxyp - kxypi
      kypi = kypi + 1
      kxypi = kxypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x fft tasks
         do 10 i = 1, nmt
         kyzip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXX,margs,f,ss(1,1,i),isign,
     1mixup,sct,indx,indy,indz,kstrt,kyzip(i),kypp,nxvh,kyp,kzp,kypd,kzp
     2d,kblok,lblok,ndim,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x fft
         call PFFT32RNXX(f,ss(1,1,nmtt),isign,mixup,sct,indx,indy,indz,k
     1strt,kypi,kypl,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,ndim,nxhyzd,nxyz
     2hd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PNTPOS3A(f,g,bs,br,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,k
     1xypd,kypd,kzpd,jblok,kblok,lblok,ndim)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 30 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXY,nargs,g,isign,mixup,sct,
     1indx,indy,indz,kstrt,kyzip(i),kxypp,nyv,kxyp,kzp,kxypd,kzpd,jblok,
     2lblok,ndim,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call PFFT32RNXY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1ypl,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,ndim,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to h
         call PWTIMERA(-1,tp,dtime)
         call PNTPOS3B(g,h,bs,br,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,k
     1xypd,kyzpd,kzpd,jblok,mblok,lblok,ndim)
         call PWTIMERA(1,tp,dtime)
c start z fft tasks
         do 50 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXZ,nargs,h,isign,mixup,sct,
     1indx,indy,indz,kstrt,kyzip(i),kxypp,nzv,kxyp,kyzp,kxypd,kyzpd,jblo
     2k,mblok,ndim,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish z fft
         call PFFT32RNXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1ypl,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,ndim,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
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
c start z fft tasks
         do 70 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXZ,nargs,h,isign,mixup,sct,
     1indx,indy,indz,kstrt,kyzip(i),kxypp,nzv,kxyp,kyzp,kxypd,kyzpd,jblo
     2k,mblok,ndim,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish z fft
         call PFFT32RNXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1ypl,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,ndim,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
c transpose h array to g
         call PWTIMERA(-1,tp,dtime)
         call PNTPOS3B(h,g,br,bs,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,k
     1xypd,kzpd,kyzpd,jblok,lblok,mblok,ndim)
         call PWTIMERA(1,tp,dtime)
c start y fft tasks
         do 90 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXY,nargs,g,isign,mixup,sct,
     1indx,indy,indz,kstrt,kyzip(i),kxypp,nyv,kxyp,kzp,kxypd,kzpd,jblok,
     2lblok,ndim,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   90    continue
c finish y fft
         call PFFT32RNXY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1ypl,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,ndim,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 100 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
  100    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PNTPOS3A(g,f,br,bs,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,k
     1ypd,kxypd,kzpd,kblok,jblok,lblok,ndim)
         call PWTIMERA(1,ttp,dtime)
c start x fft tasks
         do 110 i = 1, nmt
         kyzip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXX,margs,f,ss(1,1,i),isign,
     1mixup,sct,indx,indy,indz,kstrt,kyzip(i),kypp,nxvh,kyp,kzp,kypd,kzp
     2d,kblok,lblok,ndim,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  110    continue
c finish x fft
         call PFFT32RNXX(f,ss(1,1,nmtt),isign,mixup,sct,indx,indy,indz,k
     1strt,kypi,kypl,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,ndim,nxhyzd,nxyz
     2hd)
c wait for tasks to complete
         do 120 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  120    continue
      endif
      ttp = ttp + tp
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MPFFT32RXN(f,g,h,ss,isign,ntpose,mixup,sct,ttp,indx,ind
     1y,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxypd,kypd,kyzpd,kzpd,
     2jblok,kblok,lblok,mblok,ndim,nxhyzd,nxyzhd,kyzip,iftask,nmt,ierr)
c multi-tasking real to complex fft
c nyzip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd
      integer jblok, kblok, lblok, mblok, ndim, nxhyzd, nxyzhd
      integer mixup, kyzip, iftask, nmt, ierr
      real ttp
      complex f, g, h, ss, sct
      dimension f(ndim,nxvh,kypd,kzpd,kblok*lblok)
      dimension g(ndim,nyv,kxypd,kzpd,jblok*lblok)
      dimension h(ndim,nzv,kxypd,kyzpd,jblok*mblok)
      dimension ss(ndim,nxvh,nmt+1)
      dimension mixup(nxhyzd), sct(nxyzhd)
      dimension kyzip(nmt), iftask(nmt)
c local data
      integer i, nargs, margs, nx, ny, nz, nxh, nmtt, kypi, kxypi, kypp
      integer kxypp, kypl, kxypl
      real tp, tf
      double precision dtime
      external PFFT32RNXX, PFFT32RNXY, PFFT32RNXZ
      data nargs, margs /20,21/
      nx = 2**indx
      ny = 2**indy
      nz = 2**indz
      nxh = nx/2
      nmtt = nmt + 1
      kypp = kyp/nmtt
      kxypp = kxyp/nmtt
      kypi = kypp*nmt
      kxypi = kxypp*nmt
      kypl = kyp - kypi
      kxypl = kxyp - kxypi
      kypi = kypi + 1
      kxypi = kxypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start x fft tasks
         do 10 i = 1, nmt
         kyzip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXX,margs,f,ss(1,1,i),isign,
     1mixup,sct,indx,indy,indz,kstrt,kyzip(i),kypp,nxvh,kyp,kzp,kypd,kzp
     2d,kblok,lblok,ndim,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish x fft
         call PFFT32RNXX(f,ss(1,1,nmtt),isign,mixup,sct,indx,indy,indz,k
     1strt,kypi,kypl,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,ndim,nxhyzd,nxyz
     2hd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c transpose f array to g
         call PWTIMERA(-1,ttp,dtime)
         call PNTPOS3AX(f,g,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,kyp,kzp,kxypd,
     1kypd,kzpd,jblok,kblok,lblok,ndim)
         call PWTIMERA(1,ttp,dtime)
c start y fft tasks
         do 30 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXY,nargs,g,isign,mixup,sct,
     1indx,indy,indz,kstrt,kyzip(i),kxypp,nyv,kxyp,kzp,kxypd,kzpd,jblok,
     2lblok,ndim,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish y fft
         call PFFT32RNXY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1ypl,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,ndim,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   40    continue
c transpose g array to h
         call PWTIMERA(-1,tp,dtime)
         call PNTPOS3BX(g,h,nxh,ny,nz,kstrt,nyv,nzv,kxyp,kyzp,kzp,kxypd,
     1kyzpd,kzpd,jblok,mblok,lblok,ndim)
         call PWTIMERA(1,tp,dtime)
c start z fft tasks
         do 50 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXZ,nargs,h,isign,mixup,sct,
     1indx,indy,indz,kstrt,kyzip(i),kxypp,nzv,kxyp,kyzp,kxypd,kyzpd,jblo
     2k,mblok,ndim,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish z fft
         call PFFT32RNXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1ypl,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,ndim,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
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
c start z fft tasks
         do 70 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXZ,nargs,h,isign,mixup,sct,
     1indx,indy,indz,kstrt,kyzip(i),kxypp,nzv,kxyp,kyzp,kxypd,kyzpd,jblo
     2k,mblok,ndim,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish z fft
         call PFFT32RNXZ(h,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1ypl,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,ndim,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
c transpose h array to g
         call PWTIMERA(-1,tp,dtime)
         call PNTPOS3BX(h,g,nxh,nz,ny,kstrt,nzv,nyv,kxyp,kzp,kyzp,kxypd,
     1kzpd,kyzpd,jblok,lblok,mblok,ndim)
         call PWTIMERA(1,tp,dtime)
c start y fft tasks
         do 90 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXY,nargs,g,isign,mixup,sct,
     1indx,indy,indz,kstrt,kyzip(i),kxypp,nyv,kxyp,kzp,kxypd,kzpd,jblok,
     2lblok,ndim,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   90    continue
c finish y fft
         call PFFT32RNXY(g,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,kx
     1ypl,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,ndim,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 100 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
  100    continue
c transpose g array to f
         call PWTIMERA(-1,ttp,dtime)
         call PNTPOS3AX(g,f,ny,nxh,nz,kstrt,nyv,nxvh,kyp,kxyp,kzp,kypd,k
     1xypd,kzpd,kblok,jblok,lblok,ndim)
         call PWTIMERA(1,ttp,dtime)
c start x fft tasks
         do 110 i = 1, nmt
         kyzip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXX,margs,f,ss(1,1,i),isign,
     1mixup,sct,indx,indy,indz,kstrt,kyzip(i),kypp,nxvh,kyp,kzp,kypd,kzp
     2d,kblok,lblok,ndim,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  110    continue
c finish x fft
         call PFFT32RNXX(f,ss(1,1,nmtt),isign,mixup,sct,indx,indy,indz,k
     1strt,kypi,kypl,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,ndim,nxhyzd,nxyz
     2hd)
c wait for tasks to complete
         do 120 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  120    continue
      endif
      ttp = ttp + tp
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
c-----------------------------------------------------------------------
      subroutine MP2FFT32RN(f1,f2,g1,g2,h1,h2,bs,br,ss,isign,ntpose,mixu
     1p,sct,ttp,indx,indy,indz,kstrt,nxvh,nyv,nzv,kxyp,kyp,kyzp,kzp,kxyp
     2d,kypd,kyzpd,kzpd,kzyp,jblok,kblok,jkblok,lblok,mblok,mlblok,ndim1
     3,ndim2,nxhyzd,nxyzhd,kyzip,iftask,nmt,ierr)
c multi-tasking two real to complex ffts
c nyzip = initial index arrays for tasks
c iftask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      integer isign, ntpose, indx, indy, indz, kstrt, nxvh, nyv, nzv
      integer kxyp, kyp, kyzp, kzp, kxypd, kypd, kyzpd, kzpd, kzyp
      integer jblok, kblok, jkblok, lblok, mblok, mlblok, ndim1, ndim2
      integer nxhyzd, nxyzhd
      integer mixup, kyzip, iftask, nmt, ierr
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
      dimension ss(ndim1+ndim2,nxvh,nmt+1)
      dimension mixup(nxhyzd), sct(nxyzhd)
      dimension kyzip(nmt), iftask(nmt)
c local data
      integer i, nargs, margs, nx, ny, nz, nxh, nmtt, kypi, kxypi, kypp
      integer kxypp, kypl, kxypl
      real tp, tf
      double precision dtime
      external PFFT32RNXX, PFFT32RNXY, PFFT32RNXZ
      data nargs, margs /20,21/
      nx = 2**indx
      ny = 2**indy
      nz = 2**indz
      nxh = nx/2
      nmtt = nmt + 1
      kypp = kyp/nmtt
      kxypp = kxyp/nmtt
      kypi = kypp*nmt
      kxypi = kxypp*nmt
      kypl = kyp - kypi
      kxypl = kxyp - kxypi
      kypi = kypi + 1
      kxypi = kxypi + 1
      ierr = 0
c inverse fourier transform
      if (isign.lt.0) then
c start first x fft tasks
         do 10 i = 1, nmt
         kyzip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXX,margs,f1,ss(1,1,i),isign
     1,mixup,sct,indx,indy,indz,kstrt,kyzip(i),kypp,nxvh,kyp,kzp,kypd,kz
     2pd,kblok,lblok,ndim1,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   10    continue
c finish first x fft
         call PFFT32RNXX(f1,ss(1,1,nmtt),isign,mixup,sct,indx,indy,indz,
     1kstrt,kypi,kypl,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,ndim1,nxhyzd,nx
     2yzhd)
c wait for tasks to complete
         do 20 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   20    continue
c start second x fft tasks
         do 30 i = 1, nmt
         kyzip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXX,margs,f2,ss(1,1,i),isign
     1,mixup,sct,indx,indy,indz,kstrt,kyzip(i),kypp,nxvh,kyp,kzp,kypd,kz
     2pd,kblok,lblok,ndim2,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   30    continue
c finish second x fft
         call PFFT32RNXX(f2,ss(1,1,nmtt),isign,mixup,sct,indx,indy,indz,
     1kstrt,kypi,kypl,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,ndim2,nxhyzd,nx
     2yzhd)
c wait for tasks to complete
         do 40 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
   40    continue
c transpose f arrays to g
         call PWTIMERA(-1,ttp,dtime)
         call PN2TPOS3A(f1,f2,g1,g2,bs,br,nxh,ny,nz,kstrt,nxvh,nyv,kxyp,
     1kyp,kzp,kxypd,kypd,kzpd,jblok,kblok,lblok,ndim1,ndim2)
         call PWTIMERA(1,ttp,dtime)
c start first y fft tasks
         do 50 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXY,nargs,g1,isign,mixup,sct
     1,indx,indy,indz,kstrt,kyzip(i),kxypp,nyv,kxyp,kzp,kxypd,kzpd,jblok
     2,lblok,ndim1,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   50    continue
c finish first y fft
         call PFFT32RNXY(g1,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,k
     1xypl,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,ndim1,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 60 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   60    continue
c start second y fft tasks
         do 70 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXY,nargs,g2,isign,mixup,sct
     1,indx,indy,indz,kstrt,kyzip(i),kxypp,nyv,kxyp,kzp,kxypd,kzpd,jblok
     2,lblok,ndim2,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   70    continue
c finish second y fft
         call PFFT32RNXY(g2,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,k
     1xypl,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,ndim2,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 80 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
   80    continue
c transpose g arrays to h
         call PWTIMERA(-1,tp,dtime)
         call PN2TPOS3B(g1,g2,h1,h2,bs,br,nxh,ny,nz,kstrt,nyv,nzv,kxyp,k
     1yzp,kzp,kxypd,kyzpd,kzpd,jblok,mblok,lblok,ndim1,ndim2)
         call PWTIMERA(1,tp,dtime)
c start first z fft tasks
         do 90 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXZ,nargs,h1,isign,mixup,sct
     1,indx,indy,indz,kstrt,kyzip(i),kxypp,nzv,kxyp,kyzp,kxypd,kyzpd,jbl
     2ok,mblok,ndim1,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
   90    continue
c finish first z fft
         call PFFT32RNXZ(h1,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,k
     1xypl,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,ndim1,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 100 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
  100    continue
c start second z fft tasks
         do 110 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXZ,nargs,h2,isign,mixup,sct
     1,indx,indy,indz,kstrt,kyzip(i),kxypp,nzv,kxyp,kyzp,kxypd,kyzpd,jbl
     2ok,mblok,ndim2,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  110    continue
c finish second z fft
         call PFFT32RNXZ(h2,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,k
     1xypl,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,ndim2,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 120 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
  120    continue
c transpose h arrays to f
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PN2TPOS3B(h1,h2,g1,g2,br,bs,nxh,nz,ny,kstrt,nzv,nyv,kxy
     1p,kzp,kyzp,kxypd,kzpd,kyzpd,jblok,lblok,mblok,ndim1,ndim2)
            call PNTPOS3A(g1,g2,f1,f2,br,bs,ny,nxh,nz,kstrt,nyv,nxvh,kyp
     1,kxyp,kzp,kypd,kxypd,kzpd,kblok,jblok,lblok,ndim1,ndim2)
            call PWTIMERA(1,tf,dtime)
         endif
c forward fourier transform
      else if (isign.gt.0) then
c transpose f arrays to h
         if (ntpose.eq.0) then
            call PWTIMERA(-1,tf,dtime)
            call PN2TPOS3A(f1,f2,g1,g2,bs,br,nxh,ny,nz,kstrt,nxvh,nyv,kx
     1yp,kyp,kzp,kxypd,kypd,kzpd,jblok,kblok,lblok,ndim1,ndim2)
            call PNTPOS3B(g1,g2,h1,h2,bs,br,nxh,ny,nz,kstrt,nyv,nzv,kxyp
     1,kyzp,kzp,kxypd,kyzpd,kzpd,jblok,mblok,lblok,ndim1,ndim2)
            call PWTIMERA(1,tf,dtime)
         endif
c start first z fft tasks
         do 130 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXZ,nargs,h1,isign,mixup,sct
     1,indx,indy,indz,kstrt,kyzip(i),kxypp,nzv,kxyp,kyzp,kxypd,kyzpd,jbl
     2ok,mblok,ndim1,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  130    continue
c finish first z fft
         call PFFT32RNXZ(h1,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,k
     1xypl,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,ndim1,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 140 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
  140    continue
c start second z fft tasks
         do 150 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXZ,nargs,h2,isign,mixup,sct
     1,indx,indy,indz,kstrt,kyzip(i),kxypp,nzv,kxyp,kyzp,kxypd,kyzpd,jbl
     2ok,mblok,ndim2,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  150    continue
c finish second z fft
         call PFFT32RNXZ(h2,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,k
     1xypl,nzv,kxyp,kyzp,kxypd,kyzpd,jblok,mblok,ndim2,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 160 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
  160    continue
c transpose h arrays to g
         call PWTIMERA(-1,tp,dtime)
         call PN2TPOS3B(h1,h2,g1,g2,br,bs,nxh,nz,ny,kstrt,nzv,nyv,kxyp,k
     1zp,kyzp,kxypd,kzpd,kyzpd,jblok,lblok,mblok,ndim1,ndim2)
         call PWTIMERA(1,tp,dtime)
c start first y fft tasks
         do 170 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXY,nargs,g1,isign,mixup,sct
     1,indx,indy,indz,kstrt,kyzip(i),kxypp,nyv,kxyp,kzp,kxypd,kzpd,jblok
     2,lblok,ndim1,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  170    continue
c finish first y fft
         call PFFT32RNXY(g1,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,k
     1xypl,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,ndim1,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 180 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
  180    continue
c start second y fft tasks
         do 190 i = 1, nmt
         kyzip(i) = kxypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXY,nargs,g2,isign,mixup,sct
     1,indx,indy,indz,kstrt,kyzip(i),kxypp,nyv,kxyp,kzp,kxypd,kzpd,jblok
     2,lblok,ndim2,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  190    continue
c finish second y fft
         call PFFT32RNXY(g2,isign,mixup,sct,indx,indy,indz,kstrt,kxypi,k
     1xypl,nyv,kxyp,kzp,kxypd,kzpd,jblok,lblok,ndim2,nxhyzd,nxyzhd)
c wait for tasks to complete
         do 200 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) ierr = -2
  200    continue
c transpose g arrays to f
         call PWTIMERA(-1,ttp,dtime)
         call PN2TPOS3A(g1,g2,f1,f2,br,bs,ny,nxh,nz,kstrt,nyv,nxvh,kyp,k
     1xyp,kzp,kypd,kxypd,kzpd,kblok,jblok,lblok,ndim1,ndim2)
         call PWTIMERA(1,ttp,dtime)
c start first x fft tasks
         do 210 i = 1, nmt
         kyzip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXX,margs,f1,ss(1,1,i),isign
     1,mixup,sct,indx,indy,indz,kstrt,kyzip(i),kypp,nxvh,kyp,kzp,kypd,kz
     2pd,kblok,lblok,ndim1,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  210    continue
c finish first x fft
         call PFFT32RNXX(f1,ss(1,1,nmtt),isign,mixup,sct,indx,indy,indz,
     1kstrt,kypi,kypl,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,ndim1,nxhyzd,nx
     2yzhd)
c wait for tasks to complete
         do 220 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  220    continue
c start second x fft tasks
         do 230 i = 1, nmt
         kyzip(i) = kypp*(i - 1) + 1
         call MP_TASKSTART(iftask(i),PFFT32RNXX,margs,f2,ss(1,1,i),isign
     1,mixup,sct,indx,indy,indz,kstrt,kyzip(i),kypp,nxvh,kyp,kzp,kypd,kz
     2pd,kblok,lblok,ndim2,nxhyzd,nxyzhd)
c check for errors
         if (iftask(i).eq.0) then
           ierr = -1
           return
         endif
  230    continue
c finish second x fft
         call PFFT32RNXX(f2,ss(1,1,nmtt),isign,mixup,sct,indx,indy,indz,
     1kstrt,kypi,kypl,nxvh,kyp,kzp,kypd,kzpd,kblok,lblok,ndim2,nxhyzd,nx
     2yzhd)
c wait for tasks to complete
         do 240 i = 1, nmt
         call MP_TASKWAIT(iftask(i))
c check for errors
         if (iftask(i).ne.0) then
            ierr = -2
            return
         endif
  240    continue
      endif
      ttp = ttp + tp
      if (ntpose.eq.0) ttp = ttp + tf
      return
      end
