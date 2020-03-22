c 3d PIC parallel multi-tasking library for pushing particles
c with magnetic field and depositing current
c written by viktor k. decyk, ucla
c copyright 2000, regents of the university of california
c update: september 28, 2002
c-----------------------------------------------------------------------
      subroutine MPJDOST32(part,cux,cuy,cuz,npp,nps,noff,qm,dt,nx,idimp,
     1npmax,mnblok,nxv,nypmx,nzpmx,idds,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cux, cuy, cuz, qm, dt, cup
      integer npp, nps, noff
      integer nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), cux(nxv,nypmx,nzpmx,mnblok)
      dimension cuy(nxv,nypmx,nzpmx,mnblok), cuz(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension cup(nxv,nypmx,nzpmx,3,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m, n
      external PJDOST32
      data nargs /16/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 70 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 60 m = 1, mnblok
      do 50 n = 1, 3
      nps(m) = npt
      do 40 l = 1, nzpmx
      do 30 k = 1, nypmx
      do 20 j = 1, nx
      cup(j,k,l,n,m,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
   60 continue
      call MP_TASKSTART(idtask(i),PJDOST32,nargs,part(1,npo,1),cup(1,1,1
     1,1,1,i),cup(1,1,1,2,1,i),cup(1,1,1,3,1,i),nps,noff,qm,dt,nx,idimp,
     2npmax,mnblok,nxv,nypmx,nzpmx,idds)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   70 continue
c deposit remaining current
      npo = npl + 1
      do 80 m = 1, mnblok
      npp(m) = npp(m) - npl
   80 continue
      call PJDOST32(part(1,npo,1),cux,cuy,cuz,npp,noff,qm,dt,nx,idimp,np
     1max,mnblok,nxv,nypmx,nzpmx,idds)
      do 90 m = 1, mnblok
      npp(m) = npp(m) + npl
   90 continue
c wait for tasks to complete
      do 140 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 130 m = 1, mnblok
      do 120 l = 1, nzpmx
      do 110 k = 1, nypmx
      do 100 j = 1, nx
      cux(j,k,l,m) = cux(j,k,l,m) + cup(j,k,l,1,m,i)
      cuy(j,k,l,m) = cuy(j,k,l,m) + cup(j,k,l,2,m,i)
      cuz(j,k,l,m) = cuz(j,k,l,m) + cup(j,k,l,3,m,i)
  100 continue
  110 continue
  120 continue
  130 continue
  140 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGJPOST32(part,cu,npp,nps,noff,qm,dt,nx,ny,nz,idimp,np
     1max,mnblok,nxv,nypmx,nzpmx,idds,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer npp, nps, noff
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), cu(3,nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension cup(3,nxv,nypmx,nzpmx,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m, n
      external PGJPOST32
      data nargs /17/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 70 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 60 m = 1, mnblok
      nps(m) = npt
      do 50 l = 1, nzpmx
      do 40 k = 1, nypmx
      do 30 j = 1, nxv
      do 20 n = 1, 3
      cup(n,j,k,l,m,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
   60 continue
      call MP_TASKSTART(idtask(i),PGJPOST32,nargs,part(1,npo,1),cup(1,1,
     11,1,1,i),nps,noff,qm,dt,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpm
     2x,idds,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   70 continue
c deposit remaining current
      npo = npl + 1
      do 80 m = 1, mnblok
      npp(m) = npp(m) - npl
   80 continue
      call PGJPOST32(part(1,npo,1),cu,npp,noff,qm,dt,nx,ny,nz,idimp,npma
     1x,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
      do 90 m = 1, mnblok
      npp(m) = npp(m) + npl
   90 continue
c wait for tasks to complete
      do 150 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 140 m = 1, mnblok
      do 130 l = 1, nzpmx
      do 120 k = 1, nypmx
      do 110 j = 1, nxv
      do 100 n = 1, 3
      cu(n,j,k,l,m) = cu(n,j,k,l,m) + cup(n,j,k,l,m,i)
  100 continue
  110 continue
  120 continue
  130 continue
  140 continue
  150 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSJPOST32(part,cu,npp,nps,noff,qm,dt,nx,ny,nz,idimp,n
     1pmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer npp, nps, noff
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp, idds
      integer ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), cu(3,nxyzp,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension cup(3,nxyzp,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, m, n
      external PGSJPOST32
      data nargs /17/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 50 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 40 m = 1, mnblok
      nps(m) = npt
      do 30 j = 1, nxyzp
      do 20 n = 1, 3
      cup(n,j,m,i) = 0.
   20 continue
   30 continue
   40 continue
      call MP_TASKSTART(idtask(i),PGSJPOST32,nargs,part(1,npo,1),cup(1,1
     1,1,i),nps,noff,qm,dt,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,i
     2dds,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   50 continue
c deposit remaining current
      npo = npl + 1
      do 60 m = 1, mnblok
      npp(m) = npp(m) - npl
   60 continue
      call PGSJPOST32(part(1,npo,1),cu,npp,noff,qm,dt,nx,ny,nz,idimp,npm
     1ax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
      do 70 m = 1, mnblok
      npp(m) = npp(m) + npl
   70 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 100 m = 1, mnblok
      do 90 j = 1, nxyzp
      do 80 n = 1, 3
      cu(n,j,m) = cu(n,j,m) + cup(n,j,m,i)
   80 continue
   90 continue
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPSJOST32X(part,cu,npp,nps,noff,nn,amxyz,qm,dt,nx,idimp
     1,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n81,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, amxyz, qm, dt, cup
      integer npp, nps, noff, nn
      integer nx, idimp, npmax, mnblok, nxv, nypmx, nxvyzp, idds
      integer npd, n81
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), cu(3*nxvyzp,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension nn(n81,npd,mnblok,nmt+1), amxyz(n81,npd,mnblok,nmt+1)
      dimension cup(3*nxvyzp,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, m
      external PSJOST32X
      data nargs /18/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 30 m = 1, mnblok
      nps(m) = npt
      do 20 j = 1, 3*nxvyzp
      cup(j,m,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),PSJOST32X,nargs,part(1,npo,1),cup(1,1,
     1i),nps,noff,nn(1,1,1,i+1),amxyz(1,1,1,i+1),qm,dt,nx,idimp,npmax,mn
     2blok,nxv,nypmx,nxvyzp,idds,npd,n81)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      do 50 m = 1, mnblok
      npp(m) = npp(m) - npl
   50 continue
      call PSJOST32X(part(1,npo,1),cu,npp,noff,nn,amxyz,qm,dt,nx,idimp,n
     1pmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n81)
      do 60 m = 1, mnblok
      npp(m) = npp(m) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 80 m = 1, mnblok
      do 70 j = 1, 3*nxvyzp
      cu(j,m) = cu(j,m) + cup(j,m,i)
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSJOST32X(part,cu,npp,nps,noff,nn,amxyz,qm,dt,nx,ny,n
     1z,idimp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n81,ipbc,cup,idtask
     2,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, amxyz, qm, dt, cup
      integer npp, nps, noff, nn
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxvyzp, ipbc
      integer idds, npd, n81, idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), cu(3*nxvyzp,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension nn(n81,npd,mnblok,nmt+1), amxyz(n81,npd,mnblok,nmt+1)
      dimension cup(3*nxvyzp,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, m
      external PGSJOST32X
      data nargs /21/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 30 m = 1, mnblok
      nps(m) = npt
      do 20 j = 1, 3*nxvyzp
      cup(j,m,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),PGSJOST32X,nargs,part(1,npo,1),cup(1,1
     1,i),nps,noff,nn(1,1,1,i+1),amxyz(1,1,1,i+1),qm,dt,nx,ny,nz,idimp,n
     2pmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n81,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      do 50 m = 1, mnblok
      npp(m) = npp(m) - npl
   50 continue
      call PGSJOST32X(part(1,npo,1),cu,npp,noff,nn,amxyz,qm,dt,nx,ny,nz,
     1idimp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n81,ipbc)
      do 60 m = 1, mnblok
      npp(m) = npp(m) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 80 m = 1, mnblok
      do 70 j = 1, 3*nxvyzp
      cu(j,m) = cu(j,m) + cup(j,m,i)
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPJDOST32L(part,cux,cuy,cuz,npp,nps,noff,qm,dt,nx,idimp
     1,npmax,mnblok,nxv,nypmx,nzpmx,idds,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cux, cuy, cuz, qm, dt, cup
      integer npp, nps, noff
      integer nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), cux(nxv,nypmx,nzpmx,mnblok)
      dimension cuy(nxv,nypmx,nzpmx,mnblok), cuz(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension cup(nxv,nypmx,nzpmx,3,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m, n
      external PJDOST32L
      data nargs /16/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 70 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 60 m = 1, mnblok
      do 50 n = 1, 3
      nps(m) = npt
      do 40 l = 1, nzpmx
      do 30 k = 1, nypmx
      do 20 j = 1, nx
      cup(j,k,l,n,m,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
   60 continue
      call MP_TASKSTART(idtask(i),PJDOST32L,nargs,part(1,npo,1),cup(1,1,
     11,1,1,i),cup(1,1,1,2,1,i),cup(1,1,1,3,1,i),nps,noff,qm,dt,nx,idimp
     2,npmax,mnblok,nxv,nypmx,nzpmx,idds)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   70 continue
c deposit remaining current
      npo = npl + 1
      do 80 m = 1, mnblok
      npp(m) = npp(m) - npl
   80 continue
      call PJDOST32L(part(1,npo,1),cux,cuy,cuz,npp,noff,qm,dt,nx,idimp,n
     1pmax,mnblok,nxv,nypmx,nzpmx,idds)
      do 90 m = 1, mnblok
      npp(m) = npp(m) + npl
   90 continue
c wait for tasks to complete
      do 140 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 130 m = 1, mnblok
      do 120 l = 1, nzpmx
      do 110 k = 1, nypmx
      do 100 j = 1, nx
      cux(j,k,l,m) = cux(j,k,l,m) + cup(j,k,l,1,m,i)
      cuy(j,k,l,m) = cuy(j,k,l,m) + cup(j,k,l,2,m,i)
      cuz(j,k,l,m) = cuz(j,k,l,m) + cup(j,k,l,3,m,i)
  100 continue
  110 continue
  120 continue
  130 continue
  140 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGJPOST32L(part,cu,npp,nps,noff,qm,dt,nx,ny,nz,idimp,n
     1pmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer npp, nps, noff
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), cu(3,nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension cup(3,nxv,nypmx,nzpmx,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m, n
      external PGJPOST32L
      data nargs /17/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 70 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 60 m = 1, mnblok
      nps(m) = npt
      do 50 l = 1, nzpmx
      do 40 k = 1, nypmx
      do 30 j = 1, nxv
      do 20 n = 1, 3
      cup(n,j,k,l,m,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
   60 continue
      call MP_TASKSTART(idtask(i),PGJPOST32L,nargs,part(1,npo,1),cup(1,1
     1,1,1,1,i),nps,noff,qm,dt,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzp
     2mx,idds,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   70 continue
c deposit remaining current
      npo = npl + 1
      do 80 m = 1, mnblok
      npp(m) = npp(m) - npl
   80 continue
      call PGJPOST32L(part(1,npo,1),cu,npp,noff,qm,dt,nx,ny,nz,idimp,npm
     1ax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
      do 90 m = 1, mnblok
      npp(m) = npp(m) + npl
   90 continue
c wait for tasks to complete
      do 150 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 140 m = 1, mnblok
      do 130 l = 1, nzpmx
      do 120 k = 1, nypmx
      do 110 j = 1, nxv
      do 100 n = 1, 3
      cu(n,j,k,l,m) = cu(n,j,k,l,m) + cup(n,j,k,l,m,i)
  100 continue
  110 continue
  120 continue
  130 continue
  140 continue
  150 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSJPOST32L(part,cu,npp,nps,noff,qm,dt,nx,ny,nz,idimp,
     1npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, cup
      integer npp, nps, noff
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp, idds
      integer ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), cu(3,nxyzp,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension cup(3,nxyzp,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, m, n
      external PGSJPOST32L
      data nargs /17/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 50 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 40 m = 1, mnblok
      nps(m) = npt
      do 30 j = 1, nxyzp
      do 20 n = 1, 3
      cup(n,j,m,i) = 0.
   20 continue
   30 continue
   40 continue
      call MP_TASKSTART(idtask(i),PGSJPOST32L,nargs,part(1,npo,1),cup(1,
     11,1,i),nps,noff,qm,dt,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,
     2idds,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   50 continue
c deposit remaining current
      npo = npl + 1
      do 60 m = 1, mnblok
      npp(m) = npp(m) - npl
   60 continue
      call PGSJPOST32L(part(1,npo,1),cu,npp,noff,qm,dt,nx,ny,nz,idimp,np
     1max,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
      do 70 m = 1, mnblok
      npp(m) = npp(m) + npl
   70 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 100 m = 1, mnblok
      do 90 j = 1, nxyzp
      do 80 n = 1, 3
      cu(n,j,m) = cu(n,j,m) + cup(n,j,m,i)
   80 continue
   90 continue
  100 continue
  110 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPSJOST32XL(part,cu,npp,nps,noff,nn,amxyz,qm,dt,nx,idim
     1p,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n24,cup,idtask,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, amxyz, qm, dt, cup
      integer npp, nps, noff, nn
      integer nx, idimp, npmax, mnblok, nxv, nypmx, nxvyzp, idds
      integer npd, n24
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), cu(3*nxvyzp,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension nn(n24,npd,mnblok,nmt+1), amxyz(n24,npd,mnblok,nmt+1)
      dimension cup(3*nxvyzp,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, m
      external PSJOST32XL
      data nargs /18/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 30 m = 1, mnblok
      nps(m) = npt
      do 20 j = 1, 3*nxvyzp
      cup(j,m,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),PSJOST32XL,nargs,part(1,npo,1),cup(1,1
     1,i),nps,noff,nn(1,1,1,i+1),amxyz(1,1,1,i+1),qm,dt,nx,idimp,npmax,m
     2nblok,nxv,nypmx,nxvyzp,idds,npd,n24)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      do 50 m = 1, mnblok
      npp(m) = npp(m) - npl
   50 continue
      call PSJOST32XL(part(1,npo,1),cu,npp,noff,nn,amxyz,qm,dt,nx,idimp,
     1npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n24)
      do 60 m = 1, mnblok
      npp(m) = npp(m) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 80 m = 1, mnblok
      do 70 j = 1, 3*nxvyzp
      cu(j,m) = cu(j,m) + cup(j,m,i)
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSJOST32XL(part,cu,npp,nps,noff,nn,amxyz,qm,dt,nx,ny,
     1nz,idimp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n24,ipbc,cup,idtas
     2k,nmt,ierr)
c parallel multitasking current deposition
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, amxyz, qm, dt, cup
      integer npp, nps, noff, nn
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxvyzp, ipbc
      integer idds, npd, n24, idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), cu(3*nxvyzp,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension nn(n24,npd,mnblok,nmt+1), amxyz(n24,npd,mnblok,nmt+1)
      dimension cup(3*nxvyzp,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, m
      external PGSJOST32XL
      data nargs /21/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start current deposit tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear current arrays
      do 30 m = 1, mnblok
      nps(m) = npt
      do 20 j = 1, 3*nxvyzp
      cup(j,m,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),PGSJOST32XL,nargs,part(1,npo,1),cup(1,
     11,i),nps,noff,nn(1,1,1,i+1),amxyz(1,1,1,i+1),qm,dt,nx,ny,nz,idimp,
     2npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n24,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining current
      npo = npl + 1
      do 50 m = 1, mnblok
      npp(m) = npp(m) - npl
   50 continue
      call PGSJOST32XL(part(1,npo,1),cu,npp,noff,nn,amxyz,qm,dt,nx,ny,nz
     1,idimp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n24,ipbc)
      do 60 m = 1, mnblok
      npp(m) = npp(m) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum current arrays
      do 80 m = 1, mnblok
      do 70 j = 1, 3*nxvyzp
      cu(j,m) = cu(j,m) + cup(j,m,i)
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPBPUSH32(part,fx,fy,fz,bx,by,bz,npp,nps,noff,qbm,dt,ek
     1,nx,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ekp,idtask,nmt,ierr)
c parallel multitasking particle push with magnetic field
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, fy, fz, bx, by, bz, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok)
      dimension fx(nxv,nypmx,nzpmx,mnblok), fy(nxv,nypmx,nzpmx,mnblok)
      dimension fz(nxv,nypmx,nzpmx,mnblok), bx(nxv,nypmx,nzpmx,mnblok)
      dimension by(nxv,nypmx,nzpmx,mnblok), bz(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, m
      external PBPUSH32
      data nargs /20/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 m = 1, mnblok
      nps(m) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PBPUSH32,nargs,part(1,npo,1),fx,fy,fz,
     1bx,by,bz,nps,noff,qbm,dt,ekp(i),nx,idimp,npmax,mnblok,nxv,nypmx,nz
     2pmx,idds)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 m = 1, mnblok
      npp(m) = npp(m) - npl
   40 continue
      call PBPUSH32(part(1,npo,1),fx,fy,fz,bx,by,bz,npp,noff,qbm,dt,ek,n
     1x,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds)
      do 50 m = 1, mnblok
      npp(m) = npp(m) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGBPUSH32(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,ek,nx
     1,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask,nmt
     2,ierr)
c parallel multitasking particle push with magnetic field
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxyz, bxyz, qbm, dt, dtc, ek, ekp
      integer npp, nps, noff
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok)
      dimension fxyz(3,nxv,nypmx,nzpmx,mnblok)
      dimension bxyz(3,nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, m
      external PGBPUSH32
      data nargs /20/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 m = 1, mnblok
      nps(m) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGBPUSH32,nargs,part(1,npo,1),fxyz,bxy
     1z,nps,noff,qbm,dt,dtc,ekp(i),nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx
     2,nzpmx,idds,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 m = 1, mnblok
      npp(m) = npp(m) - npl
   40 continue
      call PGBPUSH32(part(1,npo,1),fxyz,bxyz,npp,noff,qbm,dt,dtc,ek,nx,n
     1y,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
      do 50 m = 1, mnblok
      npp(m) = npp(m) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSBPUSH32(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,ek,n
     1x,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,idtask,nm
     2t,ierr)
c parallel multitasking particle push with magnetic field
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxyz, bxyz, qbm, dt, dtc, ek, ekp
      integer npp, nps, noff
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp, idds
      integer ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok)
      dimension fxyz(3,nxyzp,mnblok), bxyz(3,nxyzp,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, m
      external PGSBPUSH32
      data nargs /20/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 m = 1, mnblok
      nps(m) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGSBPUSH32,nargs,part(1,npo,1),fxyz,bx
     1yz,nps,noff,qbm,dt,dtc,ekp(i),nx,ny,nz,idimp,npmax,mnblok,nxv,nypm
     2x,nxyzp,idds,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 m = 1, mnblok
      npp(m) = npp(m) - npl
   40 continue
      call PGSBPUSH32(part(1,npo,1),fxyz,bxyz,npp,noff,qbm,dt,dtc,ek,nx,
     1ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
      do 50 m = 1, mnblok
      npp(m) = npp(m) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPBPUSH32L(part,fx,fy,fz,bx,by,bz,npp,nps,noff,qbm,dt,e
     1k,nx,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ekp,idtask,nmt,ierr)
c parallel multitasking particle push with magnetic field
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, fy, fz, bx, by, bz, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok)
      dimension fx(nxv,nypmx,nzpmx,mnblok), fy(nxv,nypmx,nzpmx,mnblok)
      dimension fz(nxv,nypmx,nzpmx,mnblok), bx(nxv,nypmx,nzpmx,mnblok)
      dimension by(nxv,nypmx,nzpmx,mnblok), bz(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, m
      external PBPUSH32L
      data nargs /20/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 m = 1, mnblok
      nps(m) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PBPUSH32L,nargs,part(1,npo,1),fx,fy,fz
     1,bx,by,bz,nps,noff,qbm,dt,ekp(i),nx,idimp,npmax,mnblok,nxv,nypmx,n
     2zpmx,idds)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 m = 1, mnblok
      npp(m) = npp(m) - npl
   40 continue
      call PBPUSH32L(part(1,npo,1),fx,fy,fz,bx,by,bz,npp,noff,qbm,dt,ek,
     1nx,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds)
      do 50 m = 1, mnblok
      npp(m) = npp(m) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGBPUSH32L(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,ek,n
     1x,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask,nm
     2t,ierr)
c parallel multitasking particle push with magnetic field
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxyz, bxyz, qbm, dt, dtc, ek, ekp
      integer npp, nps, noff
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok)
      dimension fxyz(3,nxv,nypmx,nzpmx,mnblok)
      dimension bxyz(3,nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, m
      external PGBPUSH32L
      data nargs /20/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 m = 1, mnblok
      nps(m) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGBPUSH32L,nargs,part(1,npo,1),fxyz,bx
     1yz,nps,noff,qbm,dt,dtc,ekp(i),nx,ny,nz,idimp,npmax,mnblok,nxv,nypm
     2x,nzpmx,idds,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 m = 1, mnblok
      npp(m) = npp(m) - npl
   40 continue
      call PGBPUSH32L(part(1,npo,1),fxyz,bxyz,npp,noff,qbm,dt,dtc,ek,nx,
     1ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
      do 50 m = 1, mnblok
      npp(m) = npp(m) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSBPUSH32L(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,ek,
     1nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,idtask,n
     2mt,ierr)
c parallel multitasking particle push with magnetic field
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxyz, bxyz, qbm, dt, dtc, ek, ekp
      integer npp, nps, noff
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp, idds
      integer ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok)
      dimension fxyz(3,nxyzp,mnblok), bxyz(3,nxyzp,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, m
      external PGSBPUSH32L
      data nargs /20/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 m = 1, mnblok
      nps(m) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGSBPUSH32L,nargs,part(1,npo,1),fxyz,b
     1xyz,nps,noff,qbm,dt,dtc,ekp(i),nx,ny,nz,idimp,npmax,mnblok,nxv,nyp
     2mx,nxyzp,idds,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 m = 1, mnblok
      npp(m) = npp(m) - npl
   40 continue
      call PGSBPUSH32L(part(1,npo,1),fxyz,bxyz,npp,noff,qbm,dt,dtc,ek,nx
     1,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
      do 50 m = 1, mnblok
      npp(m) = npp(m) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPBPUSH32C(part,fx,fy,fz,bx,by,bz,npp,nps,noff,qbm,dt,e
     1k,nx,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ekp,idtask,nmt,ierr)
c parallel multitasking particle push with magnetic field
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, fy, fz, bx, by, bz, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok)
      dimension fx(nxv,nypmx,nzpmx,mnblok), fy(nxv,nypmx,nzpmx,mnblok)
      dimension fz(nxv,nypmx,nzpmx,mnblok), bx(nxv,nypmx,nzpmx,mnblok)
      dimension by(nxv,nypmx,nzpmx,mnblok), bz(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, m
      external PBPUSH32C
      data nargs /20/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 m = 1, mnblok
      nps(m) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PBPUSH32C,nargs,part(1,npo,1),fx,fy,fz
     1,bx,by,bz,nps,noff,qbm,dt,ekp(i),nx,idimp,npmax,mnblok,nxv,nypmx,n
     2zpmx,idds)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 m = 1, mnblok
      npp(m) = npp(m) - npl
   40 continue
      call PBPUSH32C(part(1,npo,1),fx,fy,fz,bx,by,bz,npp,noff,qbm,dt,ek,
     1nx,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds)
      do 50 m = 1, mnblok
      npp(m) = npp(m) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGBPUSH32C(part,fxyz,bxyz,npp,nps,noff,qbm,dt,ek,nx,ny
     1,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask,nmt,ie
     2rr)
c parallel multitasking particle push with magnetic field
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxyz, bxyz, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok)
      dimension fxyz(3,nxv,nypmx,nzpmx,mnblok)
      dimension bxyz(3,nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, m
      external PGBPUSH32C
      data nargs /19/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 m = 1, mnblok
      nps(m) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGBPUSH32C,nargs,part(1,npo,1),fxyz,bx
     1yz,nps,noff,qbm,dt,ekp(i),nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nz
     2pmx,idds,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 m = 1, mnblok
      npp(m) = npp(m) - npl
   40 continue
      call PGBPUSH32C(part(1,npo,1),fxyz,bxyz,npp,noff,qbm,dt,ek,nx,ny,n
     1z,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
      do 50 m = 1, mnblok
      npp(m) = npp(m) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPBPUSH32CL(part,fx,fy,fz,bx,by,bz,npp,nps,noff,qbm,dt,
     1ek,nx,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ekp,idtask,nmt,ierr)
c parallel multitasking particle push with magnetic field
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, fy, fz, bx, by, bz, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok)
      dimension fx(nxv,nypmx,nzpmx,mnblok), fy(nxv,nypmx,nzpmx,mnblok)
      dimension fz(nxv,nypmx,nzpmx,mnblok), bx(nxv,nypmx,nzpmx,mnblok)
      dimension by(nxv,nypmx,nzpmx,mnblok), bz(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, m
      external PBPUSH32CL
      data nargs /20/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 m = 1, mnblok
      nps(m) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PBPUSH32CL,nargs,part(1,npo,1),fx,fy,f
     1z,bx,by,bz,nps,noff,qbm,dt,ekp(i),nx,idimp,npmax,mnblok,nxv,nypmx,
     2nzpmx,idds)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 m = 1, mnblok
      npp(m) = npp(m) - npl
   40 continue
      call PBPUSH32CL(part(1,npo,1),fx,fy,fz,bx,by,bz,npp,noff,qbm,dt,ek
     1,nx,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds)
      do 50 m = 1, mnblok
      npp(m) = npp(m) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGBPUSH32CL(part,fxyz,bxyz,npp,nps,noff,qbm,dt,ek,nx,n
     1y,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask,nmt,i
     2err)
c parallel multitasking particle push with magnetic field
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxyz, bxyz, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok)
      dimension fxyz(3,nxv,nypmx,nzpmx,mnblok)
      dimension bxyz(3,nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, m
      external PGBPUSH32CL
      data nargs /19/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle push tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 m = 1, mnblok
      nps(m) = npt
   20 continue
      ekp(i) = 0.
      call MP_TASKSTART(idtask(i),PGBPUSH32CL,nargs,part(1,npo,1),fxyz,b
     1xyz,nps,noff,qbm,dt,ekp(i),nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,n
     2zpmx,idds,ipbc)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c push remaining particles
      npo = npl + 1
      do 40 m = 1, mnblok
      npp(m) = npp(m) - npl
   40 continue
      call PGBPUSH32CL(part(1,npo,1),fxyz,bxyz,npp,noff,qbm,dt,ek,nx,ny,
     1nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
      do 50 m = 1, mnblok
      npp(m) = npp(m) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
      ek = ek + ekp(i)
   60 continue
      return
      end
