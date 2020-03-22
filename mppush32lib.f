c 3d PIC parallel multi-tasking library for pushing particles
c and depositing charge
c written by viktor k. decyk, ucla
c copyright 2000, regents of the university of california
c update: june 11, 2005
c-----------------------------------------------------------------------
      subroutine MPDOST32(part,q,npp,nps,noff,qm,nx,idimp,npmax,mnblok,n
     1xv,nypmx,nzpmx,idds,qp,idtask,nmt,ierr)
c parallel multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer npp, nps, noff
      integer nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), q(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension qp(nxv,nypmx,nzpmx,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PDOST32
      data nargs /13/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start charge deposit tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear charge arrays
      do 50 m = 1, mnblok
      nps(m) = npt
      do 40 l = 1, nzpmx
      do 30 k = 1, nypmx
      do 20 j = 1, nx
      qp(j,k,l,m,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PDOST32,nargs,part(1,npo,1),qp(1,1,1,1
     1,i),nps,noff,qm,nx,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c deposit remaining charge
      npo = npl + 1
      do 70 m = 1, mnblok
      npp(m) = npp(m) - npl
   70 continue
      call PDOST32(part(1,npo,1),q,npp,noff,qm,nx,idimp,npmax,mnblok,nxv
     1,nypmx,nzpmx,idds)
      do 80 m = 1, mnblok
      npp(m) = npp(m) + npl
   80 continue
c wait for tasks to complete
      do 130 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 120 m = 1, mnblok
      do 110 l = 1, nzpmx
      do 100 k = 1, nypmx
      do 90 j = 1, nx
      q(j,k,l,m) = q(j,k,l,m) + qp(j,k,l,m,i)
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGPOST32(part,q,npp,nps,noff,qm,idimp,npmax,mnblok,nxv
     1,nypmx,nzpmx,idds,qp,idtask,nmt,ierr)
c parallel multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer npp, nps, noff
      integer idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), q(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension qp(nxv,nypmx,nzpmx,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGPOST32
      data nargs /12/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start charge deposit tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear charge arrays
      do 50 m = 1, mnblok
      nps(m) = npt
      do 40 l = 1, nzpmx
      do 30 k = 1, nypmx
      do 20 j = 1, nxv
      qp(j,k,l,m,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PGPOST32,nargs,part(1,npo,1),qp(1,1,1,
     11,i),nps,noff,qm,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c deposit remaining charge
      npo = npl + 1
      do 70 m = 1, mnblok
      npp(m) = npp(m) - npl
   70 continue
      call PGPOST32(part(1,npo,1),q,npp,noff,qm,idimp,npmax,mnblok,nxv,n
     1ypmx,nzpmx,idds)
      do 80 m = 1, mnblok
      npp(m) = npp(m) + npl
   80 continue
c wait for tasks to complete
      do 130 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 120 m = 1, mnblok
      do 110 l = 1, nzpmx
      do 100 k = 1, nypmx
      do 90 j = 1, nxv
      q(j,k,l,m) = q(j,k,l,m) + qp(j,k,l,m,i)
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSPOST32(part,q,npp,nps,noff,qm,idimp,npmax,mnblok,nx
     1v,nypmx,nxyzp,idds,qp,idtask,nmt,ierr)
c parallel multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer npp, nps, noff
      integer idimp, npmax, mnblok, nxv, nypmx, nxyzp, idds
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), q(nxyzp,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension qp(nxyzp,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, m
      external PGSPOST32
      data nargs /12/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start charge deposit tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear charge arrays
      do 30 m = 1, mnblok
      nps(m) = npt
      do 20 j = 1, nxyzp
      qp(j,m,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),PGSPOST32,nargs,part(1,npo,1),qp(1,1,i
     1),nps,noff,qm,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining charge
      npo = npl + 1
      do 50 m = 1, mnblok
      npp(m) = npp(m) - npl
   50 continue
      call PGSPOST32(part(1,npo,1),q,npp,noff,qm,idimp,npmax,mnblok,nxv,
     1nypmx,nxyzp,idds)
      do 60 m = 1, mnblok
      npp(m) = npp(m) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 80 m = 1, mnblok
      do 70 j = 1, nxyzp
      q(j,m) = q(j,m) + qp(j,m,i)
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPSOST32X(part,q,npp,nps,noff,nn,amxyz,qm,nx,idimp,npma
     1x,mnblok,nxv,nypmx,nxvyzp,idds,npd,n27,qp,idtask,nmt,ierr)
c parallel multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, amxyz, qm, qp
      integer npp, nps, noff, nn
      integer nx, idimp, npmax, mnblok, nxv, nypmx, nxvyzp, npd, n27
      integer idds
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), q(nxvyzp,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension nn(n27,npd,mnblok,nmt+1), amxyz(n27,npd,mnblok,nmt+1)   
      dimension qp(nxvyzp,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, m
      external PSOST32X
      data nargs /17/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start charge deposit tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear charge arrays
      do 30 m = 1, mnblok
      nps(m) = npt
      do 20 j = 1, nxvyzp
      qp(j,m,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),PSOST32X,nargs,part(1,npo,1),qp(1,1,i)
     1,nps,noff,nn(1,1,1,i+1),amxyz(1,1,1,i+1),qm,nx,idimp,npmax,mnblok,
     2nxv,nypmx,nxvyzp,idds,npd,n27)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining charge
      npo = npl + 1
      do 50 m = 1, mnblok
      npp(m) = npp(m) - npl
   50 continue
      call PSOST32X(part(1,npo,1),q,npp,noff,nn,amxyz,qm,nx,idimp,npmax,
     1mnblok,nxv,nypmx,nxvyzp,idds,npd,n27)
      do 60 m = 1, mnblok
      npp(m) = npp(m) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 80 m = 1, mnblok
      do 70 j = 1, nxvyzp
      q(j,m) = q(j,m) + qp(j,m,i)
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSOST32X(part,q,npp,nps,noff,nn,amxyz,qm,idimp,npmax,
     1mnblok,nxv,nypmx,nxvyzp,idds,npd,n27,qp,idtask,nmt,ierr)
c parallel multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, amxyz, qm, qp
      integer npp, nps, noff, nn
      integer idimp, npmax, mnblok, nxv, nypmx, nxvyzp, npd, n27, idds
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), q(nxvyzp,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension nn(n27,npd,mnblok,nmt+1), amxyz(n27,npd,mnblok,nmt+1)
      dimension qp(nxvyzp,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, m
      external PGSOST32X
      data nargs /16/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start charge deposit tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear charge arrays
      do 30 m = 1, mnblok
      nps(m) = npt
      do 20 j = 1, nxvyzp
      qp(j,m,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),PGSOST32X,nargs,part(1,npo,1),qp(1,1,i
     1),nps,noff,nn(1,1,1,i+1),amxyz(1,1,1,i+1),qm,idimp,npmax,mnblok,nx
     2v,nypmx,nxvyzp,idds,npd,n27)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining charge
      npo = npl + 1
      do 50 m = 1, mnblok
      npp(m) = npp(m) - npl
   50 continue
      call PGSOST32X(part(1,npo,1),q,npp,noff,nn,amxyz,qm,idimp,npmax,mn
     1blok,nxv,nypmx,nxvyzp,idds,npd,n27)
      do 60 m = 1, mnblok
      npp(m) = npp(m) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 80 m = 1, mnblok
      do 70 j = 1, nxvyzp
      q(j,m) = q(j,m) + qp(j,m,i)
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPDOST32L(part,q,npp,nps,noff,qm,nx,idimp,npmax,mnblok,
     1nxv,nypmx,nzpmx,idds,qp,idtask,nmt,ierr)
c parallel multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer npp, nps, noff
      integer nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), q(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension qp(nxv,nypmx,nzpmx,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PDOST32L
      data nargs /13/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start charge deposit tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear charge arrays
      do 50 m = 1, mnblok
      nps(m) = npt
      do 40 l = 1, nzpmx
      do 30 k = 1, nypmx
      do 20 j = 1, nx
      qp(j,k,l,m,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PDOST32L,nargs,part(1,npo,1),qp(1,1,1,
     11,i),nps,noff,qm,nx,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c deposit remaining charge
      npo = npl + 1
      do 70 m = 1, mnblok
      npp(m) = npp(m) - npl
   70 continue
      call PDOST32L(part(1,npo,1),q,npp,noff,qm,nx,idimp,npmax,mnblok,nx
     1v,nypmx,nzpmx,idds)
      do 80 m = 1, mnblok
      npp(m) = npp(m) + npl
   80 continue
c wait for tasks to complete
      do 130 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 120 m = 1, mnblok
      do 110 l = 1, nzpmx
      do 100 k = 1, nypmx
      do 90 j = 1, nx
      q(j,k,l,m) = q(j,k,l,m) + qp(j,k,l,m,i)
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGPOST32L(part,q,npp,nps,noff,qm,idimp,npmax,mnblok,nx
     1v,nypmx,nzpmx,idds,qp,idtask,nmt,ierr)
c parallel multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer npp, nps, noff
      integer idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), q(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension qp(nxv,nypmx,nzpmx,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m
      external PGPOST32L
      data nargs /12/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start charge deposit tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear charge arrays
      do 50 m = 1, mnblok
      nps(m) = npt
      do 40 l = 1, nzpmx
      do 30 k = 1, nypmx
      do 20 j = 1, nxv
      qp(j,k,l,m,i) = 0.
   20 continue
   30 continue
   40 continue
   50 continue
      call MP_TASKSTART(idtask(i),PGPOST32L,nargs,part(1,npo,1),qp(1,1,1
     1,1,i),nps,noff,qm,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c deposit remaining charge
      npo = npl + 1
      do 70 m = 1, mnblok
      npp(m) = npp(m) - npl
   70 continue
      call PGPOST32L(part(1,npo,1),q,npp,noff,qm,idimp,npmax,mnblok,nxv,
     1nypmx,nzpmx,idds)
      do 80 m = 1, mnblok
      npp(m) = npp(m) + npl
   80 continue
c wait for tasks to complete
      do 130 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 120 m = 1, mnblok
      do 110 l = 1, nzpmx
      do 100 k = 1, nypmx
      do 90 j = 1, nxv
      q(j,k,l,m) = q(j,k,l,m) + qp(j,k,l,m,i)
   90 continue
  100 continue
  110 continue
  120 continue
  130 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSPOST32L(part,q,npp,nps,noff,qm,idimp,npmax,mnblok,n
     1xv,nypmx,nxyzp,idds,qp,idtask,nmt,ierr)
c parallel multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, qm, qp
      integer npp, nps, noff
      integer idimp, npmax, mnblok, nxv, nypmx, nxyzp, idds
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), q(nxyzp,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension qp(nxyzp,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, m
      external PGSPOST32L
      data nargs /12/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start charge deposit tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear charge arrays
      do 30 m = 1, mnblok
      nps(m) = npt
      do 20 j = 1, nxyzp
      qp(j,m,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),PGSPOST32L,nargs,part(1,npo,1),qp(1,1,
     1i),nps,noff,qm,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining charge
      npo = npl + 1
      do 50 m = 1, mnblok
      npp(m) = npp(m) - npl
   50 continue
      call PGSPOST32L(part(1,npo,1),q,npp,noff,qm,idimp,npmax,mnblok,nxv
     1,nypmx,nxyzp,idds)
      do 60 m = 1, mnblok
      npp(m) = npp(m) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 80 m = 1, mnblok
      do 70 j = 1, nxyzp
      q(j,m) = q(j,m) + qp(j,m,i)
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPSOST32XL(part,q,npp,nps,noff,nn,amxyz,qm,nx,idimp,npm
     1ax,mnblok,nxv,nypmx,nxvyzp,idds,npd,ieight,qp,idtask,nmt,ierr)
c parallel multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, amxyz, qm, qp
      integer npp, nps, noff, nn
      integer nx, idimp, npmax, mnblok, nxv, nypmx, nxvyzp, npd, ieight
      integer idds
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), q(nxvyzp,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension nn(ieight,npd,mnblok,nmt+1)
      dimension amxyz(ieight,npd,mnblok,nmt+1)
      dimension qp(nxvyzp,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, m
      external PSOST32XL
      data nargs /17/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start charge deposit tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear charge arrays
      do 30 m = 1, mnblok
      nps(m) = npt
      do 20 j = 1, nxvyzp
      qp(j,m,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),PSOST32XL,nargs,part(1,npo,1),qp(1,1,i
     1),nps,noff,nn(1,1,1,i+1),amxyz(1,1,1,i+1),qm,nx,idimp,npmax,mnblok
     2,nxv,nypmx,nxvyzp,idds,npd,ieight)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining charge
      npo = npl + 1
      do 50 m = 1, mnblok
      npp(m) = npp(m) - npl
   50 continue
      call PSOST32XL(part(1,npo,1),q,npp,noff,nn,amxyz,qm,nx,idimp,npmax
     1,mnblok,nxv,nypmx,nxvyzp,idds,npd,ieight)
      do 60 m = 1, mnblok
      npp(m) = npp(m) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 80 m = 1, mnblok
      do 70 j = 1, nxvyzp
      q(j,m) = q(j,m) + qp(j,m,i)
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPGSOST32XL(part,q,npp,nps,noff,nn,amxyz,qm,idimp,npmax
     1,mnblok,nxv,nypmx,nxvyzp,idds,npd,ieight,qp,idtask,nmt,ierr)
c parallel multitasking charge deposition
c qp = charge density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, q, amxyz, qm, qp
      integer npp, nps, noff, nn
      integer idimp, npmax, mnblok, nxv, nypmx, nxvyzp, npd, ieight
      integer idds
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), q(nxvyzp,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension nn(ieight,npd,mnblok,nmt+1)
      dimension amxyz(ieight,npd,mnblok,nmt+1)
      dimension qp(nxvyzp,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, m
      external PGSOST32XL
      data nargs /16/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start charge deposit tasks
      do 40 i = 1, nmt
      npo = npt*(i - 1) + 1
c clear charge arrays
      do 30 m = 1, mnblok
      nps(m) = npt
      do 20 j = 1, nxvyzp
      qp(j,m,i) = 0.
   20 continue
   30 continue
      call MP_TASKSTART(idtask(i),PGSOST32XL,nargs,part(1,npo,1),qp(1,1,
     1i),nps,noff,nn(1,1,1,i+1),amxyz(1,1,1,i+1),qm,idimp,npmax,mnblok,n
     2xv,nypmx,nxvyzp,idds,npd,ieight)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   40 continue
c deposit remaining charge
      npo = npl + 1
      do 50 m = 1, mnblok
      npp(m) = npp(m) - npl
   50 continue
      call PGSOST32XL(part(1,npo,1),q,npp,noff,nn,amxyz,qm,idimp,npmax,m
     1nblok,nxv,nypmx,nxvyzp,idds,npd,ieight)
      do 60 m = 1, mnblok
      npp(m) = npp(m) + npl
   60 continue
c wait for tasks to complete
      do 90 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c sum charge arrays
      do 80 m = 1, mnblok
      do 70 j = 1, nxvyzp
      q(j,m) = q(j,m) + qp(j,m,i)
   70 continue
   80 continue
   90 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPPUSH32(part,fx,fy,fz,npp,nps,noff,qbm,dt,ek,nx,idimp,
     1npmax,mnblok,nxv,nypmx,nzpmx,idds,ekp,idtask,nmt,ierr)
c parallel multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, fy, fz, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok)
      dimension fx(nxv,nypmx,nzpmx,mnblok), fy(nxv,nypmx,nzpmx,mnblok)
      dimension fz(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, m
      external PPUSH32
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PPUSH32,nargs,part(1,npo,1),fx,fy,fz,n
     1ps,noff,qbm,dt,ekp(i),nx,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds)
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
      call PPUSH32(part(1,npo,1),fx,fy,fz,npp,noff,qbm,dt,ek,nx,idimp,np
     1max,mnblok,nxv,nypmx,nzpmx,idds)
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
      subroutine MPGPUSH32(part,fxyz,npp,nps,noff,qbm,dt,ek,nx,ny,nz,idi
     1mp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxyz, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok)
      dimension fxyz(3,nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, m
      external PGPUSH32
      data nargs /18/
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
      call MP_TASKSTART(idtask(i),PGPUSH32,nargs,part(1,npo,1),fxyz,nps,
     1noff,qbm,dt,ekp(i),nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idd
     2s,ipbc)
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
      call PGPUSH32(part(1,npo,1),fxyz,npp,noff,qbm,dt,ek,nx,ny,nz,idimp
     1,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
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
      subroutine MPGSPUSH32(part,fxyz,npp,nps,noff,qbm,dt,ek,nx,ny,nz,id
     1imp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxyz, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp, idds
      integer ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok)
      dimension fxyz(3,nxyzp,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, m
      external PGSPUSH32
      data nargs /18/
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
      call MP_TASKSTART(idtask(i),PGSPUSH32,nargs,part(1,npo,1),fxyz,nps
     1,noff,qbm,dt,ekp(i),nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,id
     2ds,ipbc)
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
      call PGSPUSH32(part(1,npo,1),fxyz,npp,noff,qbm,dt,ek,nx,ny,nz,idim
     1p,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
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
      subroutine MPPUSH32L(part,fx,fy,fz,npp,nps,noff,qbm,dt,ek,nx,idimp
     1,npmax,mnblok,nxv,nypmx,nzpmx,idds,ekp,idtask,nmt,ierr)
c parallel multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fx, fy, fz, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok)
      dimension fx(nxv,nypmx,nzpmx,mnblok), fy(nxv,nypmx,nzpmx,mnblok)
      dimension fz(nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, m
      external PPUSH32L
      data nargs /17/
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
      call MP_TASKSTART(idtask(i),PPUSH32L,nargs,part(1,npo,1),fx,fy,fz,
     1nps,noff,qbm,dt,ekp(i),nx,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds)
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
      call PPUSH32L(part(1,npo,1),fx,fy,fz,npp,noff,qbm,dt,ek,nx,idimp,n
     1pmax,mnblok,nxv,nypmx,nzpmx,idds)
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
      subroutine MPGPUSH32L(part,fxyz,npp,nps,noff,qbm,dt,ek,nx,ny,nz,id
     1imp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxyz, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok)
      dimension fxyz(3,nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, m
      external PGPUSH32L
      data nargs /18/
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
      call MP_TASKSTART(idtask(i),PGPUSH32L,nargs,part(1,npo,1),fxyz,nps
     1,noff,qbm,dt,ekp(i),nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,id
     2ds,ipbc)
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
      call PGPUSH32L(part(1,npo,1),fxyz,npp,noff,qbm,dt,ek,nx,ny,nz,idim
     1p,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
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
      subroutine MPGSPUSH32L(part,fxyz,npp,nps,noff,qbm,dt,ek,nx,ny,nz,i
     1dimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxyz, qbm, dt, ek, ekp
      integer npp, nps, noff
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp, idds
      integer ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok)
      dimension fxyz(3,nxyzp,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension ekp(nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, m
      external PGSPUSH32L
      data nargs /18/
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
      call MP_TASKSTART(idtask(i),PGSPUSH32L,nargs,part(1,npo,1),fxyz,np
     1s,noff,qbm,dt,ekp(i),nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,i
     2dds,ipbc)
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
      call PGSPUSH32L(part(1,npo,1),fxyz,npp,noff,qbm,dt,ek,nx,ny,nz,idi
     1mp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
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
      subroutine MPSORTP32YZ(part,pt,ip,npic,npp,nps,noff,nyzp,idimp,npm
     1ax,mnblok,nyzpm1,idds,npicp,idtask,nmt,ierr)
c multitasking particle sorting
c npicp = address offset arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, pt
      integer ip, npic, npp, nps, noff, nyzp
      integer idimp, npmax, mnblok, nyzpm1, idds
      integer npicp, idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), pt(npmax,mnblok)
      dimension ip(npmax,mnblok), npic(nyzpm1,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(mnblok)
      dimension nyzp(idds,mnblok)
      dimension npicp(nyzpm1,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, m
      external PSORTP32YZ
      data nargs /12/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle sorting tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 m = 1, mnblok
      nps(m) = npt
   20 continue
      call MP_TASKSTART(idtask(i),PSORTP32YZ,nargs,part(1,npo,1),pt(npo,
     11),ip(npo,1),npicp(1,1,i),nps,noff,nyzp,idimp,npmax,mnblok,nyzpm1,
     2idds)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c sort remaining particles
      npo = npl + 1
      do 40 m = 1, mnblok
      npp(m) = npp(m) - npl
   40 continue
      call PSORTP32YZ(part(1,npo,1),pt(npo,1),ip(npo,1),npic,npp,noff,ny
     1zp,idimp,npmax,mnblok,nyzpm1,idds)
      do 50 m = 1, mnblok
      npp(m) = npp(m) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPSORTP32YZL(part,pt,ip,npic,npp,nps,noff,nyzp,idimp,np
     1max,mnblok,nyzpm1,idds,npicp,idtask,nmt,ierr)
c multitasking particle sorting
c npicp = address offset arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, pt
      integer ip, npic, npp, nps, noff, nyzp
      integer idimp, npmax, mnblok, nyzpm1, idds
      integer npicp, idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), pt(npmax,mnblok)
      dimension ip(npmax,mnblok), npic(nyzpm1,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(mnblok)
      dimension nyzp(idds,mnblok)
      dimension npicp(nyzpm1,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, m
      external PSORTP32YZL
      data nargs /12/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle sorting tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 m = 1, mnblok
      nps(m) = npt
   20 continue
      call MP_TASKSTART(idtask(i),PSORTP32YZL,nargs,part(1,npo,1),pt(npo
     1,1),ip(npo,1),npicp(1,1,i),nps,noff,nyzp,idimp,npmax,mnblok,nyzpm1
     2,idds)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c sort remaining particles
      npo = npl + 1
      do 40 m = 1, mnblok
      npp(m) = npp(m) - npl
   40 continue
      call PSORTP32YZL(part(1,npo,1),pt(npo,1),ip(npo,1),npic,npp,noff,n
     1yzp,idimp,npmax,mnblok,nyzpm1,idds)
      do 50 m = 1, mnblok
      npp(m) = npp(m) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPDSORTP32YZ(parta,partb,npic,npp,nps,noff,nyzp,idimp,n
     1pmax,mnblok,nyzpm1,idds,npicp,idtask,nmt,ierr)
c multitasking particle sorting
c npicp = address offset arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real parta, partb
      integer npic, npp, nps, noff, nyzp
      integer idimp, npmax, mnblok, nyzpm1, idds
      integer npicp, idtask, nmt, ierr
      dimension parta(idimp,npmax,mnblok), partb(idimp,npmax,mnblok)
      dimension npic(nyzpm1,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(mnblok)
      dimension nyzp(idds,mnblok)
      dimension npicp(nyzpm1,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, m
      external PDSORTP32YZ
      data nargs /11/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle sorting tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 m = 1, mnblok
      nps(m) = npt
   20 continue
      call MP_TASKSTART(idtask(i),PDSORTP32YZ,nargs,parta(1,npo,1),partb
     1(1,npo,1),npicp(1,1,i),nps,noff,nyzp,idimp,npmax,mnblok,nyzpm1,idd
     2s)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c sort remaining particles
      npo = npl + 1
      do 40 m = 1, mnblok
      npp(m) = npp(m) - npl
   40 continue
      call PDSORTP32YZ(parta(1,npo,1),partb(1,npo,1),npic,npp,noff,nyzp,
     1idimp,npmax,mnblok,nyzpm1,idds)
      do 50 m = 1, mnblok
      npp(m) = npp(m) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   60 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPDSORTP32YZL(parta,partb,npic,npp,nps,noff,nyzp,idimp,
     1npmax,mnblok,nyzpm1,idds,npicp,idtask,nmt,ierr)
c multitasking particle sorting
c npicp = address offset arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real parta, partb
      integer npic, npp, nps, noff, nyzp
      integer idimp, npmax, mnblok, nyzpm1, idds
      integer npicp, idtask, nmt, ierr
      dimension parta(idimp,npmax,mnblok), partb(idimp,npmax,mnblok)
      dimension npic(nyzpm1,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(mnblok)
      dimension nyzp(idds,mnblok)
      dimension npicp(nyzpm1,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, m
      external PDSORTP32YZL
      data nargs /11/
c find minimum of npp
      npo = npp(1)
      do 10 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   10 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
      ierr = 0
c start particle sorting tasks
      do 30 i = 1, nmt
      npo = npt*(i - 1) + 1
      do 20 m = 1, mnblok
      nps(m) = npt
   20 continue
      call MP_TASKSTART(idtask(i),PDSORTP32YZL,nargs,parta(1,npo,1),part
     1b(1,npo,1),npicp(1,1,i),nps,noff,nyzp,idimp,npmax,mnblok,nyzpm1,id
     2ds)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   30 continue
c sort remaining particles
      npo = npl + 1
      do 40 m = 1, mnblok
      npp(m) = npp(m) - npl
   40 continue
      call PDSORTP32YZL(parta(1,npo,1),partb(1,npo,1),npic,npp,noff,nyzp
     1,idimp,npmax,mnblok,nyzpm1,idds)
      do 50 m = 1, mnblok
      npp(m) = npp(m) + npl
   50 continue
c wait for tasks to complete
      do 60 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
   60 continue
      return
      end

