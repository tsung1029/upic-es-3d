c 3d PIC parallel multi-tasking library for pushing relativistic
c particles with magnetic field and depositing current
c with 2D domain decomposition
c written by viktor k. decyk, ucla
c copyright 2000, regents of the university of california
c update: december 20, 2003
c-----------------------------------------------------------------------
      subroutine MPGRJPOST32(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,nz,idim
     1p,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer npp, nps, noff
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), cu(3,nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension cup(3,nxv,nypmx,nzpmx,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m, n
      external PGRJPOST32
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
      call MP_TASKSTART(idtask(i),PGRJPOST32,nargs,part(1,npo,1),cup(1,1
     1,1,1,1,i),nps,noff,qm,dt,ci,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,
     2nzpmx,idds,ipbc)
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
      call PGRJPOST32(part(1,npo,1),cu,npp,noff,qm,dt,ci,nx,ny,nz,idimp,
     1npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
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
      subroutine MPGSRJPOST32(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,nz,idi
     1mp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer npp, nps, noff
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp, idds
      integer ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), cu(3,nxyzp,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension cup(3,nxyzp,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, m, n
      external PGSRJPOST32
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
      call MP_TASKSTART(idtask(i),PGSRJPOST32,nargs,part(1,npo,1),cup(1,
     11,1,i),nps,noff,qm,dt,ci,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxy
     2zp,idds,ipbc)
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
      call PGSRJPOST32(part(1,npo,1),cu,npp,noff,qm,dt,ci,nx,ny,nz,idimp
     1,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
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
      subroutine MPGSRJOST32X(part,cu,npp,nps,noff,nn,amxyz,qm,dt,ci,nx,
     1ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n81,ipbc,cup,id
     2task,nmt,ierr)
c parallel multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, amxyz, qm, dt, ci, cup
      integer npp, nps, noff, nn
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxvyzp, ipbc
      integer idds, npd, n81, idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), cu(3*nxvyzp,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension nn(n81,npd,mnblok,nmt+1), amxyz(n81,npd,mnblok,nmt+1)
      dimension cup(3*nxvyzp,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, m
      external PGSRJOST32X
      data nargs /22/
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
      call MP_TASKSTART(idtask(i),PGSRJOST32X,nargs,part(1,npo,1),cup(1,
     11,i),nps,noff,nn(1,1,1,i+1),amxyz(1,1,1,i+1),qm,dt,ci,nx,ny,nz,idi
     2mp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n81,ipbc)
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
      call PGSRJOST32X(part(1,npo,1),cu,npp,noff,nn,amxyz,qm,dt,ci,nx,ny
     1,nz,idimp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n81,ipbc)
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
      subroutine MPGRJPOST32L(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,nz,idi
     1mp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer npp, nps, noff
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nzpmx, idds
      integer ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), cu(3,nxv,nypmx,nzpmx,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension cup(3,nxv,nypmx,nzpmx,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, k, l, m, n
      external PGRJPOST32L
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
      call MP_TASKSTART(idtask(i),PGRJPOST32L,nargs,part(1,npo,1),cup(1,
     11,1,1,1,i),nps,noff,qm,dt,ci,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx
     2,nzpmx,idds,ipbc)
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
      call PGRJPOST32L(part(1,npo,1),cu,npp,noff,qm,dt,ci,nx,ny,nz,idimp
     1,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
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
      subroutine MPGSRJPOST32L(part,cu,npp,nps,noff,qm,dt,ci,nx,ny,nz,id
     1imp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,cup,idtask,nmt,ierr)
c parallel multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, qm, dt, ci, cup
      integer npp, nps, noff
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxyzp, idds
      integer ipbc
      integer idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), cu(3,nxyzp,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension cup(3,nxyzp,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, m, n
      external PGSRJPOST32L
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
      call MP_TASKSTART(idtask(i),PGSRJPOST32L,nargs,part(1,npo,1),cup(1
     1,1,1,i),nps,noff,qm,dt,ci,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nx
     2yzp,idds,ipbc)
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
      call PGSRJPOST32L(part(1,npo,1),cu,npp,noff,qm,dt,ci,nx,ny,nz,idim
     1p,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
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
      subroutine MPGSRJOST32XL(part,cu,npp,nps,noff,nn,amxyz,qm,dt,ci,nx
     1,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n24,ipbc,cup,i
     2dtask,nmt,ierr)
c parallel multitasking current deposition for relativistic particles
c cup = current density arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, cu, amxyz, qm, dt, ci, cup
      integer npp, nps, noff, nn
      integer nx, ny, nz, idimp, npmax, mnblok, nxv, nypmx, nxvyzp, ipbc
      integer idds, npd, n24, idtask, nmt, ierr
      dimension part(idimp,npmax,mnblok), cu(3*nxvyzp,mnblok)
      dimension npp(mnblok), nps(mnblok), noff(idds,mnblok)
      dimension nn(n24,npd,mnblok,nmt+1), amxyz(n24,npd,mnblok,nmt+1)
      dimension cup(3*nxvyzp,mnblok,nmt), idtask(nmt)
c local data
      integer nargs, npt, npl, npo, i, j, m
      external PGSRJOST32XL
      data nargs /22/
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
      call MP_TASKSTART(idtask(i),PGSRJOST32XL,nargs,part(1,npo,1),cup(1
     1,1,i),nps,noff,nn(1,1,1,i+1),amxyz(1,1,1,i+1),qm,dt,ci,nx,ny,nz,id
     2imp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n24,ipbc)
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
      call PGSRJOST32XL(part(1,npo,1),cu,npp,noff,nn,amxyz,qm,dt,ci,nx,n
     1y,nz,idimp,npmax,mnblok,nxv,nypmx,nxvyzp,idds,npd,n24,ipbc)
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
      subroutine MPGRPUSH32(part,fxyz,npp,nps,noff,qbm,dt,ci,ek,nx,ny,nz
     1,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask,nmt,ierr)
c parallel multitasking relativistic particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxyz, qbm, dt, ci, ek, ekp
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
      external PGRPUSH32
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
      call MP_TASKSTART(idtask(i),PGRPUSH32,nargs,part(1,npo,1),fxyz,nps
     1,noff,qbm,dt,ci,ekp(i),nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx
     2,idds,ipbc)
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
      call PGRPUSH32(part(1,npo,1),fxyz,npp,noff,qbm,dt,ci,ek,nx,ny,nz,i
     1dimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
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
      subroutine MPGSRPUSH32(part,fxyz,npp,nps,noff,qbm,dt,ci,ek,nx,ny,n
     1z,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,idtask,nmt,ierr
     2)
c parallel multitasking relativistic particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxyz, qbm, dt, ci, ek, ekp
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
      external PGSRPUSH32
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
      call MP_TASKSTART(idtask(i),PGSRPUSH32,nargs,part(1,npo,1),fxyz,np
     1s,noff,qbm,dt,ci,ekp(i),nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyz
     2p,idds,ipbc)
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
      call PGSRPUSH32(part(1,npo,1),fxyz,npp,noff,qbm,dt,ci,ek,nx,ny,nz,
     1idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
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
      subroutine MPGRPUSH32L(part,fxyz,npp,nps,noff,qbm,dt,ci,ek,nx,ny,n
     1z,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask,nmt,ierr
     2)
c parallel multitasking relativistic particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxyz, qbm, dt, ci, ek, ekp
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
      external PGRPUSH32L
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
      call MP_TASKSTART(idtask(i),PGRPUSH32L,nargs,part(1,npo,1),fxyz,np
     1s,noff,qbm,dt,ci,ekp(i),nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpm
     2x,idds,ipbc)
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
      call PGRPUSH32L(part(1,npo,1),fxyz,npp,noff,qbm,dt,ci,ek,nx,ny,nz,
     1idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
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
      subroutine MPGSRPUSH32L(part,fxyz,npp,nps,noff,qbm,dt,ci,ek,nx,ny,
     1nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,idtask,nmt,ier
     2r)
c parallel multitasking relativistic particle push
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxyz, qbm, dt, ci, ek, ekp
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
      external PGSRPUSH32L
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
      call MP_TASKSTART(idtask(i),PGSRPUSH32L,nargs,part(1,npo,1),fxyz,n
     1ps,noff,qbm,dt,ci,ekp(i),nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxy
     2zp,idds,ipbc)
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
      call PGSRPUSH32L(part(1,npo,1),fxyz,npp,noff,qbm,dt,ci,ek,nx,ny,nz
     1,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
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
      subroutine MPGRBPUSH32(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,ci,e
     1k,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtask
     2,nmt,ierr)
c parallel multitasking relativistic particle push with magnetic field
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxyz, bxyz, qbm, dt, dtc, ci, ek, ekp
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
      external PGRBPUSH32
      data nargs /21/
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
      call MP_TASKSTART(idtask(i),PGRBPUSH32,nargs,part(1,npo,1),fxyz,bx
     1yz,nps,noff,qbm,dt,dtc,ci,ekp(i),nx,ny,nz,idimp,npmax,mnblok,nxv,n
     2ypmx,nzpmx,idds,ipbc)
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
      call PGRBPUSH32(part(1,npo,1),fxyz,bxyz,npp,noff,qbm,dt,dtc,ci,ek,
     1nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
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
      subroutine MPGSRBPUSH32(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,ci,
     1ek,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,idtas
     2k,nmt,ierr)
c parallel multitasking relativistic particle push with magnetic field
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxyz, bxyz, qbm, dt, dtc, ci, ek, ekp
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
      external PGSRBPUSH32
      data nargs /21/
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
      call MP_TASKSTART(idtask(i),PGSRBPUSH32,nargs,part(1,npo,1),fxyz,b
     1xyz,nps,noff,qbm,dt,dtc,ci,ekp(i),nx,ny,nz,idimp,npmax,mnblok,nxv,
     2nypmx,nxyzp,idds,ipbc)
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
      call PGSRBPUSH32(part(1,npo,1),fxyz,bxyz,npp,noff,qbm,dt,dtc,ci,ek
     1,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
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
      subroutine MPGRBPUSH32L(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,ci,
     1ek,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc,ekp,idtas
     2k,nmt,ierr)
c parallel multitasking relativistic particle push with magnetic field
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxyz, bxyz, qbm, dt, dtc, ci, ek, ekp
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
      external PGRBPUSH32L
      data nargs /21/
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
      call MP_TASKSTART(idtask(i),PGRBPUSH32L,nargs,part(1,npo,1),fxyz,b
     1xyz,nps,noff,qbm,dt,dtc,ci,ekp(i),nx,ny,nz,idimp,npmax,mnblok,nxv,
     2nypmx,nzpmx,idds,ipbc)
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
      call PGRBPUSH32L(part(1,npo,1),fxyz,bxyz,npp,noff,qbm,dt,dtc,ci,ek
     1,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nzpmx,idds,ipbc)
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
      subroutine MPGSRBPUSH32L(part,fxyz,bxyz,npp,nps,noff,qbm,dt,dtc,ci
     1,ek,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc,ekp,idta
     2sk,nmt,ierr)
c parallel multitasking relativistic particle push with magnetic field
c ekp = kinetic energy arrays for tasks
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, fxyz, bxyz, qbm, dt, dtc, ci, ek, ekp
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
      external PGSRBPUSH32L
      data nargs /21/
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
      call MP_TASKSTART(idtask(i),PGSRBPUSH32L,nargs,part(1,npo,1),fxyz,
     1bxyz,nps,noff,qbm,dt,dtc,ci,ekp(i),nx,ny,nz,idimp,npmax,mnblok,nxv
     2,nypmx,nxyzp,idds,ipbc)
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
      call PGSRBPUSH32L(part(1,npo,1),fxyz,bxyz,npp,noff,qbm,dt,dtc,ci,e
     1k,nx,ny,nz,idimp,npmax,mnblok,nxv,nypmx,nxyzp,idds,ipbc)
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
