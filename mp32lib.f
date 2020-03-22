c 3d parallel PIC multi-tasking library for MPI communications
c with 2D domain decomposition
c written by viktor k. decyk, ucla
c copyright 1995, regents of the university of california
c update: march 4, 2008
c-----------------------------------------------------------------------
      subroutine MPMOVE32(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,iho
     1le,jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,id
     2ps,nbmax,idds,ntmax,info,jssp,idtask,nmt,ierr)
c multi-tasking particle manager
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl, th
      integer npp, npq, ihole, jsr, jsl, jss, info
      integer ny, nz, kstrt, nvpy, nvpz, idimp, npmax, mblok, nblok
      integer idps, nbmax, idds, ntmax
      integer jssp, idtask, nmt, ierr
      dimension part(idimp,npmax,mblok*nblok)
      dimension edges(idps,mblok*nblok), npp(mblok*nblok)
      dimension npq(mblok*nblok)
      dimension sbufl(idimp,nbmax,mblok*nblok)
      dimension sbufr(idimp,nbmax,mblok*nblok)
      dimension rbufl(idimp,nbmax,mblok*nblok)
      dimension rbufr(idimp,nbmax,mblok*nblok)
      dimension jsl(idds,mblok*nblok), jsr(idds,mblok*nblok)
      dimension jss(idds,mblok*nblok)
      dimension ihole(ntmax,mblok*nblok)
      dimension info(9)
      dimension jssp(idds,mblok*nblok,nmt), idtask(nmt)
c local data
      integer mnblok, i, j, m, n, npr, nps, nter, nargs, npt, npl, npo
      integer mpo, mpt, mpl, mps
      real tf
      double precision dtime
      double precision bflg, work
      dimension bflg(2), work(2)
      external PMOVEH32
      data nargs /13/
      mnblok = mblok*nblok
      do 10 j = 1, 9
      info(j) = 0
   10 continue
c debugging section: count total number of particles before move
      npr = 0
      do 20 m = 1, mnblok
      npr = npr + npp(m)
   20 continue
c find minimum of npp
      npo = npp(1)
      do 30 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   30 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
c find ihole partitions
      mpt = ntmax/(nmt + 1)
      mpl = mpt*nmt
      ierr = 0
c find outgoing particles, first in y then in z direction
      do 140 n = 1, 2
   40 nter = info(n+3)
      call PWTIMERA(-1,tf,dtime)
c start tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
      mpo = mpt*(i - 1) + 1
      do 50 m = 1, nblok
      npq(m) = npt
   50 continue
      call MP_TASKSTART(idtask(i),PMOVEH32,nargs,part(1,npo,1),edges,npq
     1,ihole(mpo,1),jssp(1,1,i),idimp,npt,mblok,nblok,idps,idds,mpt,n)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c finish finding outgoing particles
      npo = npl + 1
      do 70 m = 1, nblok
      npp(m) = npp(m) - npl
   70 continue
      npl = npmax - npl
      mpo = mpl + 1
      mps = ntmax - mpl
      call PMOVEH32(part(1,npo,1),edges,npp,ihole(mpo,1),jss,idimp,npl,m
     1blok,nblok,idps,idds,mps,n)
      npl = npmax - npl
      do 80 m = 1, nblok
      npp(m) = npp(m) + npl
   80 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c compress ihole array
      if (i.eq.1) go to 110
      npo = npt*(i - 1)
      mpo = mpt*(i - 1)
      do 100 m = 1, mnblok
      do 90 j = 1, jssp(1,m,i)
      ihole(jssp(1,m,i-1)+j,m) = ihole(j+mpo,m) + npo
   90 continue
      jssp(1,m,i) = jssp(1,m,i) + jssp(1,m,i-1)
      jssp(2,m,i) = max0(jssp(2,m,i),jssp(2,m,i-1))
  100 continue
  110 continue
      if (nmt.gt.0) then
         do 130 m = 1, mnblok
         do 120 j = 1, jss(1,m)
         ihole(jssp(1,m,nmt)+j,m) = ihole(j+mpl,m) + npl
  120    continue
         jss(1,m) = jss(1,m) + jssp(1,m,nmt)
         jss(2,m) = max0(jss(2,m),jssp(2,m,nmt))
  130    continue
      endif
      call PWTIMERA(1,tf,dtime)
      th = th + tf
c send outgoing particles
      call PMOVES32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,jsl
     1,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbmax,idds
     2,ntmax,info,n)
c particle overflow error
      if (info(1).gt.0) return
c buffer overflowed and more particles remain to be checked
      if (info(n+3).gt.nter) go to 40
c iteration overflow
      if (info(1).lt.0) go to 160
  140 continue
c debugging section: count total number of particles after move
      nps = 0
      do 150 m = 1, mnblok
      nps = nps + npp(m)
  150 continue
      bflg(2) = dble(nps)
      bflg(1) = dble(npr)
      call PDSUM(bflg,work,2,1)
      info(8) = bflg(1)
      info(9) = bflg(2) - bflg(1)
      if (bflg(1).ne.bflg(2)) then
         if (kstrt.eq.1) then
           write (2,*) 'particle number error, old/new=',bflg(1),bflg(2)
         endif
         info(1) = 1
      endif
      return
c information
  160 nter = info(n+3)
      if (nter.gt.0) then
         write (2,*) 'Info: ', nter, 'Abuffer overflows, n, nbmax=', n, 
     1nbmax
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MPXMOV32(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,iho
     1le,jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,id
     2ps,nbmax,idds,ntmax,maskp,info,jssp,idtask,nmt,ierr)
c multi-tasking particle manager
c optimized for vector processor
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl, th
      integer npp, npq, ihole, jsr, jsl, jss, maskp, info
      integer ny, nz, kstrt, nvpy, nvpz, idimp, npmax, mblok, nblok
      integer idps, nbmax, idds, ntmax
      integer jssp, idtask, nmt, ierr
      dimension part(idimp,npmax,mblok*nblok), maskp(npmax,mblok*nblok)
      dimension edges(idps,mblok*nblok), npp(mblok*nblok)
      dimension npq(mblok*nblok)
      dimension sbufl(idimp,nbmax,mblok*nblok)
      dimension sbufr(idimp,nbmax,mblok*nblok)
      dimension rbufl(idimp,nbmax,mblok*nblok)
      dimension rbufr(idimp,nbmax,mblok*nblok)
      dimension jsl(idds,mblok*nblok), jsr(idds,mblok*nblok)
      dimension jss(idds,mblok*nblok)
      dimension ihole(ntmax,mblok*nblok)
      dimension info(9)
      dimension jssp(idds,mblok*nblok,nmt), idtask(nmt)
c local data
      integer mnblok, i, j, m, n, npr, nps, nter, nargs, npt, npl, npo
      integer mpo, mpt, mpl, mps
      real tf
      double precision dtime
      double precision bflg, work
      dimension bflg(2), work(2)
      external PMOVEHX32
      data nargs /14/
      mnblok = mblok*nblok
      do 10 j = 1, 9
      info(j) = 0
   10 continue
c debugging section: count total number of particles before move
      npr = 0
      do 20 m = 1, mnblok
      npr = npr + npp(m)
   20 continue
c find minimum of npp
      npo = npp(1)
      do 30 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   30 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
c find ihole partitions
      mpt = ntmax/(nmt + 1)
      mpl = mpt*nmt
      ierr = 0
c find outgoing particles, first in y then in z direction
      do 140 n = 1, 2
   40 nter = info(n+3)
      call PWTIMERA(-1,tf,dtime)
c start tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
      mpo = mpt*(i - 1) + 1
      do 50 m = 1, nblok
      npq(m) = npt
   50 continue
      call MP_TASKSTART(idtask(i),PMOVEHX32,nargs,part(1,npo,1),edges,np
     1q,ihole(mpo,1),jssp(1,1,i),idimp,npt,mblok,nblok,idps,idds,mpt,mas
     2kp(npo,1),n)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c finish finding outgoing particles
      npo = npl + 1
      do 70 m = 1, nblok
      npp(m) = npp(m) - npl
   70 continue
      npl = npmax - npl
      mpo = mpl + 1
      mps = ntmax - mpl
      call PMOVEHX32(part(1,npo,1),edges,npp,ihole(mpo,1),jss,idimp,npl,
     1mblok,nblok,idps,idds,mps,maskp(npo,1),n)
      npl = npmax - npl
      do 80 m = 1, nblok
      npp(m) = npp(m) + npl
   80 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c compress ihole array
      if (i.eq.1) go to 110
      npo = npt*(i - 1)
      mpo = mpt*(i - 1)
      do 100 m = 1, mnblok
      do 90 j = 1, jssp(1,m,i)
      ihole(jssp(1,m,i-1)+j,m) = ihole(j+mpo,m) + npo
   90 continue
      jssp(1,m,i) = jssp(1,m,i) + jssp(1,m,i-1)
      jssp(2,m,i) = max0(jssp(2,m,i),jssp(2,m,i-1))
  100 continue
  110 continue
      if (nmt.gt.0) then
         do 130 m = 1, mnblok
         do 120 j = 1, jss(1,m)
         ihole(jssp(1,m,nmt)+j,m) = ihole(j+mpl,m) + npl
  120    continue
         jss(1,m) = jss(1,m) + jssp(1,m,nmt)
         jss(2,m) = max0(jss(2,m),jssp(2,m,nmt))
  130    continue
      endif
      call PWTIMERA(1,tf,dtime)
      th = th + tf
c send outgoing particles
      call PMOVES32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,jsl
     1,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbmax,idds
     2,ntmax,info,n)
c particle overflow error
      if (info(1).gt.0) return
c buffer overflowed and more particles remain to be checked
      if (info(n+3).gt.nter) go to 40
c iteration overflow
      if (info(1).lt.0) go to 160
  140 continue
c debugging section: count total number of particles after move
      nps = 0
      do 150 m = 1, mnblok
      nps = nps + npp(m)
  150 continue
      bflg(2) = dble(nps)
      bflg(1) = dble(npr)
      call PDSUM(bflg,work,2,1)
      info(8) = bflg(1)
      info(9) = bflg(2) - bflg(1)
      if (bflg(1).ne.bflg(2)) then
         if (kstrt.eq.1) then
           write (2,*) 'particle number error, old/new=',bflg(1),bflg(2)
         endif
         info(1) = 1
      endif
      return
c information
  160 nter = info(n+3)
      if (nter.gt.0) then
         write (2,*) 'Info: ', nter, 'Bbuffer overflows, n, nbmax=', n, 
     1nbmax
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MPMOVES32(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,ih
     1ole,jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,i
     2dps,nbmax,idds,ntmax,info,jssp,idtask,nmt,ierr)
c multi-tasking particle manager
c info(6) = maximum number of particle passes required in y
c info(6) must be set on entry
c info(7) = maximum number of particle passes required in z
c info(7) must be set on entry
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl, th
      integer npp, npq, ihole, jsr, jsl, jss, info
      integer ny, nz, kstrt, nvpy, nvpz, idimp, npmax, mblok, nblok
      integer idps, nbmax, idds, ntmax
      integer jssp, idtask, nmt, ierr
      dimension part(idimp,npmax,mblok*nblok)
      dimension edges(idps,mblok*nblok), npp(mblok*nblok)
      dimension npq(mblok*nblok)
      dimension sbufl(idimp,nbmax,mblok*nblok)
      dimension sbufr(idimp,nbmax,mblok*nblok)
      dimension rbufl(idimp,nbmax,mblok*nblok)
      dimension rbufr(idimp,nbmax,mblok*nblok)
      dimension jsl(idds,mblok*nblok), jsr(idds,mblok*nblok)
      dimension jss(idds,mblok*nblok)
      dimension ihole(ntmax,mblok*nblok)
      dimension info(9)
      dimension jssp(idds,mblok*nblok,nmt), idtask(nmt)
c local data
      integer mnblok, i, j, m, n, nargs, npt, npl, npo
      integer mpo, mpt, mpl, mps
      real tf
      double precision dtime
      external PMOVEH32
      data nargs /13/
      mnblok = mblok*nblok
      do 10 j = 1, 5
      info(j) = 0
   10 continue
c find minimum of npp
      npo = npp(1)
      do 30 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   30 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
c find ihole partitions
      mpt = ntmax/(nmt + 1)
      mpl = mpt*nmt
      ierr = 0
c find outgoing particles, first in y then in z direction
      do 140 n = 1, 2
      call PWTIMERA(-1,tf,dtime)
c start tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
      mpo = mpt*(i - 1) + 1
      do 50 m = 1, nblok
      npq(m) = npt
   50 continue
      call MP_TASKSTART(idtask(i),PMOVEH32,nargs,part(1,npo,1),edges,npq
     1,ihole(mpo,1),jssp(1,1,i),idimp,npt,mblok,nblok,idps,idds,mpt,n)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c finish finding outgoing particles
      npo = npl + 1
      do 70 m = 1, nblok
      npp(m) = npp(m) - npl
   70 continue
      npl = npmax - npl
      mpo = mpl + 1
      mps = ntmax - mpl
      call PMOVEH32(part(1,npo,1),edges,npp,ihole(mpo,1),jss,idimp,npl,m
     1blok,nblok,idps,idds,mps,n)
      npl = npmax - npl
      do 80 m = 1, nblok
      npp(m) = npp(m) + npl
   80 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c compress ihole array
      if (i.eq.1) go to 110
      npo = npt*(i - 1)
      mpo = mpt*(i - 1)
      do 100 m = 1, mnblok
      do 90 j = 1, jssp(1,m,i)
      ihole(jssp(1,m,i-1)+j,m) = ihole(j+mpo,m) + npo
   90 continue
      jssp(1,m,i) = jssp(1,m,i) + jssp(1,m,i-1)
      jssp(2,m,i) = max0(jssp(2,m,i),jssp(2,m,i-1))
  100 continue
  110 continue
      if (nmt.gt.0) then
         do 130 m = 1, mnblok
         do 120 j = 1, jss(1,m)
         ihole(jssp(1,m,nmt)+j,m) = ihole(j+mpl,m) + npl
  120    continue
         jss(1,m) = jss(1,m) + jssp(1,m,nmt)
         jss(2,m) = max0(jss(2,m),jssp(2,m,nmt))
  130    continue
      endif
      call PWTIMERA(1,tf,dtime)
      th = th + tf
c send outgoing particles
      call PMOVESS32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,js
     1l,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbmax,idd
     2s,ntmax,info,n)
c particle or iteration overflow error
      if (info(1).ne.0) return
  140 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPXMOVS32(part,edges,npp,npq,sbufr,sbufl,rbufr,rbufl,ih
     1ole,jsr,jsl,jss,th,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,i
     2dps,nbmax,idds,ntmax,maskp,info,jssp,idtask,nmt,ierr)
c multi-tasking particle manager
c info(6) = maximum number of particle passes required in y
c info(6) must be set on entry
c info(7) = maximum number of particle passes required in z
c info(7) must be set on entry
c optimized for vector processor
c idtask = index to notify queue for task (0 if error)
c nmt = number of tasks
c ierr = ierror indicator (0 = no error)
      implicit none
      real part, edges, sbufr, sbufl, rbufr, rbufl, th
      integer npp, npq, ihole, jsr, jsl, jss, maskp, info
      integer ny, nz, kstrt, nvpy, nvpz, idimp, npmax, mblok, nblok
      integer idps, nbmax, idds, ntmax
      integer jssp, idtask, nmt, ierr
      dimension part(idimp,npmax,mblok*nblok), maskp(npmax,mblok*nblok)
      dimension edges(idps,mblok*nblok), npp(mblok*nblok)
      dimension npq(mblok*nblok)
      dimension sbufl(idimp,nbmax,mblok*nblok)
      dimension sbufr(idimp,nbmax,mblok*nblok)
      dimension rbufl(idimp,nbmax,mblok*nblok)
      dimension rbufr(idimp,nbmax,mblok*nblok)
      dimension jsl(idds,mblok*nblok), jsr(idds,mblok*nblok)
      dimension jss(idds,mblok*nblok)
      dimension ihole(ntmax,mblok*nblok)
      dimension info(9)
      dimension jssp(idds,mblok*nblok,nmt), idtask(nmt)
c local data
      integer mnblok, i, j, m, n, nargs, npt, npl, npo
      integer mpo, mpt, mpl, mps
      real tf
      double precision dtime
      external PMOVEHX32
      data nargs /14/
      mnblok = mblok*nblok
      do 10 j = 1, 5
      info(j) = 0
   10 continue
c find minimum of npp
      npo = npp(1)
      do 30 m = 1, mnblok
      if (npo.gt.npp(m)) npo = npp(m)
   30 continue
      npt = npo/(nmt + 1)
      npl = npt*nmt
c find ihole partitions
      mpt = ntmax/(nmt + 1)
      mpl = mpt*nmt
      ierr = 0
c find outgoing particles, first in y then in z direction
      do 140 n = 1, 2
      call PWTIMERA(-1,tf,dtime)
c start tasks
      do 60 i = 1, nmt
      npo = npt*(i - 1) + 1
      mpo = mpt*(i - 1) + 1
      do 50 m = 1, nblok
      npq(m) = npt
   50 continue
      call MP_TASKSTART(idtask(i),PMOVEHX32,nargs,part(1,npo,1),edges,np
     1q,ihole(mpo,1),jssp(1,1,i),idimp,npt,mblok,nblok,idps,idds,mpt,mas
     2kp(npo,1),n)
c check for errors
      if (idtask(i).eq.0) then
        ierr = -1
        return
      endif
   60 continue
c finish finding outgoing particles
      npo = npl + 1
      do 70 m = 1, nblok
      npp(m) = npp(m) - npl
   70 continue
      npl = npmax - npl
      mpo = mpl + 1
      mps = ntmax - mpl
      call PMOVEHX32(part(1,npo,1),edges,npp,ihole(mpo,1),jss,idimp,npl,
     1mblok,nblok,idps,idds,mps,maskp(npo,1),n)
      npl = npmax - npl
      do 80 m = 1, nblok
      npp(m) = npp(m) + npl
   80 continue
c wait for tasks to complete
      do 110 i = 1, nmt
      call MP_TASKWAIT(idtask(i))
      if (idtask(i).ne.0) ierr = -2
c compress ihole array
      if (i.eq.1) go to 110
      npo = npt*(i - 1)
      mpo = mpt*(i - 1)
      do 100 m = 1, mnblok
      do 90 j = 1, jssp(1,m,i)
      ihole(jssp(1,m,i-1)+j,m) = ihole(j+mpo,m) + npo
   90 continue
      jssp(1,m,i) = jssp(1,m,i) + jssp(1,m,i-1)
      jssp(2,m,i) = max0(jssp(2,m,i),jssp(2,m,i-1))
  100 continue
  110 continue
      if (nmt.gt.0) then
         do 130 m = 1, mnblok
         do 120 j = 1, jss(1,m)
         ihole(jssp(1,m,nmt)+j,m) = ihole(j+mpl,m) + npl
  120    continue
         jss(1,m) = jss(1,m) + jssp(1,m,nmt)
         jss(2,m) = max0(jss(2,m),jssp(2,m,nmt))
  130    continue
      endif
      call PWTIMERA(1,tf,dtime)
      th = th + tf
c send outgoing particles
      call PMOVESS32(part,edges,npp,sbufr,sbufl,rbufr,rbufl,ihole,jsr,js
     1l,jss,ny,nz,kstrt,nvpy,nvpz,idimp,npmax,mblok,nblok,idps,nbmax,idd
     2s,ntmax,info,n)
c particle or iteration overflow error
      if (info(1).ne.0) return
  140 continue
      return
      end
