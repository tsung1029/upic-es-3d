!-----------------------------------------------------------------------
! * * * periodic 3d electrostatic particle simulation kernel code * * *
! this is a simple 3d skeleton particle-in-cell code designed for
! exploring new computer architectures.  it contains the critical pieces
! needed for depositing charge and current, advancing particles, and
! solving the fields.  the code moves only electrons, with periodic
! electromagnetic forces obtained by solving maxwells's equation with
! fast fourier transforms.
! the only diagnostic is particle and field energy.
! portable gcpic kernel code, using algorithm described in:
! p. c. liewer and v. k. decyk, j. computational phys. 85, 302 (1989).
! written by viktor k. decyk, ucla
! for mpi distributed memory computers, with 2D domain decomposition
! update: february 18, 2008
      program pbbeps32
      use pinit32d
      use pempush32d
      use pfield32d
      use pdiag32d
      use mp0d, only: mpinit, ncpus
      implicit none
! idps = number of particle partition boundaries = 4
! idds = dimensionality of domain decomposition = 2
! idimp = dimension of phase space = 6
! mshare = (0,1) = (no,yes) architecture is shared memory
      integer :: idps =    4, idds =    2, idimp =   6, mshare =   0
! nmv = number of segments in v for velocity distribution
! ipbc = particle boundary condition = (0,1,2,3) =
! (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
      integer :: nmv = 40, vect = 0, ipbc = 1
      integer :: npxyz, npxyzb, np, npxyzi, npxyzbi, npi
      integer :: nx, ny, nz, nxh, nyh, nzh, nyv, nzv, nxe, nxeh
      integer :: nloop, nvpy, nvpz, nvp, iblok, nblok, inblok
      integer :: npmax, npimax = 0, kyp, kzp, nypmx, nzpmx, kxyp, kyzp
      integer :: kyb, kxb, kzb, kyzb, kxyb, kzyb, kbmin, kblok
      integer :: jbmin, jblok, lbmin, lblok, mbmin, mblok
      integer :: ngds, nxyzh, nxhyz, nx1, nyzpm1, nbmax
      integer :: idproc, id0, kstrt, itime, isign, irc, ierr
      integer :: it, modesz2
      integer :: ntasks
      real :: zero = 0.0, wki = 0.0
      real :: tpush = 0.0, tdpost = 0.0, tdjpost = 0.0, tsort = 0.0
      real :: tpushi = 0.0, tdposti = 0.0, tdjposti = 0.0, tsorti = 0.0
      real :: totpush = 0.0, totpushi = 0.0
      real :: qbme, qbmi, affp, dth, qi0, etx, we, wf, wm, wef, wke
      real :: vtxi, vtyi, vtzi, vtdxi, vtdyi, vtdzi
      double precision :: dtime
      real, dimension(:,:,:), pointer :: part, part2, parti, parti2
      real, dimension(:,:,:,:), pointer :: qe, qi
      real, dimension(:,:,:,:,:), pointer :: cu, fxyze, bxyze
      complex, dimension(:,:,:,:), pointer :: qt, qs
      complex, dimension(:,:,:,:,:), pointer :: cut, fxyzt, bxyzt
      complex, dimension(:,:,:,:,:), pointer :: fxyzs
      complex, dimension(:,:,:,:,:), pointer :: exyz, bxyz
      complex, dimension(:,:,:,:), pointer :: ffc
      integer, dimension(:), pointer :: mixup
      complex, dimension(:), pointer :: sct
      real, dimension(:,:), pointer :: edges
      integer, dimension(:,:), pointer :: nyzp, noff
      integer, dimension(:), pointer :: npp, nppi, nps
      real, dimension(:,:), pointer :: pt
      integer, dimension(:,:), pointer :: ip, npic
      real, dimension(:,:,:,:), pointer :: sfield
      complex, dimension(:,:,:,:), pointer :: sfieldt, pott
      real, dimension(:,:,:,:,:), pointer :: vfield
      real, dimension(:,:,:), pointer :: fv, fvm, fvi, fvmi
      real, dimension(:,:), pointer :: wt
      integer, dimension(2) :: ktime
! wtot = total energy
      real, dimension(7) :: wtot
! time = timing array
      real, dimension(2) :: tfft = 0.0, time = 0.0
      real, dimension(2) :: tmove = 0.0, tmovi = 0.0
      real, dimension(2) :: tjmove = 0.0, tjmovi = 0.0
! msg = heartbeat array
      double precision, dimension(12) :: msg
      character(len=10) :: cdrun
      character(len=32) :: fname
      character(len=12) :: label
  991 format (' T = ',i7)
  992 format (' field, kinetic, total energies = ',3e14.7)
  993 format (' electric(l,t), magnetic energies = ',3e14.7)
! nvp = number of real or virtual processors
! initialize for parallel processing
      call PPINIT(idproc,id0,nvp)
      kstrt = idproc + 1
      if (id0==0) then
         open(unit=8,file='pinput3',form='formatted',status='old')
         read (8,pinput3)
! create string from idrun
         write (cdrun,'(i10)') idrun
         cdrun = adjustl(cdrun)
! text output file
         fname = 'poutput3.'//cdrun
         open(unit=18,file=trim(fname),form='formatted',status='replace'&
     &)
         write (18,pinput3,iostat=irc)
         if (irc /= 0) write (18,*) 'pinput3 namelist not written'
! open initial diagnostic metafile
         fname = 'pdiag3.init.'//cdrun
         open(unit=19,file=trim(fname),form='formatted',status='replace'&
     &)
      endif
! broadcast namelist to other nodes
      call sendnml()
! np = total number of electrons in simulation
      npxyz = npx*npy*npz; npxyzb = npxb*npyb*npzb; np = npxyz + npxyzb
! npi = total number of ions in simulation
      npxyzi = npxi*npyi*npzi; npxyzbi = npxbi*npybi*npzbi
      npi = npxyzi + npxyzbi
      nx = 2**indx; ny = 2**indy; nz = 2**indz; nxh = nx/2; nzh = nz/2
      nyv = ny + 2; nzv = nz + 2; nxe = nx + 4
! nvpy/nvpz = number of real or virtual processors in y/z
      call fcomp(nvp,nx,ny,nz,nvpy,nvpz,ierr)
      nvp = nvpy*nvpz
      if (id0==0) write (18,*) 'processor partition used: nvpy, nvpz = '&
     &, nvpy, nvpz
! kyp = number of complex grids in each field partition in y direction
      kyp = (ny - 1)/nvpy + 1
! kzp = number of complex grids in each field partition in z direction
      kzp = (nz - 1)/nvpz + 1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
      nypmx = kyp + 3; nzpmx = kzp + 3
! ngds = number of guard cells
      ngds = 3*((idds - 1)/2 + 1)
      if (inorder==LINEAR) then
         ax = .912871; ay = .912871; az = .912871
         nxe = nx + 2; nypmx = kyp + 1; nzpmx = kzp + 1
         ngds = (idds - 1)/2 + 1
      endif
      nxeh = nxe/2
! initialize for multiprocessing
      ntasks = mpinit(sntasks)
      call MP_SETSTACK(262144)
      if (popt==VECTOR) vect = 1
! nloop = number of time steps in simulation
      nloop = tend/dt + .0001
! iblok/nblok = number of particle partitions in y/z
      iblok = 1 + mshare*(nvpy - 1); nblok = 1 + mshare*(nvpz - 1)
      inblok = iblok*nblok
! npmax = maximum number of particles in each partition
      npmax = (np/nvp)*1.05 + 60000
      if (movion==1) npimax = (npi/nvp)*1.05 + 60000
! kxyp = number of complex grids in each field partition in x direction
! kyzp = number of complex grids in each field partition in y direction,
! in transposed data
      kxyp = (nxh - 1)/nvpy + 1; kyzp = (ny - 1)/nvpz + 1
! kyb = number of processors in y
! kxb = number of processors in x
! kzb = number of processors in z
! kyzb = number of processors in y, in transposed data
      kyb = ny/kyp; kxb = nxh/kxyp; kzb = nz/kzp; kyzb = ny/kyzp
! kxyb = maximum(kxb,kyb)
! kzyb = maximum(kyzb,kzb)
      kxyb = max(kxb,kyb); kzyb = max(kyzb,kzb)
! kblok = number of field partitions in y direction
      kbmin = 1 + (1 - mshare)*(kxyb/kxb - 1)
      kblok = 1 + mshare*(ny/kyp - 1)
! jblok = number of field partitions in x direction
      jbmin = 1 + (1 - mshare)*(kxyb/kyb - 1)
      jblok = 1 + mshare*(nxh/kxyp - 1)
! lblok = number of field partitions in z direction
      lbmin = 1 + (1 - mshare)*(kzyb/kyzb - 1)
      lblok = 1 + mshare*(nz/kzp - 1)
! mblok = number of field partitions in x direction
      mbmin = 1 + (1 - mshare)*(kzyb/kzb - 1)
      mblok = 1 + mshare*(ny/kyzp - 1)
! nxyzh = maximum(nx,ny,nz)/2
      nxyzh = max(nx,ny,nz)/2
! nxhyz = maximum(nx/2,ny,nz)
      nxhyz = max(nxh,ny,nz)
! dimensions for index and sorting arrays
      nx1 = nx + 1; nyzpm1 = (kyp + 1)*(kzp + 1)
! nbmax = size of buffer for passing particles between processors
      nbmax = 1 + (2*(npxyz*vtz + npxyzb*vtdz) + 1.4*npxyzb*abs(vdz))*dt&
     &/nz
      if (movion==1) then
         vtxi = vtx/sqrt(rmass*rtempxi)
         vtyi = vty/sqrt(rmass*rtempyi)
         vtzi = vtz/sqrt(rmass*rtempzi)
      endif
! diagnostic information needed by diagnostic nodes
! velocity diagnostics
      if (ntv > 0) then
         allocate(fv(2*nmv+2,3,nblok),fvm(2,3,nblok))
         if (id0==0) then
            fname = 'fv3.'//cdrun
            open(unit=10,file=trim(fname),form='formatted',status='unkno&
     &wn')
! write captions
         write (10,*) 'it vdx vdy vdz vtx vty vtz'
         endif
         if (movion==1) then
            allocate(fvi(2*nmv+2,3,nblok),fvmi(2,3,nblok))
            if (id0==0) then
               fname = 'fvi3.'//cdrun
               open(unit=20,file=trim(fname),form='formatted',status='un&
     &known')
! write captions
            write (20,*) 'it vdxi vdyi vdzi vtxi vtyi vtzi'
            endif
         endif
      endif
! density or potential diagnostics
      if ((ntp > 0) .or. (ntd > 0)) then
         allocate(sfield(nxe,nypmx,nzpmx*kbmin,kblok*lblok))
      endif
! potential diagnostics
      if (ntp > 0) then
         if (modesx > nxh) modesx = nxh
         if (modesy > nyh) modesy = nyh
         if (modesz > nzh) modesz = nzh
         modesz2 = 2*modesz - 1
         allocate(pott(modesz2,min(modesx,kxyp),min(modesy,kyzp),jblok*m&
     &blok))
         if (id0==0) then
            write (19,ppot3d,iostat=irc)
            if (irc /= 0) write (18,*) 'ppot3d namelist not written'
         endif
      endif
! vector potential diagnostics
      if (nta > 0) then
         allocate(vfield(3,nxe,nypmx,nzpmx*kbmin,kblok*lblok))
      endif
! energy diagnostics
      if (ntw > 0) allocate(wt((nloop-1)/ntw+1,4))
! open restart files
      if (nustrt /= 1) then
         if (id0==0) then
            fname = 'rstrt1.'//cdrun
            open(unit=16,file=trim(fname),form='unformatted',status='old&
     &')
            fname = 'rstrt2.'//cdrun
            open(unit=17,file=trim(fname),form='unformatted',status='old&
     &')
         endif
      else if (ntr > 0) then
         if (id0==0) then
            fname = 'rstrt1.'//cdrun
            open(unit=16,file=trim(fname),form='unformatted',status='unk&
     &nown')
            fname = 'rstrt2.'//cdrun
            open(unit=17,file=trim(fname),form='unformatted',status='unk&
     &nown')
         endif
      endif
!
! diagnostic nodes have special processing
      if (idproc < 0) call diag32nodes
!
! part(1,n,m) = position x of particle n in partition m
! part(2,n,m) = position y of particle n in partition m
! part(3,n,m) = position z of particle n in partition m
! part(4,n,m) = velocity vx of particle n in partition m
! part(5,n,m) = velocity vy of particle n in partition m
! part(6,n,m) = velocity vz of particle n in partition m
      allocate(part(idimp,npmax,inblok))
! maskp = scratch array for particle addresses
!     allocate(maskp(npmax,inblok))
! in real space, qe(j+1,k,l,m) = charge density at grid point (j,kk,ll)
! in real space, cu(i,j+1,k,l,m) = i component of current density
! at grid point (j,kk,ll),
! where kk = k + noff(1,m) - 1, ll = l + noff(2,m) - 1
      allocate(qe(nxe,nypmx,nzpmx*kbmin,kblok*lblok))
      allocate(cu(3,nxe,nypmx,nzpmx*kbmin,kblok*lblok))
! in real space, fxyze(i,j+1,k,l) = i component of force/charge at 
! grid point (j,kk,ll)
! in other words, fxyze are the convolutions of the electric field
! over the particle shape, where kk = k + noff(1,m) - 1, and
! ll = l + noff(2,m) - 1
      allocate(fxyze(3,nxe,nypmx,nzpmx*kbmin,kblok*lblok))
! bxyze(i,j+1,k,l) = i component of magnetic field at 
! grid point (j,kk,ll), that is, bxyze is the convolution of the
! magnetic field over the particle shape, where kk = k + noff(1,m) - 1,
! and ll = l + noff(2,m) - 1
      allocate(bxyze(3,nxe,nypmx,nzpmx*kbmin,kblok*lblok))
! qt(l,k,j,m) = complex charge density for fourier mode jj-1,kk-1,l-1
! cut(i,l,k,j,m) = complex curent density for fourier mode jj-1,kk-1,l-1
      allocate(qt(nzv,kxyp,kyzp*mbmin,jblok*mblok))
      allocate(cut(3,nzv,kxyp,kyzp*mbmin,jblok*mblok))
! fxyzt(i,l,k,j,m) = i component of force/charge
! fxyzt(i,l,k,j,m) = i component of magnetic field
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*(my - 1), and
! kk = k + kyzp*(mz - 1), m = my + nvpy*(mz - 1)
      allocate(fxyzt(3,nzv,kxyp,kyzp*mbmin,jblok*mblok))
      allocate(bxyzt(3,nzv,kxyp,kyzp*mbmin,jblok*mblok))
! in fourier space, exyz = transverse electric field
! in fourier space, bxyz = magnetic field
      allocate(exyz(3,nzv,kxyp,kyzp*mbmin,jblok*mblok))
      allocate(bxyz(3,nzv,kxyp,kyzp*mbmin,jblok*mblok))
      allocate(qs(nyv,kxyp,nzpmx*jbmin*lbmin,jblok*lblok))
      allocate(fxyzs(3,nyv,kxyp,nzpmx*jbmin*lbmin,jblok*lblok))
! ffc = form factor array for poisson solver
      allocate(ffc(nzh,kxyp,kyzp,jblok*mblok))
! mixup, sct = arrays for fft
      allocate(mixup(nxhyz),sct(nxyzh))
! edges(1,m) = lower boundary in y of particle partition m
! edges(2,m) = upper boundary in y of particle partition m
! edges(3,m) = back boundary in z of particle partition m
! edges(4,m) = front boundary in z of particle partition m
      allocate(edges(idps,inblok))
! nyzp(1,m) = number of primary gridpoints in y in particle partition m
! nyzp(2,m) = number of primary gridpoints in z in particle partition m
! noff(1,m) = lowermost global gridpoint in y in particle partition m
! noff(2,m) = backmost global gridpoint in z in particle partition m
      allocate(nyzp(idds,inblok),noff(idds,inblok))
! npp(m) = number of particles in partition m
! nps(m) = starting address of particles in partition m
      allocate(npp(inblok),nps(inblok))
! sorting arrays
      allocate(pt(max(npmax,npimax),inblok))
      allocate(ip(max(npmax,npimax),inblok),npic(nyzpm1,inblok))
      if (sortime > 0) then
         allocate(part2(idimp,npmax,inblok))
      else
         allocate(part2(0,0,0))
      endif
! initialize parallel timer
      call pwtimer(time,dtime,-1)
! initialize constants
      itime = 0
      qbme = qme
      affp = float(nx*ny*nz)/float(np)
      dth = .5*dt
      if (movion==1) then
         qbmi = qmi/rmass
         vtdxi = vtx/sqrt(rmass*rtempdxi)
         vtdyi = vty/sqrt(rmass*rtempdyi)
         vtdzi = vtz/sqrt(rmass*rtempdzi)
      endif
! diagnostics
! velocity diagnostics
      if (ntv > 0) then
         fv(1,:,:) = 8.*max(vtx,vty,vtz)
         if (movion==1) fvi(1,:,:) = 8.*max(vtxi,vtyi,vtzi)
      endif
! density or potential diagnostics
      if ((ntp > 0) .or. (ntd > 0)) then
         allocate(sfieldt(nzv,kxyp,kyzp*mbmin,jblok*mblok))
      endif
! calculate partition variables
      call dcomp(edges,nyzp,noff,ny,nz,kstrt,nvpy,nvpz,iblok,inorder)
! initialize electromagnetic magnetic fields
      bxyze = 0.
      bxyz = cmplx(0.,0.)
      exyz = cmplx(0.,0.)
      call cguard(bxyze,nyzp,nx,inorder)
! prepare fft tables
      call fft_init(mixup,sct,indx,indy,indz)
! calculate form factors
      call pois_init(ffc,ax,ay,az,affp,nx,ny,nz,kstrt,jblok)
! allocate background charge density
      if (movion==0) allocate(qi(nxe,nypmx,nzpmx*kbmin,kblok*lblok))
! allocate ion data
      if (movion==1) then
         allocate(parti(idimp,npimax,inblok),nppi(inblok))
         allocate(parti2(0,0,0))
      endif
! new start
      if (nustrt==1) then
! initialize electrons
         nps = 1
         npp = 0
! background electrons
!        if (npxyz > 0) call distr(part,edges,npp,nps,vtx,vty,vtz,vx0,vy&
!    &0,vz0,npx,npy,npz,nx,ny,nz,ipbc,iblok)
         if (npxyz > 0) then
            call fdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,shi&
     &ftdy,ampdz,scaledz,shiftdz,npx,npy,npz,nx,ny,nz,kstrt,nvp,ipbc,ndp&
     &rof,nsrand,iblok)
            call vdistr(part,npp,nps,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,npz&
     &,kstrt,nvp,iblok)
         endif
! beam electrons
         nps = npp + 1
!        if (npxyzb > 0) call distr(part,edges,npp,nps,vtdx,vtdy,vtdz,vd&
!    &x,vdy,vdz,npxb,npyb,npzb,nx,ny,nz,ipbc,iblok)
         if (npxyzb > 0) then
            call fdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,shi&
     &ftdy,ampdz,scaledz,shiftdz,npxb,npyb,npzb,nx,ny,nz,kstrt,nvp,ipbc,&
     &ndprof,nsrand,iblok)
            call vdistr(part,npp,nps,vtdx,vtdy,vtdz,vdx,vdy,vdz,npxb,npy&
     &b,npzb,kstrt,nvp,iblok)
         endif
! calculate actual coordinates from guiding centers
!        if (relativity==1) then
!           call distr(part,bxyze,npp,noff,qbme,ci,nx,ny,nz,ipbc,iblok,i&
!    &norder)
!        else
!           call distr(part,bxyze,npp,noff,qbme,nx,ny,nz,ipbc,iblok,inor&
!    &der)
!        endif
! move electrons into appropriate spatial regions
         call pmove(part,edges,npp,tmove,ny,nz,kstrt,nvpy,nvpz,nbmax,idd&
     &s,iblok,vect,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! initialize ions
         if (movion==1) then
            nps = 1
            nppi = 0
! background ions
!           if (npxyzi > 0) call distr(parti,edges,nppi,nps,vtxi,vtyi,vt&
!    &zi,vxi0,vyi0,vzi0,npxi,npyi,npzi,nx,ny,nz,ipbc,iblok)
            if (npxyzi > 0) then
               call fdistr(parti,nps,ampdxi,scaledxi,shiftdxi,ampdyi,sca&
     &ledyi,shiftdyi,ampdzi,scaledzi,shiftdzi,npxi,npyi,npzi,nx,ny,nz,ks&
     &trt,nvp,ipbc,ndprofi,nsrandi,iblok)
               call vdistr(parti,nppi,nps,vtxi,vtyi,vtzi,vxi0,vyi0,vzi0,&
     &npxi,npyi,npzi,kstrt,nvp,iblok)
            endif
! beam ions
            nps = nppi + 1
!           if (npxyzbi > 0) call distr(parti,edges,nppi,nps,vtdxi,vtdyi&
!    &,vtdzi,vdxi,vdyi,vdzi,npxbi,npybi,npzbi,nx,ny,nz,ipbc,iblok)
            if (npxyzbi > 0) then
               call fdistr(parti,nps,ampdxi,scaledxi,shiftdxi,ampdyi,sca&
     &ledyi,shiftdyi,ampdzi,scaledzi,shiftdzi,npxbi,npybi,npzbi,nx,ny,nz&
     &,kstrt,nvp,ipbc,ndprofi,nsrandi,iblok)
               call vdistr(parti,nppi,nps,vtdxi,vtdyi,vtdzi,vdxi,vdyi,vd&
     &zi,npxbi,npybi,npzbi,kstrt,nvp,iblok)
            endif
! calculate actual coordinates from guiding centers
!           if (relativity==1) then
!              call distr(parti,bxyze,nppi,noff,qbmi,ci,nx,ny,nz,ipbc,ib&
!    &lok,inorder)
!           else
!              call distr(parti,bxyze,nppi,noff,qbmi,nx,ny,nz,ipbc,iblok&
!    &,inorder)
!           endif
! move ions into appropriate spatial regions
            call pmove(parti,edges,nppi,tmovi,ny,nz,kstrt,nvpy,nvpz,nbma&
     &x,idds,iblok,vect,ierr)
            if (ierr /= 0) then
               call MP_END
               call PPEXIT
               stop
            endif
         endif
! initialize background charge density
         if (movion==0) then
            qi0 = -qme/affp
!           call sguard(qi,nyzp,qi0,nx,inorder)
            call sguard(qi,nyzp,zero,nx,inorder)
            call dpost(part,qi,-qme,npp,noff,tdpost,inorder,dopt)
! freeze the ions now
         else if ((movion==1).and.(itime==ionoff)) then
            allocate(qi(nxe,nypmx,nzpmx*kbmin,kblok*lblok))
! initialize ion charge density to zero
            call sguard(qi,nyzp,zero,nx,inorder)
! deposit ion charge
            call dpost(parti,qi,qmi,nppi,noff,tdposti,inorder,dopt)
! delete ions
            deallocate(parti,nppi)
            movion = 0
         endif
! retard electron velocities to deposit current
         if (relativity==1) then
            call retard(part,npp,dth,ci,nx,ny,nz,ipbc)
         else
            call retard(part,npp,dth,nx,ny,nz,ipbc)
         endif
! move electrons into appropriate spatial regions
         call pmove(part,edges,npp,tmove,ny,nz,kstrt,nvpy,nvpz,nbmax,idd&
     &s,iblok,vect,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! retard ion velocities to deposit current
         if (movion==1) then
            if (relativity==1) then
               call retard(parti,nppi,dth,ci,nx,ny,nz,ipbc)
            else
               call retard(parti,nppi,dth,nx,ny,nz,ipbc)
            endif
! move ions into appropriate spatial regions
            call pmove(parti,edges,nppi,tmovi,ny,nz,kstrt,nvpy,nvpz,nbma&
     &x,idds,iblok,vect,ierr)
            if (ierr /= 0) then
               call MP_END
               call PPEXIT
               stop
            endif
         endif
! restart
      else
! determine most recent restart file
         if (id0==0) then
            read (16,iostat=ierr) ktime(1)
            if (ierr /= 0) ktime(1) = -1
            read (17,iostat=ierr) ktime(2)
            if (ierr /= 0) ktime(2) = -1
            if (ktime(1) > ktime(2)) then
               ktime(2) = 16
            else
               ktime(1) = ktime(2)
               ktime(2) = 17
            endif
         endif
         call plbcast(ktime)
         itime = ktime(1)
         if (itime < 0) go to 400
! read restart file
         it = ktime(2)
         call rddata(part,npp,it,ierr)
         if (ierr /= 0) go to 400
         call rddata(exyz,nvp,it,ierr)
         if (ierr /= 0) go to 400
         call rddata(bxyz,nvp,it,ierr)
         if (ierr /= 0) go to 400
         if (movion==1) then
            call rddata(parti,nppi,it,ierr)
            if (ierr /= 0) go to 400
         endif
         if (movion==0) then
            call rddata(qi,nvp,it,ierr)
            if (ierr /= 0) go to 400
         endif
         if (ntw > 0) then
            call rddata(wt,1,it,ierr)
            if (ierr /= 0) go to 400
            call plbcast(wt)
         endif
         if (ntp > 0) then
            if (id0==0) then
               read (it,iostat=ierr) ktime(1)
               if (ierr /= 0) ktime(1) = -1
               irc = 0
               fname = 'ppotk3.'//cdrun
               call writef(pott,modesx,modesy,modesz2,kxyp,kyzp,jblok,11&
     &,irc,trim(fname))
            endif
            call plbcast(ktime)
            nprec = ktime(1)
            if (nprec< 0) go to 400
         endif
         if (id0==0) then
            read (it,iostat=ierr) ktime(1)
            if (ierr /= 0) ktime(1) = -1
            rewind it
         endif
         call plbcast(ktime)
         irc = ktime(1)
         if (irc==itime) go to 490
! handle error
  400    if (id0==0) write (18,*) 'Restart Error'
         go to 3000
      endif
! record time
  490 call pwtimer(time,dtime)
! send initial CPU Time to diagnostic nodes
      msg(1) = time(1); msg(2) = time(2)
      call HARTBEAT(msg,2)
      if (id0==0) then
         write (18,*) 'init max/min real time=', time(1), time(2), 'sec'
      endif
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= itime) go to 2000
! send time step to diagnostic nodes
      msg(1) = itime
      call HARTBEAT(msg,1)
      if (id0==0) write (18,991) itime
      write (label,991) itime
      call LOGNAME(label)
! initialize current density to background
      call sguard(cu,nyzp,zero,zero,zero,nx,inorder)
! deposit current
      if (relativity==1) then
         call rdjpost(part,cu,npp,noff,qme,dth,ci,tdjpost,nx,ny,nz,ipbc,&
     &inorder,djopt)
      else
         call djpost(part,cu,npp,noff,qme,dth,tdjpost,nx,ny,nz,ipbc,inor&
     &der,djopt)
      endif
! move electrons into appropriate spatial regions
      call pmoves(part,edges,npp,tjmove,ny,nz,kstrt,nvpy,nvpz,nbmax,idds&
     &,iblok,vect,ierr)
      if (ierr /= 0) then
         call MP_END
         call PPABORT
         stop
      endif
! initialize charge density to background
      call sguard(qe,nyzp,zero,nx,inorder)
! deposit charge
      call dpost(part,qe,qme,npp,noff,tdpost,inorder,dopt)
! density diagnostic
      if (ntp > 0) then
         it = itime/ntp
         if (itime==ntp*it) then
            sfield = -qe
! add guard cells for density in x
            call aguard(sfield,nyzp,nx,inorder)
! add guard cells for density in y and z
            call paguard(sfield,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,iblok,in&
     &order)
! transform density to fourier space
            isign = -1
            call fft(sfield,qs,qt,isign,mixup,sct,tfft,indx,indy,indz,ks&
     &trt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! calculate smoothing in fourier space
            call pois(qt,sfieldt,ffc,nx,ny,nz,kstrt,jblok)
! transform density to real space
            isign = 1
            call fft(sfield,qs,sfieldt,isign,mixup,sct,tfft,indx,indy,in&
     &dz,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! copy to guard cells
            call pcguard(sfield,kstrt,nvpy,nvpz,kyp,kzp,ngds,iblok,inord&
     &er)
            call cguard(sfield,nyzp,nx,inorder)
         endif
      endif
! deposit ion current
      if (movion==1) then
         if (relativity==1) then
            call rdjpost(parti,cu,nppi,noff,qmi,dth,ci,tdjposti,nx,ny,nz&
     &,ipbc,inorder,djopt)
         else
            call djpost(parti,cu,nppi,noff,qmi,dth,tdjposti,nx,ny,nz,ipb&
     &c,inorder,djopt)
         endif
! move ions into appropriate spatial regions
         call pmoves(parti,edges,nppi,tjmovi,ny,nz,kstrt,nvpy,nvpz,nbmax&
     &,idds,iblok,vect,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPABORT
            stop
         endif
! add ion charge density
         call dpost(parti,qe,qmi,nppi,noff,tdposti,inorder,dopt)
      else
         qe = qe + qi
      endif
! add guard cells for current in x
      call aguard(cu,nyzp,nx,inorder)
! add up guard cells for current in y and z
      call paguard(cu,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,iblok,inorder)
! add guard cells for density in x
      call aguard(qe,nyzp,nx,inorder)
! add guard cells for density in y and z
      call paguard(qe,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,iblok,inorder)
! freeze the ions
      if ((movion==1).and.(itime==ionoff)) then
         allocate(qi(nxe,nypmx,nzpmx*kbmin,kblok*lblok))
! initialize ion charge density to zero
         call sguard(qi,nyzp,zero,nx,inorder)
! deposit ion charge
         call dpost(parti,qi,qmi,nppi,noff,tdposti,inorder,dopt)
! delete ions
         deallocate(parti,nppi)
         movion = 0
      endif
! velocity diagnostic
      if (ntv > 0) then
         it = itime/ntv
         if (itime==ntv*it) then
! calculate electron distribution function and moments
            call vdist(part,fv,fvm,npp,nmv)
            call plsum(fv(:,:,1))
            fv(1,:,:) = 8.*max(vtx,vty,vtz)
! print out velocity moments
            if (id0==0) write (10,*) it, fvm(1,:,:), fvm(2,:,:)
            if (movion==1) then
! calculate ion distribution function and moments
               call vdist(parti,fvi,fvmi,nppi,nmv)
               call plsum(fvi(:,:,1))
               fvi(1,:,:) = 8.*max(vtxi,vtyi,vtzi)
! print out velocity moments
               if (id0==0) write (20,*) it, fvmi(1,:,:), fvmi(2,:,:)
            endif
         endif
      endif
! transform charge to fourier space
      isign = -1
      call fft(qe,qs,qt,isign,mixup,sct,tfft,indx,indy,indz,kstrt,kxyp,k&
     &yp,kyzp,kzp,kblok,mblok,inorder)
! potential diagnostic
      if (ntp > 0) then
         it = itime/ntp
         if (itime==ntp*it) then
! calculate potential in fourier space
            call pois(qt,sfieldt,ffc,we,nx,ny,nz,kstrt,jblok)
! store selected fourier modes
!           call gtmodes(sfieldt,pott,nx,ny,nz,modesx,modesy,modesz,kstr&
!    &t,jblok)
! transform potential to real space
            isign = 1
            call fft(sfield,qs,sfieldt,isign,mixup,sct,tfft,indx,indy,in&
     &dz,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! copy to guard cells
            call pcguard(sfield,kstrt,nvpy,nvpz,kyp,kzp,ngds,iblok,inord&
     &er)
            call cguard(sfield,nyzp,nx,inorder)
! write diagnostic output
!           if (nprec==0) then
!              nprec = -1
!              fname = 'ppotk3.'//cdrun
!              call writef(pott,modesx,modesy,modesz2,kxyp,kyzp,jblok,11&
!    &,nprec,trim(fname))
!           else
!              call writef(pott,modesx,modesy,modesz2,kxyp,kyzp,jblok,11&
!    &,nprec)
!           endif
         endif
      endif
! transform current to fourier space
      isign = -1
      call fft(cu,fxyzs,cut,isign,mixup,sct,tfft,indx,indy,indz,kstrt,kx&
     &yp,kyp,kyzp,kzp,kblok,mblok,inorder)
!     call fftn(cu,fxyzs,cut,isign,mixup,sct,tfft,indx,indy,indz,kstrt,k&
!    &xyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! take transverse part of current
      call cuperp(cut,nx,ny,nz,kstrt,jblok)
! calculate magnetic field in fourier space
      if (itime==0) then
         call ibpois(cut,bxyz,ffc,ci,wm,nx,ny,nz,kstrt,jblok)
         wf = 0.
      else
         call maxwel(exyz,bxyz,cut,ffc,affp,ci,dt,wf,wm,nx,ny,nz,kstrt,j&
     &blok)
      endif
! vector potential diagnostic
      if (nta > 0) then
         it = itime/nta
         if (itime==nta*it) then
! calculate vector potential in fourier space
!           call bpois(cut,bxyzt,ffc,ci,wm,nx,ny,nz,kstrt,jblok)
            call avpot(bxyz,bxyzt,nx,ny,nz,kstrt,jblok)
! transform vector potential to real space
            isign = 1
            call fft(vfield,fxyzs,bxyzt,isign,mixup,sct,tfft,indx,indy,i&
     &ndz,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
            call pcguard(vfield,kstrt,nvpy,nvpz,kyp,kzp,ngds,iblok,inord&
     &er)
            call cguard(vfield,nyzp,nx,inorder)
         endif
      endif
! calculate force/charge in fourier space
      call pois(qt,fxyzt,ffc,we,nx,ny,nz,kstrt,jblok)
! add longitudinal and transverse electric fields
      isign = 1
      call emfield(fxyzt,exyz,ffc,isign,nx,ny,nz,kstrt,jblok)
! copy magnetic field
      isign = -1
      call emfield(bxyzt,bxyz,ffc,isign,nx,ny,nz,kstrt,jblok)
! transform force/charge to real space
      isign = 1
      call fft(fxyze,fxyzs,fxyzt,isign,mixup,sct,tfft,indx,indy,indz,kst&
     &rt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
!     call fftn(fxyze,fxyzs,fxyzt,isign,mixup,sct,tfft,indx,indy,indz,ks&
!    &trt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! copy data from field to particle partition, and copy to guard cells
      call pcguard(fxyze,kstrt,nvpy,nvpz,kyp,kzp,ngds,iblok,inorder)
      call cguard(fxyze,nyzp,nx,inorder)
! transform magnetic field to real space
      isign = 1
      call fft(bxyze,fxyzs,bxyzt,isign,mixup,sct,tfft,indx,indy,indz,kst&
     &rt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
!     call fftn(bxyze,fxyzs,bxyzt,isign,mixup,sct,tfft,indx,indy,indz,ks&
!    &trt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! copy data from field to particle partition, and copy to guard cells
      call pcguard(bxyze,kstrt,nvpy,nvpz,kyp,kzp,ngds,iblok,inorder)
      call cguard(bxyze,nyzp,nx,inorder)
! external pump
      if ((itpon > 0).and.(itime >= itpon)) then
         etx = (v0*vtx)*w0*cos(w0*dt*(itime - itpon))
         fxyze(1,:,:,:,:) = fxyze(1,:,:,:,:) + etx
      endif
! particle push and charge density update
      wke = 0.
! push electrons
      if (relativity==1) then
!        call rpush(part,fxyze,npp,noff,qbme,dt,ci,wke,tpush,nx,ny,nz,ip&
!    &bc,inorder,popt)
!        call rpush(part,fxyze,bxyze,npp,noff,qbme,dt,dt,ci,wke,tpush,nx&
!    &,ny,nz,ipbc,inorder,popt)
         call rpush(part,fxyze,bxyze,npp,noff,qbme,dt,dth,ci,wke,tpush,n&
     &x,ny,nz,ipbc,inorder,popt)
      else
!        call push(part,fxyze,npp,noff,qbme,dt,wke,tpush,nx,ny,nz,ipbc,i&
!    &norder,popt)
!        call push(part,fxyze,bxyze,npp,noff,qbme,dt,dt,wke,tpush,nx,ny,&
!    &nz,ipbc,inorder,popt)
         call push(part,fxyze,bxyze,npp,noff,qbme,dt,dth,wke,tpush,nx,ny&
     &,nz,ipbc,inorder,popt)
      endif
! move electrons into appropriate spatial regions
!     call pmove(part,edges,npp,tmove,ny,nz,kstrt,nvpy,nvpz,nbmax,idds,i&
!    &blok,vect,ierr)
      call pmoves(part,edges,npp,tmove,ny,nz,kstrt,nvpy,nvpz,nbmax,idds,&
     &iblok,vect,ierr)
      if (ierr /= 0) then
         call MP_END
         call PPABORT
         stop
      endif
! push ions
      if (movion==1) then
         wki = 0.
         if (relativity==1) then
!           call rpush(parti,fxyze,nppi,noff,qbmi,dt,ci,wke,tpushi,nx,ny&
!    &,nz,ipbc,inorder,popt)
!           call rpush(parti,fxyze,bxyze,nppi,noff,qbmi,dt,dt,ci,wke,tpu&
!    &shi,nx,ny,nz,ipbc,inorder,popt)
            call rpush(parti,fxyze,bxyze,nppi,noff,qbmi,dt,dth,ci,wki,tp&
     &ushi,nx,ny,nz,ipbc,inorder,popt)
         else
!           call push(parti,fxyze,nppi,noff,qbmi,dt,wke,tpushi,nx,ny,nz,&
!    &ipbc,inorder,popt)
!           call push(parti,fxyze,bxyze,nppi,noff,qbmi,dt,dt,wke,tpushi,&
!    &nx,ny,nz,ipbc,inorder,popt)
            call push(parti,fxyze,bxyze,nppi,noff,qbmi,dt,dth,wki,tpushi&
     &,nx,ny,nz,ipbc,inorder,popt)
         endif
         wki = wki*rmass
! move ions into appropriate spatial regions
!        call pmove(parti,edges,nppi,tmovi,ny,nz,kstrt,nvpy,nvpz,nbmax,i&
!    &dds,iblok,vect,ierr)
         call pmoves(parti,edges,nppi,tmovi,ny,nz,kstrt,nvpy,nvpz,nbmax,&
     &idds,iblok,vect,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPABORT
            stop
         endif
      endif
! sort electrons
      if (sortime > 0) then
         if (mod(itime,sortime)==0) then
!           call sortp(part,pt,ip,npp,noff,nyzp,npic,tsort,inorder)
            call sortp(part,part2,npp,noff,nyzp,npic,tsort,inorder)
         endif
      endif
! sort ions
      if ((movion==1) .and. (sortimi > 0)) then
         if (mod(itime,sortimi)==0) then
            call sortp(parti,pt,ip,nppi,noff,nyzp,npic,tsorti,inorder)
!           call sortp(parti,parti2,nppi,noff,nyzp,npic,tsorti,inorder)
         endif
      endif
! energy diagnostic
      if (ntw > 0) then
         it = itime/ntw
         if (itime==ntw*it) then
            wef = we + wf + wm
            wtot(1) = wef
            wtot(2) = wke
            wtot(3) = wki
            wtot(4) = wef + wke + wki
            wtot(5) = we
            wtot(6) = wf
            wtot(7) = wm
            call plsum(wtot)
! send energy values to diagnostic node
            msg(1:7) = wtot
            call HARTBEAT(msg,7)
            if (id0==0) write (18,992) wtot(1), wtot(2), wtot(4)
            if (id0==0) write (18,993) wtot(5), wtot(6), wtot(7)
         endif
      endif
      itime = itime + 1
! restart file
      if (ntr > 0) then
         it = itime/ntr
         if (itime==ntr*it) then
            it = 16 + mod(it-1,2)
            if (id0==0) write (it) itime
            call wrdata(part,npp,it)
            call wrdata(exyz,nvp,it)
            call wrdata(bxyz,nvp,it)
            if (movion==1) call wrdata(parti,nppi,it)
            if (movion==0) call wrdata(qi,nvp,it)
            if (ntw > 0) call wrdata(wt,1,it)
            if ((ntp > 0) .and. (id0==0)) write (it) nprec
            if (id0==0) then
               write (it) itime
               end file it
               rewind it
            endif
         endif
      endif
      go to 500
 2000 continue
!
! * * * end main iteration loop * * *
!
! send QUIT message to diagnostic nodes
      msg = -1.
      call HARTBEAT(msg,1)
      call pwtimer(time,dtime)
! send main CPU Time to diagnostic nodes
      msg(1:2) = time; msg(3) = tpush; msg(4) = tdpost; msg(5) = tdjpost
      msg(6) = tsort; msg(7:8) = tmove; msg(9:10) = tjmove
      msg(11:12) = tfft
      call HARTBEAT(msg,12)
      if (id0==0) then
         write (18,*) ncpus, ' processors found, ', ntasks+1, ' used'
         write (18,*) 'main max/min real time=', time(1), time(2), 'sec'
         totpush = tpush + tdpost + tdjpost
         write (18,*) 'electron push time=', tpush, 'sec'
         write (18,*) 'electron charge deposit time=', tdpost, 'sec'
         write (18,*) 'electron current deposit time = ', tdjpost, 'sec'
         write (18,*) 'total electron push time=', totpush, 'sec'
         write (18,*) 'electron sort time=', tsort, 'sec'
         tmove = tmove + tjmove
         write (18,*) 'electron move time=', tmove, 'sec'
         totpush = totpush + tsort + tmove(1)
         write (18,*) 'total electron time=', totpush, 'sec'
      endif
      if (movion==1) then
         msg(1) = tpushi; msg(2) = tdposti; msg(3) = tsorti
         msg(4:5) = tmovi; msg(6:7) = tjmovi
         call HARTBEAT(msg,7)
         if (id0==0) then
            totpushi = tpushi + tdposti + tdjposti
            write (18,*) 'ion push time=', tpushi, 'sec'
            write (18,*) 'ion charge deposit time=', tdposti, 'sec'
            write (18,*) 'ion current deposit time = ', tdjposti, 'sec'
            write (18,*) 'total ion push time=', totpushi, 'sec'
            write (18,*) 'ion sort time=', tsorti
            tmovi = tmovi + tjmovi
            write (18,*) 'ion move time=', tmovi
            totpushi = totpushi + tsorti + tmovi(1)
            write (18,*) 'total ion time=', totpushi, 'sec'
         endif
      endif
      if (id0==0) then
         write (18,*) 'total fft time=', tfft, 'sec'
         time(1) = time(1) - (totpush + totpushi + tfft(1))
         write (18,*) 'other time=', time(1), 'sec'
! write final diagnostic metafile
         close(unit=19)
         fname = 'pdiag3.final.'//cdrun
         open(unit=19,file=trim(fname),form='formatted',status='replace'&
     &)
! potential diagnostics
         if (ntp > 0) then
            nprec = nprec - 1
            write (19,ppot3d,iostat=irc)
         endif
! done
         write (18,*) '* * * q.e.d. * * *'
      endif
 3000 continue
      call MP_END
      call PPEXIT
      stop
!
      contains
!
         subroutine diag32nodes
         implicit none
! diagnostic nodes have special processing
  991    format (' T = ',i7)
  992    format (' field, kinetic, total energies = ',3e14.7)
  993    format (' electric(l,t), magnetic energies = ',3e14.7)
! allocate data for restart and/or phase space diagnostic
         if ((nts > 0) .or. (nustrt /= 1) .or. (ntr > 0)) then
            allocate(part(idimp,max(npmax,npimax),inblok))
            allocate(npp(inblok))
            if ((nustrt /= 1).or.(ntr > 0)) then
               allocate(exyz(3,nzv,kxyp,kyzp*mbmin,jblok*mblok))
               allocate(bxyz(3,nzv,kxyp,kyzp*mbmin,jblok*mblok))
            endif
         endif
         if (movion==0) allocate(qi(nxe,nypmx,nzpmx*kbmin,kblok*lblok))
! restart
         if (nustrt /= 1) then
! determine most recent restart file
            if (id0==0) then
               read (16,iostat=ierr) ktime(1)
               if (ierr /= 0) ktime(1) = -1
               read (17,iostat=ierr) ktime(2)
               if (ierr /= 0) ktime(2) = -1
               if (ktime(1) > ktime(2)) then
                  ktime(2) = 16
               else
                  ktime(1) = ktime(2)
                  ktime(2) = 17
               endif
            endif
            call plbcast(ktime)
            itime = ktime(1)
            if (itime < 0) go to 30
! read restart file
            it = ktime(2)
            call rddata(part,npp,it,ierr)
            if (ierr /= 0) go to 30
            call rddata(exyz,nvp,it,ierr)
            if (ierr /= 0) go to 30
            call rddata(bxyz,nvp,it,ierr)
            if (ierr /= 0) go to 30
            if (movion==1) then
               call rddata(part,npp,it,ierr)
               if (ierr /= 0) go to 30
            endif
            if (movion==0) then
               call rddata(qi,nvp,it,ierr)
               if (ierr /= 0) go to 30
            endif
            if (ntw > 0) then
               call rddata(wt,1,it,ierr)
               if (ierr /= 0) go to 30
               call plbcast(wt)
            endif
            if (ntp > 0) then
               if (id0==0) then
                  read (it,iostat=ierr) ktime(1)
                  if (ierr /= 0) ktime(1) = -1
                  irc = 0
                  fname = 'ppotk3.'//cdrun
                  call writef(pott,modesx,modesy,modesz2,kxyp,kyzp,jblok&
     &,11,irc,trim(fname))
               endif
               call plbcast(ktime)
               nprec = ktime(1)
               if (nprec< 0) go to 30
            endif
            if (id0==0) then
               read (it,iostat=ierr) ktime(1)
            if (ierr /= 0) ktime(1) = -1
               rewind it
            endif
            call plbcast(ktime)
            irc = ktime(1)
            if (irc==itime) go to 40
! handle error
   30       if (id0==0) write (18,*) 'Restart Error'
            call MP_END
            call PPEXIT
            stop
         endif
! get initial CPU Time
   40    call HARTBEAT(msg,2)
         time(1) = msg(1); time(2) = msg(2)
         if (id0==0) then
         write (18,*) 'init max/min real time=', time(1), time(2), 'sec'
         endif
! get time step
   10    call HARTBEAT(msg,1)
         it = msg(1)
         if (it < 0) then
! get main CPU Time
            call HARTBEAT(msg,12)
            time = msg(1:2); tpush = msg(3); tdpost = msg(4)
            tdjpost = msg(5); tsort = msg(6); tmove = msg(7:8)
            tjmove = msg(9:10); tfft = msg(11:12)
            if (id0==0) then
               write (18,*) ncpus, ' processors found, ', ntasks+1, ' us&
     &ed'
               write (18,*) 'main max/min real time=', time(1), time(2),&
     &'sec'
               totpush = tpush + tdpost + tdjpost
               write (18,*) 'electron push time=', tpush, 'sec'
               write (18,*) 'electron charge deposit time=', tdpost, 'se&
     &c'
               write (18,*) 'electron current deposit time = ', tdjpost,&
     &'sec'
               write (18,*) 'total electron push time=', totpush, 'sec'
               write (18,*) 'electron sort time=', tsort, 'sec'
               tmove = tmove + tjmove
               write (18,*) 'electron move time=', tmove, 'sec'
               totpush = totpush + tsort + tmove(1)
               write (18,*) 'total electron time=', totpush, 'sec'
            endif
            if (movion==1) then
               call HARTBEAT(msg,7)
               tpushi = msg(1); tdposti = msg(2); tsorti = msg(3)
               tmovi = msg(4:5); tjmovi = msg(6:7)
               if (id0==0) then
                  totpushi = tpushi + tdposti + tdjposti
                  write (18,*) 'ion push time=', tpushi, 'sec'
                  write (18,*) 'ion charge deposit time=', tdposti, 'sec&
     &'
                  write (18,*) 'ion current deposit time = ', tdjposti, &
     &'sec'
                  write (18,*) 'total ion push time=', totpushi, 'sec'
                  write (18,*) 'ion sort time=', tsorti
                  tmovi = tmovi + tjmovi
                  write (18,*) 'ion move time=', tmovi
                  totpushi = totpushi + tsorti + tmovi(1)
                  write (18,*) 'total ion time=', totpushi, 'sec'
               endif
            endif
            if (id0==0) then
               write (18,*) 'total fft time=', tfft, 'sec'
               time(1) = time(1) - (totpush + totpushi + tfft(1))
               write (18,*) 'other time=', time(1), 'sec'
! write final diagnostic metafile
               close(unit=19)
               fname = 'pdiag3.final.'//cdrun
               open(unit=19,file=trim(fname),form='formatted',status='re&
     &place')
! potential diagnostics
               if (ntp > 0) then
                  nprec = nprec - 1
                  write (19,ppot3d,iostat=irc)
               endif
! done
               write (18,*) '* * * q.e.d. * * *'
            endif
            call MP_END
            call PPEXIT
            stop
         else
            itime = it
         endif
         if (id0==0) write (18,991) itime
         write (label,991) itime
         call LOGNAME(label)
! velocity diagnostic
         if (ntv > 0) then
            it = itime/ntv
            if (itime==ntv*it) then
! print out velocity moments
               if (id0==0) write (10,*) it, fvm(1,:,:), fvm(2,:,:)
               if (movion==1) then
! print out velocity moments
                  if (id0==0) write (20,*) it, fvmi(1,:,:), fvmi(2,:,:)
              endif
            endif
         endif
! potential diagnostic
         if (ntp > 0) then
            it = itime/ntp
            if (itime==ntp*it) then
! write diagnostic output
!              if (nprec==0) then
!                 nprec = -1
!                 fname = 'ppotk3.'//cdrun
!                 call writef(pott,modesx,modesy,modesz2,kxyp,kyzp,jblok&
!    &,11,nprec,trim(fname))
!              else
!                 call writef(pott,modesx,modesy,modesz2,kxyp,kyzp,jblok&
!    &,11,nprec)
!              endif
            endif
         endif
! vector potential diagnostic
         if (nta > 0) then
            it = itime/nta
            if (itime==nta*it) then
               continue
            endif
         endif
! energy diagnostic
         if (ntw > 0) then
            it = itime/ntw
            if (itime==ntw*it) then
! get energy values
               call HARTBEAT(msg,7)
               wtot = msg(1:7)
               if (id0==0) write (18,992) wtot(1), wtot(2), wtot(4)
               if (id0==0) write (18,993) wtot(5), wtot(6), wtot(7)
            endif
         endif
! restart file
         itime = itime + 1
         if (ntr > 0) then
            it = itime/ntr
            if (itime==ntr*it) then
               it = 16 + mod(it-1,2)
               if (id0==0) write (it) itime
               call wrdata(part,npp,it)
               call wrdata(exyz,nvp,it)
               call wrdata(bxyz,nvp,it)
               if (movion==1) call wrdata(part,npp,it)
               if (movion==0) call wrdata(qi,nvp,it)
               if (ntw > 0) call wrdata(wt,1,it)
               if ((ntp > 0) .and. (id0==0)) write (it) nprec
               if (id0==0) then
                  write (it) itime
                  end file it
                  rewind it
               endif
            endif
         endif
         go to 10
         end subroutine
!
      end program pbbeps32

