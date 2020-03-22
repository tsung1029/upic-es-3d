!-----------------------------------------------------------------------
! * * * periodic 3d electrostatic particle simulation kernel code * * *
! this is a simple 3d skeleton particle-in-cell code designed for
! exploring new computer architectures.  it contains the critical pieces
! needed for depositing charge, advancing particles, and solving the
! field.  the code moves only electrons, with periodic electrostatic
! forces obtained by solving poisson's equation with fast fourier
! transforms.  the only diagnostic is particle and field energy.
! portable gcpic kernel code, using algorithm described in:
! p. c. liewer and v. k. decyk, j. computational phys. 85, 302 (1989).
! written by viktor k. decyk, ucla
! for mpi distributed memory computers, with 2D domain decomposition
! update: march 10, 2008
      program init_pbeps32
      use pinit32d
!     use pespush32d
!     use pfield32d
!     use pdiag32d
!     use psimul32d
      use p32d, only : dcomp, fcomp
!     use mp0d, only: mpinit, ncpus
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
! default unit numbers
      integer :: iuin = 8, iuot = 18
!     integer :: iuin = 8, iuot = 18, iudm = 19, iud = 12, iuv = 10
!     integer :: iuvi = 20, iup = 11, iuer = 2
!     integer :: npxyz, npxyzb, np, npxyzi, npxyzbi, npi
      double precision :: npxyz, npxyzb, np, npxyzi, npxyzbi, npi
      integer :: nx, ny, nz, nxh, nyh, nzh, nyv, nzv, nxe, nxeh
      integer :: nloop, nvpy, nvpz, nvp, iblok, nblok, inblok
      integer :: npmax, npimax = 0, kyp, kzp, nypmx, nzpmx, kxyp, kyzp
      integer :: kyb, kxb, kzb, kyzb, kxyb, kzyb, kbmin, kblok
      integer :: jbmin, jblok, lbmin, lblok, mbmin, mblok
      integer :: ngds, nxyzh, nxhyz, nx1, nyzpm1, nbmax
      integer :: idproc, id0, kstrt, itime, itime0, ntime, isign, irc
      integer :: j, ierr
!     integer :: j, it, itw, iur1, iur2, ierr, modesz2
!     integer :: ntasks
      real :: zero = 0.0, wki = 0.0
      real :: tpush = 0.0, tdpost = 0.0, tsort = 0.0
      real :: tpushi = 0.0, tdposti = 0.0, tsorti = 0.0
      real :: totpush = 0.0, totpushi = 0.0
      real :: qbme, qbmi, affp
!     real :: qbme, qbmi, affp, qi0, etx, we, wke
      real :: vtxi, vtyi, vtzi, vtdxi, vtdyi, vtdzi
!     double precision :: dtime
      real, dimension(:,:,:), pointer :: part, parti
!     real, dimension(:,:,:), pointer :: part, part2, parti, parti2
!     real, dimension(:,:,:,:), pointer :: qe, qi
!     real, dimension(:,:,:,:,:), pointer :: fxyze
!     real, dimension(:,:,:,:,:), pointer :: bxyze
!     complex, dimension(:,:,:,:), pointer :: qt, qs
!     complex, dimension(:,:,:,:,:), pointer :: fxyzt, fxyzs
!     complex, dimension(:,:,:,:), pointer :: ffc
!     integer, dimension(:), pointer :: mixup
!     complex, dimension(:), pointer :: sct
      real, dimension(:,:), pointer :: edges
      integer, dimension(:,:), pointer :: nyzp, noff
      integer, dimension(:), pointer :: npp, nppi, nps
!     real, dimension(:,:), pointer :: pt
!     integer, dimension(:,:), pointer :: ip, npic
!     real, dimension(:,:,:,:), pointer :: sfield
!     complex, dimension(:,:,:,:), pointer :: sfieldt, pott
!     real, dimension(:,:,:), pointer :: fv, fvm, fvi, fvmi
!     real, dimension(:,:), pointer :: wt
! wtot = total energy
!     real, dimension(4) :: wtot
! time = timing array
!     real, dimension(2) :: tfft = 0.0, time = 0.0
!     real, dimension(2) :: tmove = 0.0, tmovi = 0.0
! msg = heartbeat array
!     double precision, dimension(9) :: msg
      character(len=10) :: cdrun
      character(len=32) :: fname
!     character(len=12) :: label
      integer, external :: NDIAN, NDPREC, IDPREC
  991 format (' T = ',i7)
  992 format (' field, kinetic, total energies = ',3e14.7)
! get unit number for error file
!     iuer = get_funit(iuer)
! nvp = number of real or virtual processors
! initialize for parallel processing
!     call PPINIT(idproc,id0,nvp)
      idproc = 0; id0 = 0; nvp = 2048
!     kstrt = idproc + 1
      kstrt = nvp
! read namelist
      if (id0==0) then
!        iuin = get_funit(iuin)
         open(unit=iuin,file='pinput3',form='formatted',status='old')
         read (iuin,pinput3)
      endif
! broadcast namelist to other nodes
!     call sendnml()
! override input data
      psolve = 1
! set monitor flag
!     call SET_MON(mpimon)
! create string from idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
! text output file
      if (id0==0) then
!        iuot = get_funit(iuot)
         fname = 'poutput3.'//cdrun
         open(unit=iuot,file=trim(fname),form='formatted',status='replac&
     &e')
      endif
! np = total number of electrons in simulation
!     npxyz = npx*npy*npz; npxyzb = npxb*npyb*npzb; np = npxyz + npxyzb
      npxyz = dble(npx)*dble(npy)*dble(npz)
      npxyzb = dble(npxb)*dble(npyb)*dble(npzb)
      np = npxyz + npxyzb
! npi = total number of ions in simulation
!     npxyzi = npxi*npyi*npzi; npxyzbi = npxbi*npybi*npzbi
!     npi = npxyzi + npxyzbi
      npxyzi = dble(npxi)*dble(npyi)*dble(npzi)
      npxyzbi = dble(npxbi)*dble(npybi)*dble(npzbi)
      npi = npxyzi + npxyzbi
      nx = 2**indx; ny = 2**indy; nz = 2**indz
      nxh = nx/2; nyh = ny/2; nzh = nz/2
      nyv = ny + 2; nzv = nz + 2; nxe = nx + 4
! nvpy/nvpz = number of real or virtual processors in y/z
      call fcomp(nvp,nx,ny,nz,nvpy,nvpz,ierr)
      nvp = nvpy*nvpz
! kyp = number of complex grids in each field partition in y direction
      kyp = (ny - 1)/nvpy + 1
! kzp = number of complex grids in each field partition in z direction
      kzp = (nz - 1)/nvpz + 1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
      nypmx = kyp + 3; nzpmx = kzp + 3
! ngds = number of guard cells
      ngds = 3*((idds - 1)/2 + 1)
!     ax = .866025; ay = .866025; az = .866025
      if (inorder==LINEAR) then
         ax = .912871; ay = .912871; az = .912871
         nxe = nx + 2; nypmx = kyp + 1; nzpmx = kzp + 1
         ngds = (idds - 1)/2 + 1
      endif
      nxeh = nxe/2
! initialize for multiprocessing
!     ntasks = mpinit(sntasks)
      if (popt==VECTOR) vect = 1
! nloop = number of time steps in simulation
      nloop = tend/dt + .0001
! iblok/nblok = number of particle partitions in y/z
      iblok = 1 + mshare*(nvpy - 1); nblok = 1 + mshare*(nvpz - 1)
      inblok = iblok*nblok
! npmax = maximum number of particles in each partition
      npmax = (np/dble(nvp))*1.25
      if (movion==1) npimax = (npi/dble(nvp))*1.25
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
! debug
      write (iuot,*) 'np,npmax,nbmax=',np,npmax,nbmax
! end debug
      nbmax = 2*nbmax
      if (movion==1) then
         vtxi = vtx/sqrt(rmass*rtempxi)
         vtyi = vty/sqrt(rmass*rtempyi)
         vtzi = vtz/sqrt(rmass*rtempzi)
      endif
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
! where kk = k + noff(1,m) - 1, ll = l + noff(2,m) - 1
!     allocate(qe(nxe,nypmx,nzpmx*kbmin,kblok*lblok))
! in real space, fxyze(i,j+1,k,l) = i component of force/charge at 
! grid point (j,kk,ll)
! in other words, fxyze are the convolutions of the electric field
! over the particle shape, where kk = k + noff(1,m) - 1, and
! ll = l + noff(2,m) - 1
!     allocate(fxyze(3,nxe,nypmx,nzpmx*kbmin,kblok*lblok))
! bxyze(i,j+1,k,l) = i component of magnetic field at 
! grid point (j,kk,ll), that is, bxyze is the convolution of the
! magnetic field over the particle shape, where kk = k + noff(1,m) - 1,
! and ll = l + noff(2,m) - 1
!     allocate(bxyze(3,nxe,nypmx,nzpmx*kbmin,kblok*lblok))
! qt(l,k,j,m) = complex charge density for fourier mode jj-1,kk-1,l-1
! fxyzt(1,l,k,j,m) = x component of force/charge
! fxyzt(2,l,k,j,m) = y component of force/charge
! fxyzt(3,l,k,j,m) = z component of force/charge
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*(my - 1), and
! kk = k + kyzp*(mz - 1), m = my + nvpy*(mz - 1)
!     allocate(qt(nzv,kxyp,kyzp*mbmin,jblok*mblok))
!     allocate(fxyzt(3,nzv,kxyp,kyzp*mbmin,jblok*mblok))
!     allocate(qs(nyv,kxyp,nzpmx*jbmin*lbmin,jblok*lblok))
!     allocate(fxyzs(3,nyv,kxyp,nzpmx*jbmin*lbmin,jblok*lblok))
! ffc = form factor array for poisson solver
!     allocate(ffc(nzh,kxyp,kyzp,jblok*mblok))
! mixup, sct = arrays for fft
!     allocate(mixup(nxhyz),sct(nxyzh))
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
!     allocate(pt(max(npmax,npimax),inblok))
!     allocate(ip(max(npmax,npimax),inblok),npic(nyzpm1,inblok))
!     if (sortime > 0) then
!        allocate(part2(idimp,npmax,inblok))
!     else
!        allocate(part2(0,0,0))
!     endif
! initialize parallel timer
!     call pwtimer(time,dtime,-1)
! initialize constants
      itime0 = 0
      itime = itime0
      ntime = itime + itime0
      qbme = qme
      affp = dble(nx)*dble(ny)*dble(nz)/np
      if (movion==1) then
         qbmi = qmi/rmass
         vtdxi = vtx/sqrt(rmass*rtempdxi)
         vtdyi = vty/sqrt(rmass*rtempdyi)
         vtdzi = vtz/sqrt(rmass*rtempdzi)
      endif
! set initial time
      t0 = dt*real(itime0)
! determine number format and default precisions
      indian = NDIAN()
      rlprec = NDPREC()
      inprec = IDPREC()
! calculate partition variables
      call dcomp(edges,nyzp,noff,ny,nz,kstrt,nvpy,nvpz,iblok,inorder)
!     bxyze(1,:,:,:,:) = omx
!     bxyze(2,:,:,:,:) = omy
!     bxyze(3,:,:,:,:) = omz
!     call cguard(bxyze,nyzp,nx,inorder)
! prepare fft tables
!     call fft_init(mixup,sct,indx,indy,indz)
! calculate form factors
!     call pois_init(ffc,ax,ay,az,affp,nx,ny,nz,kstrt,jblok)
! allocate background charge density
!     if (movion==0) allocate(qi(nxe,nypmx,nzpmx*kbmin,kblok*lblok))
! allocate ion data
      if (movion==1) then
         allocate(parti(idimp,npimax,inblok),nppi(inblok))
!        allocate(parti2(0,0,0))
      endif
! new start
      if (nustrt==1) then
! initialize electrons
         nps = 1
         npp = 0
! background electrons
!        if (npxyz > 0) call distr(part,edges,npp,nps,vtx,vty,vtz,vx0,vy&
!    &0,vz0,npx,npy,npz,nx,ny,nz,ipbc,iblok)
!        if (npxyz > 0) call pvdistr(part,npp,nps,vtx,vty,vtz,vx0,vy0,vz&
!    &0,npx,npy,npz,nx,ny,nz,ipbc,kstrt,nvp,iblok)
         if (npxyz > 0) then
            call fdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,shi&
     &ftdy,ampdz,scaledz,shiftdz,npx,npy,npz,nx,ny,nz,kstrt,nvp,ipbc,ndp&
     &rof,nsrand,iblok)
!           call vdistr(part,npp,nps,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,npz&
!    &,kstrt,nvp,iblok)
            call vvdistr(part,npp,nps,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,np&
     &z,kstrt,nvp,iblok)
         endif
! beam electrons
         nps = npp + 1
!        if (npxyzb > 0) call distr(part,edges,npp,nps,vtdx,vtdy,vtdz,vd&
!    &x,vdy,vdz,npxb,npyb,npzb,nx,ny,nz,ipbc,iblok)
!        if (npxyzb > 0) call pvdistr(part,npp,nps,vtdx,vtdy,vtdz,vdx,vd&
!    &y,vdz,npxb,npyb,npzb,nx,ny,nz,ipbc,kstrt,nvp,iblok)
         if (npxyzb > 0) then
            call fdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,shi&
     &ftdy,ampdz,scaledz,shiftdz,npxb,npyb,npzb,nx,ny,nz,kstrt,nvp,ipbc,&
     &ndprof,nsrand,iblok)
!           call vdistr(part,npp,nps,vtdx,vtdy,vtdz,vdx,vdy,vdz,npxb,npy&
!    &b,npzb,kstrt,nvp,iblok)
            call vvdistr(part,npp,nps,vtdx,vtdy,vtdz,vdx,vdy,vdz,npxb,np&
     &yb,npzb,kstrt,nvp,iblok)
         endif
! calculate actual coordinates from guiding centers
!        if (relativity==1) then
!           call distr(part,bxyze,npp,noff,qbme,ci,nx,ny,nz,ipbc,iblok,i&
!    &norder)
!        else
!           call distr(part,bxyze,npp,noff,qbme,nx,ny,nz,ipbc,iblok,inor&
!    &der)
!        endif
! debug
      write (iuot,*) 'done electron initialization, kstrt = ', kstrt
! end debug
! move electrons into appropriate spatial regions
!        call pmove(part,edges,npp,tmove,ny,nz,kstrt,nvpy,nvpz,nbmax,idd&
!    &s,iblok,vect,ierr)
!        if (ierr /= 0) then
!           call MP_END
!           call PPEXIT
!           stop
!        endif
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
!              call vdistr(parti,nppi,nps,vtxi,vtyi,vtzi,vxi0,vyi0,vzi0,&
!    &npxi,npyi,npzi,kstrt,nvp,iblok)
               call vvdistr(parti,nppi,nps,vtxi,vtyi,vtzi,vxi0,vyi0,vzi0&
     &,npxi,npyi,npzi,kstrt,nvp,iblok)
            endif
! beam ions
            nps = nppi + 1
!           if (npxyzbi > 0) call distr(parti,edges,nppi,nps,vtdxi,vtdyi&
!    &,vtdzi,vdxi,vdyi,vdzi,npxbi,npybi,npzbi,nx,ny,nz,ipbc,iblok)
            if (npxyzbi > 0) then
               call fdistr(parti,nps,ampdxi,scaledxi,shiftdxi,ampdyi,sca&
     &ledyi,shiftdyi,ampdzi,scaledzi,shiftdzi,npxbi,npybi,npzbi,nx,ny,nz&
     &,kstrt,nvp,ipbc,ndprofi,nsrandi,iblok)
!              call vdistr(parti,nppi,nps,vtdxi,vtdyi,vtdzi,vdxi,vdyi,vd&
!    &zi,npxbi,npybi,npzbi,kstrt,nvp,iblok)
               call vvdistr(parti,nppi,nps,vtdxi,vtdyi,vtdzi,vdxi,vdyi,v&
     &dzi,npxbi,npybi,npzbi,kstrt,nvp,iblok)
            endif
! calculate actual coordinates from guiding centers
!           if (relativity==1) then
!              call distr(parti,bxyze,nppi,noff,qbmi,ci,nx,ny,nz,ipbc,ib&
!    &lok,inorder)
!           else
!              call distr(parti,bxyze,nppi,noff,qbmi,nx,ny,nz,ipbc,iblok&
!    &,inorder)
!           endif
! debug
!     write (iuot,*) 'done ion initialization, kstrt = ', kstrt
! end debug
! move ions into appropriate spatial regions
!           call pmove(parti,edges,nppi,tmovi,ny,nz,kstrt,nvpy,nvpz,nbma&
!    &x,idds,iblok,vect,ierr)
!           if (ierr /= 0) then
!              call MP_END
!              call PPEXIT
!              stop
!           endif
         endif
! initialize background charge density
!        if (movion==0) then
!           qi0 = -qme/affp
!           call sguard(qi,nyzp,zero,nx,inorder)
!           call dpost(part,qi,-qme,npp,noff,tdpost,inorder,dopt)
! debug
!           call sguard(qi,nyzp,qi0,nx,inorder)
! freeze the ions now
!        else if ((movion==1).and.(ntime==ionoff)) then
!           allocate(qi(nxe,nypmx,nzpmx*kbmin,kblok*lblok))
! initialize ion charge density to zero
!           call sguard(qi,nyzp,zero,nx,inorder)
! deposit ion charge
!           call dpost(parti,qi,qmi,nppi,noff,tdposti,inorder,dopt)
! delete ions
!           deallocate(parti,nppi)
!           movion = 0
!        endif
      endif
! debug
      do j = 1, npp(1)
         write (62,*) j,part(:,j,1)
      enddo
! end debug
      write (iuot,*) '* * * q.e.d. * * *'
      call MP_END
      call PPEXIT
      stop
!
      end program init_pbeps32

