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
! Fortran 90 for Macintosh G3
! copyright 1999, regents of the university of california
! update: may 2, 2008
      program pbeps32
      use pinit32d
      use pespush32d
      use pfield32d
      use pdiag32d
      use psimul32d
      use ext_driver32_jf
      use pinit32d_jf
      use diag32_jf
      use hdf_write32_jf
!      use p32d_jf
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
      integer :: npimax = 0
! default unit numbers
      integer :: iuin = 8, iuot = 18, iudm = 19, iud = 12, iuv = 10
      integer :: iuvi = 20, iup = 11, ium = 21, iuer = 2
!     integer :: npxyz, npxyzb, np, npxyzi, npxyzbi, npi
      double precision :: npxyz, npxyzb, np, npxyzi, npxyzbi, npi
      integer :: nx, ny, nz, nxh, nyh, nzh, nyv, nzv, nxe, nxeh
      integer :: nloop, nvpy, nvpz, nvp, iblok, nblok, inblok
      integer :: npmax, kyp, kzp, nypmx, nzpmx, kxyp, kyzp
      integer :: kyb, kxb, kzb, kyzb, kxyb, kzyb, kbmin, kblok
      integer :: jbmin, jblok, lbmin, lblok, mbmin, mblok
      integer :: ngds, nxyzh, nxhyz, nx1, nyzpm1, nbmax
      integer :: idproc, id0, kstrt, itime, itime0, ntime, isign, irc
      integer :: j, it, itw, iur1, iur2, ierr, modesz2
      integer :: ntasks
      real :: zero = 0.0, wki = 0.0, ltime = 0.0, tloop = 0.0
      real :: tpush = 0.0, tdpost = 0.0, tsort = 0.0
      real :: tpushi = 0.0, tdposti = 0.0, tsorti = 0.0
      real :: totpush = 0.0, totpushi = 0.0
      real :: qbme, qbmi, affp, qi0, etx, we, wke
      real :: sx = 0.0, sy = 0.0, sz = 0.0
      real :: pxi = 0.0, pyi = 0.0, pzi = 0.0
      real :: pxe, pye, pze, wx, wy, wz
      real :: vtxi, vtyi, vtzi, vtdxi, vtdyi, vtdzi
      double precision :: dtime, ldtime
      real, dimension(:,:,:), pointer :: part, part2, parti, parti2
      real, dimension(:,:,:,:), pointer :: qe, qi
      real, dimension(:,:,:,:,:), pointer :: fxyze
      
      real, dimension(:,:), allocatable :: U_sumover_x
      real, dimension(:,:,:,:,:), pointer :: grad_phi
			complex,dimension(:,:,:,:,:),pointer :: grad_phis, grad_phit !For field diagnostic			
! for line diagnostics
			real,dimension(:,:), pointer :: lines
			integer :: i
      
!     real, dimension(:,:,:,:,:), pointer :: bxyze
      complex, dimension(:,:,:,:), pointer :: qt, qs
      complex, dimension(:,:,:,:,:), pointer :: fxyzt, fxyzs
      complex, dimension(:,:,:,:), pointer :: ffc
      integer, dimension(:), pointer :: mixup
      complex, dimension(:), pointer :: sct
      real, dimension(:,:), pointer :: edges
      integer, dimension(:,:), pointer :: nyzp, noff
      integer, dimension(:), pointer :: npp, nppi, nps
      real, dimension(:,:), pointer :: pt
      integer, dimension(:,:), pointer :: ip, npic
      real, dimension(:,:,:,:), pointer :: sfield
      complex, dimension(:,:,:,:), pointer :: sfieldt, dent, pott
      real, dimension(:,:,:), pointer :: fv, fvm, fvi, fvmi
      real, dimension(:,:), pointer :: wt
! wtot = total energy
      real, dimension(4) :: wtot
! time = timing array
      real, dimension(2) :: tfft = 0.0, time = 0.0
      real, dimension(2) :: tmove = 0.0, tmovi = 0.0
! msg = heartbeat array
      double precision, dimension(10) :: msg
      character(len=10) :: cdrun
      character(len=32) :: fname
      character(len=12) :: label
      integer, external :: NDIAN, NDPREC, IDPREC
  991 format (' T = ',i7)
  992 format (' field, kinetic, total energies = ',3e14.7)
  996 format (' total momentum = ',3e14.7)
! get unit number for error file
      iuer = get_funit(iuer)
! nvp = number of real or virtual processors
! initialize for parallel processing
      call PPINIT(idproc,id0,nvp)
      kstrt = idproc + 1
! read namelist
      if (id0==0) then
         iuin = get_funit(iuin)
         open(unit=iuin,file='pinput3',form='formatted',status='old')
         read (iuin,pinput3)
         read (iuin,pinput3_jf)
      endif
! broadcast namelist to other nodes
      call sendnml()
      call sendnml_jf()
! override input data
      psolve = 1
! set monitor flag
      call SET_MON(mpimon)
! create string from idrun
      write (cdrun,'(i10)') idrun
      cdrun = adjustl(cdrun)
! text output file
      if (id0==0) then
         iuot = get_funit(iuot)
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
      ntasks = mpinit(sntasks)
      if (popt==VECTOR) vect = 1
! nloop = number of time steps in simulation
      nloop = tend/dt + .0001
! iblok/nblok = number of particle partitions in y/z
      iblok = 1 + mshare*(nvpy - 1); nblok = 1 + mshare*(nvpz - 1)
      inblok = iblok*nblok
! npmax = maximum number of particles in each partition
!      npmax = (np/dble(nvp))*1.75
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
      nbmax = 2*nbmax
      if (movion==1) then
         vtxi = vtx/sqrt(rmass*rtempxi)
         vtyi = vty/sqrt(rmass*rtempyi)
         vtzi = vtz/sqrt(rmass*rtempzi)
      endif
! initialize time constants
      itime0 = 0
      itime = itime0
      ntime = itime + itime0
! diagnostic information needed by diagnostic nodes
! set default diagnostic file names
      if (ntd > 0) fdname = 'denk32.'//cdrun
      if (ntp > 0) fpname = 'potk32.'//cdrun
! energy time history
      if (ndw > 0) then
         allocate(wt((nloop-1)/ndw-(itime0/ndw)+1,4))
         itw = 0
      endif
      call mkdir_structure(idproc)
! open restart files
      if (id0==0) then
         if (nustrt==0) then
            call restart_open(nustrt,ntr,idrun0,iur1,iur2,iuer)
         else
            call restart_open(nustrt,ntr,idrun,iur1,iur2,iuer)
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
! where kk = k + noff(1,m) - 1, ll = l + noff(2,m) - 1
      allocate(qe(nxe,nypmx,nzpmx*kbmin,kblok*lblok))
! in real space, qi(j+1,k,l,m) = ion charge charge density
      allocate(qi(nxe,nypmx,nzpmx*kbmin,kblok*lblok))
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
!     allocate(bxyze(3,nxe,nypmx,nzpmx*kbmin,kblok*lblok))
! qt(l,k,j,m) = complex charge density for fourier mode jj-1,kk-1,l-1
! fxyzt(1,l,k,j,m) = x component of force/charge
! fxyzt(2,l,k,j,m) = y component of force/charge
! fxyzt(3,l,k,j,m) = z component of force/charge
! all for fourier mode jj-1,kk-1,l-1, where jj = j + kxyp*(my - 1), and
! kk = k + kyzp*(mz - 1), m = my + nvpy*(mz - 1)
      allocate(qt(nzv,kxyp,kyzp*mbmin,jblok*mblok))
      allocate(fxyzt(3,nzv,kxyp,kyzp*mbmin,jblok*mblok))
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
! initialize parallel timer
      call pwtimer(time,dtime,-1)
! initialize constants
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
      call fft_init(mixup,sct,indx,indy,indz)
! calculate form factors
      call pois_init(ffc,ax,ay,az,affp,nx,ny,nz,kstrt,jblok)
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
!            call fdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,shi&
!     &ftdy,ampdz,scaledz,shiftdz,npx,npy,npz,nx,ny,nz,kstrt,nvp,ipbc,ndp&
!     &rof,nsrand,iblok)
!            call vdistr(part,npp,nps,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,npz&
!     &,kstrt,nvp,iblok)
!					 call uniform_density_init(part,nx,ny,nz,npx,npy,npz,nvpy,nvpz,idproc)
!print*,"before init"
!				call uniform_random_init(part,npp,nps,nx,ny,nz,npx,npy,npz,&
!				&	kstrt,nvp,mblok,nvpy,nvpz)
           call vfdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,sh&
    &iftdy,ampdz,scaledz,shiftdz,npx,npy,npz,nx,ny,nz,kstrt,nvp,ipbc,nd&
    &prof,nsrand,iblok)
           call vvdistr(part,npp,nps,vtx,vty,vtz,vx0,vy0,vz0,npx,npy,np&
    &z,kstrt,nvp,iblok)
!print*,"after init"
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
!           call vfdistr(part,nps,ampdx,scaledx,shiftdx,ampdy,scaledy,sh&
!    &iftdy,ampdz,scaledz,shiftdz,npxb,npyb,npzb,nx,ny,nz,kstrt,nvp,ipbc&
!    &,ndprof,nsrand,iblok)
!           call vvdistr(part,npp,nps,vtdx,vtdy,vtdz,vdx,vdy,vdz,npxb,np&
!    &yb,npzb,kstrt,nvp,iblok)
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
!print*,"before pmove"
         nbmax = nbmax / 8
         call pmove(part,edges,npp,tmove,ny,nz,kstrt,nvpy,nvpz,nbmax,idd&
     &s,iblok,vect,ierr)
!         nbmax = nbmax * 8
!print*,"after pmove"
! calculate initial electron momentum
         if (ntm > 0) call initmomt3(part,npp,pxe,pye,pze)
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
!              call vfdistr(parti,nps,ampdxi,scaledxi,shiftdxi,ampdyi,sc&
!    &aledyi,shiftdyi,ampdzi,scaledzi,shiftdzi,npxi,npyi,npzi,nx,ny,nz,k&
!    &strt,nvp,ipbc,ndprofi,nsrandi,iblok)
!              call vvdistr(parti,nppi,nps,vtxi,vtyi,vtzi,vxi0,vyi0,vzi0&
!    &,npxi,npyi,npzi,kstrt,nvp,iblok)
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
!              call vfdistr(parti,nps,ampdxi,scaledxi,shiftdxi,ampdyi,sc&
!    &aledyi,shiftdyi,ampdzi,scaledzi,shiftdzi,npxbi,npybi,npzbi,nx,ny,n&
!    &z,kstrt,nvp,ipbc,ndprofi,nsrandi,iblok)
!              call vvdistr(parti,nppi,nps,vtdxi,vtdyi,vtdzi,vdxi,vdyi,v&
!    &dzi,npxbi,npybi,npzbi,kstrt,nvp,iblok)
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
! calculate initial ion momentum
            if (ntm > 0) call initmomt3(parti,nppi,pxi,pyi,pzi)
         endif
! initialize background charge density

         if (movion==0) then
            qi0 = -qme/affp
            call sguard(qi,nyzp,zero,nx,inorder)
            call dpost(part,qi,-qme,npp,noff,tdpost,inorder,dopt)
            call aguard(qi,nyzp,nx,inorder)
            call paguard(qi,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,iblok,inorde&
     &r)
! debug
!           call sguard(qi,nyzp,qi0,nx,inorder)
! freeze the ions now
         else if ((movion==1).and.(ntime==ionoff)) then
! initialize ion charge density to zero
            call sguard(qi,nyzp,zero,nx,inorder)
! deposit ion charge
            call dpost(parti,qi,qmi,nppi,noff,tdposti,inorder,dopt)
! add guard cells for density in x
            call aguard(qi,nyzp,nx,inorder)
! add guard cells for density in y and z
            call paguard(qi,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,iblok,inorde&
     &r)
! delete ions
            deallocate(parti,nppi)
            movion = 0
         endif
! restart
      else
! read restart files
         it = 0
         call restart_bread(iur1,iur2,id0,it,itime,itime0,nvp,npp,part, &
     &movion,nppi,parti,qi,irc,iuer)
         if (irc /= 0) go to 400
! extend run
         if (nustrt==0) then
            itime0 = itime + itime0
            t0 = dt*real(itime0)
            itime = 0
            ntime = itime + itime0
            if (id0==0) then
               if (iur1 >= 0) close (unit=iur1)
               if (iur2 >= 0) close (unit=iur2)
               call restart_open(1,ntr,idrun,iur1,iur2,iuer)
            endif
            go to 490
         endif
! read diagnostics
         call restart_dread(it,id0,itime,itw,wt,iud,ndrec,fdname,iup,   &
     &nprec,fpname,irc,iuer)
         if (irc /= 0) go to 400
         ntime = itime + itime0
         t0 = dt*real(itime0)
         if (id0==0) rewind it
         go to 490
! handle error'
  400    if (id0==0) then
            write (iuer,*) 'Restart Error, irc = ', irc
         endif
         go to 3000
      endif
!
! sorting arrays
  490 allocate(pt(max(npmax,npimax),inblok))
      allocate(ip(max(npmax,npimax),inblok),npic(nyzpm1,inblok))
!Commented by JF to reduce memory usage from sort
!      if (sortime > 0) then
!         allocate(part2(idimp,npmax,inblok))
!      else
!         allocate(part2(0,0,0))
!      endif
! reduce size of particle manager buffers
!     nbmax = nbmax/2
! initialize diagnostics
! open initial diagnostic metafile
      if (id0==0) then
         iudm = get_funit(iudm)
         fname = 'pdiag32.init.'//cdrun
         open(unit=iudm,file=trim(fname),form='formatted',status='replac&
     &e')
      endif
      
      call write_jf_diag_file(idproc,cdrun)

      if (ntlines .ne. 0) then
				do i = 1, 5
					if (linepos_3d(i,1) .ne. -1) then
						nlines = nlines + 1
					endif
				enddo
      	allocate(nlinerec(nlines))
      	allocate(lines(nlines,nx))
      	lines=0.
      	nlinerec=0.
      endif
      
      if (nt_write_U_sumover_x .ne. 0) then
      	allocate(grad_phi(3,nxe,nypmx,nzpmx,1))
	      allocate(grad_phis(3,nyv,kxyp,nzpmx*jbmin*lbmin,jblok*lblok))
      	allocate(grad_phit(3,nzv,kxyp,kyzp*mbmin,jblok*mblok))
      	allocate(U_sumover_x(ny,nz))
      endif
			if ((nt_write_U_sumover_x .ne. 0)) then
				if (idproc==0) then
					fname = './DIAG/U_sumover_x_file'
					nt_write_U_sumover_x_funit = get_funit(20)
					open(unit=nt_write_U_sumover_x_funit,file=trim(fname),form='unformatted',status='replace')
					write (nt_write_U_sumover_x_funit) ny,nz,nt_write_U_sumover_x,dt,(timerise+timeflat+timefall),tend
				endif
			endif
! ion density or potential diagnostics
      if ((ntp > 0) .or. (ndp > 0) .or. (ntd > 0) .or. (ndd > 0) .or. (ntfield > 0)) then
         allocate(sfield(nxe,nypmx,nzpmx*kbmin,kblok*lblok))
         allocate(sfieldt(nzv,kxyp,kyzp*mbmin,jblok*mblok))
      endif
! ion density diagnostics
      call initmodediag(dent,ntd,id0,nxh,nyh,nzh,kxyp,kyzp,modesxd,     &
     &modesyd,modeszd,jblok,mblok,iud,ndrec,fdname)
      if (ntd > 0) then
         if (id0==0) then
            ceng = zero
            write (iudm,pden32d,iostat=irc)
         endif
      endif
! velocity diagnostics
      fname = 'fv3.'//cdrun
      call initveldiag(fv,fvm,vtx,vty,vtz,ntv,ndv,id0,nmv,nblok,iuv,fnam&
     &e)
      if (movion==1) then
         fname = 'fvi3.'//cdrun
         call initveldiag(fvi,fvmi,vtxi,vtyi,vtzi,ntv,ndv,id0,nmv,nblok,&
     &iuvi,fname)
      endif
! potential diagnostics
      call initmodediag(pott,ntp,id0,nxh,nyh,nzh,kxyp,kyzp,modesxp,     &
     &modesyp,modeszp,jblok,mblok,iup,nprec,fpname)
      if (ntp > 0) then
         if (id0==0) then
            ceng = zero
            write (iudm,ppot32d,iostat=irc)
         endif
      endif
! momentum diagnostics
      fname = 'pmomentum3.'//cdrun
      if (ntm > 0) then
         if (id0==0) then
            ium = get_funit(ium)
            open(unit=ium,file=trim(fname),form='formatted',status='unkn&
     &own')
         endif
      endif
! write out and close input file
      if (id0==0) then
         write (iudm,pinput3,iostat=irc)
         close(unit=iudm)
      endif
! record time
      call pwtimer(time,dtime)
! send initial CPU Time to diagnostic nodes
      msg(1) = time(1); msg(2) = time(2)
      call HARTBEAT(msg,2)
      call wtimer(ltime,ldtime,-1)
      if (id0==0) then
         write (iuot,*) 'init max/min real time=',time(1),time(2), 'sec'
      endif
!
! * * * start main iteration loop * * *
!
  500 if (nloop <= ntime) go to 2000
! send time step to diagnostic nodes
      msg(1) = ntime
      call HARTBEAT(msg,1)
      if (id0==0) write (iuot,991) ntime
      write (label,991) ntime
      call LOGNAME(label)
! deposit electron charge density
! initialize electron charge density to zero
      call sguard(qe,nyzp,zero,nx,inorder)
! deposit electron charge
      call dpost(part,qe,qme,npp,noff,tdpost,inorder,dopt)
! add guard cells for density in x
      call aguard(qe,nyzp,nx,inorder)
! add guard cells for density in y and z
      call paguard(qe,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,iblok,inorder)
! deposit ion charge density
      if (movion==1) then
! initialize ion charge density to zero
         call sguard(qi,nyzp,zero,nx,inorder)
! deposit ion charge
         call dpost(parti,qi,qmi,nppi,noff,tdposti,inorder,dopt)
! add guard cells for density in x
         call aguard(qi,nyzp,nx,inorder)
! add guard cells for density in y and z
         call paguard(qi,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,iblok,inorder)
! freeze the ions
         if (ntime==ionoff) then
            deallocate(parti,nppi)
            movion = 0
         endif
      endif
! add electron and ion densities
      call addqei(qe,qi,nyzp,nx,inorder)
! ion density diagnostic
      call dendiag(qt,qs,qi,sfield,dent,sfieldt,ffc,nyzp,mixup,sct,tfft,&
     &ntd,ndd,nx,ny,nz,modesxd,modesyd,modeszd,iud,ndrec,indx,indy,indz,&
     &ntime,nvpy,nvpz,kstrt,kxyp,kyp,kyzp,kzp,ngds,iblok,jblok,kblok,   &
     &mblok,inorder)
! velocity diagnostic
      call veldiag(part,fv,fvm,npp,ntv,ndv,id0,nmv,iuv,ntime)
      if (movion==1) then
         call veldiag(parti,fvi,fvmi,nppi,ntv,ndv,id0,nmv,iuvi,ntime)
      endif
! transform charge to fourier space
      isign = -1
      call fft(qe,qs,qt,isign,mixup,sct,tfft,indx,indy,indz,kstrt,kxyp,k&
     &yp,kyzp,kzp,kblok,mblok,inorder)
! potential diagnostic
      call potdiag(qt,qs,sfield,pott,sfieldt,ffc,nyzp,mixup,sct,tfft,ntp&
     &,ndp,nx,ny,nz,modesxp,modesyp,modeszp,iup,nprec,indx,indy,indz,   &
     &ntime,nvpy,nvpz,kstrt,kxyp,kyp,kyzp,kzp,ngds,iblok,jblok,kblok,   &
     &mblok,inorder)
! calculate force/charge in fourier space
      call pois(qt,fxyzt,ffc,we,nx,ny,nz,kstrt,jblok)
! transform force/charge to real space
      isign = 1
      call fft(fxyze,fxyzs,fxyzt,isign,mixup,sct,tfft,indx,indy,indz,kst&
     &rt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
!     call fftn(fxyze,fxyzs,fxyzt,isign,mixup,sct,tfft,indx,indy,indz,ks&
!    &trt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! copy data from field to particle partition, and copy to guard cells
      call pcguard(fxyze,kstrt,nvpy,nvpz,kyp,kzp,ngds,iblok,inorder)
      call cguard(fxyze,nyzp,nx,inorder)
! external pump
!     if ((itpon > 0).and.(ntime >= itpon)) then
!        etx = (v0*vtx)*w0*cos(w0*dt*(ntime - itpon))
!        fxyze(1,:,:,:,:) = fxyze(1,:,:,:,:) + etx
!     endif
! particle push and charge density update

!External driver
			select case (driver_select)
			case (0) 
				call do_nothing()
			case (1)
				call plane_wave(fxyze,real(itime)*dt,nx,nxe,nvp,idproc)
			case (7)
				call gauss_tran_per_wavelen(fxyze,real(itime)*dt,nx,nxe,ny,nypmx,nz,nzpmx,nvpy,nvpz,idproc)
			case (12)
				call supergauss_tran_per_wavelen(fxyze,real(itime)*dt,nx,nxe,ny,nypmx,nz,nzpmx,nvpy,nvpz,idproc)
			case default
				print*,"driver_select = ",driver_select," is not yet implemented.  Exiting..."
				call MP_END
				call PPEXIT
				stop
			end select

			if ((dump_start < 0.) .OR. &
				&((real(itime)*dt >= dump_start) .and. (real(itime)*dt <= dump_end))) then
			if (ntfield > 0) then
				it = itime/ntfield
				if (itime==ntfield*it) then
				
					fname = './DIAG/Ex/'//'pfieldx.'//cdrun
					sfield = fxyze(1,:,:,:,:)
					call writef(sfield,nxe,nypmx,nzpmx,nvpy,nvpz,itime,itime*dt,&
						&EX,trim(fname),inorder,field_write_stride)
					fname = './DIAG/Ey/'//'pfieldy.'//cdrun
					sfield = fxyze(2,:,:,:,:)
					call writef(sfield,nxe,nypmx,nzpmx,nvpy,nvpz,itime,itime*dt,&
						&EY,trim(fname),inorder,field_write_stride)
					fname = './DIAG/Ez/'//'pfieldz.'//cdrun
					sfield = fxyze(3,:,:,:,:)
					call writef(sfield,nxe,nypmx,nzpmx,nvpy,nvpz,itime,itime*dt,&
						&EZ,trim(fname),inorder,field_write_stride)
				endif
			endif
			endif

			if (nt_write_U_sumover_x > 0) then
				it = itime/nt_write_U_sumover_x
				if (itime==nt_write_U_sumover_x*it) then
					
		      call pois(qt,sfieldt,ffc,we,nx,ny,nz,kstrt,jblok)
					call ipgradf32(sfieldt,grad_phit,nx,ny,nz,kstrt,jblok)
					isign = 1
					call fft(grad_phi,grad_phis,grad_phit,isign,mixup,sct,tfft,indx,indy,indz,&
						&kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
					call pcguard(grad_phi,kstrt,nvpy,nvpz,kyp,kzp,ngds,iblok,inorder)
					call cguard(grad_phi,nyzp,nx,inorder)
					
					U_sumover_x = 0.
					
					call get_3D_sumover_x(0.5*(grad_phi**2),U_sumover_x,nx,nxe,ny,nypmx,nz,nzpmx,nvpy,&
							& nvpz,idproc,inorder)
					if (idproc == 0) then
						write (nt_write_U_sumover_x_funit) U_sumover_x
					endif

				endif
			endif

	! field lines diagnostic added by JF
			if (ntlines > 0) then
				it = itime / ntlines
				if (itime == ntlines*it) then
					call write_Eoft_aty(lines,fxyze,nx,ny,nz,nvpy,nvpz,idproc)
				endif
			endif

      wke = 0.
! push electrons
      if (relativity==1) then
         call rpush(part,fxyze,npp,noff,qbme,dt,ci,wke,tpush,nx,ny,nz,ip&
     &bc,inorder,popt)
!        call rpush(part,fxyze,bxyze,npp,noff,qbme,dt,dt,ci,wke,tpush,nx&
!    &,ny,nz,ipbc,inorder,popt)
      else
         call push(part,fxyze,npp,noff,qbme,dt,wke,tpush,nx,ny,nz,ipbc,i&
     &norder,popt)
!        call push(part,fxyze,bxyze,npp,noff,qbme,dt,dt,wke,tpush,nx,ny,&
!    &nz,ipbc,inorder,popt)
      endif
! move electrons into appropriate spatial regions
     call pmove(part,edges,npp,tmove,ny,nz,kstrt,nvpy,nvpz,nbmax,idds,i&
    &blok,vect,ierr)
!      call pmoves(part,edges,npp,tmove,ny,nz,kstrt,nvpy,nvpz,nbmax,idds,&
!     &iblok,vect,ierr)
! push ions
      if (movion==1) then
         wki = 0.
         if (relativity==1) then
            call rpush(parti,fxyze,nppi,noff,qbmi,dt,ci,wke,tpushi,nx,ny&
     &,nz,ipbc,inorder,popt)
!           call rpush(parti,fxyze,bxyze,nppi,noff,qbmi,dt,dt,ci,wke,tpu&
!    &shi,nx,ny,nz,ipbc,inorder,popt)
         else
            call push(parti,fxyze,nppi,noff,qbmi,dt,wke,tpushi,nx,ny,nz,&
     &ipbc,inorder,popt)
!           call push(parti,fxyze,bxyze,nppi,noff,qbmi,dt,dt,wke,tpushi,&
!    &nx,ny,nz,ipbc,inorder,popt)
         endif
! move ions into appropriate spatial regions
        call pmove(parti,edges,nppi,tmovi,ny,nz,kstrt,nvpy,nvpz,nbmax,i&
    &dds,iblok,vect,ierr)
!         call pmoves(parti,edges,nppi,tmovi,ny,nz,kstrt,nvpy,nvpz,nbmax,&
!     &idds,iblok,vect,ierr)
      endif
! momentum diagnostic
      if (ntm > 0) then
         it = itime/ntm
         it = itime - ntm*it + 1
         if (it > 1) it = it - ntm
! calculate electron momentum
         if (it >= 0) then
            call premoment3(part,itime,npp,id0,ium,pxe,pye,pze,sx,sy,sz,&
     &wx,wy,wz,nprint=it)
! calculate ion momentum
            if (movion==0) then
               if (it==1) then
                  call imoment(qi,fxyze,nyzp,id0,ium,pxi,pyi,pzi,dt,wx,w&
     &y,wz,nx,inorder)
               endif
            else if (movion==1) then
               call primoment3(parti,nppi,id0,ium,rmass,pxi,pyi,pzi,wx,w&
     &y,wz,nprint=it)
            endif
! print total momentum
            if (it==1) then
               if (id0==0) write (ium,996) wx, wy, wz
            endif
         endif
      endif
! sort electrons
      if (sortime > 0) then
         if (mod(ntime,sortime)==0) then
           call sortp(part,pt,ip,npp,noff,nyzp,npic,tsort,inorder)
!            call sortp(part,part2,npp,noff,nyzp,npic,tsort,inorder)
         endif
      endif
! sort ions
      if ((movion==1) .and. (sortimi > 0)) then
         if (mod(ntime,sortimi)==0) then
            call sortp(parti,pt,ip,nppi,noff,nyzp,npic,tsorti,inorder)
!           call sortp(parti,parti2,nppi,noff,nyzp,npic,tsorti,inorder)
         endif
      endif
! energy diagnostic
      call esenergy(wt,wtot,msg,we,wke,wki,ntw,ndw,id0,itw,iuot,ntime)
      itime = itime + 1
      ntime = itime + itime0
! restart file
      if (ntr > 0) then
         it = ntime/ntr
         if (ntime==ntr*it) then
            it = iur1 + mod(it-1,2)*(iur2 - iur1)
! write file
            call restart_bwrite(it,id0,itime,itime0,nvp,npp,part,movion,&
     &nppi,parti,qi)
            call restart_dwrite(it,id0,itime,itw,wt,ndrec,fdname,nprec, &
     &fpname)
         endif
      endif
      call wtimer(tloop,ldtime)
      ltime = ltime + tloop
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
      msg(1:2) = time; msg(3) = tpush; msg(4) = tdpost; msg(5) = tsort
      msg(6:7) = tmove; msg(8:9) = tfft; msg(10) = ltime
      call HARTBEAT(msg,10)
      if (id0==0) then
         write (iuot,*) 'processor partition used: nvpy, nvpz = ', nvpy,&
     &nvpz
         write (iuot,*) ncpus, ' processors found, ', ntasks+1, ' used'
         write (iuot,*) 'main max/min real time=',time(1),time(2), 'sec'
         totpush = tpush + tdpost
         write (iuot,*) 'electron push time=', tpush, 'sec'
         write (iuot,*) 'electron charge deposit time=', tdpost, 'sec'
         write (iuot,*) 'total electron push time=', totpush, 'sec'
         write (iuot,*) 'electron sort time=', tsort, 'sec'
         write (iuot,*) 'electron move time=', tmove, 'sec'
         totpush = totpush + tsort + tmove(1)
         write (iuot,*) 'total electron time=', totpush, 'sec'
      endif
      if (movion==1) then
         msg(1) = tpushi; msg(2) = tdposti; msg(3) = tsorti
         msg(4:5) = tmovi
         call HARTBEAT(msg,5)
         if (id0==0) then
            totpushi = tpushi + tdposti
            write (iuot,*) 'ion push time=', tpushi, 'sec'
            write (iuot,*) 'ion charge deposit time=', tdposti, 'sec'
            write (iuot,*) 'total ion push time=', totpushi, 'sec'
            write (iuot,*) 'ion sort time=', tsorti
            write (iuot,*) 'ion move time=', tmovi
            totpushi = totpushi + tsorti + tmovi(1)
            write (iuot,*) 'total ion time=', totpushi, 'sec'
         endif
      endif
      if (id0==0) then
         write (iuot,*) 'total fft time=', tfft, 'sec'
         time(1) = time(1) - (totpush + totpushi + tfft(1))
         write (iuot,*) 'other time=', time(1), ltime, 'sec'
! write final diagnostic metafile
         fname = 'pdiag32.'//cdrun
         open(unit=iudm,file=trim(fname),form='formatted',status='replac&
     &e')
! ion density diagnostics
         if (ntd > 0) then
            if (id0==0) then
               ndrec = ndrec - 1
               ceng = zero
               write (iudm,pden32d,iostat=irc)
               if (irc /= 0) then
                  write (iuer,*) 'pden32d namelist not written'
               endif
            endif
         endif
! potential diagnostics
         if (ntp > 0) then
            if (id0==0) then
               nprec = nprec - 1
               ceng = affp
               write (iudm,ppot32d,iostat=irc)
               if (irc /= 0) then
                  write (iuer,*) 'ppot32d namelist not written'
               endif
            endif
         endif
! write out input file
         write (iudm,pinput3,iostat=irc)
         if (irc /= 0) write (iuer,*) 'pinput3 namelist not written'
! done
         write (iuot,*) '* * * q.e.d. * * *'
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
! allocate data for restart and/or phase space diagnostic
         if ((nts > 0) .or. (nustrt /= 1) .or. (ntr > 0)) then
            allocate(part(idimp,max(npmax,npimax),inblok))
            allocate(npp(inblok))
            parti => part
            nppi => npp
         endif
         if (movion==0) then
            allocate(qi(nxe,nypmx,nzpmx*kbmin,kblok*lblok))
         else
            allocate(qi(0,0,0,0))
         endif
! restart
         if (nustrt /= 1) then
            it = 0
            call restart_bread(iur1,iur2,id0,it,itime,itime0,nvp,npp,   &
     &part,movion,nppi,parti,qi,irc,iuer)
            if (irc /= 0) go to 30
! extend run
            if (nustrt==0) then
               itime0 = itime + itime0
               t0 = dt*real(itime0)
               itime = 0
               if (id0==0) then
                  if (iur1 >= 0) close (unit=iur1)
                  if (iur2 >= 0) close (unit=iur2)
                  call restart_open(1,ntr,idrun,iur1,iur2,iuer)
               endif
               go to 40
            endif
! read diagnostics
            call restart_dread(it,id0,itime,itw,wt,iud,ndrec,fdname,iup,&
     &nprec,fpname,irc,iuer)
            if (irc /= 0) go to 30
            t0 = dt*real(itime0)
            if (id0==0) rewind it
            go to 40
! handle error'
   30       if (id0==0) then
               write (iuer,*) 'Restart Error, irc = ', irc
            endif
            call MP_END
            call PPEXIT
            stop
         endif
! initialize diagnostics
   40    if (id0==0) then
            iudm = get_funit(iudm)
            fname = 'pdiag32.init.'//cdrun
            open(unit=iudm,file=trim(fname),form='formatted',status='rep&
     &lace')
         endif
! ion density diagnostics
         call initmodediag(dent,ntd,id0,nxh,nyh,nzh,kxyp,kyzp,modesxd,  &
     &modesyd,modeszd,jblok,mblok,iud,ndrec,fdname)
         if (ntd > 0) then
            if (id0==0) then
               ceng = zero
               write (iudm,pden32d,iostat=irc)
            endif
         endif
! velocity diagnostics
         fname = 'fv3.'//cdrun
         call initveldiag(fv,fvm,vtx,vty,vtz,ntv,ndv,id0,nmv,nblok,iuv,f&
     &name)
         if (movion==1) then
            fname = 'fvi3.'//cdrun
            call initveldiag(fvi,fvmi,vtxi,vtyi,vtzi,ntv,ndv,id0,nmv,nbl&
     &ok,iuvi,fname)
         endif
! potential diagnostics
         call initmodediag(pott,ntp,id0,nxh,nyh,nzh,kxyp,kyzp,modesxp,  &
     &modesyp,modeszp,jblok,mblok,iup,nprec,fpname)
         if (ntp > 0) then
            if (id0==0) then
               ceng = affp
               write (iudm,ppot32d,iostat=irc)
            endif
         endif
! write out and close input file
         if (id0==0) then
            write (iudm,pinput3,iostat=irc)
            close(unit=iudm)
         endif
! get initial CPU Time
         call HARTBEAT(msg,2)
         time(1) = msg(1); time(2) = msg(2)
         if (id0==0) then
            write (iuot,*) 'init max/min real time=',time(1),time(2), 's&
     &ec'
         endif
! get time step
   10    call HARTBEAT(msg,1)
         it = msg(1)
         if (it < 0) then
! get main CPU Time
            call HARTBEAT(msg,10)
            time = msg(1:2); tpush = msg(3); tdpost = msg(4)
            tsort = msg(5); tmove = msg(6:7); tfft = msg(8:9)
            ltime = msg(10)
            if (id0==0) then
               write (iuot,*) 'processor partition used: nvpy, nvpz = ',&
     & nvpy, nvpz
               write (iuot,*) ncpus, ' processors found, ', ntasks+1, ' &
     &used'
               write (iuot,*) 'main max/min real time=',time(1),time(2),&
     &'sec'
               totpush = tpush + tdpost
               write (iuot,*) 'electron push time=', tpush, 'sec'
               write (iuot,*) 'electron charge deposit time=', tdpost, '&
     &sec'
               write (iuot,*) 'total electron push time=',totpush, 'sec'
               write (iuot,*) 'electron sort time=', tsort, 'sec'
               write (iuot,*) 'electron move time=', tmove, 'sec'
               totpush = totpush + tsort + tmove(1)
               write (iuot,*) 'total electron time=', totpush, 'sec'
            endif
            if (movion==1) then
               call HARTBEAT(msg,5)
               tpushi = msg(1); tdposti = msg(2); tsorti = msg(3)
               tmovi = msg(4:5)
               if (id0==0) then
                  totpushi = tpushi + tdposti
                  write (iuot,*) 'ion push time=', tpushi, 'sec'
                  write (iuot,*) 'ion charge deposit time=', tdposti, 's&
     &ec'
                  write (iuot,*) 'total ion push time=', totpushi, 'sec'
                  write (iuot,*) 'ion sort time=', tsorti
                  write (iuot,*) 'ion move time=', tmovi
                  totpushi = totpushi + tsorti + tmovi(1)
                  write (iuot,*) 'total ion time=', totpushi, 'sec'
               endif
            endif
            if (id0==0) then
               write (iuot,*) 'total fft time=', tfft, 'sec'
               time(1) = time(1) - (totpush + totpushi + tfft(1))
               write (iuot,*) 'other time=', time(1), ltime, 'sec'
! write final diagnostic metafile
               fname = 'pdiag32.'//cdrun
               open(unit=iudm,file=trim(fname),form='formatted',status='&
     &replace')
! ion density diagnostics
               if (ntp > 0) then
                  ndrec = ndrec - 1
                  ceng = zero
                  write (iudm,pden32d,iostat=irc)
                  if (irc /= 0) then
                     write (iuer,*) 'pden32d namelist not written'
                  endif
               endif
! potential diagnostics
               if (ntp > 0) then
                  nprec = nprec - 1
                  ceng = affp
                  write (iudm,ppot32d,iostat=irc)
                  if (irc /= 0) then
                     write (iuer,*) 'ppot32d namelist not written'
                  endif
               endif
! write out input file
               write (iudm,pinput3,iostat=irc)
               if (irc /= 0) then
                  write (iuer,*) 'pinput3 namelist not written'
               endif
! done
               write (iuot,*) '* * * q.e.d. * * *'
            endif
            call MP_END
            call PPEXIT
            stop
         else
            ntime = it
         endif
         if (id0==0) write (iuot,991) ntime
         write (label,991) ntime
         call LOGNAME(label)
! ion density diagnostic
         if (ntd > 0) then
            it = ntime/ntd
            if (ntime==ntd*it) then
! write diagnostic output
               modesz2 = 2*modeszd - 1
               call writebf(dent,modesxd,modesyd,modesz2,kxyp,kyzp,jblok&
     &,iud,ndrec)
            endif
         endif
! velocity diagnostic
         if (ntv > 0) then
            it = ntime/ntv
            if (ntime==ntv*it) then
! print out velocity moments
               if (id0==0) write (iuv,*) it, fvm(1,:,:), fvm(2,:,:)
               if (movion==1) then
! print out velocity moments
                  if (id0==0) write (iuvi,*) it, fvmi(1,:,:),fvmi(2,:,:)
              endif
            endif
         endif
! potential diagnostic
         if (ntp > 0) then
            it = ntime/ntp
            if (ntime==ntp*it) then
! write diagnostic output
               modesz2 = 2*modeszp - 1
               call writebf(pott,modesxp,modesyp,modesz2,kxyp,kyzp,jblok&
     &,iup,nprec)
            endif
         endif
! energy diagnostic
         if (ntw > 0) then
            it = ntime/ntw
            if (ntime==ntw*it) then
! get energy values
               call HARTBEAT(msg,4)
               wtot = msg(1:4)
               if (id0==0) write (iuot,992) wtot(1), wtot(2), wtot(4)
!
!              it = ntime/ndw
!              if (ntime==ndw*it) then
!                 itw = itw + 1
!                 wt(itw,:) = wtot
!              endif
!
            endif
         endif
         itime = itime + 1
         ntime = itime + itime0
! restart file
         if (ntr > 0) then
            it = ntime/ntr
            if (ntime==ntr*it) then
               it = iur1 + mod(it-1,2)*(iur2 - iur1)
! write file
               call restart_bwrite(it,id0,itime,itime0,nvp,npp,part,    &
     &movion,nppi,parti,qi)
               call restart_dwrite(it,id0,itime,itw,wt,ndrec,fdname,    &
     &nprec,fpname)
            endif
         endif
         go to 10
      end subroutine
!
      end program pbeps32

