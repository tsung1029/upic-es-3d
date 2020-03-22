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
! update: june 14, 2008
      program pbeps32
      use pinit32d
      use pespush32d
      use pnpfield32d
      use pdiag32d
      use mp0d, only: mpinit, ncpus
      implicit none
! idps = number of particle partition boundaries = 4
! idds = dimensionality of domain decomposition = 2
! idimp = dimension of phase space = 6
! mshare = (0,1) = (no,yes) architecture is shared memory
      integer :: idps =    4, idds =    2, idimp =   6, mshare =   0
! nmv = number of segments in v for velocity distribution
      integer :: nmv = 40, vect = 0
      integer :: npxyz, npxyzb, np, npxyzi, npxyzbi, npi
      integer :: nx, ny, nz, nxh, nyh, nzh, nyv, nzv, nxe, nxeh, ipbc
      integer :: nloop, nvpy, nvpz, nvp, iblok, nblok, inblok
      integer :: npav, npmax, npimax = 0, kyp, kzp, nypmx, nzpmx
      integer :: kxyp, kyzp, kyb, kxb, kzb, kyzb, kxyb, kzyb
      integer :: kbmin, kblok, jbmin, jblok, lbmin, lblok, mbmin, mblok
      integer :: ngds, nxyzh, nxhyz, nx1, nypm1, nzpm1, nyzpm1, myzpm1
      integer :: nbmax, idproc, id0, kstrt, itime, isign, irc, ierr
      integer :: mterf, mterg, nterf, nterg, it, modesz2
      integer :: ntasks
      real :: zero = 0.0, wki = 0.0, anpav = 0.0, pibal = 0.0
      real :: tpush = 0.0, tdpost = 0.0, tsort = 0.0
      real :: tpushi = 0.0, tdposti = 0.0, tsorti = 0.0
      real :: totpush = 0.0, totpushi = 0.0
      real :: trepart = 0.0, tfmove = 0.0, ts = 0.0
      real :: qbme, qbmi, affp, qi0, etx, we, wke
      real :: vtxi, vtyi, vtzi, vtdxi, vtdyi, vtdzi
      double precision :: dtime, etime
      real, dimension(:,:,:), pointer :: part, part2, parti, parti2
      real, dimension(:,:,:,:), pointer :: qe, qi
      real, dimension(:,:,:,:,:), pointer :: fxyze
!     real, dimension(:,:,:,:,:), pointer :: bxyze
      complex, dimension(:,:,:,:), pointer :: qt, qs
      complex, dimension(:,:,:,:,:), pointer :: fxyzt, fxyzs
      complex, dimension(:,:,:,:), pointer :: ffc
      integer, dimension(:), pointer :: mixup
      complex, dimension(:), pointer :: sct
      real, dimension(:,:), pointer :: edges
      integer, dimension(:,:), pointer :: nyzp, noff, nyzpu, noffu
      integer, dimension(:), pointer :: npp, nppi, nps
      real, dimension(:,:), pointer :: pt
      integer, dimension(:,:), pointer :: ip, npic
      integer, dimension(:,:), pointer :: nyzps, noffs
      integer, dimension(:,:,:), pointer :: npicyz
      real, dimension(:,:,:,:), pointer :: sfield
      complex, dimension(:,:,:,:), pointer :: sfieldt, pott
      real, dimension(:,:,:), pointer :: fv, fvm, fvi, fvmi
      real, dimension(:,:), pointer :: wt
! dirichlet or vacuum boundary conditions
      integer :: indx1, indy1, indz1, nx2, ny2, nz2, nx2e
      integer :: kyp2, kzp2, kxyp2, kyzp2, kyb2, kxb2, kzb2, kyzb2
      integer :: kxyb2, kzyb2, k2bmin, k2blok, j2bmin, j2blok
      integer :: l2bmin, l2blok, m2bmin, m2blok
      real, dimension(:,:,:,:), pointer :: q3, sfield3
      real, dimension(:,:,:,:,:), pointer :: fxyz3
      complex, dimension(:,:,:,:), pointer :: qt3, qs3, sfieldt3
      complex, dimension(:,:,:,:,:), pointer :: fxyzt3, fxyzs3
      complex, dimension(:,:,:,:), pointer :: ffd
      integer, dimension(:), pointer :: mixup3
      complex, dimension(:), pointer :: sct3
! semi-periodic boundary conditions
      integer :: nxhy2z, nxy2zh, kzyb1, l1bmin, m1bmin
      real, dimension(:,:,:,:), pointer :: q2, sfield2
      real, dimension(:,:,:,:,:), pointer :: fxyz2
      complex, dimension(:,:,:,:), pointer :: qt2, qs2, sfieldt2
      complex, dimension(:,:,:,:,:), pointer :: fxyzt2, fxyzs2
      complex, dimension(:,:,:,:), pointer :: ffb
      integer, dimension(:), pointer :: mixup2
      complex, dimension(:), pointer :: sct2
! vacuum boundary conditions
      real, dimension(:,:,:,:,:), pointer :: ffg
!
      integer, dimension(2) :: ktime
! wtot = total energy
      real, dimension(4) :: wtot
! time = timing array
      real, dimension(2) :: tfft = 0.0, time = 0.0
      real, dimension(2) :: tmove = 0.0, tmovi = 0.0
! msg = heartbeat array
      double precision, dimension(11) :: msg
      character(len=10) :: cdrun
      character(len=32) :: fname
      character(len=12) :: label
  991 format (' T = ',i7)
  992 format (' field, kinetic, total energies = ',3e14.7)
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
! open initial diagnostic metafile
         fname = 'pdiag3.init.'//cdrun
         open(unit=19,file=trim(fname),form='formatted',status='replace'&
     &)
      endif
! broadcast namelist to other nodes
      call sendnml()
! set mpi monitor value
      call SET_MON(mpimon)
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
      kyp = (ny - 1)/nvpy + 1; nypm1 = 3*kyp + 1
! kzp = number of complex grids in each field partition in z direction
      kzp = (nz - 1)/nvpz + 1; nzpm1 = 3*kzp + 1
! nypmx = maximum size of particle partition in y, including guard cells
! nzpmx = maximum size of particle partition in z, including guard cells
      nypmx = nypm1 + 2; nzpmx = nzpm1 + 2
! ngds = number of guard cells
      ngds = 3*((idds - 1)/2 + 1)
      if (inorder==LINEAR) then
         ax = .912871; ay = .912871; az = .912871
         nxe = nx + 2; nypmx = nypm1; nzpmx = nzpm1
         ngds = (idds - 1)/2 + 1
      endif
      nxeh = nxe/2
! boundary conditions
      ipbc = psolve
      if (psolve==VACUUM_3D) ipbc = 2
! initialize for multiprocessing
      ntasks = mpinit(sntasks)
      if (popt==VECTOR) vect = 1
! nloop = number of time steps in simulation
      nloop = tend/dt + .0001
! iblok/nblok = number of particle partitions in y/z
      iblok = 1 + mshare*(nvpy - 1); nblok = 1 + mshare*(nvpz - 1)
      inblok = iblok*nblok
! npav = average number of particles per processor
! npmax = maximum number of particles in each partition
      npav = np/nvp; npmax = npav*1.2 + 60000
      if (movion==1) npimax = (npi/nvp)*1.2 + 60000
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
      nx1 = nx + 1; nyzpm1 = nypm1*nzpm1; myzpm1 = max(nypm1,nzpm1)
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
         allocate(fv(2*nmv+2,3,nblok),fvm(3,3,nblok))
         if (id0==0) then
            fname = 'fv3.'//cdrun
            open(unit=10,file=trim(fname),form='formatted',status='unkno&
     &wn')
! write captions
         write (10,*) 'it vdx vdy vdz vtx vty vtz'
         endif
         if (movion==1) then
            allocate(fvi(2*nmv+2,3,nblok),fvmi(3,3,nblok))
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
         if (modesxp > nxh) modesxp = nxh
         if (modesyp > nyh) modesyp = nyh
         if (modeszp > nzh) modeszp = nzh
         modesz2 = 2*modeszp - 1
         allocate(pott(modesz2,min(modesxp,kxyp),min(modesyp,kyzp),jblok&
     &*mblok))
         if (id0==0) then
            write (19,ppot32d,iostat=irc)
         endif
      endif
! write out and close input file
      if (id0==0) then
         write (19,pinput3,iostat=irc)
         close(unit=19)
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
! where kk = k + noff(1,m) - 1, ll = l + noff(2,m) - 1
      allocate(qe(nxe,nypmx,nzpmx*kbmin,kblok*lblok))
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
! nyzpu(1,m) = number of primary gridpoints in y in uniform partition m
! nyzpu(2,m) = number of primary gridpoints in z in uniform partition m
! noffu(1,m) = lowermost global gridpoint in y in uniform partition m
! noffu(2,m) = backmost global gridpoint in z in uniform partition m
      allocate(nyzp(idds,inblok),noff(idds,inblok))
      allocate(nyzpu(idds,inblok),noffu(idds,inblok))
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
! data for moving field partitions
      allocate(nyzps(idds,inblok),noffs(idds,inblok))
! data for repartitioning
      allocate(npicyz(myzpm1,idds,inblok))
! non-periodic boundary conditions
      if (psolve /= PERIODIC_3D) then
         indx1 = indx + 1; indy1 = indy + 1; indz1 = indz + 1
         nx2 = 2*nx; ny2 = 2*ny; nz2 = 2*nz; nx2e = 2*nxe
         kyp2 = (ny2 - 1)/nvpy + 1; kzp2 = (nz2 - 1)/nvpz + 1
         kxyp2 = (nx - 1)/nvpy + 1; kyzp2 = (ny2 - 1)/nvpz + 1
         kyb2 = ny2/kyp2; kxb2 = nx/kxyp2; kzb2 = nz2/kzp2
         kyzb2 = ny2/kyzp2
         kxyb2 = max(kxb2,kyb2); kzyb2 = max(kyzb2,kzb2)
         k2bmin = 1 + (1 - mshare)*(kxyb2/kxb2 - 1)
         k2blok = 1 + mshare*(ny2/kyp2 - 1)
         j2bmin = 1 + (1 - mshare)*(kxyb2/kyb2 - 1)
         j2blok = 1 + mshare*(nx/kxyp2 - 1)
         m2bmin = 1 + (1 - mshare)*(kzyb2/kzb2 - 1)
         m2blok = 1 + mshare*(ny2/kyzp2 - 1)
      endif
! dirichlet boundary conditions
      if (psolve==DIRICHLET_3D) then
         nxy2zh = 2*nxyzh; nxhy2z = 2*nxhyz
         l2bmin = 1 + (1 - mshare)*(kzyb2/kyzb2 - 1)
         l2blok = 1 + mshare*(nz2/kzp2 - 1)
         allocate(q3(nx2e,2*nypmx,kzp2*k2bmin,k2blok*l2blok))
         allocate(fxyz3(3,nx2e,2*nypmx,kzp2*k2bmin,k2blok*l2blok))
         allocate(sfield3(nx2e,2*nypmx,kzp2*k2bmin,k2blok*l2blok))
         allocate(qt3(nz2,kxyp2,kyzp2*m2bmin,j2blok*m2blok))
         allocate(fxyzt3(3,nz2,kxyp2,kyzp2*m2bmin,j2blok*m2blok))
         allocate(qs3(ny2,kxyp2,kzp2*j2bmin*l2bmin,j2blok*l2blok))
         allocate(fxyzs3(3,ny2,kxyp2,kzp2*j2bmin*l2bmin,j2blok*l2blok))
         allocate(sfieldt3(nz2,kxyp2,kyzp2*m2bmin,j2blok*m2blok))
         allocate(ffd(nz,kxyp2,kyzp2,j2blok*m2blok))
         allocate(mixup3(nxhy2z), sct3(nxy2zh))
! semi-periodic boundary conditions
      else if (psolve==MIXED_3D) then
         nxy2zh = max(nx2,ny2,nz)/2
         nxhy2z = max(nx,ny2,nz)
         kzyb1 = max(kyzb2,kzb)
         l1bmin=1 + (1 - mshare)*(kzyb1/kyzb2 - 1)
         m1bmin=1 + (1 - mshare)*(kzyb1/kzb - 1)
         allocate(q2(nx2e,2*nypmx,kzp*kbmin,k2blok*lblok))
         allocate(fxyz2(3,nx2e,2*nypmx,kzp*kbmin,k2blok*lblok))
         allocate(sfield2(nx2e,2*nypmx,kzp*kbmin,k2blok*lblok))
         allocate(qt2(nz,kxyp2,kyzp2*m1bmin,j2blok*m2blok))
         allocate(fxyzt2(3,nz,kxyp2,kyzp2*m1bmin,j2blok*m2blok))
         allocate(qs2(ny2,kxyp2,kzp*jbmin*l1bmin,j2blok*lblok))
         allocate(fxyzs2(3,ny2,kxyp2,kzp*jbmin*l1bmin,j2blok*lblok))
         allocate(sfieldt2(nz,kxyp2,kyzp2*m1bmin,j2blok*m2blok))
         allocate(ffb(nzh,kxyp2,kyzp2,j2blok*m2blok))
         allocate(mixup2(nxhy2z), sct2(nxy2zh))
! vacuum boundary conditions
      else if (psolve==VACUUM_3D) then
         nxy2zh = 2*nxyzh; nxhy2z = 2*nxhyz
         l2bmin = 1 + (1 - mshare)*(kzyb2/kyzb2 - 1)
         l2blok = 1 + mshare*(nz2/kzp2 - 1)
         allocate(q3(nx2e,2*nypmx,kzp2*k2bmin,k2blok*l2blok))
         allocate(fxyz3(3,nx2e,2*nypmx,kzp2*k2bmin,k2blok*l2blok))
         allocate(sfield3(nx2e,2*nypmx,kzp2*k2bmin,k2blok*l2blok))
         allocate(qt3(nz2,kxyp2,kyzp2*m2bmin,j2blok*m2blok))
         allocate(fxyzt3(3,nz2,kxyp2,kyzp2*m2bmin,j2blok*m2blok))
         allocate(qs3(ny2,kxyp2,kzp2*j2bmin*l2bmin,j2blok*l2blok))
         allocate(fxyzs3(3,ny2,kxyp2,kzp2*j2bmin*l2bmin,j2blok*l2blok))
         allocate(sfieldt3(nz2,kxyp2,kyzp2*m2bmin,j2blok*m2blok))
         allocate(ffg(5,nz+1,kxyp2+1,kyzp2+1,j2blok))
         allocate(mixup3(nxhy2z), sct3(nxy2zh))
      endif
!
! initialize parallel timer
      call pwtimer(time,dtime,-1)
! initialize constants
      itime = 0
      mterf = 0
      mterg = kyp - 1
      nterf = 0
      nterg = kzp - 1
      qbme = qme
      if (ipbc==1) then
         affp = float(nx*ny*nz)/float(np)
      else if (ipbc==2) then
         affp = float((nx-2)*(ny-2)*(nz-2))/float(np)
      else if (ipbc==3) then
         affp = float((nx-2)*(ny-2)*nz)/float(np)
      endif
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
      call dcomp(edges,nyzpu,noffu,ny,nz,kstrt,nvpy,nvpz,iblok,inorder)
! prepare fft tables
      call fft_init(mixup,sct,indx,indy,indz)
! calculate form factors
      call pois_init(ffc,ax,ay,az,affp,nx,ny,nz,kstrt,jblok)
! dirichlet boundary conditions
      if (psolve==DIRICHLET_3D) then
         call fft_init(mixup3,sct3,indx1,indy1,indz1)
         call poisd_init(ffd,ax,ay,az,affp,nx,ny,nz,kstrt,j2blok)
! semi-periodic boundary conditions
      else if (psolve==MIXED_3D) then
         call fft_init(mixup2,sct2,indx1,indy1,indz)
         call poism_init(ffb,ax,ay,az,affp,nx,ny,nz,kstrt,j2blok)
! vacuum boundary conditions
      else if (psolve==VACUUM_3D) then
         call fft_init(mixup3,sct3,indx1,indy1,indz1)
         call poisc3_init(ffg,q3,qs3,qt3,mixup3,sct3,ax,affp,indx,indy,i&
     &ndz,kstrt,kyp2,j2blok,k2blok)
      endif
!
! allocate background charge density
      if (movion==0) allocate(qi(nxe,nypmx,nzpmx*kbmin,kblok*lblok))
! allocate ion data
      if (movion==1) then
         allocate(parti(idimp,npimax,inblok),nppi(inblok))
         allocate(parti2(0,0,0))
      endif
! debug
      if (ipbc==2) vdz = 0.
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
! find new partition analytically
         if (imbalance >= 0.0) then
            call fedges(edges,noff,nyzp,ampdy,scaledy,shiftdy,ampdz,scal&
     &edz,shiftdz,ny,nz,kstrt,nvpy,nvpz,nypmx,nzpmx,ipbc,ndprof,iblok,mt&
     &erg,nterg,ierr,inorder)
            if (ierr /= 0) then
               call MP_END
               call PPEXIT
               stop
            endif
! use uniform partition
         else
            noff = noffu; nyzp = nyzpu
         endif
! move electrons into appropriate spatial regions
         call pmove(part,edges,npp,tmove,ny,nz,kstrt,nvpy,nvpz,nbmax,idd&
     &s,iblok,vect,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! calculate actual coordinates from guiding centers
!        if (relativity==1) then
!           call distr(part,bxyze,npp,noff,qbme,ci,nx,ny,nz,ipbc,iblok,i&
!    &norder)
!        else
!           call distr(part,bxyze,npp,noff,qbme,nx,ny,nz,ipbc,iblok,inor&
!    &der)
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
! use electron partition for ions
! move ions into appropriate spatial regions
            call pmove(parti,edges,nppi,tmovi,ny,nz,kstrt,nvpy,nvpz,nbma&
     &x,idds,iblok,vect,ierr)
            if (ierr /= 0) then
               call MP_END
               call PPEXIT
               stop
            endif
! calculate actual coordinates from guiding centers
!           if (relativity==1) then
!              call distr(parti,bxyze,nppi,noff,qbmi,ci,nx,ny,nz,ipbc,ib&
!    &lok,inorder)
!           else
!              call distr(parti,bxyze,nppi,noff,qbmi,nx,ny,nz,ipbc,iblok&
!    &,inorder)
!           endif
         endif
! initialize background charge density
         if (movion==0) then
            qi0 = -qme/affp
!           call sguardp(qi,kstrt,nvpy,nvpz,noff,nyzp,qi0,nx,ny,nz,kblok&
!    &,ipbc,inorder)
            call sguardp(qi,kstrt,nvpy,nvpz,noff,nyzp,zero,nx,ny,nz,kblo&
     &k,ipbc,inorder)
            call dpost(part,qi,-qme,npp,noff,tdpost,inorder,dopt)
! debug:
!           call sguardp(qi,kstrt,nvpy,nvpz,noff,nyzp,qi0,nx,ny,nz,kblok&
!    &,ipbc,inorder)
! freeze the ions now
         else if ((movion==1).and.(itime==ionoff)) then
            allocate(qi(nxe,nypmx,nzpmx*kbmin,kblok*lblok))
! initialize ion charge density to zero
            call sguardp(qi,kstrt,nvpy,nvpz,noff,nyzp,zero,nx,ny,nz,kblo&
     &k,ipbc,inorder)
! deposit ion charge
            call dpost(parti,qi,qmi,nppi,noff,tdposti,inorder,dopt)
! delete ions
            deallocate(parti,nppi)
            movion = 0
         endif
! add guard cells for ion density in x
         if (movion==0) then
            call aguardp(qi,nyzp,nx,ipbc,inorder)
! add guard cells for ion density in y and z
            call pnaguardp(qi,nyzp,kstrt,nvpy,nvpz,nx,mterg,nterg,ngds,i&
     &pbc,kblok,inorder)
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
         call rddata(edges,nvp,it,ierr)
         if (ierr /= 0) go to 400
         call fnoff(edges,noff,nyzp,nypmx,nzpmx,mterg,nterg,ierr,inorder&
     &)
         if (ierr /= 0) go to 400
         if (movion==1) then
            call rddata(parti,nppi,it,ierr)
            if (ierr /= 0) go to 400
         endif
         if (movion==0) then
            call rddata(qi,nvp,it,ierr)
            if (ierr /= 0) go to 400
            isign = 1
            call pfmove(qi,noff,nyzp,isign,tfmove,kyp,kzp,kstrt,nvpy,nvp&
     &z,kblok,mterf,nterf,ierr,inorder)
            if (ierr /= 0) then
               call MP_END
               call PPEXIT
               stop
            endif
            call zguard(qi,nyzp,nx,inorder)
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
               call writebf(pott,modesxp,modesyp,modesz2,kxyp,kyzp,jblok&
     &,11,irc,trim(fname))
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
! initialize charge density to background
      call sguardp(qe,kstrt,nvpy,nvpz,noff,nyzp,zero,nx,ny,nz,kblok,ipbc&
     &,inorder)
! deposit electron charge
      call dpost(part,qe,qme,npp,noff,tdpost,inorder,dopt)
! density diagnostic
      if (ntp > 0) then
         it = itime/ntp
         if (itime==ntp*it) then
            sfield = -qe
! add guard cells for density in x
            call aguardp(sfield,nyzp,nx,ipbc,inorder)
! add guard cells for density in y and z
            call pnaguardp(sfield,nyzp,kstrt,nvpy,nvpz,nx,mterg,nterg,ng&
     &ds,ipbc,kblok,inorder)
! move density to uniform field partition
            isign = -1
            call pfmove(sfield,noff,nyzp,isign,tfmove,kyp,kzp,kstrt,nvpy&
     &,nvpz,kblok,mterf,nterf,ierr,inorder)
            if (ierr /= 0) then
               call MP_END
               call PPEXIT
               stop
            endif
         endif
      endif
! add ion charge density
      if (movion==1) then
         call dpost(parti,qe,qmi,nppi,noff,tdposti,inorder,dopt)
      else
         qe = qe + qi
      endif
! add guard cells for density in x
      call aguardp(qe,nyzp,nx,ipbc,inorder)
! add guard cells for density in y and z
      call pnaguardp(qe,nyzp,kstrt,nvpy,nvpz,nx,mterg,nterg,ngds,ipbc,kb&
     &lok,inorder)
! freeze the ions
      if ((movion==1).and.(itime==ionoff)) then
         allocate(qi(nxe,nypmx,nzpmx*kbmin,kblok*lblok))
! initialize ion charge density to zero
         call sguardp(qi,kstrt,nvpy,nvpz,noff,nyzp,zero,nx,ny,nz,kblok,i&
     &pbc,inorder)
! deposit ion charge
         call dpost(parti,qi,qmi,nppi,noff,tdposti,inorder,dopt)
! add guard cells for ion density in x
         call aguardp(qi,nyzp,nx,ipbc,inorder)
! add guard cells for ion density in y and z
         call pnaguardp(qi,nyzp,kstrt,nvpy,nvpz,nx,mterg,nterg,ngds,ipbc&
     &,kblok,inorder)
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
! move charge density to uniform field partition
      isign = -1
      call pfmove(qe,noff,nyzp,isign,tfmove,kyp,kzp,kstrt,nvpy,nvpz,kblo&
     &k,mterf,nterf,ierr,inorder)
      if (ierr /= 0) then
         call MP_END
         call PPEXIT
         stop
      endif
!
! dirichlet boundary conditions
!
      if (psolve==DIRICHLET_3D) then
! density diagnostic
      if (ntd > 0) then
         it = itime/ntd
         if (itime==ntd*it) then
! transform electron density to fourier space
            call trpsin(sfield,q3,nx,ny,nz,kstrt,nvpy,nvpz,kyp,kyp2,kzp,&
     &kzp2,kblok,k2blok,inorder)
            isign = -1
            call fft(q3,qs3,qt3,isign,mixup3,sct3,tfft,indx1,indy1,indz1&
     &,kstrt,kxyp2,kyp2,kyzp2,kzp2,k2blok,m2blok,LINEAR)
! calculate smoothing in fourier space
            call poisd(qt3,sfieldt3,ffd,nx,ny,nz,kstrt,j2blok)
! transform electron density to real space
            isign = 1
            call fft(sfield3,qs3,sfieldt3,isign,mixup3,sct3,tfft,indx1,i&
     &ndy1,indz1,kstrt,kxyp2,kyp2,kyzp2,kzp2,k2blok,m2blok,LINEAR)
            call haftrp(sfield,sfield3,nx,ny,nz,kstrt,nvpy,kyp,kyp2,kzp,&
     &kzp2,kblok,k2blok,inorder)
            call pcguardp(sfield,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,ipbc,kb&
     &lok,inorder)
            call cguardp(sfield,nyzpu,nx,ipbc,inorder)
         endif
      endif
! transform charge to fourier space
      call trpsin(qe,q3,nx,ny,nz,kstrt,nvpy,nvpz,kyp,kyp2,kzp,kzp2,kblok&
     &,k2blok,inorder)
      isign = -1
      call fft(q3,qs3,qt3,isign,mixup3,sct3,tfft,indx1,indy1,indz1,kstrt&
     &,kxyp2,kyp2,kyzp2,kzp2,k2blok,m2blok,LINEAR)
! potential diagnostic
      if (ntp > 0) then
         it = itime/ntp
         if (itime==ntp*it) then
! solve for potential
            call poisd(qt3,sfieldt3,ffd,we,nx,ny,nz,kstrt,j2blok)
            isign = 1
            call fft(sfield3,qs3,sfieldt3,isign,mixup3,sct3,tfft,indx1,i&
     &ndy1,indz1,kstrt,kxyp2,kyp2,kyzp2,kzp2,k2blok,m2blok,LINEAR)
            call haftrp(sfield,sfield3,nx,ny,nz,kstrt,nvpy,kyp,kyp2,kzp,&
     &kzp2,kblok,k2blok,inorder)
            call pcguardp(sfield,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,ipbc,kb&
     &lok,inorder)
            call cguardp(sfield,nyzpu,nx,ipbc,inorder)
! write diagnostic output
!           write (nlabel,'(i4)') it
!           fname = trim(potname)//'_'//trim(adjustl(nlabel))
!           nrec = -lprec
!           call writef(sfield,nx1,ny1,kzp,11,nrec,trim(fname),inorder)
         endif
      endif
! calculate force/charge in fourier space
      call poisd(qt3,fxyzt3,ffd,we,nx,ny,nz,kstrt,j2blok)
      isign = 1
      call fft(fxyz3,fxyzs3,fxyzt3,isign,mixup3,sct3,tfft,indx1,indy1,in&
     &dz1,kstrt,kxyp2,kyp2,kyzp2,kzp2,k2blok,m2blok,LINEAR)
      call haftrp(fxyze,fxyz3,nx,ny,nz,kstrt,nvpy,kyp,kyp2,kzp,kzp2,kblo&
     &k,k2blok,inorder)
!
! semi-periodic boundary conditions
!
      else if (psolve==MIXED_3D) then
! density diagnostic
      if (ntd > 0) then
         it = itime/ntd
         if (itime==ntd*it) then
! transform electron density to fourier space
            call dblsin(sfield,q2,nx,ny,kstrt,nvpy,kyp,kyp2,kzp,kblok,k2&
     &blok,inorder)
            isign = -1
            call fft(q2,qs2,qt2,isign,mixup2,sct2,tfft,indx1,indy1,indz,&
     &kstrt,kxyp2,kyp2,kyzp2,kzp,k2blok,m2blok,LINEAR)
! calculate smoothing in fourier space
            call poism(qt2,sfieldt2,ffb,nx,ny,nz,kstrt,j2blok)
! transform electron density to real space
            isign = 1
            call fft(sfield2,qs2,sfieldt2,isign,mixup2,sct2,tfft,indx1,i&
     &ndy1,indz,kstrt,kxyp2,kyp2,kyzp2,kzp,k2blok,m2blok,LINEAR)
            call hafdbl(sfield,sfield2,nx,ny,kstrt,nvpy,nvpz,kyp,kyp2,kz&
     &p,kblok,k2blok,inorder)
            call pcguardp(sfield,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,ipbc,kb&
     &lok,inorder)
            call cguardp(sfield,nyzpu,nx,ipbc,inorder)
         endif
      endif
! transform charge to fourier space
      call dblsin(qe,q2,nx,ny,kstrt,nvpy,kyp,kyp2,kzp,kblok,k2blok,inord&
     &er)
      isign = -1
      call fft(q2,qs2,qt2,isign,mixup2,sct2,tfft,indx1,indy1,indz,kstrt,&
     &kxyp2,kyp2,kyzp2,kzp,k2blok,m2blok,LINEAR)
! potential diagnostic
      if (ntp > 0) then
         it = itime/ntp
         if (itime==ntp*it) then
! solve for potential
            call poism(qt2,sfieldt2,ffb,we,nx,ny,nz,kstrt,j2blok)
            isign = 1
            call fft(sfield2,qs2,sfieldt2,isign,mixup2,sct2,tfft,indx1,i&
     &ndy1,indz,kstrt,kxyp2,kyp2,kyzp2,kzp,k2blok,m2blok,LINEAR)
            call hafdbl(sfield,sfield2,nx,ny,kstrt,nvpy,nvpz,kyp,kyp2,kz&
     &p,kblok,k2blok,inorder)
            call pcguardp(sfield,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,ipbc,kb&
     &lok,inorder)
            call cguardp(sfield,nyzpu,nx,ipbc,inorder)
! write diagnostic output
!           write (nlabel,'(i4)') it
!           fname = trim(potname)//'_'//trim(adjustl(nlabel))
!           nrec = -lprec
!           call writef(sfield,nx1,ny1,kzp,11,nrec,trim(fname),inorder)
         endif
      endif
! calculate force/charge in fourier space
      call poism(qt2,fxyzt2,ffb,we,nx,ny,nz,kstrt,j2blok)
      isign = 1
      call fft(fxyz2,fxyzs2,fxyzt2,isign,mixup2,sct2,tfft,indx1,indy1,in&
     &dz,kstrt,kxyp2,kyp2,kyzp2,kzp,k2blok,m2blok,LINEAR)
      call hafdbl(fxyze,fxyz2,nx,ny,kstrt,nvpy,nvpz,kyp,kyp2,kzp,kblok,k&
     &2blok,inorder)
!
! vacuum boundary conditions
!
      else if (psolve==VACUUM_3D) then
! density diagnostic
      if (ntd > 0) then
         it = itime/ntd
         if (itime==ntd*it) then
! transform electron density to fourier space
            call ztrp(sfield,q3,nx,ny,nz,kstrt,nvpy,kyp,kyp2,kzp,kzp2,kb&
     &lok,k2blok,inorder)
            isign = -1
            call fft(q3,qs3,qt3,isign,mixup3,sct3,tfft,indx1,indy1,indz1&
     &,kstrt,kxyp2,kyp2,kyzp2,kzp2,k2blok,m2blok,LINEAR)
! calculate smoothing in fourier space
            call poisc(qt3,sfieldt3,ffg,nx,ny,nz,kstrt,j2blok)
! transform electron density to real space
            isign = 1
            call fft(sfield3,qs3,sfieldt3,isign,mixup3,sct3,tfft,indx1,i&
     &ndy1,indz1,kstrt,kxyp2,kyp2,kyzp2,kzp2,k2blok,m2blok,LINEAR)
            call haftrp(sfield,sfield3,nx,ny,nz,kstrt,nvpy,kyp,kyp2,kzp,&
     &kzp2,kblok,k2blok,inorder)
            call pcguardp(sfield,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,ipbc,kb&
     &lok,inorder)
            call cguardp(sfield,nyzpu,nx,ipbc,inorder)
         endif
      endif
! transform charge to fourier space
      call ztrp(qe,q3,nx,ny,nz,kstrt,nvpy,kyp,kyp2,kzp,kzp2,kblok,k2blok&
     &,inorder)
      isign = -1
      call fft(q3,qs3,qt3,isign,mixup3,sct3,tfft,indx1,indy1,indz1,kstrt&
     &,kxyp2,kyp2,kyzp2,kzp2,k2blok,m2blok,LINEAR)
! potential diagnostic
      if (ntp > 0) then
         it = itime/ntp
         if (itime==ntp*it) then
! solve for potential
            call poisc(qt3,sfieldt3,ffg,we,nx,ny,nz,kstrt,j2blok)
            isign = 1
            call fft(sfield3,qs3,sfieldt3,isign,mixup3,sct3,tfft,indx1,i&
     &ndy1,indz1,kstrt,kxyp2,kyp2,kyzp2,kzp2,k2blok,m2blok,LINEAR)
            call haftrp(sfield,sfield3,nx,ny,nz,kstrt,nvpy,kyp,kyp2,kzp,&
     &kzp2,kblok,k2blok,inorder)
            call pcguardp(sfield,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,ipbc,kb&
     &lok,inorder)
            call cguardp(sfield,nyzpu,nx,ipbc,inorder)
! write diagnostic output
!           write (nlabel,'(i4)') it
!           fname = trim(potname)//'_'//trim(adjustl(nlabel))
!           nrec = -lprec
!           call writef(sfield,nx1,ny1,kzp,11,nrec,trim(fname),inorder)
         endif
      endif
! calculate force/charge in fourier space
      call poisc(qt3,fxyzt3,ffg,we,nx,ny,nz,kstrt,j2blok)
      isign = 1
      call fft(fxyz3,fxyzs3,fxyzt3,isign,mixup3,sct3,tfft,indx1,indy1,in&
     &dz1,kstrt,kxyp2,kyp2,kyzp2,kzp2,k2blok,m2blok,LINEAR)
      call haftrp(fxyze,fxyz3,nx,ny,nz,kstrt,nvpy,kyp,kyp2,kzp,kzp2,kblo&
     &k,k2blok,inorder)
!
! periodic boundary conditions
!
      else if (psolve==PERIODIC_3D) then
! density diagnostic
      if (ntp > 0) then
         it = itime/ntp
         if (itime==ntp*it) then
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
            call pcguardp(sfield,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,ipbc,kb&
     &lok,inorder)
            call cguardp(sfield,nyzpu,nx,ipbc,inorder)
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
            call pcguardp(sfield,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,ipbc,kb&
     &lok,inorder)
            call cguardp(sfield,nyzpu,nx,ipbc,inorder)
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
! calculate force/charge in fourier space
      call pois(qt,fxyzt,ffc,we,nx,ny,nz,kstrt,jblok)
! transform force/charge to real space
      isign = 1
      call fft(fxyze,fxyzs,fxyzt,isign,mixup,sct,tfft,indx,indy,indz,kst&
     &rt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
!
      endif
!
! move force/charge to non-uniform field partition
      isign = 1
      call pfmove(fxyze,noff,nyzp,isign,tfmove,kyp,kzp,kstrt,nvpy,nvpz,k&
     &blok,mterf,nterf,ierr,inorder)
      if (ierr /= 0) then
         call MP_END
         call PPEXIT
         stop
      endif
! copy to guard cells
      call pncguardp(fxyze,nyzp,kstrt,nvpy,nvpz,nx,mterg,nterg,ngds,ipbc&
     &,kblok,inorder)
      call cguardp(fxyze,nyzp,nx,ipbc,inorder)
! external pump
      if ((itpon > 0).and.(itime >= itpon)) then
         etx = (v0*vtx)*w0*cos(w0*dt*(itime - itpon))
         fxyze(1,:,:,:,:) = fxyze(1,:,:,:,:) + etx
      endif
! particle push and charge density update
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
      call pmove(part,edges,npp,anpav,pibal,tmove,ny,nz,kstrt,nvpy,nvpz,&
     &nbmax,idds,iblok,vect,ierr)
!     call pmove(part,edges,npp,tmove,ny,nz,kstrt,nvpy,nvpz,nbmax,idds,i&
!    &blok,info,tinfo,ierr)
!     totinfo = totinfo + tinfo
!     npav = info(9)/(nvpy*nvpz)
!     if (npav > 0) then
!        pibal = real(max0(info(2)-npav,npav-info(3)))/real(npav)
!     endif
      if (ierr /= 0) then
         call MP_END
         call PPEXIT
         stop
      endif
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
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
      endif
!
! begin repartitioning
      if (imbalance >= 0.0) then
      if (pibal > imbalance) then
! initialize timer
         call wtimer(ts,etime,-1)
! count the number of electrons per cell
         call countp(part,npicyz,npp,noff,nyzp,kstrt,nvpy,nvpz,iblok)
! save old repartitioning boundaries
         noffs = noff
         nyzps = nyzp
! determine new repartitioning boundaries
         call repart(edges,npicyz,noff,nyzp,anpav,kstrt,nvpy,nvpz,nypmx,&
     &nzpmx,iblok,mterg,nterg,ierr,inorder)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         mterf = 0
         nterf = 0
         if (movion==0) then
! move background ion density to new field partition.
            isign = 2
            call pfmove(qi,noff,nyzp,noffs,nyzps,isign,tfmove,kstrt,nvpy&
     &,nvpz,kblok,ierr,inorder)
            if (ierr /= 0) then
               call MP_END
               call PPEXIT
               stop
            endif
! zero out guard cells
            call zguard(qi,nyzp,nx,inorder)
         endif
! move electrons into new spatial regions
         call pmove(part,edges,npp,tmove,ny,nz,kstrt,nvpy,nvpz,nbmax,idd&
     &s,iblok,vect,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
! use electron partition for ions
         if (movion==1) then
! move ions into new spatial regions
            call pmove(parti,edges,nppi,tmovi,ny,nz,kstrt,nvpy,nvpz,nbma&
     &x,idds,iblok,vect,ierr)
            if (ierr /= 0) then
               call MP_END
               call PPEXIT
               stop
            endif
         endif
! record time
         call wtimer(ts,etime)
         trepart = trepart + ts
         if (id0==0) then
            write (18,*) 'repartitioning complete, imbalance = ', pibal
         endif
      endif
      endif
! end repartitioning
!
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
            wtot(1) = we
            wtot(2) = wke
            wtot(3) = wki
            wtot(4) = we + wke + wki
            call plsum(wtot)
! send energy values to diagnostic node
            msg(1:4) = wtot
            call HARTBEAT(msg,4)
            if (id0==0) write (18,992) wtot(1), wtot(2), wtot(4)
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
            call wrdata(edges,nvp,it)
            if (movion==1) call wrdata(parti,nppi,it)
            if (movion==0) then
               isign = -1
               call pfmove(qi,noff,nyzp,isign,tfmove,kyp,kzp,kstrt,nvpy,&
     &nvpz,kblok,mterf,nterf,ierr,inorder)
               if (ierr /= 0) then
                  call MP_END
                  call PPEXIT
                  stop
               endif
               call wrdata(qi,nvp,it)
               isign = 1
               call pfmove(qi,noff,nyzp,isign,tfmove,kyp,kzp,kstrt,nvpy,&
     &nvpz,kblok,mterf,nterf,ierr,inorder)
               if (ierr /= 0) then
                  call MP_END
                  call PPEXIT
                  stop
               endif
               call zguard(qi,nyzp,nx,inorder)
            endif
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
      msg(1:2) = time; msg(3) = tpush; msg(4) = tdpost; msg(5) = tsort
      msg(6:7) = tmove; msg(8:9) = tfft; msg(10) = tfmove
      msg(11) = trepart
      call HARTBEAT(msg,11)
      if (id0==0) then
         write (18,*) ncpus, ' processors found, ', ntasks+1, ' used'
         write (18,*) 'main max/min real time=', time(1), time(2), 'sec'
         totpush = tpush + tdpost
         write (18,*) 'electron push time=', tpush, 'sec'
         write (18,*) 'electron charge deposit time=', tdpost, 'sec'
         write (18,*) 'total electron push time=', totpush, 'sec'
         write (18,*) 'electron sort time=', tsort, 'sec'
         write (18,*) 'electron move time=', tmove, 'sec'
         totpush = totpush + tsort + tmove(1)
         write (18,*) 'total electron time=', totpush, 'sec'
      endif
      if (movion==1) then
         msg(1) = tpushi; msg(2) = tdposti; msg(3) = tsorti
         msg(4:5) = tmovi
         call HARTBEAT(msg,5)
         if (id0==0) then
            totpushi = tpushi + tdposti
            write (18,*) 'ion push time=', tpushi, 'sec'
            write (18,*) 'ion charge deposit time=', tdposti, 'sec'
            write (18,*) 'total ion push time=', totpushi, 'sec'
            write (18,*) 'ion sort time=', tsorti
            write (18,*) 'ion move time=', tmovi
            totpushi = totpushi + tsorti + tmovi(1)
            write (18,*) 'total ion time=', totpushi, 'sec'
         endif
      endif
      if (id0==0) then
         write (18,*) 'total fft time=', tfft, 'sec'
         write (18,*) 'total partition time=', tfmove, 'sec'
         write (18,*) 'total repartition time=', trepart, 'sec'
         tfmove = tfmove + trepart
         time(1) = time(1) - (totpush + totpushi + tfft(1) + tfmove)
         write (18,*) 'other time=', time(1), 'sec'
! write final diagnostic metafile
         fname = 'pdiag3.'//cdrun
         open(unit=19,file=trim(fname),form='formatted',status='replace'&
     &)
! potential diagnostics
         if (ntp > 0) then
            nprec = nprec - 1
            write (19,ppot32d,iostat=irc)
         endif
         write (19,pinput3,iostat=irc)
         if (irc /= 0) write (18,*) 'pinput3 namelist not written'
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
! allocate data for restart and/or phase space diagnostic
         if ((nts > 0) .or. (nustrt /= 1) .or. (ntr > 0)) then
            allocate(part(idimp,max(npmax,npimax),inblok))
            allocate(npp(inblok))
            allocate(edges(idps,inblok))
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
            call rddata(edges,nvp,it,ierr)
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
                  call writebf(pott,modesxp,modesyp,modesz2,kxyp,kyzp,jb&
     &lok,11,irc,trim(fname))
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
            call HARTBEAT(msg,11)
            time = msg(1:2); tpush = msg(3); tdpost = msg(4)
            tsort = msg(5); tmove = msg(6:7); tfft = msg(8:9)
            tfmove = msg(10); trepart = msg(11)
            if (id0==0) then
               write (18,*) ncpus, ' processors found, ', ntasks+1, ' us&
     &ed'
               write (18,*) 'main max/min real time=', time(1), time(2),&
     &'sec'
               totpush = tpush + tdpost
               write (18,*) 'electron push time=', tpush, 'sec'
               write (18,*) 'electron charge deposit time=', tdpost, 'se&
     &c'
               write (18,*) 'total electron push time=', totpush, 'sec'
               write (18,*) 'electron sort time=', tsort, 'sec'
               write (18,*) 'electron move time=', tmove, 'sec'
               totpush = totpush + tsort + tmove(1)
               write (18,*) 'total electron time=', totpush, 'sec'
            endif
            if (movion==1) then
               call HARTBEAT(msg,5)
               tpushi = msg(1); tdposti = msg(2); tsorti = msg(3)
               tmovi = msg(4:5)
               if (id0==0) then
                  totpushi = tpushi + tdposti
                  write (18,*) 'ion push time=', tpushi, 'sec'
                  write (18,*) 'ion charge deposit time=', tdposti, 'sec&
     &'
                  write (18,*) 'total ion push time=', totpushi, 'sec'
                  write (18,*) 'ion sort time=', tsorti
                  write (18,*) 'ion move time=', tmovi
                  totpushi = totpushi + tsorti + tmovi(1)
                  write (18,*) 'total ion time=', totpushi, 'sec'
               endif
            endif
            if (id0==0) then
               write (18,*) 'total fft time=', tfft, 'sec'
               write (18,*) 'total partition time=', tfmove, 'sec'
               write (18,*) 'total repartition time=', trepart, 'sec'
               tfmove = tfmove + trepart
               time(1) = time(1) - (totpush + totpushi + tfft(1) + tfmov&
     &e)
               write (18,*) 'other time=', time(1), 'sec'
! write final diagnostic metafile
               fname = 'pdiag3.'//cdrun
               open(unit=19,file=trim(fname),form='formatted',status='re&
     &place')
! potential diagnostics
               if (ntp > 0) then
                  nprec = nprec - 1
                  write (19,ppot32d,iostat=irc)
               endif
               write (19,pinput3,iostat=irc)
               if (irc /= 0) write (18,*) 'pinput3 namelist not written'
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
! energy diagnostic
         if (ntw > 0) then
            it = itime/ntw
            if (itime==ntw*it) then
! get energy values
               call HARTBEAT(msg,4)
               wtot = msg(1:4)
               if (id0==0) write (18,992) wtot(1), wtot(2), wtot(4)
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
               call wrdata(edges,nvp,it)
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
      end program pbeps32

