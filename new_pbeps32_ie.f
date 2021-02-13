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
! update: February 1, 2013 by Ian N. Ellis, UCLA and LLNL



program pbeps32
    use pinit32d
    use pinit32d_ie
    use pinit32d_jf
    use pespush32d_ie
    use pfield32d
    use pdiag32d
    use psimul32d_ie
    use ext_driver32_jf
    use diag32_ie
    use hdf5_write32_ie
    use sqn_ie
    use par_track_ie
    use m_h5_diagnostic_utilities, only: start_hdf5, stop_hdf5
    use m_pdiagnostic_utilities, only: create_io_comm, destroy_io_comm

    use m_system

    use mp0d, only: mpinit, ncpus
    implicit none
    include 'mpif.h'

    ! idps = number of particle partition boundaries = 4
    ! idds = dimensionality of domain decomposition = 2
    ! idimp = dimension of phase space = 6
    ! mshare = (0,1) = (no,yes) architecture is shared memory
    integer :: idps = 4, idds = 2, idimp = 6, mshare = 0
    ! nmv = number of segments in v for velocity distribution
    ! ipbc = particle boundary condition = (0,1,2,3) =
    ! (none,3d periodic,3d reflecting,mixed 2d reflecting/1d periodic)
    integer :: nmv = 40, vect = 0, ipbc = 1
    integer :: npimax = 0
    ! default unit numbers
    integer :: iuin = 8, iuot = 18, iudm = 19, iud = 12, iuv = 10, iusqn = 25
    integer :: iuvi = 20, iup = 11, ium = 21, iuer = 2
    logical :: fexist
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
    real :: qbme, qbmi, qbmsqn, affp, qi0, etx, we, wke
    real :: sx = 0.0, sy = 0.0, sz = 0.0
    real :: pxi = 0.0, pyi = 0.0, pzi = 0.0
    real :: pxe, pye, pze, wx, wy, wz
    real :: vtxi, vtyi, vtzi, vtdxi, vtdyi, vtdzi
    double precision :: dtime, ldtime
    real, dimension(:,:,:), pointer :: part, part2, parti, parti2
    real, dimension(:,:,:,:), pointer :: qe, qi, qiSQN
    real, dimension(:,:,:,:,:), pointer :: fxyze

    ! For Ian's test charge stuff
    real, dimension(:,:,:), pointer :: parteTest, partiTest
    integer, dimension(:), pointer :: nppeTest, nppiTest
    real :: qbmetest, qbmitest
    real :: keteste, ketesti

    ! for particle tracking
	 type (t_track_set) :: tracks
			! junk is an array used to pass the boundary conditions and moving window info
			! to the track create_file subroutine.  Since this code is always periodic
			! and doesn't have moving window, they are junk arrays.
	 logical, dimension(2) :: junk
    real :: write_tracks_time = 0., write_dtime = 0., temp_time=0.
    real :: add_tracks_time = 0., add_dtime = 0.


    ! stuff tracking time remaining in the run
    integer :: maxtime = 0.0, average_period_time
    character(len=20) :: inputbuf

    ! for the Strange Quark Nugget
    real, parameter :: pi = 3.14159265
    real, dimension(:,:,:), pointer :: partiSQN, parteSQN
    integer, dimension(:), pointer :: nppiSQN, nppeSQN
    real :: dxSQN, qmSQNpart, keSQN
    integer :: npiSQNmax, numSQNparti, numSQNparte

    real, dimension(:,:), allocatable :: U_sumover_x
    real, dimension(:,:,:,:,:), pointer :: grad_phi
    complex, dimension(:,:,:,:,:), pointer :: grad_phis, grad_phit !For field diagnostic
    ! for line diagnostics
    real, dimension(:,:), pointer :: lines
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
!!!!!!!!!!!!!!!!!!
! for diagnostics Ben
	real, dimension(:,:,:), pointer :: fvxx, fvyx, fvzx
	real, dimension(:,:,:), pointer :: fvxy, fvyy, fvzy
	real, dimension(:,:,:), pointer :: fvxz, fvyz, fvzz
!!!!!!!!!!!!!!!!!!!
    real, dimension(:,:), pointer :: wt
    ! wtot = total energy
    real, dimension(4) :: wtot
    ! time = timing array
    real, dimension(2) :: tfft = 0.0, time = 0.0
    real, dimension(2) :: tmove = 0.0, tmovi = 0.0
    ! msg = heartbeat array
    double precision, dimension(10) :: msg
    character(len = 10) :: cdrun
    character(len = 32) :: fname
    character(len = 12) :: label
    integer, external :: NDIAN, NDPREC, IDPREC
    991 format (' T = ', i7)
    992 format (' field, kinetic, total energies = ', 3e14.7)
    996 format (' total momentum = ', 3e14.7)
    ! get unit number for error file
    iuer = get_funit(iuer)
    ! nvp = number of real or virtual processors
    ! initialize for parallel processing
    call PPINIT(idproc, id0, nvp)
    kstrt = idproc + 1
    ! read namelist
    if (id0 == 0) then
        iuin = get_funit(iuin)
        open(unit = iuin, file = 'pinput3', form = 'formatted', status = 'old')
        read (iuin, pinput3)
        read (iuin, pinput3_jf)
        read (iuin, pinput3_ie)
    endif
    ! broadcast namelist to other nodes
    call sendnml()
    call sendnml_jf()
    call sendnml_ie()
    ! override input data
    psolve = 1
    ! set monitor flag
    call SET_MON(mpimon)
    ! create string from idrun
    write (cdrun, '(i10)') idrun
    cdrun = adjustl(cdrun)
    ! text output file
    if (id0 == 0) then
        iuot = get_funit(iuot)
        fname = 'poutput3.' // cdrun
        inquire(file=trim(fname),exist=fexist)

        if (fexist) then
            open(unit = iuot, file = trim(fname), form = 'formatted', access='sequential', status = 'unknown', position='append')
        else
            open(unit = iuot, file = trim(fname), form = 'formatted', access='sequential', status = 'unknown')
        endif

        ! for the Strange Quark Nugget
        if (rsqn > 0.0 .and. qmsqn > 0.0 .and. rmsqn > 0.0) then
            iusqn = get_funit(iusqn)
            fname = 'sqninfo.' // cdrun
            inquire(file=trim(fname),exist=fexist)

            if (fexist) then
                open(unit = iusqn, file = trim(fname),&
               &form = 'formatted', access='sequential', status = 'unknown', position='append')
            else
                open(unit = iusqn, file = trim(fname), form = 'formatted', access='sequential', status = 'unknown')
            endif

            write(iusqn,'(A7, 7A23)') 'time','x','y','z','vx','vy','vz','kinetic'
            flush(iusqn)
        endif
    endif
    ! np = total number of electrons in simulation
    !     npxyz = npx*npy*npz; npxyzb = npxb*npyb*npzb; np = npxyz + npxyzb
    npxyz = dble(npx) * dble(npy) * dble(npz)
    npxyzb = dble(npxb) * dble(npyb) * dble(npzb)
    np = npxyz + npxyzb
    ! npi = total number of ions in simulation
    !     npxyzi = npxi*npyi*npzi; npxyzbi = npxbi*npybi*npzbi
    !     npi = npxyzi + npxyzbi
    npxyzi = dble(npxi) * dble(npyi) * dble(npzi)
    npxyzbi = dble(npxbi) * dble(npybi) * dble(npzbi)
    npi = npxyzi + npxyzbi
    nx = 2**indx; ny = 2**indy; nz = 2**indz
    nxh = nx/2; nyh = ny/2; nzh = nz/2
    nyv = ny + 2; nzv = nz + 2; nxe = nx + 4
    ! nvpy/nvpz = number of real or virtual processors in y/z
    call fcomp(nvp, nx, ny, nz, nvpy, nvpz, ierr)
    nvp = nvpy * nvpz
    ! kyp = number of complex grids in each field partition in y direction
    kyp = (ny - 1)/nvpy + 1
    ! kzp = number of complex grids in each field partition in z direction
    kzp = (nz - 1)/nvpz + 1
    ! nypmx = maximum size of particle partition in y, including guard cells
    ! nzpmx = maximum size of particle partition in z, including guard cells
    nypmx = kyp + 3; nzpmx = kzp + 3
    ! ngds = number of guard cells
    ngds = 3 * ((idds - 1)/2 + 1)
    !     ax = .866025; ay = .866025; az = .866025
    if (inorder == LINEAR) then
        ax = .912871; ay = .912871; az = .912871
        nxe = nx + 2; nypmx = kyp + 1; nzpmx = kzp + 1
        ngds = (idds - 1)/2 + 1
    endif
    nxeh = nxe/2
    ! initialize for multiprocessing
    ntasks = mpinit(sntasks)
    if (popt == VECTOR) vect = 1
    ! nloop = number of time steps in simulation
    nloop = tend/dt + .0001
    ! iblok/nblok = number of particle partitions in y/z
    iblok = 1 + mshare * (nvpy - 1); nblok = 1 + mshare * (nvpz - 1)
    inblok = iblok * nblok
    ! npmax = maximum number of particles in each partition
    ! npmax = (np/dble(nvp)) * 1.5
    npmax = (np/dble(nvp)) * 2.0
    npmax = max(npmax, 2000)
    !      npmax = (np/dble(nvp))*1.25+2
    if (movion == 1) then
       npimax = (npi/dble(nvp)) * 2.0
       npimax = max(npimax, 300)
    endif
    ! if (movion == 1) npimax = (npi/dble(nvp)) * 1.25
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
    kxyb = max(kxb, kyb); kzyb = max(kyzb, kzb)
    ! kblok = number of field partitions in y direction
    kbmin = 1 + (1 - mshare)*(kxyb/kxb - 1)
    kblok = 1 + mshare * (ny/kyp - 1)
    ! jblok = number of field partitions in x direction
    jbmin = 1 + (1 - mshare)*(kxyb/kyb - 1)
    jblok = 1 + mshare * (nxh/kxyp - 1)
    ! lblok = number of field partitions in z direction
    lbmin = 1 + (1 - mshare)*(kzyb/kyzb - 1)
    lblok = 1 + mshare * (nz/kzp - 1)
    ! mblok = number of field partitions in x direction
    mbmin = 1 + (1 - mshare)*(kzyb/kzb - 1)
    mblok = 1 + mshare * (ny/kyzp - 1)
    ! nxyzh = maximum(nx,ny,nz)/2
    nxyzh = max(nx, ny, nz)/2
    ! nxhyz = maximum(nx/2,ny,nz)
    nxhyz = max(nxh, ny, nz)
    ! dimensions for index and sorting arrays
    nx1 = nx + 1; nyzpm1 = (kyp + 1)*(kzp + 1)
    ! nbmax = size of buffer for passing particles between processors
!    nbmax = 4 + (2*(npxyz*vtz + npxyzb*vtdz) + 1.4*npxyzb*abs(vdz))*dt/nz
    nbmax = 0.5*np/dble(nvp)
    nbmax = max(nbmax,2000)
    ! nbmax = 2*nbmax
    if (movion == 1) then
        vtxi = vtx/sqrt(rmass * rtempxi)
        vtyi = vty/sqrt(rmass * rtempyi)
        vtzi = vtz/sqrt(rmass * rtempzi)
    endif
    ! initialize time constants
    itime0 = 0
    itime = itime0
    ntime = itime + itime0
    ! diagnostic information needed by diagnostic nodes
    ! set default diagnostic file names
    !if (ntd > 0) fdname = 'denk32.' // cdrun
    !if (ntp > 0) fpname = 'potk32.' // cdrun
    ! energy time history
    if (ndw > 0) then
        allocate(wt((nloop - 1)/ndw - (itime0/ndw) + 1, 4))
        itw = 0
    endif

    call mkdir_structure(idproc)

    if ((nttrack < 0) .or. (ntraw < 0) .or. (ntp < 0) .or. (ndp < 0) .or. (ntd < 0) .or. (ndd < 0) &
        & .or. (ntfield < 0) .or. (nt_dump_track < 0) .or. (nphxx < 0) .or. (nphxy < 0).or. (nphxz < 0) &
        & .or. (nphyx < 0).or. (nphyy < 0).or. (nphyz < 0).or. (nphzx < 0).or. (nphzy < 0).or. (nphzz < 0)) then
       if (idproc==0) then
          write(*,'(///)')
          print*, "You think you're really clever, setting an iteration step check negative."
          print*, "But being as this simulation starts from zero and counts forward,"
          print*, 'you''ve got to ask yourself one question: "Do I feel lucky?" Well, do ya, punk?'
          write(*,'(///)')
          call MPI_ABORT( MPI_COMM_WORLD, iur1, ierr )
       endif
       goto 3000
    endif

    ! tracking stuff
    if ((nttrack .ne. 0) .or. (ntraw .ne. 0)) then
       if (track_teste .ne. 0 .and. track_testi .ne. 0) then
          if (idproc==0) then
             print*, "Error: can only have one test particle species when tracking"
          endif
          goto 3000
       endif
         
       ! must use double-precision when tracking due to the way Jay adds the tags
       ! to the particle array
       if ((track_teste .ne. 0 .or. track_testi .ne. 0 ) .and. .not. is_real_double() ) then
          if (idproc==0) then
             print*, "Error: code must be compiled with double-precision reals to use particle tracking"
          endif
          goto 3000
       endif

     	addtag = 1
    else
    	nt_dump_track = 0
    endif
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ben phasespace diagnostic grabbed from 2d beps
	if ((nphxx .ne. 0) .or. (nphxy .ne. 0) .or. (nphxz .ne. 0) .or. (ntphsl_x .ne. 0)) then
		nphbx = nphbx*2+1
		if ((nphxx .ne. 0) .or. (ntphsl_x .ne. 0)) then
			allocate(fvxx(nx,nphbx,1))
		endif
		if (nphxy .ne. 0) then
			allocate(fvxy(ny,nphbx,1))
		endif
		if (nphxz .ne. 0) then
			allocate(fvxz(nz,nphbx,1))
		endif
	endif
	
	if ((nphyx .ne. 0) .or. (nphyy .ne. 0) .or. (nphyz .ne. 0) .or. (ntphsl_x .ne. 0)) then
		nphby = nphby*2+1
		if (nphyx .ne. 0 .or. (ntphsl_x .ne. 0)) then
			allocate(fvyx(nx,nphby,1))
		endif
		if (nphyy .ne. 0) then
			allocate(fvyy(ny,nphby,1))
		endif
		if (nphyz .ne. 0) then
			allocate(fvyz(nz,nphby,1))
		endif
	endif

	if ((nphzx .ne. 0) .or. (nphzy .ne. 0) .or. (nphzz .ne. 0)) then
		nphbz = nphbz*2+1
		if (nphzx .ne. 0 .or. (ntphsl_x .ne. 0)) then
			allocate(fvzx(nx,nphbz,1))
		endif
		if (nphzy .ne. 0) then
			allocate(fvzy(ny,nphbz,1))
		endif
		if (nphzz .ne. 0) then
			allocate(fvzz(nz,nphbz,1))
		endif
	endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    


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
    allocate(part(idimp, npmax, inblok))

    ! maskp = scratch array for particle addresses
    !     allocate(maskp(npmax,inblok))
    ! in real space, qe(j+1,k,l,m) = charge density at grid point (j,kk,ll)
    ! where kk = k + noff(1,m) - 1, ll = l + noff(2,m) - 1
    allocate(qe(nxe, nypmx, nzpmx * kbmin, kblok * lblok))
    ! in real space, qi(j+1,k,l,m) = ion charge charge density
    allocate(qi(nxe, nypmx, nzpmx * kbmin, kblok * lblok))
    ! in real space, fxyze(i,j+1,k,l) = i component of force/charge at
    ! grid point (j,kk,ll)
    ! in other words, fxyze are the convolutions of the electric field
    ! over the particle shape, where kk = k + noff(1,m) - 1, and
    ! ll = l + noff(2,m) - 1
    allocate(fxyze(3, nxe, nypmx, nzpmx * kbmin, kblok * lblok))
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
    allocate(qt(nzv, kxyp, kyzp * mbmin, jblok * mblok))
    allocate(fxyzt(3, nzv, kxyp, kyzp * mbmin, jblok * mblok))
    allocate(qs(nyv, kxyp, nzpmx * jbmin * lbmin, jblok * lblok))
    allocate(fxyzs(3, nyv, kxyp, nzpmx * jbmin * lbmin, jblok * lblok))
    ! ffc = form factor array for poisson solver
    allocate(ffc(nzh, kxyp, kyzp, jblok * mblok))
    ! mixup, sct = arrays for fft
    allocate(mixup(nxhyz), sct(nxyzh))
    ! edges(1,m) = lower boundary in y of particle partition m
    ! edges(2,m) = upper boundary in y of particle partition m
    ! edges(3,m) = back boundary in z of particle partition m
    ! edges(4,m) = front boundary in z of particle partition m
    allocate(edges(idps, inblok))
    ! nyzp(1,m) = number of primary gridpoints in y in particle partition m
    ! nyzp(2,m) = number of primary gridpoints in z in particle partition m
    ! noff(1,m) = lowermost global gridpoint in y in particle partition m
    ! noff(2,m) = backmost global gridpoint in z in particle partition m
    allocate(nyzp(idds, inblok), noff(idds, inblok))
    ! npp(m) = number of particles in partition m
    ! nps(m) = starting address of particles in partition m
    allocate(npp(inblok), nps(inblok))

    ! initialize parallel timer
    call pwtimer(time, dtime, -1)
    ! initialize constants
    qbme = qme
    affp = dble(nx) * dble(ny) * dble(nz)/np

    if (movion == 1) then
        qbmi = qmi/rmass
        vtdxi = vtx/sqrt(rmass * rtempdxi)
        vtdyi = vty/sqrt(rmass * rtempdyi)
        vtdzi = vtz/sqrt(rmass * rtempdzi)
    endif
    ! set initial time
    t0 = dt*real(itime0)
    ! determine number format and default precisions
    indian = NDIAN()
    rlprec = NDPREC()
    inprec = IDPREC()
    ! calculate partition variables
    call dcomp(edges, nyzp, noff, ny, nz, kstrt, nvpy, nvpz, iblok, inorder)
    !     bxyze(1,:,:,:,:) = omx
    !     bxyze(2,:,:,:,:) = omy
    !     bxyze(3,:,:,:,:) = omz
    !     call cguard(bxyze,nyzp,nx,inorder)
    ! prepare fft tables
    call fft_init(mixup, sct, indx, indy, indz)
    ! calculate form factors
    call pois_init(ffc, ax, ay, az, affp, nx, ny, nz, kstrt, jblok)
    ! allocate ion data, we need it if we're reading in the ions during a restart
    if (movion == 1) then
        allocate(parti(idimp, npimax, inblok), nppi(inblok))
        allocate(parti2(0, 0, 0))
    endif
    
    ! for Ian's test charge stuff
    if ((track_teste>0  .or. track_testi>0) .and. nttrack > 0) then
       allocate(parteTest(idimp+1+idimp/2, nteste, inblok))
       allocate(partiTest(idimp+1+idimp/2, ntesti, inblok))
       addtag_loc = idimp + 1
    else
       allocate(parteTest(idimp, nteste, inblok))
       allocate(partiTest(idimp, ntesti, inblok))
    endif
    parteTest = 0
    partiTest = 0

    allocate(nppeTest(inblok), nppiTest(inblok))
    nppeTest = 0
    nppiTest = 0
    qbmetest = qmetest/rmteste
    qbmitest = qmitest/rmtesti

    ! for the Strange Quark Nugget
    if (rsqn > 0.0 .and. qmsqn > 0 .and. rmsqn > 0) then
        ! simple checks
        if (rsqn>nxh .or. rsqn>nyh .or. rsqn>nzh) then
            if(idproc==0) print*,'Error: Strange Quark Nugget too big'
            call MP_END
            call PPEXIT
            stop
        endif

        if (inorder/=QUADRATIC) then
            if(idproc==0) print*,'Error: Only quadratic splines are supported with the Strange Quark Nugget.'
            call MP_END
            call PPEXIT
            stop
        endif

        if (popt/=LOOKAHEAD) then
            if(idproc==0) print*,'Error: Only the lookahead particle optimization scheme is supported with the SQN.'
            call MP_END
            call PPEXIT
            stop
        endif

        if (relativity == 1) then
            if(idproc==0) print*,'Error: Special relativity is not supported with the Strange Quark Nugget.'
            call MP_END
            call PPEXIT
            stop
        endif

        if (ipbc /= 1) then
            if(idproc==0) print*,'Error: Only periodic boundary conditions are supported with the Strange Quark Nugget.'
            call MP_END
            call PPEXIT
            stop
        endif

        allocate(nppiSQN(inblok), nppeSQN(inblok))
        nppiSQN = 0

        ! The SQN is implemented as a uniform blob of positive charges with
        ! a number density of the number of electrons per cell, so we need the
        ! particle separation.
        dxSQN = (affp)**(1.0/3.0)
        numSQNparti = 4.0/3.0 * pi * rsqn**3/affp
        
        qbmsqn = qmsqn/rmsqn

        ! estimate maximum number of particles in a partition--this is very approximate
        if ( (nyzp(1,1) + 1) >= rsqn .and. (nyzp(2,1) + 1) >= rsqn) then
            npiSQNmax = numSQNparti*1.1
        else
            npiSQNmax = ( 2 * rsqn * (nyzp(1,1) + 1) * (nyzp(2,1) + 1) ) / affp + 1
        endif

        if(idproc==0) then
            print*, 'numSQNpart = ', numSQNparti, 'npiSQNmax = ', npiSQNmax, 'dxSQN = ', dxSQN
            print*, 'nyp = ', nyzp(1,1), 'nzp = ', nyzp(2,1)
        endif
        allocate(partiSQN(idimp,npiSQNmax,inblok), parteSQN(idimp,npiSQNmax,inblok))

        ! in real space, qiSQN(j+1,k,l,m) = SQN charge charge density
        allocate(qiSQN(nxe, nypmx, nzpmx * kbmin, kblok * lblok))

    endif

    ! start HDF5 and setup the I/O communicator
    !call start_hdf5(ierr)
    call create_io_comm()

    ! new start
    if (nustrt == 1) then
        nustrt = 0
        if (idproc==0) then
            print*, 'nustrt = 1, setting to 0.  No worries, we will automatically detect restart files.'
        endif
    endif

    if (.not. restart_exists()) then
    !if (0 /= 1) then
        ! initialize electrons
        nps = 1
        npp = 0

        ! test charge stuff
        if (id0 == 0) then
            parteTest(1:idimp, 1:nteste, inblok) = reshape(qetest, shape(parteTest(1:idimp, 1:nteste, inblok)))
            partiTest(1:idimp, 1:ntesti, inblok) = reshape(qitest, shape(partiTest(1:idimp, 1:ntesti, inblok)))

            ! scale the test charge positions and velocities
            ! makes for far less input file editing when shrinking cell widths
            if (scale_test == 1) then
               parteTest(1, 1:nteste, inblok) = parteTest(1, 1:nteste, inblok)*vtx
               parteTest(2, 1:nteste, inblok) = parteTest(2, 1:nteste, inblok)*vty
               parteTest(3, 1:nteste, inblok) = parteTest(3, 1:nteste, inblok)*vtz
               parteTest(4, 1:nteste, inblok) = parteTest(4, 1:nteste, inblok)*vtx
               parteTest(5, 1:nteste, inblok) = parteTest(5, 1:nteste, inblok)*vty
               parteTest(6, 1:nteste, inblok) = parteTest(6, 1:nteste, inblok)*vtz

               partiTest(1, 1:ntesti, inblok) = partiTest(1, 1:ntesti, inblok)*vtx
               partiTest(2, 1:ntesti, inblok) = partiTest(2, 1:ntesti, inblok)*vty
               partiTest(3, 1:ntesti, inblok) = partiTest(3, 1:ntesti, inblok)*vtz
               partiTest(4, 1:ntesti, inblok) = partiTest(4, 1:ntesti, inblok)*vtx
               partiTest(5, 1:ntesti, inblok) = partiTest(5, 1:ntesti, inblok)*vty
               partiTest(6, 1:ntesti, inblok) = partiTest(6, 1:ntesti, inblok)*vtz
            endif
 
            nppeTest = nteste;
            nppiTest = ntesti;
            if (addtag .ne. 0) then
               if (track_teste .ne. 0) then
			         call assign_tags(parteTest,nppeTest,idproc)
               elseif (track_testi .ne. 0) then
			         call assign_tags(partiTest,nppiTest,idproc)
               endif
		      endif
         endif

         if (addtag .ne. 0) then
            ! store particle tracking data
			   tracks%niter = nttrack

            if (track_teste .ne. 0) then
               call setup(tracks, nt_dump_track, .false., parteTest, nteste)
            elseif (track_testi .ne. 0) then
               call setup(tracks, nt_dump_track, .false., partiTest, ntesti)
            endif

            call create_file(tracks,'testpart',nt_dump_track,nx,ny,nz,dt,junk,junk,relativity,idproc,1)

		   endif

			
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
            call vfdistr(part, nps, ampdx, scaledx, shiftdx, ampdy, scaledy, sh&
            &iftdy, ampdz, scaledz, shiftdz, npx, npy, npz, nx, ny, nz, kstrt, nvp, ipbc, nd&
            &prof, nsrand, iblok)
            call vvdistr(part, npp, nps, vtx, vty, vtz, vx0, vy0, vz0, npx, npy, np&
            &z, kstrt, nvp, iblok)

            ! print*,"after init"
        endif

        ! beam electrons
        nps = npp + 1
        !        if (npxyzb > 0) call distr(part,edges,npp,nps,vtdx,vtdy,vtdz,vd&
        !    &x,vdy,vdz,npxb,npyb,npzb,nx,ny,nz,ipbc,iblok)

        if (npxyzb > 0) then
            call fdistr(part, nps, ampdx, scaledx, shiftdx, ampdy, scaledy, shi&
            &ftdy, ampdz, scaledz, shiftdz, npxb, npyb, npzb, nx, ny, nz, kstrt, nvp, ipbc, &
            &ndprof, nsrand, iblok)
            call vdistr(part, npp, nps, vtdx, vtdy, vtdz, vdx, vdy, vdz, npxb, npy&
            &b, npzb, kstrt, nvp, iblok)
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
        call pmove(part, edges, npp, tmove, ny, nz, kstrt, nvpy, nvpz, nbmax, idd&
        &s, iblok, vect, ierr)
        ! move Ian's test electron charges into appropriate spatial regions
        if (track_teste .ne. 0) then
            call pmove(parteTest, edges, nppeTest, tmove, ny, nz, kstrt, nvpy, nvpz, nteste, idd&
           &s, iblok, vect, ierr,tracks)
        else
           call pmove(parteTest, edges, nppeTest, tmove, ny, nz, kstrt, nvpy, nvpz, nteste, idd&
           &s, iblok, vect, ierr)
        endif

        ! calculate initial electron momentum
        if (ntm > 0) call initmomt3(part, npp, pxe, pye, pze)
        
        !SQN
        if (rsqn > 0.0 .and. qmsqn > 0 .and. rmsqn > 0) then
            ! initialize the Strange Quark Nugget
            call SQNinit(sqn, rsqn, dxSQN, edges, nppiSQN, partiSQN, npiSQNmax)
            
            ! neutralizing SQN charge
            nppeSQN(:) = nppiSQN(:)
            parteSQN(:,:,:) = partiSQN(:,:,:)

            ! check number of particles
            call MPI_REDUCE(nppiSQN(1),numSQNparti, 1, MPI_INTEGER, MPI_SUM,0,MPI_COMM_WORLD,ierr)
            if (id0==0) print*, 'Total SQN particles = ', numSQNparti
            
            ! determine the charge on each SQN particle
            qmSQNpart = qmsqn/numSQNparti
        endif

        ! initialize ions
        if (movion == 1) then
            nps = 1
            nppi = 0
            ! background ions
            !           if (npxyzi > 0) call distr(parti,edges,nppi,nps,vtxi,vtyi,vt&
            !    &zi,vxi0,vyi0,vzi0,npxi,npyi,npzi,nx,ny,nz,ipbc,iblok)
            if (npxyzi > 0) then
!                call fdistr(parti, nps, ampdxi, scaledxi, shiftdxi, ampdyi, sca&
!                &ledyi, shiftdyi, ampdzi, scaledzi, shiftdzi, npxi, npyi, npzi, nx, ny, nz, ks&
!                &trt, nvp, ipbc, ndprofi, nsrandi, iblok)
!                call vdistr(parti, nppi, nps, vtxi, vtyi, vtzi, vxi0, vyi0, vzi0, &
!                &npxi, npyi, npzi, kstrt, nvp, iblok)
                call vfdistr(parti,nps,ampdxi,scaledxi,shiftdxi,ampdyi,sc&
                &aledyi,shiftdyi,ampdzi,scaledzi,shiftdzi,npxi,npyi,npzi,nx,ny,nz,k&
                &strt,nvp,ipbc,ndprofi,nsrandi,iblok)
                call vvdistr(parti,nppi,nps,vtxi,vtyi,vtzi,vxi0,vyi0,vzi0&
                &,npxi,npyi,npzi,kstrt,nvp,iblok)
            endif
            ! beam ions
            nps = nppi + 1
            !           if (npxyzbi > 0) call distr(parti,edges,nppi,nps,vtdxi,vtdyi&
            !    &,vtdzi,vdxi,vdyi,vdzi,npxbi,npybi,npzbi,nx,ny,nz,ipbc,iblok)
            if (npxyzbi > 0) then
                call fdistr(parti, nps, ampdxi, scaledxi, shiftdxi, ampdyi, sca&
                &ledyi, shiftdyi, ampdzi, scaledzi, shiftdzi, npxbi, npybi, npzbi, nx, ny, nz&
                &, kstrt, nvp, ipbc, ndprofi, nsrandi, iblok)
                call vdistr(parti, nppi, nps, vtdxi, vtdyi, vtdzi, vdxi, vdyi, vd&
                &zi, npxbi, npybi, npzbi, kstrt, nvp, iblok)
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
            call pmove(parti, edges, nppi, tmovi, ny, nz, kstrt, nvpy, nvpz, nbma&
            &x, idds, iblok, vect, ierr)

            ! Ian's test ion stuff
            if (track_testi .ne. 0) then
               call pmove(partiTest, edges, nppiTest, tmovi, ny, nz, kstrt, nvpy, nvpz, ntesti, &
               &idds, iblok, vect, ierr, tracks)
            else
               call pmove(partiTest, edges, nppiTest, tmovi, ny, nz, kstrt, nvpy, nvpz, ntesti, &
               &idds, iblok, vect, ierr)
            endif

            ! calculate initial ion momentum
            if (ntm > 0) call initmomt3(parti, nppi, pxi, pyi, pzi)
        endif

        ! initialize background charge density
        if (movion == 0) then
            if (smoothIon==1) then ! for smooth neutralizing ion background
                qi0 = -qme/affp
                call sguard(qi, nyzp, qi0, nx, inorder)
            else
                call sguard(qi, nyzp, zero, nx, inorder)
                call dpost(part, qi, -qme, npp, noff, tdpost, inorder, d&
                &opt)

                ! dump the positions now
                if (ntr>0) call ion_write(idrun, npp, part, nppiTest, partiTest, itime, itime0)
                
            endif
            ! Ian's test charge stuff
            ! put the test ion in the appropriate spatial region
            call pmove(partiTest, edges, nppiTest, tmovi, ny, nz, kstrt, nvpy, n&
            &vpz, ntesti, idds, iblok, vect, ierr)
            ! deposit the test ion charge
            call dpost(partiTest, qi, qmitest, nppiTest, noff, tdposti, inorder, d&
            &opt)
            deallocate(partiTest, nppiTest)
            call aguard(qi, nyzp, nx, inorder)
            call paguard(qi, kstrt, nvpy, nvpz, nx, kyp, kzp, ngds, iblok, inorde&
            &r)
            
            ! freeze the ions now   !!Ian is confused about why this is here.....
        else if ((movion == 1).and.(ntime == ionoff)) then
            ! initialize ion charge density to zero
            call sguard(qi, nyzp, zero, nx, inorder)

            ! Ian's test charge stuff
            ! make sure test ion is in the right place....ntime should be 0 here, but just in case....
            if (fixedvel==1) then
               partiTest(1:3, 1:nppiTest(1), 1) = partiTest(1:3, 1:nppiTest(1), 1) + &
               &ntime*dt * partiTest(4:6, 1:nppiTest(1), 1)
            endif
            ! put the test ion in the appropriate spatial region
            call pmove(partiTest, edges, nppiTest, tmovi, ny, nz, kstrt, nvpy, n&
            &vpz, ntesti, idds, iblok, vect, ierr)
            ! deposit the test ion charge
            call dpost(partiTest, qi, qmitest, nppiTest, noff, tdposti, inorder, dopt)

            ! deposit ion charge
            call dpost(parti, qi, qmi, nppi, noff, tdposti, inorder, dopt)
            ! add guard cells for density in x
            call aguard(qi, nyzp, nx, inorder)
            ! add guard cells for density in y and z
            call paguard(qi, kstrt, nvpy, nvpz, nx, kyp, kzp, ngds, iblok, inorde&
            &r)

            ! dump the ion positions now
            if (ntr>0) call ion_write(idrun, npp, part, nppiTest, partiTest, itime, itime0)
            
            ! delete ions
            deallocate(parti, nppi, partiTest, nppiTest)
            movion = 0
        endif

    ! restart
    else
        ! Exit if freeze ions.  It's probably easy to implement, but I don't have any use for it.
        if (ionoff>0) then
            print*, 'Error, Ian''s restarts do not support ionoff>0'
            call MP_END
            call PPEXIT
            stop
        endif

        if (movion==0 .and. smoothIon==0) allocate(parti(idimp, npmax, inblok), nppi(inblok))

        if (rsqn > 0.0 .and. qmsqn > 0 .and. rmsqn > 0) then
            ! initialize the neutralizing Strange Quark Nugget
            call SQNinit(sqn, rsqn, dxSQN, edges, nppeSQN, parteSQN, npiSQNmax)
        endif

        ! read restart
        if (nustrt==0) then
            ! read and extend run

            ! particle tracks do need memory....unfortunately
            if (addtag .ne. 0) then
               ! store particle tracking data
			      tracks%niter = nttrack

               if (track_teste .ne. 0) then
                  call setup(tracks, nt_dump_track, .true., parteTest, nteste)
               elseif (track_testi .ne. 0) then
                  call setup(tracks, nt_dump_track, .true., partiTest, ntesti)
               endif

               call restart_read(idrun, npp, part, movion, nppi, parti, sqn, nppiSQN, partiSQN,&
                     & nppeTest, parteTest, nppiTest, partiTest, smoothIon, itime, itime0, itw, wt, tracks)

            else

               call restart_read(idrun, npp, part, movion, nppi, parti, sqn, nppiSQN, partiSQN,&
                     & nppeTest, parteTest, nppiTest, partiTest, smoothIon, itime, itime0, itw, wt)

            endif

            itime0 = itime + itime0
        else
            ! use restart data as new initial conditions
            call restart_read(idrun, npp, part, movion, nppi, parti, sqn, nppiSQN, partiSQN, &
                  & nppeTest, parteTest, nppiTest, partiTest, smoothIon)
            itime0 = 0
            
            ! a particularly ambitious geek could make this work with test particles

        endif

        t0 = dt*real(itime0)
        itime = 0
        ntime = itime + itime0

        !SQN checks
        if (rsqn > 0.0 .and. qmsqn > 0 .and. rmsqn > 0) then
            ! check number of particles
            call MPI_REDUCE(nppiSQN(1),numSQNparti, 1, MPI_INTEGER, MPI_SUM,0,MPI_COMM_WORLD,ierr)
            if (id0==0) print*, 'Total SQN particles = ', numSQNparti
            call MPI_REDUCE(nppeSQN(1),numSQNparte, 1, MPI_INTEGER, MPI_SUM,0,MPI_COMM_WORLD,ierr)
            if (id0==0) print*, 'Total neutralizing SQN particles = ', numSQNparte

            if(numSQNparti /= numSQNparte) then
                print*, 'FATAL Error: number old SQN particles and new neutralizing SQN particles differ!'
                call MP_END
                call PPEXIT
                stop
            endif

            ! determine the charge on each SQN particle
            qmSQNpart = qmsqn/numSQNparti
        endif

        if (movion==0) then
            if (smoothIon==1) then ! for smooth neutralizing ion background
                qi0 = -qme/affp
                call sguard(qi, nyzp, qi0, nx, inorder)
                ! put the test ions in the appropriate spatial regions b/c won't have read in from restarts
                
                if (ntesti > 0) then 
                   if(id0==0) then
                     partiTest(1:idimp, 1:ntesti, inblok) = reshape(qitest, shape(partiTest(1:idimp, 1:ntesti, inblok)))
                   endif

                   call pmove(partiTest, edges, nppiTest, tmovi, ny, nz, kstrt, nvpy, n&
                   &vpz, ntesti, idds, iblok, vect, ierr)
                endif
            else
                call sguard(qi, nyzp, zero, nx, inorder)
                call dpost(parti, qi, -qme, nppi, noff, tdpost, inorder, d&
                &opt)
                ! free up memory
                deallocate(parti, nppi)
            endif
            ! deposit the test ion charge
            call dpost(partiTest, qi, -qme, nppiTest, noff, tdposti, inorder, dopt)
            ! free up memory
            deallocate(partiTest, nppiTest)
            call aguard(qi, nyzp, nx, inorder)
            call paguard(qi, kstrt, nvpy, nvpz, nx, kyp, kzp, ngds, iblok, inorde&
            &r)
        endif
    endif

    ! sorting arrays
490 allocate(pt(max(npmax, npimax), inblok))
    allocate(ip(max(npmax, npimax), inblok), npic(nyzpm1, inblok))
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
    if (id0 == 0) then
        iudm = get_funit(iudm)
        fname = 'pdiag32.init.' // cdrun
        open(unit = iudm, file = trim(fname), form = 'formatted', status = 'replace')
    endif

    call write_jf_diag_file(idproc, cdrun)

    if (ntlines .ne. 0) then
        do i = 1, 5
            if (linepos_3d(i, 1) .ne. - 1) then
                nlines = nlines + 1
            endif
        enddo
        allocate(nlinerec(nlines))
        allocate(lines(nlines, nx))
        lines = 0.
        nlinerec = 0.
    endif

    if (nt_write_U_sumover_x .ne. 0) then
        allocate(grad_phi(3, nxe, nypmx, nzpmx, 1))
        allocate(grad_phis(3, nyv, kxyp, nzpmx * jbmin * lbmin, jblok * lblok))
        allocate(grad_phit(3, nzv, kxyp, kyzp * mbmin, jblok * mblok))
        allocate(U_sumover_x(ny, nz))
    endif
    if ((nt_write_U_sumover_x .ne. 0)) then
        if (idproc == 0) then
            fname = './DIAG/U_sumover_x_file'
            nt_write_U_sumover_x_funit = get_funit(20)
            open(unit = nt_write_U_sumover_x_funit, file = trim(fname), form = 'unformatted', status = 'replace')
            write (nt_write_U_sumover_x_funit) ny, nz, nt_write_U_sumover_x, dt, (timerise + timeflat + timefall), tend
        endif
    endif
    ! ion density or potential diagnostics
    if ((ntp > 0) .or. (ndp > 0) .or. (ntd > 0) .or. (ndd > 0) .or. (ntfield > 0)) then
        allocate(sfield(nxe, nypmx, nzpmx * kbmin, kblok * lblok))
        allocate(sfieldt(nzv, kxyp, kyzp * mbmin, jblok * mblok))
    endif
    ! ion density diagnostics
    call initmodediag(dent, ntd, id0, nxh, nyh, nzh, kxyp, kyzp, modesxd, &
    &modesyd, modeszd, jblok, mblok, iud, ndrec, fdname)
    if (ntd > 0) then
        if (id0 == 0) then
            ceng = zero
            write (iudm, pden32d, iostat = irc)
        endif
    endif
    ! velocity diagnostics
    fname = 'fv3.' // cdrun
    call initveldiag(fv, fvm, vtx, vty, vtz, ntv, ndv, id0, nmv, nblok, iuv, fnam&
    &e)
    if (movion == 1) then
        fname = 'fvi3.' // cdrun
        call initveldiag(fvi, fvmi, vtxi, vtyi, vtzi, ntv, ndv, id0, nmv, nblok, &
        &iuvi, fname)
    endif
    ! potential diagnostics
    call initmodediag(pott, ntp, id0, nxh, nyh, nzh, kxyp, kyzp, modesxp, &
    &modesyp, modeszp, jblok, mblok, iup, nprec, fpname)
    if (ntp > 0) then
        if (id0 == 0) then
            ceng = zero
            write (iudm, ppot32d, iostat = irc)
        endif
    endif
    ! momentum diagnostics
    fname = 'pmomentum3.' // cdrun
    if (ntm > 0) then
        if (id0 == 0) then
            ium = get_funit(ium)
            open(unit = ium, file = trim(fname), form = 'formatted', status = 'unkn&
            &own')
        endif
    endif
    ! write out and close input file
    if (id0 == 0) then
        write (iudm, pinput3, iostat = irc)
        close(unit = iudm)
    endif

    ! get final time of the run from the command line
    
    if (getargc()>0) then
        call getarg(1,inputbuf)
        read(inputbuf,*) maxtime
    endif



    ! record time
    call pwtimer(time, dtime)
    ! send initial CPU Time to diagnostic nodes
    msg(1) = time(1); msg(2) = time(2)
    call HARTBEAT(msg, 2)


    

    call wtimer(ltime, ldtime, -1)
    if (id0 == 0) then
        write (iuot, *) 'init max/min real time=', time(1), time(2), 'sec'
    endif
    !
    ! * * * start main iteration loop * * *
    !
500 if (nloop < ntime) go to 2000
    ! send time step to diagnostic nodes
    msg(1) = ntime
    call HARTBEAT(msg, 1)
    if (id0 == 0) then
        write (iuot, 991) ntime
        flush(iuot) ! required to ever see anything on BGL
    endif
    write (label, 991) ntime
    call LOGNAME(label)
    ! deposit electron charge density
    ! initialize electron charge density to zero
    call sguard(qe, nyzp, zero, nx, inorder)

    ! deposit electron charge
    call dpost(part, qe, qme, npp, noff, tdpost, inorder, dopt)

    ! deposit Ian's test electron charge
    call dpost(parteTest, qe, qmetest, nppeTest, noff, tdpost, inorder, dopt)

    ! add guard cells for density in x
    call aguard(qe, nyzp, nx, inorder)
    ! add guard cells for density in y and z
    call paguard(qe, kstrt, nvpy, nvpz, nx, kyp, kzp, ngds, iblok, inorder)

    if (rsqn > 0.0 .and. qmsqn > 0.0 .and. rmsqn > 0.0) then
        ! deposit SQN charge density

        ! initialize SQN charge density to zero
        call sguard(qiSQN, nyzp, zero, nx, inorder)

        ! deposit SQN charge
        call dpost(partiSQN, qiSQN, qmSQNpart, nppiSQN, noff, tdposti, inorder, dopt)
        call dpost(parteSQN, qiSQN, -qmSQNpart, nppeSQN, noff, tdpost, inorder, dopt)

        ! add guard cells for density in x
        call aguard(qiSQN, nyzp, nx, inorder)
        ! add guard cells for density in y and z
        call paguard(qiSQN, kstrt, nvpy, nvpz, nx, kyp, kzp, ngds, iblok, inorder)
    endif

    ! deposit ion charge density
    if (movion == 1) then
        ! initialize ion charge density to zero
        call sguard(qi, nyzp, zero, nx, inorder)
        ! deposit ion charge
        call dpost(parti, qi, qmi, nppi, noff, tdposti, inorder, dopt)
        !call dpost(parti, qi, -qme, nppi, noff, tdposti, inorder, dopt)
        ! deposit Ian's test ion charge
        call dpost(partiTest, qi, qmitest, nppiTest, noff, tdposti, inorder, dopt)
        ! add guard cells for density in x
        call aguard(qi, nyzp, nx, inorder)
        ! add guard cells for density in y and z
        call paguard(qi, kstrt, nvpy, nvpz, nx, kyp, kzp, ngds, iblok, inorder)
        ! freeze the ions
        if (ntime == ionoff) then
            deallocate(parti, nppi)
            movion = 0
        endif
    endif

    ! charge density diagnostics
    if (ntd > 0) then
        it = ntime/ntd
        if (ntime == ntd * it) then
            !electron charge density
            call dendiag(qt, qs, qe, sfield, dent, sfieldt, ffc, nyzp, mixup, sct, tfft, &
                &ntd, ntd, nx, ny, nz, modesxd, modesyd, modeszd, iud, ndrec, indx, indy, indz, &
                &ntime, nvpy, nvpz, kstrt, kxyp, kyp, kyzp, kzp, ngds, iblok, jblok, kblok, &
                &mblok, inorder)
            call writef(sfield, nxe, nypmx, nzpmx, nvpy, nvpz, ntime, dt, &
            & E_CHARGE, inorder, write_stride)

            !ion charge density
            if (ntime == itime0 .or. movion == 1 .or. ntesti > 0) then 
               call dendiag(qt, qs, qi, sfield, dent, sfieldt, ffc, nyzp, mixup, sct, tfft, &
                   &ntd, ntd, nx, ny, nz, modesxd, modesyd, modeszd, iud, ndrec, indx, indy, indz, &
                   &ntime, nvpy, nvpz, kstrt, kxyp, kyp, kyzp, kzp, ngds, iblok, jblok, kblok, &
                   &mblok, inorder)
               call writef(sfield, nxe, nypmx, nzpmx, nvpy, nvpz, ntime, dt, &
               & I_CHARGE, inorder, write_stride)
            endif

!            if (rsqn > 0.0 .and. qmsqn > 0.0 .and. rmsqn > 0.0) then
!                !SQN charge density
!                call dendiag(qt, qs, qiSQN, sfield, dent, sfieldt, ffc, nyzp, mixup, sct, tfft, &
!                    &ntd, ndd, nx, ny, nz, modesxd, modesyd, modeszd, iud, ndrec, indx, indy, indz, &
!                    &ntime, nvpy, nvpz, kstrt, kxyp, kyp, kyzp, kzp, ngds, iblok, jblok, kblok, &
!                    &mblok, inorder)
!                fname = 'qiSQN#id' // cdrun
!                call writef(sfield, nxe, nypmx, nzpmx, nvpy, nvpz, ntime, ntime * dt, &
!                &SQN_CHARGE, trim(fname), inorder)
!            endif
        endif
    endif

    ! add electron and ion densities
    call addqei(qe, qi, nyzp, nx, inorder)

    if (rsqn > 0.0 .and. qmsqn > 0.0 .and. rmsqn > 0.0) then
        ! add in SQN charge density
        call addqei(qe, qiSQN, nyzp, nx, inorder)
    endif

    ! ion density diagnostic
!    call dendiag(qt, qs, qi, sfield, dent, sfieldt, ffc, nyzp, mixup, sct, tfft, &
!    &ntd, ndd, nx, ny, nz, modesxd, modesyd, modeszd, iud, ndrec, indx, indy, indz, &
!    &ntime, nvpy, nvpz, kstrt, kxyp, kyp, kyzp, kzp, ngds, iblok, jblok, kblok, &
!    &mblok, inorder)

    ! velocity diagnostic
    call veldiag(part, fv, fvm, npp, ntv, ndv, id0, nmv, iuv, ntime)
    if (movion == 1) then
        call veldiag(parti, fvi, fvmi, nppi, ntv, ndv, id0, nmv, iuvi, ntime)
    endif
    ! transform charge to fourier space
    isign = -1
    call fft(qe, qs, qt, isign, mixup, sct, tfft, indx, indy, indz, kstrt, kxyp, k&
    &yp, kyzp, kzp, kblok, mblok, inorder)
    ! potential diagnostic
    if (ntp > 0) then
        it = ntime/ntp
        if (ntime == ntp * it) then
            call potdiag(qt, qs, sfield, pott, sfieldt, ffc, nyzp, mixup, sct, tfft, ntp&
            &, ntp, nx, ny, nz, modesxp, modesyp, modeszp, iup, nprec, indx, indy, indz, &
            &ntime, nvpy, nvpz, kstrt, kxyp, kyp, kyzp, kzp, ngds, iblok, jblok, kblok, &
            &mblok, inorder)
            call writef(sfield, nxe, nypmx, nzpmx, nvpy, nvpz, ntime, dt, &
            &POT, inorder, write_stride)
        endif
    endif

    ! calculate force/charge in fourier space
    call pois(qt, fxyzt, ffc, we, nx, ny, nz, kstrt, jblok)
    ! transform force/charge to real space
    isign = 1
    call fft(fxyze, fxyzs, fxyzt, isign, mixup, sct, tfft, indx, indy, indz, kst&
    &rt, kxyp, kyp, kyzp, kzp, kblok, mblok, inorder)
    !     call fftn(fxyze,fxyzs,fxyzt,isign,mixup,sct,tfft,indx,indy,indz,ks&
    !    &trt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
    ! copy data from field to particle partition, and copy to guard cells
    call pcguard(fxyze, kstrt, nvpy, nvpz, kyp, kzp, ngds, iblok, inorder)
    call cguard(fxyze, nyzp, nx, inorder)
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
        call plane_wave(fxyze, real(ntime) * dt, nx, nxe, nvp, idproc)
    case (4)
        call laguerre_gaussian_ponderomotive_force(fxyze, real(ntime) * dt, nx, nxe, ny, nypmx, nz, nzpmx, nvpy, nvpz, idproc, 1)
    case (5)
        call two_laguerre_gaussian_beams(fxyze,real(ntime) * dt, nx,nxe,ny,nypmx,nz,nzpmx,nvpy,nvpz,idproc,l_number1,l_number2)
    case (7)
        call gauss_tran_per_wavelen(fxyze, real(ntime) * dt, nx, nxe, ny, nypmx, nz, nzpmx, nvpy, nvpz, idproc)
    case (12)
        call supergauss_tran_per_wavelen(fxyze, real(ntime) * dt, nx, nxe, ny, nypmx, nz, nzpmx, nvpy, nvpz, idproc)
    case default
        print*, "driver_select = ", driver_select, " is not yet implemented.  Exiting..."
        call MP_END
        call PPEXIT
        stop
    end select

    if ((dump_start < 0.) .OR. &
        &((real(ntime) * dt >= dump_start) .and. (real(ntime) * dt <= dump_end))) then
        if (ntfield > 0) then
            it = ntime/ntfield
            if (ntime == ntfield * it) then

                sfield = fxyze(1,:,:,:,:)
                call writef(sfield,nxe,nypmx,nzpmx,nvpy,nvpz,ntime,dt,&
                	&EX,inorder,write_stride)

                sfield = fxyze(2,:,:,:,:)
                call writef(sfield,nxe,nypmx,nzpmx,nvpy,nvpz,ntime,dt,&
                	&EY,inorder,write_stride)

                sfield = fxyze(3,:,:,:,:)
                call writef(sfield,nxe,nypmx,nzpmx,nvpy,nvpz,ntime,dt,&
                	&EZ,inorder,write_stride)
            endif
        endif
    endif

    ! update tracks
    if ( nttrack > 0) then
		it = ntime / nttrack
		if ( ntime == nttrack*it) then
			 temp_time = 0.
          call PWTIMERA(-1,temp_time,add_dtime)
          if (track_teste .ne. 0) then
            call add_track_data( tracks, parteTest, ntime, real(ntime)*dt, qmetest, rmteste, relativity )
          elseif (track_testi .ne. 0) then
            call add_track_data( tracks, partiTest, ntime, real(ntime)*dt, qmitest, rmtesti, relativity )
          endif
          call PWTIMERA(1,temp_time,add_dtime)
          add_tracks_time = add_tracks_time + temp_time
		endif
	 endif

    !Write track data to disk
		if ( nt_dump_track > 0) then
			it = ntime / nt_dump_track
			if ( ntime == nt_dump_track*it) then
				temp_time = 0.
				call PWTIMERA(-1,temp_time,write_dtime)
				call write_tracks( tracks, nvp, idproc )
				call PWTIMERA(1,temp_time,write_dtime)
				write_tracks_time = write_tracks_time + temp_time
			endif
		endif


    if (nt_write_U_sumover_x > 0) then
        it = ntime/nt_write_U_sumover_x
        if (ntime == nt_write_U_sumover_x * it) then

            call pois(qt, sfieldt, ffc, we, nx, ny, nz, kstrt, jblok)
            call ipgradf32(sfieldt, grad_phit, nx, ny, nz, kstrt, jblok)
            isign = 1
            call fft(grad_phi, grad_phis, grad_phit, isign, mixup, sct, tfft, indx, indy, indz, &
            &kstrt, kxyp, kyp, kyzp, kzp, kblok, mblok, inorder)
            call pcguard(grad_phi, kstrt, nvpy, nvpz, kyp, kzp, ngds, iblok, inorder)
            call cguard(grad_phi, nyzp, nx, inorder)

            U_sumover_x = 0.

            call get_3D_sumover_x(0.5 * (grad_phi**2), U_sumover_x, nx, nxe, ny, nypmx, nz, nzpmx, nvpy, &
            & nvpz, idproc, inorder)
            if (idproc == 0) then
                write (nt_write_U_sumover_x_funit) U_sumover_x
            endif

        endif
    endif

    ! field lines diagnostic added by JF
    			!if (ntlines > 0) then
    			!	it = ntime / ntlines
    			!	if (ntime == ntlines*it) then
    			!		call write_Eoft_aty(lines,fxyze,nx,ny,nz,nvpy,nvpz,idproc)
    			!	endif
    			!endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Ben phasespace diagnostic grabbed from 2d beps
! Summed over transverse dim phase space diagnostic added by JF
		if (nphxx > 0) then
			it = ntime / nphxx
			if (ntime==nphxx*it) then
				fvxx = 0.
				call phase_space_vxyz_vs_x(part,fvxx,1,nx,npp,dt,ntime,nblok,idproc)
			endif
		endif
		if (nphyx > 0) then
			it = ntime / nphyx
			if (ntime==nphyx*it) then
				fvyx = 0.
				call phase_space_vxyz_vs_x(part,fvyx,2,nx,npp,dt,ntime,nblok,idproc)
			endif
		endif
		if (nphzx > 0) then
			it = ntime / nphzx
			if (ntime==nphzx*it) then
				fvzx = 0.
				call phase_space_vxyz_vs_x(part,fvzx,3,nx,npp,dt,ntime,nblok,idproc)
			endif
		endif
		if (nphxy > 0) then
			it = ntime / nphxy
			if (ntime==nphxy*it) then
				fvxy = 0.
				call phase_space_vxyz_vs_y(part,fvxy,1,ny,npp,dt,ntime,nblok,idproc)
			endif
		endif
		if (nphyy > 0) then
			it = ntime / nphyy
			if (ntime==nphyy*it) then
				fvyy = 0.
				call phase_space_vxyz_vs_y(part,fvyy,2,ny,npp,dt,ntime,nblok,idproc)
			endif
		endif
		if (nphzy > 0) then
			it = ntime / nphzy
			if (ntime==nphzy*it) then
				fvzy = 0.
				call phase_space_vxyz_vs_y(part,fvzy,3,ny,npp,dt,ntime,nblok,idproc)
			endif
		endif
		if (nphxz > 0) then
			it = ntime / nphxz
			if (ntime==nphxz*it) then
				fvxz = 0.
				call phase_space_vxyz_vs_z(part,fvxz,1,nz,npp,dt,ntime,nblok,idproc)
			endif
		endif
		if (nphyz > 0) then
			it = ntime / nphyz
			if (ntime==nphyz*it) then
				fvyz = 0.
				call phase_space_vxyz_vs_z(part,fvyz,2,nz,npp,dt,ntime,nblok,idproc)
			endif
		endif
		if (nphzz > 0) then
			it = ntime / nphzz
			if (ntime==nphzz*it) then
				fvzz = 0.
				call phase_space_vxyz_vs_z(part,fvzz,3,nz,npp,dt,ntime,nblok,idproc)
			endif
		endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    wke = 0.

    ! push electrons
    if (relativity == 1) then
        call rpush(part, fxyze, npp, noff, qbme, dt, ci, wke, tpush, nx, ny, nz, ip&
        &bc, inorder, popt)
        !        call rpush(part,fxyze,bxyze,npp,noff,qbme,dt,dt,ci,wke,tpush,nx&
        !    &,ny,nz,ipbc,inorder,popt)
    else
        call push(part, fxyze, npp, noff, qbme, dt, wke, tpush, nx, ny, nz, ipbc, i&
        &norder, popt)
        !        call push(part,fxyze,bxyze,npp,noff,qbme,dt,dt,wke,tpush,nx,ny,&
        !    &nz,ipbc,inorder,popt)
    endif

    ! move electrons into appropriate spatial regions
!    call pmove(part, edges, npp, tmove, ny, nz, kstrt, nvpy, nvpz, nbmax, idds, i&
!    &blok, vect, ierr)
     call pmoves(part,edges,npp,tmove,ny,nz,kstrt,nvpy,nvpz,nbmax,idds,&
     &iblok,vect,ierr)

    ! push Ian's electron test charges
    keteste = 0.
    if (fixedvel==1) then
       parteTest(1:3, 1:nppeTest(1), 1) = parteTest(1:3, 1:nppeTest(1), 1) + &
       &dt * parteTest(4:6, 1:nppeTest(1), 1)
    else
        if (relativity == 1) then
           call rpush(parteTest, fxyze, nppeTest, noff, qbmetest, dt, ci, keteste, tpush, nx, ny, nz, ip&
           &bc, inorder, popt)
        else
           call push(parteTest, fxyze, nppeTest, noff, qbmetest, dt, keteste, tpush, nx, ny, nz, ipbc, i&
           &norder, popt)
        endif
        keteste = keteste * rmteste
    endif

    ! move Ian's test electrons into appropriate spatial regions
    if (track_teste .ne. 0) then
       call pmoves(parteTest, edges, nppeTest, tmove, ny, nz, kstrt, nvpy, nvpz, nteste, idd&
       &s, iblok, vect, ierr, tracks)
    else
       call pmoves(parteTest, edges, nppeTest, tmove, ny, nz, kstrt, nvpy, nvpz, nteste, idd&
       &s, iblok, vect, ierr)
    endif
    
    if (rsqn > 0.0 .and. qmsqn > 0 .and. rmsqn > 0) then
        ! push the Strange Quark Nugget
        call ipushSQN(sqn, partiSQN, fxyze, nppiSQN, numSQNparti, noff, qbmsqn, dt, keSQN, tpushi, nx, ny, nz, ipbc)
        keSQN = keSQN*rmsqn

        ! move the Strange Quark Nugget particles to the appropriate spatial regions
        call pmoves(partiSQN, edges, nppiSQN, tmovi, ny, nz, kstrt, nvpy, nvpz, nbmax, i&
        &dds, iblok, vect, ierr)

        if(idproc==0) then
            !write position, velocity, and energy to the file
            write(iusqn,'(f11.6, 7e23.15)') &
            &ntime*dt, sqn(1,1), sqn(2,1), sqn(3,1), sqn(4,1), sqn(5,1), sqn(6,1), keSQN
            flush(iusqn)
        endif
    endif

    ! push ions
    if (movion == 1) then
        wki = 0.
        if (relativity == 1) then
            call rpush(parti, fxyze, nppi, noff, qbmi, dt, ci, wki, tpushi, nx, ny&
            &, nz, ipbc, inorder, popt)
            !           call rpush(parti,fxyze,bxyze,nppi,noff,qbmi,dt,dt,ci,wki,tpu&
            !    &shi,nx,ny,nz,ipbc,inorder,popt)
        else
            call push(parti, fxyze, nppi, noff, qbmi, dt, wki, tpushi, nx, ny, nz, &
            &ipbc, inorder, popt)
            !           call push(parti,fxyze,bxyze,nppi,noff,qbmi,dt,dt,wki,tpushi,&
            !    &nx,ny,nz,ipbc,inorder,popt)
        endif
        wki = wki * rmass

        ! move ions into appropriate spatial regions
!        call pmove(parti, edges, nppi, tmovi, ny, nz, kstrt, nvpy, nvpz, nbmax, i&
!        &dds, iblok, vect, ierr)
        call pmoves(parti,edges,nppi,tmovi,ny,nz,kstrt,nvpy,nvpz,nbmax,&
        &idds,iblok,vect,ierr)


        if (fixedvel==1) then
           ! push Ian's ion test charges
           partiTest(1:3, 1:nppiTest(1), 1) = partiTest(1:3, 1:nppiTest(1), 1) + &
           &dt * partiTest(4:6, 1:nppiTest(1), 1)
        else
           if (relativity == 1) then
               call rpush(partiTest, fxyze, nppiTest, noff, qbmitest, dt, ci, ketesti, tpushi, nx, ny&
               &, nz, ipbc, inorder, popt)
           else
               call push(partiTest, fxyze, nppiTest, noff, qbmitest, dt, ketesti, tpushi, nx, ny, nz, &
               &ipbc, inorder, popt)
           endif
           ketesti = ketesti * rmtesti

        endif

        ! move Ian's test ions into appropriate spatial regions
        if (track_testi .ne. 0) then
           call pmoves(partiTest, edges, nppiTest, tmovi, ny, nz, kstrt, nvpy, nvpz, ntesti, idd&
           &s, iblok, vect, ierr, tracks)
        else
           call pmoves(partiTest, edges, nppiTest, tmovi, ny, nz, kstrt, nvpy, nvpz, ntesti, idd&
           &s, iblok, vect, ierr)
        endif
    endif




    ! momentum diagnostic
    if (ntm > 0) then
        it = ntime/ntm
        it = ntime - ntm * it + 1
        if (it > 1) it = it - ntm
        ! calculate electron momentum
        if (it >= 0) then
            call premoment3(part, ntime, npp, id0, ium, pxe, pye, pze, sx, sy, sz, &
            &wx, wy, wz, nprint = it)
            ! calculate ion momentum
            if (movion == 0) then
                if (it == 1) then
                    call imoment(qi, fxyze, nyzp, id0, ium, pxi, pyi, pzi, dt, wx, w&
                    &y, wz, nx, inorder)
                endif
            else if (movion == 1) then
                call primoment3(parti, nppi, id0, ium, rmass, pxi, pyi, pzi, wx, w&
                &y, wz, nprint = it)
            endif
            ! print total momentum
            if (it == 1) then
                if (id0 == 0) write (ium, 996) wx, wy, wz
            endif
        endif
    endif
    ! sort electrons
    if (sortime > 0) then
        if (mod(ntime, sortime) == 0) then
            call sortp(part, pt, ip, npp, noff, nyzp, npic, tsort, inorder)
            !            call sortp(part,part2,npp,noff,nyzp,npic,tsort,inorder)
        endif
    endif
    ! sort ions
    if ((movion == 1) .and. (sortimi > 0)) then
        if (mod(ntime, sortimi) == 0) then
            call sortp(parti, pt, ip, nppi, noff, nyzp, npic, tsorti, inorder)
            !           call sortp(parti,parti2,nppi,noff,nyzp,npic,tsorti,inorder)
        endif
    endif
    ! energy diagnostic
    call esenergy(wt, wtot, msg, we, wke, wki, ntw, ndw, id0, itw, iuot, ntime)
    itime = itime + 1
    ntime = itime + itime0
    ! restart file
    if (ntr > 0) then
        it = ntime/ntr
        if (ntime == ntr * it) then
            it = mod(it - 1, 2) + 1
            ! write file

            if (addtag .ne. 0) then

               call restart_write(itime, itime0, npp, part, movion, nppi, parti, sqn, nppiSQN, partiSQN,&
                     & nppeTest, parteTest, nppiTest, partiTest, itw, wt, tracks)

            else

               call restart_write(itime, itime0, npp, part, movion, nppi, parti, sqn, nppiSQN, partiSQN,&
                     & nppeTest, parteTest, nppiTest, partiTest, itw, wt)

            endif

            ! if (nppeTest(1)==1) print*,parteTest(1:3,1,1)

        endif
    endif
    call wtimer(tloop, ldtime)
    ltime = ltime + tloop

    ! check to see if we can run for another period
    
    if (ntr>0 .and. check_period > 0 .and. maxtime > 0.0) then
        it = ntime/check_period
        !print*, 'maxtime=',maxtime,'unix_time=',unix_time()
        if (ntime == it * check_period) then

            ! calculate average time for a period so far
            average_period_time = check_period * ltime/real(itime)

            ! if (idproc==0) print*, 'systime =', time_(), 'maxtime =', maxtime, 'diff =', maxtime - time_()
            
            ! write out restart if we don't have time left for another period
            ! note the -100 is an attempt to give us enough time to write a restart
            if ( maxtime - 100 - unix_time() < average_period_time ) then
                if ( ntr > 0) then
                    it = ntime/ntr
                    if (ntime == ntr * it) goto 2000
                endif

                ! write restart file
                if (addtag .ne. 0) then
                   call restart_write(itime, itime0, npp, part, movion, nppi, parti, sqn, nppiSQN, partiSQN, &
                        & nppeTest, parteTest, nppiTest, partiTest, itw, wt, tracks)
                else
                   call restart_write(itime, itime0, npp, part, movion, nppi, parti, sqn, nppiSQN, partiSQN, &
                        & nppeTest, parteTest, nppiTest, partiTest, itw, wt)
                endif
                     
                goto 2000
                
            endif
        endif
    endif

    go to 500
2000 continue
    !
    ! * * * end main iteration loop * * *
    !
    ! send QUIT message to diagnostic nodes

    !stop HDF5 and destroy the I/O comm
    !call stop_hdf5(ierr)
    call destroy_io_comm()

    msg = -1.
    call HARTBEAT(msg, 1)
    call pwtimer(time, dtime)
    ! send main CPU Time to diagnostic nodes
    msg(1:2) = time; msg(3) = tpush; msg(4) = tdpost; msg(5) = tsort
    msg(6:7) = tmove; msg(8:9) = tfft; msg(10) = ltime
    call HARTBEAT(msg, 10)
    if (id0 == 0) then
        write (iuot, *) 'processor partition used: nvpy, nvpz = ', nvpy, &
        &nvpz
        write (iuot, *) ncpus, ' processors found, ', ntasks + 1, ' used'
        write (iuot, *) 'main max/min real time=', time(1), time(2), 'sec'
        totpush = tpush + tdpost
        write (iuot, *) 'electron push time=', tpush, 'sec'
        write (iuot, *) 'electron charge deposit time=', tdpost, 'sec'
        write (iuot, *) 'total electron push time=', totpush, 'sec'
        write (iuot, *) 'electron sort time=', tsort, 'sec'
        write (iuot, *) 'electron move time=', tmove, 'sec'
        totpush = totpush + tsort + tmove(1)
        write (iuot, *) 'total electron time=', totpush, 'sec'
        write (iuot,*) 'total add_tracks_time=',add_tracks_time
        write (iuot,*) 'total write_tracks_time=',write_tracks_time
    endif
    if (movion == 1) then
        msg(1) = tpushi; msg(2) = tdposti; msg(3) = tsorti
        msg(4:5) = tmovi
        call HARTBEAT(msg, 5)
        if (id0 == 0) then
            totpushi = tpushi + tdposti
            write (iuot, *) 'ion push time=', tpushi, 'sec'
            write (iuot, *) 'ion charge deposit time=', tdposti, 'sec'
            write (iuot, *) 'total ion push time=', totpushi, 'sec'
            write (iuot, *) 'ion sort time=', tsorti
            write (iuot, *) 'ion move time=', tmovi
            totpushi = totpushi + tsorti + tmovi(1)
            write (iuot, *) 'total ion time=', totpushi, 'sec'
        endif
    endif
    if (id0 == 0) then
        write (iuot, *) 'total fft time=', tfft, 'sec'
        time(1) = time(1) - (totpush + totpushi + tfft(1))
        write (iuot, *) 'other time=', time(1), ltime, 'sec'
        ! write final diagnostic metafile
        fname = 'pdiag32.' // cdrun
        open(unit = iudm, file = trim(fname), form = 'formatted', status = 'replac&
        &e')
        ! ion density diagnostics
        if (ntd > 0) then
            if (id0 == 0) then
                ndrec = ndrec - 1
                ceng = zero
                write (iudm, pden32d, iostat = irc)
                if (irc /= 0) then
                    write (iuer, *) 'pden32d namelist not written'
                endif
            endif
        endif
        ! potential diagnostics
        if (ntp > 0) then
            if (id0 == 0) then
                nprec = nprec - 1
                ceng = affp
                write (iudm, ppot32d, iostat = irc)
                if (irc /= 0) then
                    write (iuer, *) 'ppot32d namelist not written'
                endif
            endif
        endif
        ! write out input file
        write (iudm, pinput3, iostat = irc)
        if (irc /= 0) write (iuer, *) 'pinput3 namelist not written'
        ! done
        write (iuot, *) '* * * q.e.d. * * *'
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
991     format (' T = ', i7)
992     format (' field, kinetic, total energies = ', 3e14.7)
        ! allocate data for restart and/or phase space diagnostic
        if ((nts > 0) .or. (nustrt /= 1) .or. (ntr > 0)) then
            allocate(part(idimp, max(npmax, npimax), inblok))
            allocate(npp(inblok))
            parti => part
            nppi => npp
        endif
        if (movion == 0) then
            allocate(qi(nxe, nypmx, nzpmx * kbmin, kblok * lblok))
        else
            allocate(qi(0, 0, 0, 0))
        endif

        ! restart
        if (nustrt == 0) then
            ! read restart and extend run
            call restart_read(idrun, npp, part, movion, nppi, parti, sqn, nppiSQN, partiSQN, &
                     & nppeTest, parteTest, nppiTest, partiTest, smoothIon, itime, itime0, itw, wt)
            itime0 = itime + itime0
            t0 = dt*real(itime0)
            itime = 0
            ntime = itime + itime0
        elseif (nustrt == 2) then
            ! read restart but start everything except particle positions anew
            call restart_read(idrun, npp, part, movion, nppi, parti, sqn, nppiSQN, partiSQN, &
                   & nppeTest, parteTest, nppiTest, partiTest, smoothIon)
            itime0 = 0
            itime = itime0
            ntime = itime + itime0
        endif

        ! initialize diagnostics
40  if (id0 == 0) then
        iudm = get_funit(iudm)
        fname = 'pdiag32.init.' // cdrun
        open(unit = iudm, file = trim(fname), form = 'formatted', status = 'rep&
        &lace')
    endif
    ! ion density diagnostics
    call initmodediag(dent, ntd, id0, nxh, nyh, nzh, kxyp, kyzp, modesxd, &
    &modesyd, modeszd, jblok, mblok, iud, ndrec, fdname)
    if (ntd > 0) then
        if (id0 == 0) then
            ceng = zero
            write (iudm, pden32d, iostat = irc)
        endif
    endif
    ! velocity diagnostics
    fname = 'fv3.' // cdrun
    call initveldiag(fv, fvm, vtx, vty, vtz, ntv, ndv, id0, nmv, nblok, iuv, f&
    &name)
    if (movion == 1) then
        fname = 'fvi3.' // cdrun
        call initveldiag(fvi, fvmi, vtxi, vtyi, vtzi, ntv, ndv, id0, nmv, nbl&
        &ok, iuvi, fname)
    endif
    ! potential diagnostics
    call initmodediag(pott, ntp, id0, nxh, nyh, nzh, kxyp, kyzp, modesxp, &
    &modesyp, modeszp, jblok, mblok, iup, nprec, fpname)
    if (ntp > 0) then
        if (id0 == 0) then
            ceng = affp
            write (iudm, ppot32d, iostat = irc)
        endif
    endif
    ! write out and close input file
    if (id0 == 0) then
        write (iudm, pinput3, iostat = irc)
        close(unit = iudm)
    endif
    ! get initial CPU Time
    call HARTBEAT(msg, 2)
    time(1) = msg(1); time(2) = msg(2)
    if (id0 == 0) then
        write (iuot, *) 'init max/min real time=', time(1), time(2), 's&
        &ec'
    endif
    ! get time step
    10 call HARTBEAT(msg, 1)
    it = msg(1)
    if (it < 0) then
        ! get main CPU Time
        call HARTBEAT(msg, 10)
        time = msg(1:2); tpush = msg(3); tdpost = msg(4)
        tsort = msg(5); tmove = msg(6:7); tfft = msg(8:9)
        ltime = msg(10)
        if (id0 == 0) then
            write (iuot, *) 'processor partition used: nvpy, nvpz = ', &
            & nvpy, nvpz
            write (iuot, *) ncpus, ' processors found, ', ntasks + 1, ' &
            &used'
            write (iuot, *) 'main max/min real time=', time(1), time(2), &
            &'sec'
            totpush = tpush + tdpost
            write (iuot, *) 'electron push time=', tpush, 'sec'
            write (iuot, *) 'electron charge deposit time=', tdpost, '&
            &sec'
            write (iuot, *) 'total electron push time=', totpush, 'sec'
            write (iuot, *) 'electron sort time=', tsort, 'sec'
            write (iuot, *) 'electron move time=', tmove, 'sec'
            totpush = totpush + tsort + tmove(1)
            write (iuot, *) 'total electron time=', totpush, 'sec'
        endif
        if (movion == 1) then
            call HARTBEAT(msg, 5)
            tpushi = msg(1); tdposti = msg(2); tsorti = msg(3)
            tmovi = msg(4:5)
            if (id0 == 0) then
                totpushi = tpushi + tdposti
                write (iuot, *) 'ion push time=', tpushi, 'sec'
                write (iuot, *) 'ion charge deposit time=', tdposti, 's&
                &ec'
                write (iuot, *) 'total ion push time=', totpushi, 'sec'
                write (iuot, *) 'ion sort time=', tsorti
                write (iuot, *) 'ion move time=', tmovi
                totpushi = totpushi + tsorti + tmovi(1)
                write (iuot, *) 'total ion time=', totpushi, 'sec'
            endif
        endif
        if (id0 == 0) then
            write (iuot, *) 'total fft time=', tfft, 'sec'
            time(1) = time(1) - (totpush + totpushi + tfft(1))
            write (iuot, *) 'other time=', time(1), ltime, 'sec'
            ! write final diagnostic metafile
            fname = 'pdiag32.' // cdrun
            open(unit = iudm, file = trim(fname), form = 'formatted', status = '&
            &replace')
            ! ion density diagnostics
            if (ntp > 0) then
                ndrec = ndrec - 1
                ceng = zero
                write (iudm, pden32d, iostat = irc)
                if (irc /= 0) then
                    write (iuer, *) 'pden32d namelist not written'
                endif
            endif
            ! potential diagnostics
            if (ntp > 0) then
                nprec = nprec - 1
                ceng = affp
                write (iudm, ppot32d, iostat = irc)
                if (irc /= 0) then
                    write (iuer, *) 'ppot32d namelist not written'
                endif
            endif
            ! write out input file
            write (iudm, pinput3, iostat = irc)
            if (irc /= 0) then
                write (iuer, *) 'pinput3 namelist not written'
            endif
            ! done
            write (iuot, *) '* * * q.e.d. * * *'
        endif
        call MP_END
        call PPEXIT
        stop
    else
        ntime = it
    endif
    if (id0 == 0) write (iuot, 991) ntime
    write (label, 991) ntime
    call LOGNAME(label)
    ! ion density diagnostic
    if (ntd > 0) then
        it = ntime/ntd
        if (ntime == ntd * it) then
            ! write diagnostic output
            modesz2 = 2 * modeszd - 1
            call writebf(dent, modesxd, modesyd, modesz2, kxyp, kyzp, jblok&
            &, iud, ndrec)
        endif
    endif
    ! velocity diagnostic
    if (ntv > 0) then
        it = ntime/ntv
        if (ntime == ntv * it) then
            ! print out velocity moments
            if (id0 == 0) write (iuv, *) it, fvm(1,:,:), fvm(2,:,:)
            if (movion == 1) then
                ! print out velocity moments
                if (id0 == 0) write (iuvi, *) it, fvmi(1,:,:), fvmi(2,:,:)
            endif
        endif
    endif
    ! potential diagnostic
    if (ntp > 0) then
        it = ntime/ntp
        if (ntime == ntp * it) then
            ! write diagnostic output
            modesz2 = 2 * modeszp - 1
            call writebf(pott, modesxp, modesyp, modesz2, kxyp, kyzp, jblok&
            &, iup, nprec)
        endif
    endif
    ! energy diagnostic
    if (ntw > 0) then
        it = ntime/ntw
        if (ntime == ntw * it) then
            ! get energy values
            call HARTBEAT(msg, 4)
            wtot = msg(1:4)
            if (id0 == 0) write (iuot, 992) wtot(1), wtot(2), wtot(4)
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
        if (ntime == ntr * it) then
            it = mod(it - 1, 2) + 1
            ! write file
            call restart_write(itime, itime0, npp, part, movion, nppi, parti, sqn, nppiSQN, partiSQN, &
                  & nppeTest, parteTest, nppiTest, partiTest, itw, wt)
        endif
    endif
    go to 10
end subroutine

end program pbeps32

