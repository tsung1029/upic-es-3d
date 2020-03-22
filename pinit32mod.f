!-----------------------------------------------------------------------
!
      module pinit32d
!
! Fortran90 interface to 3d parallel PIC Fortran77 library pinit32lib.f
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: june 11, 2008
!
      use globals, only: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR,&
     & PERIODIC_3D, DIRICHLET_3D, MIXED_3D, VACUUM_3D
      implicit none
      private
      public :: LINEAR, QUADRATIC, STANDARD, LOOKAHEAD, VECTOR
      public :: PERIODIC_3D, DIRICHLET_3D, MIXED_3D, VACUUM_3D
      public :: idrun, idrun0, indx, indy, indz, npx, npy, npz
      public :: npxb, npyb, npzb, inorder
      public :: popt, dopt, djopt, nustrt, ntr
      public :: ntw, ntp, ntd, nta, ntv, nts, ntm, nte
      public :: ndw, ndp, ndd, nda, ndv, nds, ndm, nde
      public :: tend, dt, qme, vtx, vty, vtz, vx0, vy0, vz0
      public :: vdx, vdy, vdz, vtdx, vtdy, vtdz
      public :: psolve, relativity, omx, omy, omz, ci, ax, ay, az
      public :: ndc, movion, npxi, npyi, npzi, npxbi, npybi, npzbi
      public :: qmi, rmass, rtempxi, rtempyi, rtempzi, vxi0, vyi0, vzi0
      public :: vdxi, vdyi, vdzi, rtempdxi, rtempdyi, rtempdzi, v0, w0
      public :: sortime, sortimi, nplot, idpal, ndstyle, sntasks
      public :: itpon, ionoff, nsrand, ndprof, nsrandi, ndprofi
      public :: ampdx, scaledx, shiftdx, ampdy, scaledy, shiftdy
      public :: ampdz, scaledz, shiftdz
      public :: ampdxi, scaledxi, shiftdxi, ampdyi, scaledyi, shiftdyi
      public :: ampdzi, scaledzi, shiftdzi
      public :: modesxd, modesyd, modeszd, modesxp, modesyp, modeszp
      public :: modesxa, modesya, modesza, modesxe, modesye, modesze
      public :: imbalance, mpimon
      public :: pinput3, sendnml, distr, pvdistr, ldistr, fdistr, vdistr
      public :: vfdistr, vvdistr, fedges
      public :: t0, ceng, indian, rlprec, inprec, pden32d, ndrec, fdname
      public :: ppot32d, nprec, fpname, pvpot32d, narec, faname, pem32d
      public :: nerec, fename
!
! Namelist Input
      save
! idrun/idrun0 = run identifier for current/old run
      integer :: idrun = 0, idrun0 = 0
! indx/indy/indz = exponent which determines length in x/y/z direction,
! where nx=2**indx, ny=2**indy, nz=2**indz
!     integer :: indx =   4, indy =   3, indz =   5
      integer :: indx =   5, indy =   4, indz =   6
!     integer :: indx =   6, indy =   5, indz =   7
!     integer :: indx =   7, indy =   6, indz =   8
! npx/npy/npz = initial number of electrons distributed in x/y/z
! direction
!     integer :: npx =      24, npy =      12, npz =     48
      integer :: npx =      64, npy =      32, npz =    128
!     integer :: npx =     192, npy =      96, npz =    384
!     integer :: npx =     512, npy =     256, npz =   1024
! npxb/npyb/npzb = initial number of electrons in beam in x/y/z
! direction
!     integer :: npxb =   0, npyb =   0, npzb =  0
!     integer :: npxb =  12, npyb =   6, npzb =  24
      integer :: npxb =  32, npyb =  16, npzb =  64
!     integer :: npxb =  96, npyb =  48, npzb = 192
!     integer :: npxb = 256, npyb = 128, npzb = 512
! inorder = interpolation order
! popt = particle optimization scheme
! dopt = charge deposit optimization scheme
! djopt = current deposit optimization scheme
      integer :: inorder = LINEAR, popt = STANDARD, dopt = LOOKAHEAD
      integer :: djopt = STANDARD
! nustrt = (0,1,2) = this is an (old start,new start,restart)
! ntr = number of time steps between restart routine
      integer :: nustrt = 1, ntr = 0
! ntw, ndw = number of time steps between energy diagnostic
! ntp, ndp = number of time steps between potential diagnostic
! ntd, ndd = number of time steps between ion density diagnostic
! nta, nda = number of time steps between vector potential diagnostic
! ntv, ndv = number of time steps between velocity-space diagnostic
! nts, nds = number of time steps between phase space diagnostic
      integer :: ntw = 1, ntp = 0, ntd = 0, nta = 0, ntv = 0, nts = 0
      integer :: ndw = 1, ndp = 0, ndd = 0, nda = 0, ndv = 0, nds = 0
! ntm, ndm = number of time steps between momentum diagnostic
! nte, nde = number of time steps between electromagnetic diagnostic
      integer :: ntm = 0, nte = 0
      integer :: ndm = 0, nde = 0
! tend = time at end of simulation, in units of plasma frequency
! dt = time interval between successive calculations
      real :: tend =  85.000, dt = 0.2000000e+00
! qme = charge on electron, in units of e
! vtx/vty/vtz = thermal velocity of electrons in x/y/z direction
      real :: qme = -1.0, vtx = 1.0, vty = 1.0, vtz = 1.0
! vx0/vy0/vz0 = drift velocity of electrons in x/y/z direction
      real :: vx0 = 0.0, vy0 = 0.0, vz0 = 0.0
! vdx/vdy/vdz = drift velocity of beam electrons in x/y/z direction
      real :: vdx = 0.0, vdy = 0.0, vdz = 0.0
! vtdx/vtdy/vtdz = thermal velocity of beam electrons in x/y/z direction
      real :: vtdx = 1.0, vtdy = 1.0, vtdz = 1.0
! psolve = type of poisson solver = (1,2,3)
      integer :: psolve = PERIODIC_3D
! relativity = (no,yes) = (0,1) = relativity is used
      integer :: relativity = 0
! omx/omy/omz = magnetic field electron cyclotron frequency in x/y/z 
      real :: omx = 0.0, omy = 0.0, omz = 0.0
! ci = reciprical of velocity of light
      real :: ci = 0.1
! ax/ay/az = half-width of particle in x/y/z direction
!     real :: ax = .816497, ay = .816497, az = .816497
!     real :: ax = .866025, ay = .866025, az = .866025
      real :: ax = .912871, ay = .912871, az = .912871
! ndc = number of corrections in darwin iteration
      integer :: ndc = 2
! movion = (0,1) = (no,yes) move the ions
! npxi/npyi/npzi = initial number of ions distributed in x/y/z direction
      integer :: movion = 0, npxi =  64, npyi =  32, npzi = 128
! npxbi/npybi/npzbi = initial number of ions in beam in x/y/z direction
      integer :: npxbi = 32, npybi = 16, npzbi = 64
! qmi = charge on ion, in units of 3
! rmass = ion/electron mass ratio
      real :: qmi = 1.0, rmass = 16.0
! rtempxi/rtempyi/rtempzi = electron/ion temperature ratio of background
! ions in x/y/z direction
      real :: rtempxi = 1.0, rtempyi = 1.0, rtempzi = 1.0
! vxi0/vyi0/vzi0 = drift velocity of ions in x/y/z direction
      real :: vxi0 = 0.0, vyi0 = 0.0, vzi0 = 0.0
! vdxi/vdyi/vdzi = drift velocity of beam ions in x/y/z direction
      real :: vdxi = 0.0, vdyi = 0.0, vdzi = 0.0
! rtempdxi/rtempdyi/rtempdzi = electron/ion temperature ratio of beam
! ions in x/y/z direction
      real :: rtempdxi = 1.0, rtempdyi = 1.0, rtempdzi = 1.0     
! v0 = external pump strength, in units vos/vthermal
! w0 = external pump frequency, in units of wpe
      real :: v0 = 0.0, w0 = 0.0
! sortime = number of time steps between electron sorting
! sortimi = number of time steps between ion sorting
      integer :: sortime = 50, sortimi = 250
! nplot = maximum number of plots per page
! idpal = palette id number: 1 = cold/hot, 2 = color wheel, 3 = rainbow
! ndstyle = (1,2,3) = display (color map,contour plot,both)
      integer :: nplot = 4, idpal = 1, ndstyle = 1
! sntasks = (-1,n) = set maximum number of tasks (-1 = number of cpus-1)
      integer :: sntasks = -1
! itpon = time when external pump is turned on (-1=never)
! ionoff = time when ions are frozen and their charge saved (-1=never)
      integer :: itpon = -1, ionoff = -1
! nsrand = (0,1) = (no,yes) randomize spatially positions locally
! ndprof = profile type (uniform=0,linear=1,sinusoidal=2,gaussian=3,
!                        hyperbolic secant squared=4)
! nsrandi = (0,1) = (no,yes) randomize spatially ion positions locally
! ndprofi = ion profile (uniform=0,linear=1,sinusoidal=2,gaussian=3,
!                        hyperbolic secant squared=4)
      integer :: nsrand = 0, ndprof = 0, nsrandi = 0, ndprofi = 0
! ampdx/ampdx = amplitude of density compared to uniform in x/y
! scaledx/scaledx = scale length for spatial coordinate in x/y
! shiftdx/shiftdx = shift of spatial coordinate in x/y
      real :: ampdx = 0.0, scaledx = 0.0, shiftdx = 0.0
      real :: ampdy = 0.0, scaledy = 0.0, shiftdy = 0.0
      real :: ampdz = 0.0, scaledz = 0.0, shiftdz = 0.0
! ampdxi/ampdxi = amplitude of ion density compared to uniform in x/y
! scaledxi/scaledxi = scale length for spatial ion coordinate in x/y
! shiftdxi/shiftdxi = shift of spatial ion coordinate in x/y
      real :: ampdxi = 0.0, scaledxi = 0.0, shiftdxi = 0.0
      real :: ampdyi = 0.0, scaledyi = 0.0, shiftdyi = 0.0
      real :: ampdzi = 0.0, scaledzi = 0.0, shiftdzi = 0.0
! modesxd/modesyd/modeszd = number of modes in x/y/z to keep for ion
! density diagnostic
      integer :: modesxd = 11, modesyd = 11, modeszd = 11
! modesxp/modesyp/modeszp = number of modes in x/y/z to keep for 
! potential diagnostic
      integer :: modesxp = 11, modesyp = 11, modeszp = 11
! modesxa/modesya/modesza = number of modes in x/y/z to keep for vector
! potential diagnostic
      integer :: modesxa = 11, modesya = 11, modesza = 11
! modesxe/modesye/modesze = number of modes in x/y/z to keep for
! electromagnetic diagnostic
      integer :: modesxe = 11, modesye = 11, modesze = 11
! imbalance = load imbalance fraction repartition trigger
! (< 0.  to suppress repartion)
      real :: imbalance = .08
! mpimon = (0,1,2) = (suppress,display,display and log) mpi messages
      integer :: mpimon = 1
! define namelist
      namelist /pinput3/ idrun, idrun0, indx, indy, indz, npx, npy, npz,&
     &npxb, npyb, npzb, inorder, popt, dopt, djopt, nustrt, ntr, ntw,   &
     &ntp, ntd, nta, ntv, nts, ntm, nte, ndw, ndp, ndd, nda, ndv, nds,  &
     &ndm, nde, tend, dt, qme, vtx, vty, vtz, vx0, vy0, vz0, vdx, vdy,  &
     &vdz, vtdx, vtdy, vtdz, psolve, relativity, omx, omy, omz, ci, ax, &
     &ay, az, ndc, movion, npxi, npyi, npzi, npxbi, npybi, npzbi, qmi,  &
     &rmass, rtempxi, rtempyi, rtempzi, vxi0, vyi0, vzi0, vdxi, vdyi,   &
     &vdzi, rtempdxi, rtempdyi, rtempdzi, v0, w0, sortime, sortimi,     &
     &nplot, idpal, ndstyle, sntasks, itpon, ionoff, nsrand, ndprof,    &
     &nsrandi, ndprofi, ampdx, scaledx, shiftdx, ampdy, scaledy, shiftdy&
     &, ampdz, scaledz, shiftdz, ampdxi, scaledxi, shiftdxi, ampdyi,    &
     &scaledyi, shiftdyi, ampdzi, scaledzi, shiftdzi, modesxd, modesyd, &
     &modeszd, modesxp, modesyp, modeszp, modesxa, modesya, modesza,    &
     &modesxe, modesye, modesze, imbalance, mpimon
!
! t0 = initial time value
! ceng = energy normalization
      real :: t0 = 0.0, ceng = 1.0
! indian = (0,1) = architecture is (little-endian,big-endian)
! rlprec = (0,1) = default reals are (normal,double-precision)
! inprec = (0,1) = default integers are (normal,double-precision)
      integer :: indian = 1, rlprec = 1, inprec = 0
!
! Namelist output for ion density diagnostic
! ndrec = current record number for ion density writes
      integer :: ndrec = 0
! fdname = file name for potential diagnostic
      character(len=32) :: fdname
! define namelist
      namelist /pden32d/ idrun, indx, indy, ntd, modesxd, modesyd,      &
     &modeszd, psolve, ndrec, t0, tend, dt, ceng, indian, rlprec, inprec&
     &, fdname 
!
! Namelist output for potential diagnostic
! nprec = current record number for potential writes
      integer :: nprec = 0
! fpname = file name for potential diagnostic
      character(len=32) :: fpname
! define namelist
      namelist /ppot32d/ idrun, indx, indy, ntp, modesxp, modesyp,      &
     &modeszp, psolve, omx, omy, omz, nprec, t0, tend, dt, ceng, indian,&
     & rlprec, inprec, fpname
!
! Namelist output for vector potential diagnostic
! narec = current record number for vector potential writes
      integer :: narec = 0
! faname = file name for potential diagnostic
      character(len=32) :: faname
! define namelist
      namelist /pvpot32d/ idrun, indx, indy, nta, modesxa, modesya,     &
     &modesza, psolve, omx, omy, omz, ci, narec, t0, tend, dt, ceng,    &
     &indian, rlprec, inprec, faname
!
! Namelist output for electromagnetic diagnostic
! nerec = current record number for electromagnetic writes
      integer :: nerec = 0
! fename = file name for potential diagnostic
      character(len=32) :: fename
! define namelist
      namelist /pem32d/ idrun, indx, indy, nte, modesxe, modesye,       &
     &modesze, psolve, omx, omy, omz, ci, nerec, t0, tend, dt, ceng,    &
     &indian, rlprec, inprec, fename
!
! define interface to original Fortran77 procedures
      interface
         subroutine PISTR32(part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,n&
     &px,npy,npz,nx,ny,nz,idimp,npmax,mblok,nblok,idps,ipbc,ierr)
         implicit none
         integer :: npx, npy, npz, nx, ny, nz, idimp, npmax
         integer :: mblok, nblok, idps, ipbc, ierr
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,npmax,mblok*nblok) :: part
         real, dimension(idps,mblok*nblok) :: edges
         integer, dimension(mblok*nblok) :: npp, nps
         end subroutine
      end interface
      interface
         subroutine PVISTR32(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,np&
     &y,npz,nx,ny,nz,idimp,npmax,mblok,nblok,ipbc,vranx,vrany,vranz,kstr&
     &t,nvp,ndv,nvrp,ierr)
         implicit none
         integer :: npx, npy, npz, nx, ny, nz, idimp, npmax
         integer :: mblok, nblok, ipbc, kstrt, nvp, ndv, nvrp, ierr
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,npmax,mblok*nblok) :: part
         integer, dimension(mblok*nblok) :: npp, nps
         double precision, dimension(nvrp,mblok*nblok) :: vranx, vrany, &
     &vranz
         end subroutine
      end interface
      interface
         subroutine PLDISTR32(part,nps,anlx,anly,anlz,npx,npy,npz,nx,ny,&
     &nz,idimp,npmax,mblok,nblok,kstrt,nvp,ipbc,ierr)
         implicit none
         integer :: npx, npy, npz, nx, ny, nz, idimp, npmax
         integer :: mblok, nblok, kstrt, nvp, ipbc, ierr
         real :: anlx, anly, anlz
         real, dimension(idimp,npmax,mblok*nblok) :: part
         integer, dimension(mblok*nblok) :: nps
         end subroutine
      end interface
      interface
         subroutine PFDISTR32(part,nps,fnx,argx1,argx2,argx3,fny,argy1,a&
     &rgy2,argy3,fnz,argz1,argz2,argz3,npx,npy,npz,nx,ny,nz,idimp,npmax,&
     &mblok,nblok,kstrt,nvp,ipbc,ierr)
         implicit none
         integer :: npx, npy, npz, nx, ny, nz, idimp, npmax
         integer :: mblok, nblok, kstrt, nvp, ipbc, ierr
         real :: argx1, argx2, argx3, argy1, argy2, argy3
         real :: argz1, argz2, argz3
         real, dimension(idimp,npmax,mblok*nblok) :: part
         integer, dimension(mblok*nblok) :: nps
         real, external :: fnx, fny, fnz
         end subroutine
      end interface
      interface
         subroutine PRDISTR32(part,nps,fnx,argx1,argx2,argx3,fny,argy1,a&
     &rgy2,argy3,fnz,argz1,argz2,argz3,npx,npy,npz,nx,ny,nz,idimp,npmax,&
     &mnblok,kstrt,nvp,ipbc)
         implicit none
         integer :: npx, npy, npz, nx, ny, nz, idimp, npmax, mnblok
         integer :: kstrt, nvp, ipbc
         real :: argx1, argx2, argx3, argy1, argy2, argy3
         real :: argz1, argz2, argz3
         real, dimension(idimp,npmax,mnblok) :: part
         integer, dimension(mnblok) :: nps
         real, external :: fnx, fny, fnz
         end subroutine
      end interface
      interface
         subroutine PVRDISTR32(part,nps,fnx,argx1,argx2,argx3,fny,argy1,&
     &argy2,argy3,fnz,argz1,argz2,argz3,npx,npy,npz,nx,ny,nz,idimp,npmax&
     &,mnblok,vranx,vrany,vranz,kstrt,nvp,ndv,nvrp,ipbc,ierr)
         implicit none
         integer :: npx, npy, npz, nx, ny, nz, idimp, npmax, mnblok
         integer :: ipbc, kstrt, nvp, ndv, nvrp, ierr
         real :: argx1, argx2, argx3, argy1, argy2, argy3
         real :: argz1, argz2, argz3
         real, dimension(idimp,npmax,mnblok) :: part
         integer, dimension(mnblok) :: nps
         double precision, dimension(nvrp,mnblok) :: vranx, vrany, vranz
         real, external :: fnx, fny, fnz
         end subroutine
      end interface
      interface
         subroutine PVDISTR32(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,n&
     &py,npz,idimp,npmax,mblok,nblok,kstrt,nvp,ierr)
         implicit none
         integer :: npx, npy, npz, idimp, npmax
         integer :: mblok, nblok, kstrt, nvp, ierr
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,npmax,mblok*nblok) :: part
         integer, dimension(mblok*nblok) :: npp, nps
         end subroutine
      end interface
      interface
         subroutine PVVISTR32(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,n&
     &py,npz,idimp,npmax,mblok,nblok,vranx,vrany,vranz,kstrt,nvp,ndv,nvr&
     &p,ierr)
         implicit none
         integer :: npx, npy, npz, idimp, npmax
         integer :: mblok, nblok, kstrt, nvp, ndv, nvrp, ierr
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(idimp,npmax,mblok*nblok) :: part
         integer, dimension(mblok*nblok) :: npp, nps
         double precision, dimension(nvrp,mblok*nblok) :: vranx, vrany, &
     &vranz
         end subroutine
      end interface
      interface
         subroutine PBDISTR32L(part,bx,by,bz,npp,noff,qbm,nx,ny,nz,idimp&
     &,npmax,mblok,nblok,nxv,nypmx,nzpmx,idds)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mblok, nblok
         integer :: nxv, nypmx, nzpmx, idds
         real :: qbm
         real, dimension(idimp,npmax,mblok*nblok) :: part
         real, dimension(nxv,nypmx,nzpmx,mblok*nblok) :: bx, by, bz
         integer, dimension(mblok*nblok) :: npp
         integer, dimension(idds,mblok*nblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGBDISTR32L(part,bxyz,npp,noff,qbm,nx,ny,nz,idimp,np&
     &max,mblok,nblok,nxv,nypmx,nzpmx,idds,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mblok, nblok
         integer :: nxv, nypmx, nzpmx, idds, ipbc
         real :: qbm
         real, dimension(idimp,npmax,mblok*nblok) :: part
!        real, dimension(*) :: bxyz
         real :: bxyz
         integer, dimension(mblok*nblok) :: npp
         integer, dimension(idds,mblok*nblok) :: noff
         end subroutine
      end interface
      interface
         subroutine PGRBDISTR32L(part,bxyz,npp,noff,qbm,ci,nx,ny,nz,idim&
     &p,npmax,mblok,nblok,nxv,nypmx,nzpmx,idds,ipbc)
         implicit none
         integer :: nx, ny, nz, idimp, npmax, mblok, nblok
         integer :: nxv, nypmx, nzpmx, idds, ipbc
         real :: qbm, ci
         real, dimension(idimp,npmax,mblok*nblok) :: part
!        real, dimension(*) :: bxyz
         real :: bxyz
         integer, dimension(mblok*nblok) :: npp
         integer, dimension(idds,mblok*nblok) :: noff
         end subroutine
      end interface
      interface
         subroutine FEDGES32(edges,noff,nyzp,fny,argy1,argy2,argy3,fnz,a&
     &rgz1,argz2,argz3,ny,nz,nzpmin,nypmax,nypmin,nzpmax,kstrt,nvpy,nvpz&
     &,mblok,nblok,idps,idds,ipbc)
         implicit none
         integer :: ny, nz, nypmin, nypmax, nzpmin, nzpmax, kstrt
         integer :: nvpy, nvpz, mblok, nblok, idps, idds, ipbc
         real :: argy1, argy2, argy3, argz1, argz2, argz3
         real, dimension(idps,mblok*nblok) :: edges
         integer, dimension(idds,mblok*nblok) :: nyzp
         integer, dimension(idds,mblok*nblok) :: noff
         real, external :: fny, fnz
         end subroutine
      end interface
!
! define generic interfaces to Fortran90 library
!
      interface distr
         module procedure ipistr32
         module procedure ipbdistr32
         module procedure iprbdistr32
      end interface
!
      interface pvdistr
         module procedure ipvistr32
      end interface
!
      interface ldistr
         module procedure ipldistr32
      end interface
!
      interface fdistr
         module procedure ipfdistr32
      end interface
!
      interface vdistr
         module procedure ipvdistr32
      end interface
!
      interface vfdistr
         module procedure ipvfdistr32
      end interface
!
      interface vvdistr
         module procedure ipvvdistr32
      end interface
!
      interface fedges
         module procedure ifedges32
      end interface
!
      interface sendnml
         module procedure sendnml32
      end interface
!
! define Fortran90 interface functions to Fortran77 library
!
      contains
!
         subroutine sendnml32()
! this subroutine packs 3d namelist variables into a double precision
! buffer and broadcasts them to other nodes
         integer, parameter :: lenml = 125
         double precision, dimension(lenml) :: ddata
! pack data
         ddata(1) = idrun; ddata(2) = idrun0
         ddata(3) = indx; ddata(4) = indy; ddata(5) = indz
         ddata(6) = npx; ddata(7) = npy; ddata(8) = npz
         ddata(9) = npxb; ddata(10) = npyb; ddata(11) = npzb
         ddata(12) = inorder; ddata(13) = popt; ddata(14) = dopt
         ddata(15) = djopt; ddata(16) = nustrt; ddata(17) = ntr
         ddata(18) = ntw; ddata(19) = ntp; ddata(20) = ntd
         ddata(21) = nta; ddata(22) = ntv; ddata(23) = nts
         ddata(24) = ntm; ddata(25) = nte
         ddata(26) = ndw; ddata(27) = ndp; ddata(28) = ndd
         ddata(29) = nda; ddata(30) = ndv; ddata(31) = nds
         ddata(32) = ndm; ddata(33) = nde
         ddata(34) = tend; ddata(35) = dt; ddata(36) = qme
         ddata(37) = vtx; ddata(38) = vty; ddata(39) = vtz
         ddata(40) = vx0; ddata(41) = vy0; ddata(42) = vz0
         ddata(43) = vdx; ddata(44) = vdy; ddata(45) = vdz
         ddata(46) = vtdx; ddata(47) = vtdy; ddata(48) = vtdz
         ddata(49) = psolve; ddata(50) = relativity
         ddata(51) = omx; ddata(52) = omy; ddata(53) = omz
         ddata(54) = ci; ddata(55) = ax; ddata(56) = ay; ddata(57) = az
         ddata(58) = ndc; ddata(59) = movion
         ddata(60) = npxi; ddata(61) = npyi; ddata(62) = npzi
         ddata(63) = npxbi; ddata(64) = npybi; ddata(65) = npzbi
         ddata(66) = qmi; ddata(67) = rmass
         ddata(68) = rtempxi; ddata(69) = rtempyi; ddata(70) = rtempzi
         ddata(71) = vxi0; ddata(72) = vyi0; ddata(73) = vzi0
         ddata(74) = vdxi; ddata(75) = vdyi; ddata(76) = vdzi
         ddata(77) = rtempdxi; ddata(78) = rtempdyi
         ddata(79) = rtempdzi; ddata(80) = v0; ddata(81) = w0
         ddata(82) = sortime; ddata(83) = sortimi
         ddata(84) = nplot; ddata(85) = idpal; ddata(86) = ndstyle
         ddata(87) = sntasks; ddata(88) = itpon; ddata(89) = ionoff
         ddata(90) = nsrand; ddata(91) = ndprof
         ddata(92) = nsrandi; ddata(93) = ndprofi
         ddata(94) = ampdx; ddata(95) = scaledx; ddata(96) = shiftdx
         ddata(97) = ampdy; ddata(98) = scaledy; ddata(99) = shiftdy
         ddata(100) = ampdz; ddata(101) = scaledz; ddata(102) = shiftdz
         ddata(103) = ampdxi; ddata(104) = scaledxi
         ddata(105) = shiftdxi; ddata(106) = ampdyi
         ddata(107) = scaledyi; ddata(108) = shiftdyi
         ddata(109) = ampdzi; ddata(110) = scaledzi;
         ddata(111) = shiftdzi
         ddata(112) = modesxd; ddata(113) = modesyd
         ddata(114) = modeszd; ddata(115) = modesxp
         ddata(116) = modesyp; ddata(117) = modeszp
         ddata(118) = modesxa; ddata(119) = modesya
         ddata(120) = modesza; ddata(121) = modesxe
         ddata(122) = modesye; ddata(123) = modesze
         ddata(124) = imbalance; ddata(125) = mpimon
! broadcast data
         call PBDCAST(ddata,lenml)
! unpack data
         idrun = ddata(1); idrun0 = ddata(2)
         indx = ddata(3); indy = ddata(4); indz = ddata(5)
         npx = ddata(6); npy = ddata(7); npz = ddata(8)
         npxb = ddata(9); npyb = ddata(10); npzb = ddata(11)
         inorder = ddata(12); popt = ddata(13); dopt = ddata(14)
         djopt = ddata(15); nustrt = ddata(16); ntr = ddata(17)
         ntw = ddata(18); ntp = ddata(19); ntd = ddata(20)
         nta = ddata(21); ntv = ddata(22); nts = ddata(23)
         ntm = ddata(24); nte = ddata(25)
         ndw = ddata(26); ndp = ddata(27); ndd = ddata(28)
         nda = ddata(29); ndv = ddata(30); nds = ddata(31)
         ndm = ddata(32); nde = ddata(33)
         tend = ddata(34); dt = ddata(35); qme = ddata(36)
         vtx = ddata(37); vty = ddata(38); vtz = ddata(39)
         vx0 = ddata(40); vy0 = ddata(41); vz0 = ddata(42)
         vdx = ddata(43); vdy = ddata(44); vdz = ddata(45)
         vtdx = ddata(46); vtdy = ddata(47); vtdz = ddata(48)
         psolve = ddata(49); relativity = ddata(50)
         omx = ddata(51); omy = ddata(52); omy = ddata(53)
         ci = ddata(54); ax = ddata(55); ay = ddata(56); az = ddata(57)
         ndc = ddata(58); movion = ddata(59)
         npxi = ddata(60); npyi = ddata(61); npzi = ddata(62)
         npxbi = ddata(63); npybi = ddata(64); npzbi = ddata(65)
         qmi = ddata(66); rmass = ddata(67)
         rtempxi = ddata(68); rtempyi = ddata(69); rtempzi = ddata(70)
         vxi0 = ddata(71); vyi0 = ddata(72); vzi0 = ddata(73)
         vdxi = ddata(74); vdyi = ddata(75); vdzi = ddata(76)
         rtempdxi = ddata(77); rtempdyi = ddata(78)
         rtempdzi = ddata(79); v0 = ddata(80); w0 = ddata(81)
         sortime = ddata(82); sortimi = ddata(83)
         nplot = ddata(84); idpal = ddata(85); ndstyle = ddata(86)
         sntasks = ddata(87); itpon = ddata(88); ionoff = ddata(89)
         nsrand = ddata(90); ndprof = ddata(91)
         nsrandi = ddata(92); ndprofi = ddata(93)
         ampdx = ddata(94); scaledx = ddata(95); shiftdx = ddata(96)
         ampdy = ddata(97); scaledy = ddata(98); shiftdy = ddata(99)
         ampdz = ddata(100); scaledz = ddata(101); shiftdz = ddata(102)
         ampdxi = ddata(103); scaledxi = ddata(104)
         shiftdxi = ddata(105); ampdyi = ddata(106)
         scaledyi = ddata(107); shiftdyi = ddata(108)
         ampdzi = ddata(109); scaledzi = ddata(110)
         shiftdzi = ddata(111)
         modesxd = ddata(112); modesyd = ddata(113)
         modeszd = ddata(114); modesxp = ddata(115)
         modesyp = ddata(116); modeszp = ddata(117)
         modesxa = ddata(118); modesya = ddata(119)
         modesza = ddata(120); modesxe = ddata(121)
         modesye = ddata(122); modesze = ddata(123)
         imbalance = ddata(124); mpimon = ddata(125)
         end subroutine sendnml32
!
         subroutine ipistr32(part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,&
     &npx,npy,npz,nx,ny,nz,ipbc,mblok)
! calculates initial particle co-ordinates and velocities in 3d
! with uniform density and maxwellian velocity with drift
         implicit none
         integer :: npx, npy, npz, nx, ny, nz, ipbc, mblok
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:), pointer :: edges
         integer, dimension(:), pointer :: npp, nps
! local data
         integer :: idimp, npmax, nblok, idps, ierr
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)/mblok; idps = size(edges,1)
         call PISTR32(part,edges,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,npy&
     &,npz,nx,ny,nz,idimp,npmax,mblok,nblok,idps,ipbc,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         end subroutine ipistr32
!
         subroutine ipvistr32(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,n&
     &py,npz,nx,ny,nz,ipbc,kstrt,nvp,mblok)
! calculates initial particle co-ordinates and velocities in 3d
! with uniform density and maxwellian velocity with drift
! using parallel random number generator
         implicit none
         integer :: npx, npy, npz, nx, ny, nz, ipbc, kstrt, nvp, mblok
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp, nps
! local data
         integer, parameter :: ndv = 256
         integer :: idimp, npmax, nblok, nvrp, ierr
         double precision, dimension((ndv-1)/nvp+1,size(part,3)) :: vran&
     &x, vrany, vranz
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)/mblok
         nvrp = (ndv - 1)/nvp + 1
         call PVISTR32(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,npz,&
     &nx,ny,nz,idimp,npmax,mblok,nblok,ipbc,vranx,vrany,vranz,kstrt,nvp,&
     &ndv,nvrp,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         end subroutine ipvistr32
!
         subroutine ipldistr32(part,nps,anlx,anly,anlz,npx,npy,npz,nx,ny&
     &,nz,kstrt,nvp,ipbc,mblok)
! calculates initial particle co-ordinates in 3d
! with tri-linear density profile
         implicit none
         integer :: npx, npy, npz, nx, ny, nz, kstrt, nvp, ipbc, mblok
         real :: anlx, anly, anlz
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: nps
! local data
         integer :: idimp, npmax, nblok, ierr
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)/mblok
         call PLDISTR32(part,nps,anlx,anly,anlz,npx,npy,npz,nx,ny,nz,idi&
     &mp,npmax,mblok,nblok,kstrt,nvp,ipbc,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         end subroutine ipldistr32
!
         subroutine ipfdistr32(part,nps,ampx,scalex,shiftx,ampy,scaley,s&
     &hifty,ampz,scalez,shiftz,npx,npy,npz,nx,ny,nz,kstrt,nvp,ipbc,ndpro&
     &,nsran,mblok)
! calculates initial particle co-ordinates in 3d
! with various density profiles
         implicit none
         integer :: npx, npy, npz, nx, ny, nz, kstrt, nvp, ipbc, ndpro
         integer :: nsran, mblok
         real :: ampx, scalex, shiftx, ampy, scaley, shifty
         real :: ampz, scalez, shiftz
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: nps
! local data
         integer :: idimp, npmax, mnblok, nblok, ierr
         real :: sxi, syi, szi, zero
         real, external :: FLDISTR1, FSDISTR1, FGDISTR1, FHDISTR1
         idimp = size(part,1); npmax = size(part,2)
         mnblok = size(part,3); nblok = mnblok/mblok
         sxi = 0.
         if (scalex /= 0.) sxi = 1.0/scalex
         syi = 0.
         if (scaley /= 0.) syi = 1.0/scaley
         szi = 0.
         if (scalez /= 0.) szi = 1.0/scalez
         zero = 0.0
! uniform density
         if (ndpro==0) then
            call PFDISTR32(part,nps,FLDISTR1,zero,zero,zero,FLDISTR1,zer&
     &o,zero,zero,FLDISTR1,zero,zero,zero,npx,npy,npz,nx,ny,nz,idimp,npm&
     &ax,mblok,nblok,kstrt,nvp,ipbc,ierr)
            if ((ierr.eq.0).and.(nsran /= 0)) then
               call PRDISTR32(part,nps,FLDISTR1,zero,zero,zero,FLDISTR1,&
     &zero,zero,zero,FLDISTR1,zero,zero,zero,npx,npy,npz,nx,ny,nz,idimp,&
     &npmax,mnblok,kstrt,nvp,ipbc)
            endif
! linear density
         else if (ndpro==1) then
            call PFDISTR32(part,nps,FLDISTR1,ampx,sxi,shiftx,FLDISTR1,am&
     &py,syi,shifty,FLDISTR1,ampz,szi,shiftz,npx,npy,npz,nx,ny,nz,idimp,&
     &npmax,mblok,nblok,kstrt,nvp,ipbc,ierr)
            if ((ierr.eq.0).and.(nsran /= 0)) then
               call PRDISTR32(part,nps,FLDISTR1,ampx,sxi,shiftx,FLDISTR1&
     &,ampy,syi,shifty,FLDISTR1,ampz,szi,shiftz,npx,npy,npz,nx,ny,nz,idi&
     &mp,npmax,mnblok,kstrt,nvp,ipbc)
            endif
! sinusoidal density
         else if (ndpro==2) then
            call PFDISTR32(part,nps,FSDISTR1,ampx,sxi,shiftx,FSDISTR1,am&
     &py,syi,shifty,FSDISTR1,ampz,szi,shiftz,npx,npy,npz,nx,ny,nz,idimp,&
     &npmax,mblok,nblok,kstrt,nvp,ipbc,ierr)
            if ((ierr.eq.0).and.(nsran /= 0)) then
               call PRDISTR32(part,nps,FSDISTR1,ampx,sxi,shiftx,FSDISTR1&
     &,ampy,syi,shifty,FSDISTR1,ampz,szi,shiftz,npx,npy,npz,nx,ny,nz,idi&
     &mp,npmax,mnblok,kstrt,nvp,ipbc)
            endif
! gaussian density
         else if (ndpro==3) then
            call PFDISTR32(part,nps,FGDISTR1,ampx,sxi,shiftx,FGDISTR1,am&
     &py,syi,shifty,FGDISTR1,ampz,szi,shiftz,npx,npy,npz,nx,ny,nz,idimp,&
     &npmax,mblok,nblok,kstrt,nvp,ipbc,ierr)
            if ((ierr.eq.0).and.(nsran /= 0)) then
               call PRDISTR32(part,nps,FGDISTR1,ampx,sxi,shiftx,FGDISTR1&
     &,ampy,syi,shifty,FGDISTR1,ampz,szi,shiftz,npx,npy,npz,nx,ny,nz,idi&
     &mp,npmax,mnblok,kstrt,nvp,ipbc)
            endif
! hyperbolic secant squared density
         else if (ndpro==4) then
            call PFDISTR32(part,nps,FHDISTR1,ampx,sxi,shiftx,FHDISTR1,am&
     &py,syi,shifty,FHDISTR1,ampz,szi,shiftz,npx,npy,npz,nx,ny,nz,idimp,&
     &npmax,mblok,nblok,kstrt,nvp,ipbc,ierr)
            if ((ierr.eq.0).and.(nsran /= 0)) then
               call PRDISTR32(part,nps,FHDISTR1,ampx,sxi,shiftx,FHDISTR1&
     &,ampy,syi,shifty,FHDISTR1,ampz,szi,shiftz,npx,npy,npz,nx,ny,nz,idi&
     &mp,npmax,mnblok,kstrt,nvp,ipbc)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         end subroutine ipfdistr32
!
         subroutine ipvdistr32(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,&
     &npy,npz,kstrt,nvp,mblok)
! calculates initial particle velocities in 3d
! with maxwellian velocity with drift
         implicit none
         integer :: npx, npy, npz, kstrt, nvp, mblok
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp, nps
! local data
         integer :: idimp, npmax, nblok, ierr
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)/mblok
         call PVDISTR32(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,npz&
     &,idimp,npmax,mblok,nblok,kstrt,nvp,ierr)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         end subroutine ipvdistr32
!
         subroutine ipvfdistr32(part,nps,ampx,scalex,shiftx,ampy,scaley,&
     &shifty,ampz,scalez,shiftz,npx,npy,npz,nx,ny,nz,kstrt,nvp,ipbc,ndpr&
     &o,nsran,mblok)
! calculates initial particle co-ordinates in 3d
! with various density profiles
! using parallel random number generator
         implicit none
         integer :: npx, npy, npz, nx, ny, nz, kstrt, nvp, ipbc, ndpro
         integer :: nsran, mblok
         real :: ampx, scalex, shiftx, ampy, scaley, shifty
         real :: ampz, scalez, shiftz
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: nps
! local data
         integer, parameter :: ndv = 256
         integer :: idimp, npmax, mnblok, nblok, nvrp, ierr
         real :: sxi, syi, szi, zero
         real, external :: FLDISTR1, FSDISTR1, FGDISTR1, FHDISTR1
         double precision, dimension((ndv-1)/nvp+1,size(part,3)) :: vran&
     &x, vrany, vranz
         idimp = size(part,1); npmax = size(part,2)
         mnblok = size(part,3); nblok = mnblok/mblok
         nvrp = (ndv - 1)/nvp + 1
         sxi = 0.
         if (scalex /= 0.) sxi = 1.0/scalex
         syi = 0.
         if (scaley /= 0.) syi = 1.0/scaley
         szi = 0.
         if (scalez /= 0.) szi = 1.0/scalez
         zero = 0.0
! uniform density
         if (ndpro==0) then
            call PFDISTR32(part,nps,FLDISTR1,zero,zero,zero,FLDISTR1,zer&
     &o,zero,zero,FLDISTR1,zero,zero,zero,npx,npy,npz,nx,ny,nz,idimp,npm&
     &ax,mblok,nblok,kstrt,nvp,ipbc,ierr)
            if ((ierr.eq.0).and.(nsran /= 0)) then
               call PVRDISTR32(part,nps,FLDISTR1,zero,zero,zero,FLDISTR1&
     &,zero,zero,zero,FLDISTR1,zero,zero,zero,npx,npy,npz,nx,ny,nz,idimp&
     &,npmax,mnblok,vranx,vrany,vranz,kstrt,nvp,ndv,nvrp,ipbc,ierr)
            endif
! linear density
         else if (ndpro==1) then
            call PFDISTR32(part,nps,FLDISTR1,ampx,sxi,shiftx,FLDISTR1,am&
     &py,syi,shifty,FLDISTR1,ampz,szi,shiftz,npx,npy,npz,nx,ny,nz,idimp,&
     &npmax,mblok,nblok,kstrt,nvp,ipbc,ierr)
            if ((ierr.eq.0).and.(nsran /= 0)) then
               call PVRDISTR32(part,nps,FLDISTR1,ampx,sxi,shiftx,FLDISTR&
     &1,ampy,syi,shifty,FLDISTR1,ampz,szi,shiftz,npx,npy,npz,nx,ny,nz,id&
     &imp,npmax,mnblok,vranx,vrany,vranz,kstrt,nvp,ndv,nvrp,ipbc,ierr)
            endif
! sinusoidal density
         else if (ndpro==2) then
            call PFDISTR32(part,nps,FSDISTR1,ampx,sxi,shiftx,FSDISTR1,am&
     &py,syi,shifty,FSDISTR1,ampz,szi,shiftz,npx,npy,npz,nx,ny,nz,idimp,&
     &npmax,mblok,nblok,kstrt,nvp,ipbc,ierr)
            if ((ierr.eq.0).and.(nsran /= 0)) then
               call PVRDISTR32(part,nps,FSDISTR1,ampx,sxi,shiftx,FSDISTR&
     &1,ampy,syi,shifty,FSDISTR1,ampz,szi,shiftz,npx,npy,npz,nx,ny,nz,id&
     &imp,npmax,mnblok,vranx,vrany,vranz,kstrt,nvp,ndv,nvrp,ipbc,ierr)
            endif
! gaussian density
         else if (ndpro==3) then
            call PFDISTR32(part,nps,FGDISTR1,ampx,sxi,shiftx,FGDISTR1,am&
     &py,syi,shifty,FGDISTR1,ampz,szi,shiftz,npx,npy,npz,nx,ny,nz,idimp,&
     &npmax,mblok,nblok,kstrt,nvp,ipbc,ierr)
            if ((ierr.eq.0).and.(nsran /= 0)) then
               call PVRDISTR32(part,nps,FGDISTR1,ampx,sxi,shiftx,FGDISTR&
     &1,ampy,syi,shifty,FGDISTR1,ampz,szi,shiftz,npx,npy,npz,nx,ny,nz,id&
     &imp,npmax,mnblok,vranx,vrany,vranz,kstrt,nvp,ndv,nvrp,ipbc,ierr)
            endif
! hyperbolic secant squared density
         else if (ndpro==4) then
            call PFDISTR32(part,nps,FHDISTR1,ampx,sxi,shiftx,FHDISTR1,am&
     &py,syi,shifty,FHDISTR1,ampz,szi,shiftz,npx,npy,npz,nx,ny,nz,idimp,&
     &npmax,mblok,nblok,kstrt,nvp,ipbc,ierr)
            if ((ierr.eq.0).and.(nsran /= 0)) then
               call PVRDISTR32(part,nps,FHDISTR1,ampx,sxi,shiftx,FHDISTR&
     &1,ampy,syi,shifty,FHDISTR1,ampz,szi,shiftz,npx,npy,npz,nx,ny,nz,id&
     &imp,npmax,mnblok,vranx,vrany,vranz,kstrt,nvp,ndv,nvrp,ipbc,ierr)
            endif
         endif
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         end subroutine ipvfdistr32
!
         subroutine ipvvdistr32(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx&
     &,npy,npz,kstrt,nvp,mblok)
! calculates initial particle velocities in 3d
! with maxwellian velocity with drift
! using parallel random number generator
         implicit none
         integer :: npx, npy, npz, kstrt, nvp, mblok
         real :: vtx, vty, vtz, vdx, vdy, vdz
         real, dimension(:,:,:), pointer :: part
         integer, dimension(:), pointer :: npp, nps
! local data
         integer, parameter :: ndv = 256
         integer :: idimp, npmax, nblok, nvrp, ierr
         double precision, dimension((ndv-1)/nvp+1,size(part,3)) :: vran&
     &x, vrany, vranz
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)/mblok
         nvrp = (ndv - 1)/nvp + 1
         call PVVISTR32(part,npp,nps,vtx,vty,vtz,vdx,vdy,vdz,npx,npy,npz&
     &,idimp,npmax,mblok,nblok,vranx,vrany,vranz,kstrt,nvp,ndv,nvrp,ierr&
     &)
         if (ierr /= 0) then
            call MP_END
            call PPEXIT
            stop
         endif
         end subroutine ipvvdistr32
!
         subroutine ipbdistr32(part,bxyz,npp,noff,qbm,nx,ny,nz,ipbc,mblo&
     &k,inorder)
! reinterprets curent particle positions as positions of guiding centers
! and calculates the actual particle positions for 3d
         implicit none
         integer :: nx, ny, nz, ipbc, mblok
         integer, optional :: inorder
         real :: qbm
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:,:), pointer :: bxyz
         integer, dimension(:), pointer :: npp
         integer, dimension(:,:), pointer :: noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, nzpmx, idds
         integer :: order
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)/mblok; idds = size(noff,1)
         nxv = size(bxyz,2); nypmx = size(bxyz,3); nzpmx = size(bxyz,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PGBDISTR32L(part,bxyz(1,1,1,1,1),npp,noff,qbm,nx,ny,nz,&
     &idimp,npmax,mblok,nblok,nxv,nypmx,nzpmx,idds,ipbc)
         else
            call PGBDISTR32L(part,bxyz(1,2,2,2,1),npp,noff,qbm,nx,ny,nz,&
     &idimp,npmax,mblok,nblok,nxv,nypmx,nzpmx,idds,ipbc)
         endif
         end subroutine ipbdistr32
!
         subroutine iprbdistr32(part,bxyz,npp,noff,qbm,ci,nx,ny,nz,ipbc,&
     &mblok,inorder)
! reinterprets curent particle positions as positions of guiding centers
! and calculates the actual particle positions for relativistic 3d
         implicit none
         integer :: nx, ny, nz, ipbc, mblok
         integer, optional :: inorder
         real :: qbm, ci
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:,:,:), pointer :: bxyz
         integer, dimension(:), pointer :: npp
         integer, dimension(:,:), pointer :: noff
! local data
         integer :: idimp, npmax, nblok, nxv, nypmx, nzpmx, idds
         integer :: order
         idimp = size(part,1); npmax = size(part,2)
         nblok = size(part,3)/mblok; idds = size(noff,1)
         nxv = size(bxyz,2); nypmx = size(bxyz,3); nzpmx = size(bxyz,4)
         order = QUADRATIC
         if (present(inorder)) order = inorder
         if (order==LINEAR) then
            call PGRBDISTR32L(part,bxyz(1,1,1,1,1),npp,noff,qbm,ci,nx,ny&
     &,nz,idimp,npmax,mblok,nblok,nxv,nypmx,nzpmx,idds,ipbc)
         else
            call PGRBDISTR32L(part,bxyz(1,2,2,2,1),npp,noff,qbm,ci,nx,ny&
     &,nz,idimp,npmax,mblok,nblok,nxv,nypmx,nzpmx,idds,ipbc)
         endif
         end subroutine iprbdistr32
!
         subroutine ifedges32(edges,noff,nyzp,ampy,scaley,shifty,ampz,sc&
     &alez,shiftz,ny,nz,kstrt,nvpy,nvpz,nypmx,nzpmx,ipbc,ndpro,mblok,mte&
     &rg,nterg,ierr,inorder)
! finds new 3d partitions from initial analytic distribution function
         implicit none
         integer :: ny, nz, kstrt, nvpy, nvpz, nypmx, nzpmx, ipbc, ndpro
         integer :: mblok, mterg, nterg, ierr
         real :: ampy, scaley, shifty, ampz, scalez, shiftz
         integer, optional :: inorder
         real, dimension(:,:), pointer :: edges
         integer, dimension(:,:), pointer :: noff, nyzp
! local data
         integer :: idps, idds, nblok, nypmin, nypmax, nzpmin, nzpmax
         integer :: order
         real :: syi, szi, zero
         real, external :: FLDISTR1, FSDISTR1, FGDISTR1, FHDISTR1
         idps = size(edges,1); nblok = size(edges,2)/mblok
         idds = size(noff,1)
         syi = 0.
         if (scaley /= 0.) syi = 1.0/scaley
         szi = 0.
         if (scalez /= 0.) szi = 1.0/scalez
         zero = 0.0
         ierr = 0
         order = QUADRATIC
         if (present(inorder)) order = inorder
! uniform density
         if (ndpro==0) then
            call FEDGES32(edges,noff,nyzp,FLDISTR1,zero,zero,zero,FLDIST&
     &R1,zero,zero,zero,ny,nz,nypmin,nypmax,nzpmin,nzpmax,kstrt,nvpy,nvp&
     &z,mblok,nblok,idps,idds,ipbc)
! linear density
         else if (ndpro==1) then
            call FEDGES32(edges,noff,nyzp,FLDISTR1,ampy,syi,shifty,FLDIS&
     &TR1,ampz,szi,shiftz,ny,nz,nypmin,nypmax,nzpmin,nzpmax,kstrt,nvpy,n&
     &vpz,mblok,nblok,idps,idds,ipbc)
! sinusoidal density
         else if (ndpro==2) then
            call FEDGES32(edges,noff,nyzp,FSDISTR1,ampy,syi,shifty,FSDIS&
     &TR1,ampz,szi,shiftz,ny,nz,nypmin,nypmax,nzpmin,nzpmax,kstrt,nvpy,n&
     &vpz,mblok,nblok,idps,idds,ipbc)
! gaussian density
         else if (ndpro==3) then
            call FEDGES32(edges,noff,nyzp,FGDISTR1,ampy,syi,shifty,FGDIS&
     &TR1,ampz,szi,shiftz,ny,nz,nypmin,nypmax,nzpmin,nzpmax,kstrt,nvpy,n&
     &vpz,mblok,nblok,idps,idds,ipbc)
! hyperbolic secant squared density
         else if (ndpro==4) then
            call FEDGES32(edges,noff,nyzp,FHDISTR1,ampy,syi,shifty,FHDIS&
     &TR1,ampz,szi,shiftz,ny,nz,nypmin,nypmax,nzpmin,nzpmax,kstrt,nvpy,n&
     &vpz,mblok,nblok,idps,idds,ipbc)
         endif
         if (order==LINEAR) then
            nypmax = nypmax + 1
            nzpmax = nzpmax + 1
         else
            nypmax = nypmax + 3
            nzpmax = nzpmax + 3
         endif
         if ((nypmin.lt.1).or.(nypmax.gt.nypmx)) then
            write (2,*) 'Field size error: nypmin,nypmax=',nypmin,nypmax
            ierr = 1
         endif
         if ((nzpmin.lt.1).or.(nzpmax.gt.nzpmx)) then
            write (2,*) 'Field size error: nzpmin,nzpmax=',nzpmin,nzpmax
            ierr = 2
         endif
         mterg = nypmin - 1
         nterg = nzpmin - 1
         end subroutine ifedges32
!
      end module pinit32d
