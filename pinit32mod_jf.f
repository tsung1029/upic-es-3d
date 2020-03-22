!-----------------------------------------------------------------------
!
      module pinit32d_jf
!
! written by Jay Fahlen as new input parameters and namelists for
! Viktor Decyk's 2d beps code.
! written Jan. 2008

      implicit none
      private
      public :: ntfield,uniform_density_init,uniform_random_init
      public :: pinput3_jf,sendnml_jf
      public :: amp,wavemode,wavew,angle
      public :: fvxmax,fvymax,fvzmax
      public :: nphbx,nphby,nphbz
      public :: nphxx,nphyx,nphzx,nphxy,nphyy,nphzy,nphxz,nphyz,nphzz
      public :: driver_select,rise,flat,fall,timerise,timeflat,timefall,phase_offset
      public :: yrise_fall,linepos,linepos_3d,ntlines,center1,center2
      public :: phsl_x_pos,phsl_x_thick,ntphsl_x,num_phsl_x
      public :: phsl_y_pos,phsl_y_thick,ntphsl_y,num_phsl_y
      public :: ntden,ntpene,ntj,ntvdotE
      public :: ntraw,nttrack,nt_dump_track,raw_fraction,raw_range,raw_dump_tagged,track_no_h5
      public :: ntESPoynt,ntESPoynt_int
      public :: keep_init_vel,initEneRange,finalEneRange,nEneBin_init,nEneBin_final
      public :: initVxLow,initVxHigh,finalVxLow,finalVxHigh,nVxBin_init,nVxBin_final
      public :: initVyLow,initVyHigh,finalVyLow,finalVyHigh,nVyBin_init,nVyBin_final
      public :: nvdotE_part,nvdotE_int,fvx_xy_max,nphb_xy,ntphxy,fvy_xy_max
      public :: fvxy_xy_xrange,fvxy_xy_yrange,nfvxy_xy_xrange,nfvxy_xy_yrange
      public :: nvdotE_follow_part,ntfield_ene_int
      public :: dump_start,dump_end,nt_div_ESPoynt_int,nt_div_ESPoynt,nt_phase_vx_vy
      public :: phase_keep_only_tracked,nt_write_field_ene_in_driver_region,int_pos,nt_write_ene_int
      public :: nt_write_Py_vs_y,nt_write_div_P_vs_y,nt_write_jE_sumover_x,nt_write_Py_sumover_x
      public :: nt_write_U_sumover_x, npfac
      public :: nt_bin_E,bin_E,bin_E_xrange,bin_E_yrange,bin_E_Exrange,bin_E_Eyrange,raw_cycle_x_move
      public :: write_stride
      
      save
      
      integer :: ntfield = 0,wavemode=0
      ! write_stride specifiecs the stride, or number of data points to skip, for each dimension
      ! of the 3D data to reduce the file size
      integer, dimension(3) :: write_stride
      real :: amp = 0., wavew = 0.,angle = 0.
      ! these three are the ranges for the phase space velocities
      ! the phase space will be from v = -fvxmax to fvxmax
      real :: fvxmax = 10., fvymax = 10., fvzmax = 10.
      ! nphbx,nphby,nphbz = # phase space bins in vx,vy,vz
      ! nphxx,nphyx,nphxy,nphyy etc = dump factor for phase space in vx vs. x, vx vs. y, and so on
      integer :: nphbx = 10, nphby = 10, nphbz = 10
      integer :: nphxx = 0, nphyx = 0, nphzx = 0
      integer :: nphxy = 0, nphyy = 0, nphzy = 0
      integer :: nphxz = 0, nphyz = 0, nphzz = 0
      ! driver_select = 0 -> no external driver
      ! driver_select = 1 -> plane wave driver
      ! driver_select = 2 -> plane finite wavelength driver
      integer :: driver_select = 0
	  	! rise, flat, and fall times for finite wavelength driver
			! these are given in # of wavelengths!
			! so, for example, rise = 3 means the rise distance is 3*lambda
      real :: rise=1.,flat=1.,fall=1.
			! rise time flat time and fall time
      real :: timerise = 0.,timeflat = 1000000., timefall = 0.
      ! yrise_fall =  spatial rise and fall for the transverse profile
      ! phase_offset is the distance that the phase of the center of the transverse profile is 
      !  offset from the top and bottom edges, for use with the circular wavefront driver
      real :: yrise_fall, phase_offset
      !these two specify the center of the two drivers for doub_gauss_tran_finite_wavelen
      integer :: center1,center2
      ! nlines is the number of evenly spaced lines of E(x=0:nx,y=i*ny/nlines)
      ! where i goes from 0 to nlines-1
      ! ntlines specifies how often to dump the lines data
      integer,dimension(5) :: linepos !Don't use in this code!  just there so I dont have to renumber below
      integer,dimension(5,2) :: linepos_3d
      integer :: ntlines = 0
      ! phsl_x_pos is an array that lists the positions for each thin phasespace dump
      ! phsl_x_thick is the thickness of the sum over the y direction for the
      ! narrow phasespace.  If the range over y to sum falls across processors then 
      ! it will only get those particles in the lowest processor.  Only five places
      ! work for now, if values is set to -1 then don't do it.
      ! those with x in the name are a lineout in y, with y in the name a lineout in x
      integer,dimension(5) :: phsl_x_pos,phsl_y_pos
      integer :: phsl_x_thick=2,ntphsl_x = 0,num_phsl_x=0
      integer :: phsl_y_thick=2,ntphsl_y = 0,num_phsl_y=0
      ! ntden is how often to dump density, ntpene is particle energy dump,
      ! ntj is current dump
      ! ntvdotE is how often to dumb v dot E
      integer :: ntden=0,ntpene=0,ntj=0,ntvdotE=0
      ! ntraw is how often to dump raw particle data, nt_dump_track is how often to
      ! flush the track data from the buffer to disk
      ! raw_fraction is fraction of
      ! particles to dump data for, nttrack adds the tag data to track the particles
      ! raw_dump_tagged means dump to raw file only tagged particles, 0 off, 1 on
      ! it cannot be used with nttrack
      ! track_no_h5 writes track data in plain Fortran unformatted to be converted to H5 in postprocess
      integer :: ntraw = 0, nttrack = 0, nt_dump_track = 0,track_no_h5=1
      ! raw_dump_tagged means dump to raw file only tagged particles, 0 off, 1 on
      ! it cannot be used with nttrack, file input_tags must be present to read list of
      ! particles to dump.
      ! raw_cycle_x_move moves the window specified by raw_range at the specified speed 
      ! in the x direction
      real :: raw_fraction=0., raw_dump_tagged = 0, raw_cycle_x_move = -1.
      ! raw_range is the 4D cube inside of which the raw data is stored for the particles
      ! raw_range(1,*) is the lower and upper boundary in x
      ! raw_range(2,*) is the lower and upper boundary in y
      ! raw_range(3,*) is the lower and upper boundary in vx
      ! raw_range(4,*) is the lower and upper boundary in vy
      real,dimension(4,2) :: raw_range = -1.
      ! ntESPoynt is the electrostatic component of the Poynting vector 
      ! See VK Decyk Phys. Fluids 25, 1205 (1982)
      integer :: ntESPoynt = 0, ntESPoynt_int = 0
      ! keep_init_vel = 1 means to add two more arrays to particle array that 
      ! store the initial velocity of the particles so comparisons and statistics
      ! can be calculated at the end of the run
      integer :: keep_init_vel = 0
      ! the EneRanges are the maximum energies to bin for the initial and final energy diagnostic
      real :: initEneRange,finalEneRange
      ! nEneBin_init and nEneBin_final are the number of bins for the diagnostic above
      integer :: nEneBin_init,nEneBin_final
      ! same as Ene stuff above, but for vx and must give high and low values of vx to save
      real :: initVxLow,initVxHigh,finalVxLow,finalVxHigh
      real :: initVyLow,initVyHigh,finalVyLow,finalVyHigh
      integer :: nVxBin_init,nVxBin_final
      integer :: nVyBin_init,nVyBin_final
      ! These two calculate the v dot E for each particle and store it in the particle array.
      ! For diagnostics, this is then deposited just like the current would be deposited.  
      ! These are better diagnostics than the other vdotE (ntvdotE) and are meant to replace it.
      ! nvdotE_part dumps vdot E at specified interval, while nvdotE_int sums the vdotE at every
      ! timestep  and then dumps the running sum at spedified intervals.
      ! nvdotE_part MUST be specified to use nvdotE_int, but if one is nonzero, then BOTH must be.
      integer :: nvdotE_part=0,nvdotE_int=0
      real :: fvx_xy_max=10.,fvy_xy_max=10.
      ! These specifiy the spatial extent of the phase space
      real,dimension(2) :: fvxy_xy_xrange=(/-1.,-1./),fvxy_xy_yrange=(/-1.,-1./)
      ! These are the number of velocity bins, timestep skip, number of bins in x and y
      integer :: nphb_xy,ntphxy=0,nfvxy_xy_xrange=0,nfvxy_xy_yrange=0
			! if nvdotE_follow_part is set, then the particles will carry their total vdotE rather than
			! their instantaneous vdotE.  This cannot be used with either nvdotE_int or nvdotE_part
      integer :: nvdotE_follow_part=0
      integer :: ntfield_ene_int = 0		!this will write out the integrated field energy
      ! These two are to specify at what time to start and stop dumping diagnostic data
      ! So far only works for e-field and raw, negative means dump from beginning to end
      real :: dump_start=-1.,dump_end=-1.
      ! This dumps the divergence and integrated divergence of the ES Poynting vector
      integer :: nt_div_ESPoynt = 0, nt_div_ESPoynt_int = 0
      ! nt_phase_vx_vy plots vx vs. vy phase space
      integer :: nt_phase_vx_vy
      ! phase_keep_only_tracked, this deposits only the tracked particles to the phase space 
      ! diagnostics.  Therefore, must be used with the tracking stuff on, but just set them
      ! so that they are too big to ever actually be written.
      integer :: phase_keep_only_tracked
      ! nt_write_field_ene_in_driver_region writes electric field energy in the region:
      ! -yrise_fall < y - ny/2 < yrise_fall and all x to the formatted file on unit 77
      integer :: nt_write_field_ene_in_driver_region = 0
      ! int_pos are the y positions of the three boxes used in integrating the Poynting flow energy diagnostic
      ! int_pos(1) is y pos of bottom of bottom box, int_pos(2) is top of bottom box and bottom of middle box
      ! and so on
      integer,dimension(4) :: int_pos=(/-1.,-1.,-1.,-1./)
      ! nt_write_ene_int writes to file the integrated energies for the boxes indicated in int_pos
      integer :: nt_write_ene_int = 0
      ! nt_write_Py_vs_y sums Py over all x and writes it for each y position
      integer :: nt_write_Py_vs_y = 0, nt_write_div_P_vs_y = 0
      ! nt_write_jE_sumover_x sums jE over all x and writes it for each y position
      ! nt_write_Py_sumover_x sums Py over all x and writes it for each y position
      integer :: nt_write_jE_sumover_x = 0,nt_write_Py_sumover_x = 0, nt_write_U_sumover_x = 0
      ! npfac is the percentage increase to np for allocating the part array to ensure that when
      ! new particles come onto the current processor the part array isn't over filled.
      real :: npfac = 0.05
      ! nt_bin_E writes Ex and Ey as 3d data, bin_E is the number of bins in Ex and Ey
      ! bin_E_xrange and bin_E_yrange are the space ranges over which to do the binning, must be
      ! 	cells, cannot do finer or coarser grid than the simulation grid
      ! bin_E_Exrange,bin_E_Eyrange are the ranges in Ex and Ey to do it
      integer :: nt_bin_E = 0, bin_E = 100
      integer, dimension(2) :: bin_E_xrange = (/0,1/), bin_E_yrange = (/0,1/)
      real, dimension(2) :: bin_E_Exrange = (/-1.,1./), bin_E_Eyrange = (/-1.,1./)
      
      namelist /pinput3_jf/ ntfield, amp,wavemode,wavew,fvxmax,fvymax,&
      	&nphbx,nphby,nphxx,nphyx,rise,flat,fall,driver_select,timerise,&
      	&timeflat,timefall,yrise_fall,linepos,ntlines,phsl_x_pos,phsl_x_thick,ntphsl_x,&
      	&phsl_y_pos,phsl_y_thick,ntphsl_y,ntden,ntpene,ntj,ntvdotE, ntraw,nttrack, &
      	&raw_fraction, nt_dump_track,ntESPoynt,angle, keep_init_vel,initEneRange,finalEneRange,&
      	&nEneBin_init,nEneBin_final,initVxLow,initVxHigh,finalVxLow,finalVxHigh,nVxBin_init,&
      	&nVxBin_final,initVyLow,initVyHigh,finalVyLow,finalVyHigh,nVyBin_init,nVyBin_final,&
      	&nvdotE_part,nvdotE_int,fvx_xy_max,nphb_xy,ntphxy,fvy_xy_max,fvxy_xy_xrange,fvxy_xy_yrange,&
      	&nfvxy_xy_xrange,nfvxy_xy_yrange,nvdotE_follow_part,raw_range,ntESPoynt_int,ntfield_ene_int,&
      	&center1,center2,raw_dump_tagged,dump_start,dump_end,nt_div_ESPoynt_int,track_no_h5,phase_offset,&
      	&nt_phase_vx_vy,nphxy,nphyy,phase_keep_only_tracked,nt_write_field_ene_in_driver_region,int_pos,&
      	&nt_write_ene_int,nt_write_Py_vs_y,nt_div_ESPoynt,nt_write_div_P_vs_y,nt_write_jE_sumover_x,&
      	&nt_write_Py_sumover_x, nt_write_U_sumover_x, npfac, nt_bin_E,bin_E,bin_E_xrange,bin_E_yrange,&
      	&bin_E_Exrange,bin_E_Eyrange,raw_cycle_x_move,write_stride,linepos_3d,&
      	&fvzmax,nphbz,nphzx,nphzy,nphxz,nphyz,nphzz
      	      
      contains

!***************************************
! This function needs to be tested more!      
!***************************************
			subroutine uniform_random_init(part,npp,nps,nx,ny,nz,npx,npy,npz,&
				&	kstrt,nvp,mblok,nvpy,nvpz)
				! calculates initial particle positions in 3d
				! with random positions
				! using parallel random number generator
				! taken from Viktor's velocity distribution function and modified by Jay Fahlen
				implicit none
				integer :: nx,ny,nz,nvpy,nvpz
				integer :: npx, npy, npz, kstrt, nvp, mblok
				real, dimension(:,:,:), pointer :: part
				integer, dimension(:), pointer :: npp, nps
				! local data
				integer, parameter :: ndv = 256
				integer :: idimp, npmax, nblok, nvrp, ierr
				real :: procy,procz
				double precision, dimension(:,:),allocatable :: posx, posy, posz
				allocate(posx((ndv-1)/nvp+1,size(part,3)),posy((ndv-1)/nvp+1,size(part,3)))
				allocate(posz((ndv-1)/nvp+1,size(part,3)))
				if (kstrt == 1) print*,"A"
				idimp = size(part,1); npmax = size(part,2)
				nblok = size(part,3)/mblok
				nvrp = (ndv - 1)/nvp + 1
				procy = real(mod(kstrt-1,nvpy) * ny / nvpy)
				procz = real((kstrt-1) / nvpy * nz / nvpz)
				if (kstrt == 1) print*,"B",npmax
				call UNIF_RNDM_INIT(part,npp,nps,0.,procy,procz,real(nx),real(ny/nvpy),real(nz/nvpz),&
				&npx,npy,npz,idimp,npmax,mblok,nblok,posx,posy,posz,kstrt,nvp,ndv,nvrp,ierr)
				if (kstrt == 1) print*,"C",ierr
				if (ierr /= 0) then
					call MP_END
					call PPEXIT
					stop
				endif
				deallocate(posx,posy,posz)
			end subroutine uniform_random_init
      
			subroutine sendnml_jf()
				implicit none
				integer,parameter :: lenml = 145
				double precision, dimension(lenml) :: ddata
				ddata(1) = ntfield
				ddata(2) = amp
				ddata(3) = wavemode
				ddata(4) = wavew
				ddata(5) = fvxmax
				ddata(6) = fvymax
				ddata(7) = nphbx
				ddata(8) = nphby
				ddata(9) = nphxx
				ddata(10) = nphyx
				ddata(11) = driver_select
				ddata(12) = rise
				ddata(13) = flat
				ddata(14) = fall
				ddata(15) = timerise
				ddata(16) = timeflat
				ddata(17) = yrise_fall
				ddata(18:22) = linepos
				ddata(23) = ntlines
				ddata(24:28) = phsl_x_pos
				ddata(29) = phsl_x_thick
				ddata(30) = ntphsl_x
				ddata(31) = ntden
				ddata(32) = ntpene
				ddata(33) = ntj
				ddata(34) = ntvdotE
				ddata(35:39) = phsl_y_pos
				ddata(40) = phsl_y_thick
				ddata(41) = ntphsl_y
				ddata(42) = ntraw
				ddata(43) = nttrack
				ddata(44) = raw_fraction
				ddata(45) = nt_dump_track
				ddata(46) = ntESPoynt
				ddata(47) = angle
				ddata(48) = keep_init_vel
				ddata(49) = initEneRange
				ddata(50) = finalEneRange
				ddata(51) = nEneBin_init
				ddata(52) = nEneBin_final
				ddata(53) = initVxLow
				ddata(54) = initVxHigh
				ddata(55) = finalVxLow
				ddata(56) = finalVxHigh
				ddata(57) = nVxBin_init
				ddata(58) = nVxBin_final
				ddata(59) = initVyLow
				ddata(60) = initVyHigh
				ddata(61) = finalVyLow
				ddata(62) = finalVyHigh
				ddata(63) = nVyBin_init
				ddata(64) = nVyBin_final
				ddata(65) = nvdotE_part
				ddata(66) = nvdotE_int
				ddata(67) = fvx_xy_max
				ddata(68) = nphb_xy
				ddata(69) = ntphxy
				ddata(70) = fvy_xy_max
				ddata(71:72) = fvxy_xy_xrange
				ddata(73:74) = fvxy_xy_yrange
				ddata(75) = nfvxy_xy_xrange
				ddata(76) = nfvxy_xy_yrange
				ddata(77) = nvdotE_follow_part
				ddata(78:79) = raw_range(1,:)
				ddata(80:81) = raw_range(2,:)
				ddata(82:83) = raw_range(3,:)
				ddata(84:85) = raw_range(4,:)
				ddata(86) = ntESPoynt_int
				ddata(87) = ntfield_ene_int
				ddata(88) = center1
				ddata(89) = center2
				ddata(90) = raw_dump_tagged
				ddata(91) = dump_start
				ddata(92) = dump_end
				ddata(93) = nt_div_ESPoynt_int
				ddata(94) = track_no_h5
				ddata(95) = phase_offset
				ddata(96) = nt_phase_vx_vy
				ddata(97) = nphxy
				ddata(98) = nphyy
				ddata(99) = phase_keep_only_tracked
				ddata(100) = nt_write_field_ene_in_driver_region
				ddata(101:104) = int_pos
				ddata(105) = nt_write_ene_int
				ddata(106) = nt_write_Py_vs_y
				ddata(107) = nt_div_ESPoynt
				ddata(108) = nt_write_div_P_vs_y
				ddata(109) = nt_write_jE_sumover_x
				ddata(110) = nt_write_Py_sumover_x
				ddata(111) = nt_write_U_sumover_x
				ddata(112) = timefall
				ddata(113) = npfac
				ddata(114) = nt_bin_E
				ddata(115) = bin_E
				ddata(116:117) = bin_E_xrange
				ddata(118:119) = bin_E_yrange
				ddata(120:121) = bin_E_Exrange
				ddata(122:123) = bin_E_Eyrange
				ddata(124) = raw_cycle_x_move
				ddata(125:127) = write_stride
				ddata(128:132) = linepos_3d(:,1)
				ddata(133:137) = linepos_3d(:,2)
				ddata(139) = fvzmax
				ddata(140) = nphbz
				ddata(141) = nphzx
				ddata(142) = nphzy
				ddata(143) = nphxz
				ddata(144) = nphyz
				ddata(145) = nphzz
				
				call PBCAST(ddata,lenml)
				ntfield = ddata(1)
				amp = ddata(2)
				wavemode = ddata(3)
				wavew = ddata(4)
				fvxmax = ddata(5)
				fvymax = ddata(6)
				nphbx = ddata(7)
				nphby = ddata(8)
				nphxx = ddata(9)
				nphyx = ddata(10)
				driver_select = ddata(11)
				rise = ddata(12)
				flat = ddata(13)
				fall = ddata(14)
				timerise = ddata(15)
				timeflat = ddata(16)
				yrise_fall = ddata(17)
				linepos = ddata(18:22)
				ntlines = ddata(23)
				phsl_x_pos = ddata(24:28)
				phsl_x_thick = ddata(29)
				ntphsl_x = ddata(30)
				ntden = ddata(31)
				ntpene = ddata(32)
				ntj = ddata(33)
				ntvdotE = ddata(34)
				phsl_y_pos = ddata(35:39)
				phsl_y_thick = ddata(40)
				ntphsl_y = ddata(41)
				ntraw = ddata(42)
				nttrack = ddata(43)
				raw_fraction = ddata(44)
				nt_dump_track = ddata(45)
				ntESPoynt = ddata(46)
				angle = ddata(47)
				keep_init_vel = ddata(48)
				initEneRange = ddata(49)
				finalEneRange = ddata(50)
				nEneBin_init = ddata(51)
				nEneBin_final = ddata(52)
				initVxLow = ddata(53)
				initVxHigh = ddata(54)
				finalVxLow = ddata(55)
				finalVxHigh = ddata(56)
				nVxBin_init = ddata(57)
				nVxBin_final = ddata(58)
				initVyLow = ddata(59)
				initVyHigh = ddata(60)
				finalVyLow = ddata(61)
				finalVyHigh = ddata(62)
				nVyBin_init = ddata(63)
				nVyBin_final = ddata(64)
				nvdotE_part = ddata(65)
				nvdotE_int = ddata(66)
				fvx_xy_max = ddata(67)
				nphb_xy = ddata(68)
				ntphxy = ddata(69)
				fvy_xy_max = ddata(70)
				fvxy_xy_xrange = ddata(71:72)
				fvxy_xy_yrange = ddata(73:74)
				nfvxy_xy_xrange = ddata(75)
				nfvxy_xy_yrange = ddata(76)
				nvdotE_follow_part = ddata(77)
				raw_range(1,:) = ddata(78:79)
				raw_range(2,:) = ddata(80:81)
				raw_range(3,:) = ddata(82:83)
				raw_range(4,:) = ddata(84:85)
				ntESPoynt_int = ddata(86)
				ntfield_ene_int = ddata(87)
				center1 = ddata(88)
				center2 = ddata(89)
				raw_dump_tagged = ddata(90)
				dump_start = ddata(91)
				dump_end = ddata(92)
				nt_div_ESPoynt_int = ddata(93)
				track_no_h5 = ddata(94)
				phase_offset = ddata(95)
				nt_phase_vx_vy = ddata(96)
				nphxy = ddata(97)
				nphyy = ddata(98)
				phase_keep_only_tracked = ddata(99)
				nt_write_field_ene_in_driver_region = ddata(100)
				int_pos = ddata(101:104)
				nt_write_ene_int = ddata(105)
				nt_write_Py_vs_y = ddata(106)
				nt_div_ESPoynt = ddata(107)
				nt_write_div_P_vs_y = ddata(108)
				nt_write_jE_sumover_x = ddata(109)
				nt_write_Py_sumover_x = ddata(110)
				nt_write_U_sumover_x = ddata(111)
				timefall = ddata(112)
				npfac = ddata(113)
				nt_bin_E = ddata(114)
				bin_E = ddata(115)
				bin_E_xrange = ddata(116:117)
				bin_E_yrange = ddata(118:119)
				bin_E_Exrange = ddata(120:121)
				bin_E_Eyrange = ddata(122:123)
				raw_cycle_x_move = ddata(124)
				write_stride = ddata(125:127)
				linepos_3d(:,1) = ddata(128:132)
				linepos_3d(:,2) = ddata(133:137)
				fvzmax = ddata(139)
				nphbz = ddata(140)
				nphzx = ddata(141)
				nphzy = ddata(142)
				nphxz = ddata(143)
				nphyz = ddata(144)
				nphzz = ddata(145)

			end subroutine sendnml_jf
      	
			subroutine uniform_density_init(part,nx,ny,nz,npx,npy,npz,nvpy,nvpz,idproc)
				implicit none
				real,dimension(:,:,:) :: part
				integer :: idproc,nx,ny,nz,npx,npy,npz,nvpy,nvpz
				
				integer :: i,j,k, arr, procy, procz, checky, checkz
				real :: dx, dy, dz, x, y, z, lx, ly, lz
							
				!For this function, npx must be evenly divisible by nvp!
				checky = npy / nvpy
				checkz = npz / nvpz
				if (checky * nvpy .ne. npy .or. checkz * nvpz .ne. npz) then
					write (2,*) 'In uniform_density_init, npy must be evenly divisible by nvpy and the',&
					&' same for npz and nvpz.  Exiting...'
					call MP_END
					call PPEXIT
					stop
				endif
								
				lx = nx
				ly = ny / nvpy
				lz = nz / nvpz
				
				dx = real(nx) / real(npx)
				dy = real(ny) / real(npy)
				dz = real(nz) / real(npz)
				
				procy = mod(idproc,nvpy) * ny / nvpy
				procz = idproc / nvpy * nz / nvpz

				arr = 1
				x = 0
				do i = 1, lx / dx
					y = procy
					do j = 1, ly / dy
						z = procz
						do k = 1, lz / dz				
							part(1,arr,1) = x
							part(2,arr,1) = y
							part(3,arr,1) = z

							arr = arr + 1
							z = z + dz
						enddo
						y = y + dy
					enddo
					x = x + dx
				enddo
				return
			end subroutine uniform_density_init

      	
      	
      end module pinit32d_jf
