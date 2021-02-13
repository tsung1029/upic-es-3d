!-----------------------------------------------------------------------
!
module ext_driver32_jf
	use pinit32d_jf
	implicit none
	private
	public :: plane_wave,plane_finite_wavelen,gauss_tran_finite_wavelen
	public :: angle_wave,cross_wave,rect_tran_finite_wavelen
	public :: gauss_tran_per_wavelen,doub_gauss_tran_finite_wavelen
	public :: gauss_tran_per_wavelen_circle_wavefronts
	public :: gauss_tran_per_wavelen_pois,rect_tran_per_wavelen
	public :: do_nothing,supergauss_tran_per_wavelen
	public :: laguerre_gaussian_ponderomotive_force
        public :: two_laguerre_gaussian_beams

	contains
	
	subroutine do_nothing()
	return
	end subroutine do_nothing
	
	subroutine cross_wave(fxye,t,nx,nxe,ny,nypmx,&
		&nvp,idproc)
		integer :: nx,nxe,ny,nypmx,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i,j,stoppos, nyproc,ystart,yend
		integer, dimension(nvp) :: proc_pos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac
		real :: ylow, yhigh, ypos, yhalf
		real :: tempamp, tempamp2,ytemp
		real :: mfallinv,falloffset,r_f,r_f_f, phase
		real :: tan_value
		
!		wavek = real(nx)/wavemode
!		wavek = 6.283185307/wavek
!		r_f = rise + flat
!		r_f_f = rise + flat + fall
!		riseinv = 1. / rise
!		mfallinv = -1. / fall
!		falloffset = rise/fall + flat/fall + 1.
!		tfac = 1.
!		
!		nyproc = ny / nvp
!		do i = 0, nvp-1
!			proc_pos(i+1) = ny / nvp * (i)
!		enddo
!
!		tfac = 1.
!		if (t < timerise) then
!			tfac = t / timerise
!		else if (t < timerise + timeflat) then
!			tfac = 1.
!		else if (t < timerise + timeflat + timefall) then
!			tfac = 1. - (t - (timerise+timeflat))/timefall
!		endif
!		
!		if (t < timerise + timeflat + timefall) then
!				
!			ylow = real(ny/2-yrise_fall)
!			yhigh = real(ny/2+yrise_fall)
!			yhalf = real(ny/2)
!			tan_value = tan(angle)
!			do j = 1,nypmx
!		
!				tempamp = amp * tfac
!		
!				ypos = real(proc_pos(idproc+1) + j - 2)
!				if (ypos >= ylow .AND. ypos <= yhigh) then
!					if (ypos > ny/2) then
!						ytemp = 1.
!						! phase should be tan(angle) * (y - y_0)
!						phase = -1. * tan_value * (ypos - yhalf)
!					else
!						ytemp = 1.
!						phase = tan_value * (ypos - yhalf)
!					endif
!				else
!					ytemp = 0.
!				endif
!		
!				tempamp = tempamp * ytemp
!		
!				! Array position 2 corresponds to x=0 cause of gaurdcells
!				stoppos = floor(r_f_f)+3
!				do i = 2, stoppos-2
!					xpos = real(i)-2.
!					if (xpos < rise) then
!						xfac = xpos * riseinv
!					else if (xpos < r_f) then
!						xfac = 1.
!					else if (xpos < r_f_f) then
!						xfac = mfallinv*xpos + falloffset
!					endif
!					
!					tempamp2 = tempamp*xfac  
!					fxye(1,i,j,1) = fxye(1,i,j,1) + tempamp2*sin(wavek*xpos-wavew*t + phase)
!				enddo
!			enddo
!		endif
	end subroutine cross_wave
	
	subroutine angle_wave(fxye,t,nx,nxe,ny,nypmx,&
		&nvp,idproc)
		integer :: nx,nxe,ny,nypmx,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i,j,stoppos, nyproc,ystart,yend
		integer, dimension(nvp) :: proc_pos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac
		real :: ylow, yhigh, ypos
		real :: tempamp, tempamp2,ytemp
		real :: mfallinv,falloffset,r_f,r_f_f, phase,tan_value
		
!		wavek = real(nx)/wavemode
!		wavek = 6.283185307/wavek
!		r_f = rise + flat
!		r_f_f = rise + flat + fall
!		riseinv = 1. / rise
!		mfallinv = -1. / fall
!		falloffset = rise/fall + flat/fall + 1.
!		tfac = 1.
!		
!		nyproc = ny / nvp
!		do i = 0, nvp-1
!			proc_pos(i+1) = ny / nvp * (i)
!		enddo
!		
!		tfac = 1.
!		if (t < timerise) then
!			tfac = t / timerise
!		else if (t < timerise + timeflat) then
!			tfac = 1.
!		else if (t < timerise + timeflat + timefall) then
!			tfac = 1. - (t - (timerise+timeflat))/timefall
!		endif
!		
!		if (t < timerise + timeflat + timefall) then
!		
!			ylow = real(ny/2-yrise_fall)
!			yhigh = real(ny/2+yrise_fall)
!		!				tan_value = tan(1.5707963 - angle)
!			tan_value = tan(angle)
!			do j = 1,nypmx
!		
!				tempamp = amp * tfac
!		
!				ypos = real(proc_pos(idproc+1) + j - 2)
!				if (ypos >= ylow .AND. ypos <= yhigh) then
!					ytemp = 1.
!					! phase should be tan(angle) * (y - y_0)
!					phase = tan_value * (ypos - ylow)
!				else
!					ytemp = 0.
!				endif
!		
!				tempamp = tempamp * ytemp
!		
!				! Array position 2 corresponds to x=0 cause of gaurdcells
!				stoppos = floor(r_f_f)+3
!				do i = 2, stoppos-2
!					xpos = real(i)-2.
!					if (xpos < rise) then
!						xfac = xpos * riseinv
!					else if (xpos < r_f) then
!						xfac = 1.
!					else if (xpos < r_f_f) then
!						xfac = mfallinv*xpos + falloffset
!					endif
!					
!					tempamp2 = tempamp*xfac  
!					fxye(1,i,j,1) = fxye(1,i,j,1) + tempamp2*sin(wavek*xpos-wavew*t - phase)
!				enddo
!			enddo
!		endif
	end subroutine angle_wave
	
	subroutine plane_wave(fxyze,t,nx,nxe,nvp,idproc)
		integer :: nx,nxe,nvp,idproc
		real :: t
		real,dimension(:,:,:,:,:) :: fxyze
		integer :: i
		real :: xstart,xpos,wavek,tempamp,fac, tfac
		
		wavek = real(nx)/wavemode
		wavek = 6.283185307/wavek
		fac = 1.
		
		if (t < timerise) then
			tfac = t / timerise
		else if (t < timerise + timeflat) then
			tfac = 1.
		else if (t < timerise + timeflat + timefall) then
			tfac = 1. - (t - (timerise+timeflat))/timefall
		endif
		
		tfac = tfac * amp
		
		if (t < timerise + timeflat + timefall) then
		
			! Array position 2 corresponds to x=0 cause of gaurdcells
			do i = 1, nxe
				xpos = real(i)-2.
				fxyze(1,i,:,:,1) = fxyze(1,i,:,:,1) + tfac * sin(wavek*xpos-wavew*t)
			enddo
		endif
		
	end subroutine plane_wave
	
	subroutine plane_finite_wavelen(fxye,t,nx,nxe,nvp,idproc)
		integer :: nx,nxe,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i,stoppos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac
		real :: tempamp, tempamp2
		real :: mfallinv,falloffset,r_f,r_f_f
		
!		wavek = real(nx)/wavemode
!		wavek = 6.283185307/wavek
!		r_f = rise + flat
!		r_f_f = rise + flat + fall
!		riseinv = 1. / rise
!		mfallinv = -1. / fall
!		falloffset = rise/fall + flat/fall + 1.
!		tfac = 1.
!		
!		tfac = 1.
!		if (t < timerise) then
!			tfac = t / timerise
!		else if (t < timerise + timeflat) then
!			tfac = 1.
!		else if (t < timerise + timeflat + timefall) then
!			tfac = 1. - (t - (timerise+timeflat))/timefall
!		endif
!		
!		if (t < timerise + timeflat + timefall) then
!			tempamp = amp * tfac
!		
!			! Array position 2 corresponds to x=0 cause of gaurdcells      	
!			stoppos = floor(r_f_f)+3
!			do i = 2, stoppos-2
!				xpos = real(i)-2.
!				if (xpos < rise) then
!					xfac = xpos * riseinv
!				else if (xpos < r_f) then
!					xfac = 1.
!				else if (xpos < r_f_f) then
!					xfac = mfallinv*xpos + falloffset
!				endif
!				
!				tempamp2 = tempamp*xfac      		
!				fxye(1,i,:,1) = fxye(1,i,:,1) + &
!					& tempamp2*sin(wavek*xpos-wavew*t)
!			enddo
!		endif
	end subroutine plane_finite_wavelen
	
	subroutine gauss_tran_finite_wavelen(fxye,t,nx,nxe,ny,nypmx,&
		&nvp,idproc)
		integer :: nx,nxe,ny,nypmx,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i,j,stoppos, nyproc,ystart,yend
		integer, dimension(nvp) :: proc_pos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac
		real :: ylow, yhigh, ypos
		real :: tempamp, tempamp2,ytemp
		real :: mfallinv,falloffset,r_f,r_f_f
		
!		wavek = real(nx)/wavemode
!		wavek = 6.283185307/wavek
!		r_f = rise + flat
!		r_f_f = rise + flat + fall
!		riseinv = 1. / rise
!		mfallinv = -1. / fall
!		falloffset = rise/fall + flat/fall + 1.
!		tfac = 1.
!		
!		nyproc = ny / nvp
!		do i = 0, nvp-1
!			proc_pos(i+1) = ny / nvp * (i)
!		enddo
!		
!		tfac = 1.
!		if (t < timerise) then
!			tfac = t / timerise
!		else if (t < timerise + timeflat) then
!			tfac = 1.
!		else if (t < timerise + timeflat + timefall) then
!			tfac = 1. - (t - (timerise+timeflat))/timefall
!		endif
!		
!		if (t < timerise + timeflat + timefall) then
!		
!			ylow = real(ny/2-yrise_fall)
!			yhigh = real(ny/2+yrise_fall)
!			do j = 1,nypmx
!		
!				tempamp = amp * tfac
!		
!				ypos = real(proc_pos(idproc+1) + j - 2)
!				if (ypos >= ylow .AND. ypos <= yhigh) then
!					if (ypos <= ny / 2) then
!						ytemp = ( ypos - ylow ) / yrise_fall
!						ytemp = 10.*(ytemp**3) - 15.*(ytemp**4) + 6.*(ytemp**5)
!					else 
!						ytemp = ( ypos - real(ny/2) ) / yrise_fall
!						ytemp = 1. - (10.*(ytemp**3) - 15.*(ytemp**4) + 6.*(ytemp**5))
!					endif
!				else
!					ytemp = 0.
!				endif
!		
!				tempamp = tempamp * ytemp
!		
!				! Array position 2 corresponds to x=0 cause of gaurdcells
!				stoppos = floor(r_f_f)+3
!				do i = 2, stoppos-2
!					xpos = real(i)-2.
!					if (xpos < rise) then
!						xfac = xpos * riseinv
!					else if (xpos < r_f) then
!						xfac = 1.
!					else if (xpos < r_f_f) then
!						xfac = mfallinv*xpos + falloffset
!					endif
!					
!					tempamp2 = tempamp*xfac  
!		
!					fxye(1,i,j,1) = fxye(1,i,j,1) + &
!						& tempamp2*sin(wavek*xpos-wavew*t)
!				enddo
!			enddo
!		endif
	end subroutine gauss_tran_finite_wavelen
	
	subroutine doub_gauss_tran_finite_wavelen(fxye,t,nx,nxe,ny,nypmx,nvp,idproc)
		integer :: nx,nxe,ny,nypmx,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i,j,stoppos, nyproc,ystart,yend
		integer, dimension(nvp) :: proc_pos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac
		real :: ylow, yhigh, ypos
		real :: tempamp, tempamp2,ytemp
		real :: mfallinv,falloffset,r_f,r_f_f
		
!		wavek = real(nx)/wavemode
!		wavek = 6.283185307/wavek
!		r_f = rise + flat
!		r_f_f = rise + flat + fall
!		riseinv = 1. / rise
!		mfallinv = -1. / fall
!		falloffset = rise/fall + flat/fall + 1.
!		tfac = 1.
!		
!		nyproc = ny / nvp
!		do i = 0, nvp-1
!			proc_pos(i+1) = ny / nvp * (i)
!		enddo
!		
!		tfac = 1.
!		if (t < timerise) then
!			tfac = t / timerise
!		else if (t < timerise + timeflat) then
!			tfac = 1.
!		else if (t < timerise + timeflat + timefall) then
!			tfac = 1. - (t - (timerise+timeflat))/timefall
!		endif
!		
!		if (t < timerise + timeflat + timefall) then
!		
!			!Do driver centered around center1 first
!			!both driver one and two have the same yrise_fall
!			ylow = real(center1-yrise_fall)
!			yhigh = real(center1+yrise_fall)
!			do j = 1,nypmx
!		
!				tempamp = amp * tfac
!		
!				ypos = real(proc_pos(idproc+1) + j - 2)
!				if (ypos >= ylow .AND. ypos <= yhigh) then
!					if (ypos <= center1) then
!						ytemp = ( ypos - ylow ) / yrise_fall
!						ytemp = 10.*(ytemp**3) - 15.*(ytemp**4) + 6.*(ytemp**5)
!					else 
!						ytemp = ( ypos - real(center1) ) / yrise_fall
!						ytemp = 1. - (10.*(ytemp**3) - 15.*(ytemp**4) + 6.*(ytemp**5))
!					endif
!				else
!					ytemp = 0.
!				endif
!		
!				tempamp = tempamp * ytemp
!		
!				! Array position 2 corresponds to x=0 cause of gaurdcells
!				stoppos = floor(r_f_f)+3
!				do i = 2, stoppos-2
!					xpos = real(i)-2.
!					if (xpos < rise) then
!						xfac = xpos * riseinv
!					else if (xpos < r_f) then
!						xfac = 1.
!					else if (xpos < r_f_f) then
!						xfac = mfallinv*xpos + falloffset
!					endif
!					
!					tempamp2 = tempamp*xfac  
!		
!					fxye(1,i,j,1) = fxye(1,i,j,1) + &
!						& tempamp2*sin(wavek*xpos-wavew*t)
!				enddo
!			enddo
!		
!			!Now do driver centered around center2
!			!both driver one and two have the same yrise_fall
!			ylow = real(center2-yrise_fall)
!			yhigh = real(center2+yrise_fall)
!			do j = 1,nypmx
!		
!				tempamp = amp * tfac
!		
!				ypos = real(proc_pos(idproc+1) + j - 2)
!				if (ypos >= ylow .AND. ypos <= yhigh) then
!					if (ypos <= center2) then
!						ytemp = ( ypos - ylow ) / yrise_fall
!						ytemp = 10.*(ytemp**3) - 15.*(ytemp**4) + 6.*(ytemp**5)
!					else 
!						ytemp = ( ypos - real(center2) ) / yrise_fall
!						ytemp = 1. - (10.*(ytemp**3) - 15.*(ytemp**4) + 6.*(ytemp**5))
!					endif
!				else
!					ytemp = 0.
!				endif
!		
!				tempamp = tempamp * ytemp
!		
!				! Array position 2 corresponds to x=0 cause of gaurdcells
!				stoppos = floor(r_f_f)+3
!				do i = 2, stoppos-2
!					xpos = real(i)-2.
!					if (xpos < rise) then
!						xfac = xpos * riseinv
!					else if (xpos < r_f) then
!						xfac = 1.
!					else if (xpos < r_f_f) then
!						xfac = mfallinv*xpos + falloffset
!					endif
!					
!					tempamp2 = tempamp*xfac  
!		
!					fxye(1,i,j,1) = fxye(1,i,j,1) + &
!						& tempamp2*sin(wavek*xpos-wavew*t)
!				enddo
!			enddo
!		
!		endif
	end subroutine doub_gauss_tran_finite_wavelen
	
	!This makes a transverse circular wavefront with yrise_fall in both y and z
	subroutine gauss_tran_per_wavelen(fxyze,t,nx,nxe,ny,nypmx,nz,nzpmx,nvpy,nvpz,idproc)
		integer :: nx,nxe,ny,nypmx,nz,nzpmx,nvpy,nvpz,idproc
		real :: t
		real,dimension(:,:,:,:,:) :: fxyze
		integer :: i,j,k,stoppos, nyproc,ystart,yend, yproc, zproc
		integer, dimension(:,:,:),allocatable,save :: proc_pos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac, r, rmax
		real :: ny2myr,nz2myr,ny2pyr,nz2pyr
		real :: ylow, yhigh, ypos,zpos,ny2,nz2, yrise_fall_inv
		real :: tempamp, tempamp2,ytemp
		real :: mfallinv,falloffset,r_f,r_f_f
		
		wavek = real(nx)/wavemode
		wavek = 6.283185307/wavek
		r_f = rise + flat
		r_f_f = rise + flat + fall
		riseinv = 1. / rise
		mfallinv = -1. / fall
		falloffset = rise/fall + flat/fall + 1.
		tfac = 1.
		
		if (.not. allocated(proc_pos)) then
			allocate(proc_pos(nvpy,nvpz,2))
			do i = 0, nvpy-1
				do j = 0, nvpz-1
					proc_pos(i+1,j+1,1) = ny / nvpy
					proc_pos(i+1,j+1,1) = proc_pos(i+1,j+1,1) * i
					proc_pos(i+1,j+1,2) = nz / nvpz
					proc_pos(i+1,j+1,2) = proc_pos(i+1,j+1,2) * j
				enddo
			enddo
		endif
		
		if (t < timerise) then
			tfac = t / timerise
		else if (t < timerise + timeflat) then
			tfac = 1.
		else if (t < timerise + timeflat + timefall) then
			tfac = 1. - (t - (timerise+timeflat))/timefall
		endif
		
		if (t < timerise + timeflat + timefall) then
			yproc = mod(idproc,nvpy) + 1
			zproc = idproc / nvpy + 1
			
			rmax = yrise_fall
			yrise_fall_inv = 1. / real(yrise_fall)
		
			ylow = real(ny/2-yrise_fall)
			yhigh = real(ny/2+yrise_fall)
			ny2 = real(ny/2)
			nz2 = real(nz/2)
			ny2myr = ny2 - yrise_fall
			nz2myr = nz2 - yrise_fall
			ny2pyr = ny2 + yrise_fall
			nz2pyr = nz2 + yrise_fall
			
			do k = 1, nzpmx
				do j = 1,nypmx
					ypos = real(proc_pos(yproc,zproc,1) + j)
					zpos = real(proc_pos(yproc,zproc,2) + k)
					!Avoid square root for the cases where we are not within the square defined
					!by the midpoint plus and minus yrise_fall
					if ( (ypos .ge. ny2myr) .AND. (ypos .le. ny2pyr) ) then
					if ( (zpos .ge. nz2myr) .AND. (zpos .le. nz2pyr) ) then
						r = sqrt( (ypos - ny2)**2 + (zpos - nz2)**2 )
				
						if (r .le. rmax) then
							tempamp = amp * tfac
							r = r * yrise_fall_inv
							ytemp = 1. - (10.*(r**3) - 15.*(r**4) + 6.*(r**5))
							tempamp = tempamp * ytemp
						else
							tempamp = 0.
							
						endif
				
						! Array position 2 corresponds to x=0 cause of gaurdcells
						do i = 2, nx
							xpos = real(i)-2.
							
							fxyze(1,i,j,k,1) = fxyze(1,i,j,k,1) + tempamp*sin(wavek*xpos-wavew*t)
						enddo
					endif
					endif
				enddo
			enddo
		endif
	end subroutine gauss_tran_per_wavelen

	!This makes a transverse circular wavefront with yrise_fall in both y and z
	subroutine supergauss_tran_per_wavelen(fxyze,t,nx,nxe,ny,nypmx,nz,nzpmx,nvpy,nvpz,idproc)
		integer :: nx,nxe,ny,nypmx,nz,nzpmx,nvpy,nvpz,idproc
		real :: t
		real,dimension(:,:,:,:,:) :: fxyze
		integer :: i,j,k,stoppos, nyproc,ystart,yend, yproc, zproc
		integer, dimension(:,:,:),allocatable,save :: proc_pos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac, r, rmax, rmax2
		real :: ny2myr,nz2myr,ny2pyr,nz2pyr
		real :: ylow, yhigh, ypos,zpos,ny2,nz2, yrise_fall_inv,yrise_fall_inv2
		real :: tempamp, tempamp2,ytemp
		real :: mfallinv,falloffset,r_f,r_f_f
		
		wavek = real(nx)/wavemode
		wavek = 6.283185307/wavek
		r_f = rise + flat
		r_f_f = rise + flat + fall
		riseinv = 1. / rise
		mfallinv = -1. / fall
		falloffset = rise/fall + flat/fall + 1.
		tfac = 1.
		
		if (.not. allocated(proc_pos)) then
			allocate(proc_pos(nvpy,nvpz,2))
			do i = 0, nvpy-1
				do j = 0, nvpz-1
					proc_pos(i+1,j+1,1) = ny / nvpy
					proc_pos(i+1,j+1,1) = proc_pos(i+1,j+1,1) * i
					proc_pos(i+1,j+1,2) = nz / nvpz
					proc_pos(i+1,j+1,2) = proc_pos(i+1,j+1,2) * j
				enddo
			enddo
		endif
		
		if (t < timerise) then
			tfac = t / timerise
		else if (t < timerise + timeflat) then
			tfac = 1.
		else if (t < timerise + timeflat + timefall) then
			tfac = 1. - (t - (timerise+timeflat))/timefall
		endif
		
		if (t < timerise + timeflat + timefall) then
			yproc = mod(idproc,nvpy) + 1
			zproc = idproc / nvpy + 1
			
			rmax = yrise_fall
			rmax2 = rmax*rmax * 2.5 * 2.5	!The supergaussian is 3E-9 smaller by r = 2.5*yrise_fall
			yrise_fall_inv = 1. / real(yrise_fall)
			yrise_fall_inv2 = yrise_fall_inv * yrise_fall_inv
		
			ylow = real(ny/2-yrise_fall)
			yhigh = real(ny/2+yrise_fall)
			ny2 = real(ny/2)
			nz2 = real(nz/2)
			ny2myr = ny2 - yrise_fall
			nz2myr = nz2 - yrise_fall
			ny2pyr = ny2 + yrise_fall
			nz2pyr = nz2 + yrise_fall
			
			do k = 1, nzpmx
				do j = 1,nypmx
					ypos = real(proc_pos(yproc,zproc,1) + j)
					zpos = real(proc_pos(yproc,zproc,2) + k)
					!Avoid square root for the cases where we are not within the square defined
					!by the midpoint plus and minus yrise_fall
					r = (ypos - ny2)**2 + (zpos - nz2)**2
!						r = sqrt( (ypos - ny2)**2 + (zpos - nz2)**2 )
			
					if (r .le. rmax2) then
!						print*,sqrt(r)
						tempamp = amp * tfac
						r = r * yrise_fall_inv2
						ytemp = exp( -0.5 * (r)**2 )	!r is already squared above, so only square again for 4th power
						tempamp = tempamp * ytemp
						! Array position 2 corresponds to x=0 cause of gaurdcells
						do i = 2, nx
							xpos = real(i)-2.
							fxyze(1,i,j,k,1) = fxyze(1,i,j,k,1) + tempamp*sin(wavek*xpos-wavew*t)
						enddo
					endif
			
				enddo
			enddo
		endif
	end subroutine supergauss_tran_per_wavelen

	
	subroutine gauss_tran_per_wavelen_circle_wavefronts(fxye,t,nx,nxe,ny,nypmx,nvp,idproc)
		integer :: nx,nxe,ny,nypmx,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i,j,stoppos, nyproc,ystart,yend,yarr
		integer, dimension(nvp) :: proc_pos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac
		real :: ylow, yhigh, ypos
		real :: tempamp, tempamp2,ytemp
		real :: mfallinv,falloffset,r_f,r_f_f
		real,allocatable,save,dimension(:) :: profile,phase_off
		real :: r, phi
						
!		wavek = real(nx)/wavemode
!		wavek = 6.283185307/wavek
!		r_f = rise + flat
!		r_f_f = rise + flat + fall
!		riseinv = 1. / rise
!		mfallinv = -1. / fall
!		falloffset = rise/fall + flat/fall + 1.
!		tfac = 1.
!		
!		nyproc = ny / nvp
!		do i = 0, nvp-1
!			proc_pos(i+1) = ny / nvp * (i)
!		enddo
!		
!		!Initialize phase_offset array
!		if (.not. allocated(phase_off)) then
!			allocate(phase_off(0:int(yrise_fall)))
!			r = ( phase_offset**2 + yrise_fall**2 ) / 2. / phase_offset
!		!				print*,"r",r
!			phase_off(0) = -1.*phase_offset
!			do j = 1, int(yrise_fall)
!				phi = asin( (real(j)) / r)
!				phase_off(j) = -((real(j)) / tan(phi)-(r-phase_offset))
!		!					print*,"phase_off",j,phase_off(j)-(r-phase_offset)
!			enddo
!		endif
!		
!		!Initialize transverse profile array
!		if (.not. allocated(profile)) then
!			allocate(profile(0:int(yrise_fall)))
!			do j = 0,int(yrise_fall)
!				ypos = real(j)
!				ytemp = (yrise_fall-ypos) / yrise_fall
!				profile(j) = 10.*(ytemp**3) - 15.*(ytemp**4) + 6.*(ytemp**5)
!			enddo
!		endif
!					
!		tfac = 1.
!		if (t < timerise) then
!			tfac = t / timerise
!		else if (t < timerise + timeflat) then
!			tfac = 1.
!		else if (t < timerise + timeflat + timefall) then
!			tfac = 1. - (t - (timerise+timeflat))/timefall
!		endif
!		
!		if (t < timerise + timeflat + timefall) then
!		
!			ylow = real(ny/2-yrise_fall)
!			yhigh = real(ny/2+yrise_fall)
!			do j = 1,nypmx
!				
!				ypos = real(proc_pos(idproc+1) + j - 2)
!				if (ypos >= ylow .AND. ypos <= yhigh) then
!					yarr = abs( ny/2 - int(ypos))
!		
!					tempamp = amp * tfac * profile(yarr)
!		
!				! Array position 2 corresponds to x=0 cause of gaurdcells
!					do i = 2, nx
!						xpos = real(i)-2.
!						
!						fxye(1,i,j,1) = fxye(1,i,j,1) + tempamp*sin(wavek*xpos-wavew*t + phase_off( yarr ) )
!					enddo
!				endif
!			enddo
!		endif
	end subroutine gauss_tran_per_wavelen_circle_wavefronts
	
	!This driver is like gauss_tran_per_wavelen but also adds a component in Ey to make wave like one
	! that would come from a potential with gaussian profile
	subroutine gauss_tran_per_wavelen_pois(fxye,t,nx,nxe,ny,nypmx,nvp,idproc)
		integer :: nx,nxe,ny,nypmx,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i,j,stoppos, nyproc,ystart,yend,yarr
		integer, dimension(nvp) :: proc_pos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac
		real :: ylow, yhigh, ypos
		real :: tempamp, tempamp2,ytemp, kW2inv, yfac
		real :: mfallinv,falloffset,r_f,r_f_f
		real,allocatable,save,dimension(:) :: profile
		real :: r, phi
						
!		wavek = real(nx)/wavemode
!		wavek = 6.283185307/wavek
!		r_f = rise + flat
!		r_f_f = rise + flat + fall
!		riseinv = 1. / rise
!		mfallinv = -1. / fall
!		falloffset = rise/fall + flat/fall + 1.
!		tfac = 1.
!		kW2inv = 1./wavek/yrise_fall/yrise_fall
!		
!		nyproc = ny / nvp
!		do i = 0, nvp-1
!			proc_pos(i+1) = ny / nvp * (i)
!		enddo
!		
!		!Initialize transverse profile array
!		if (.not. allocated(profile)) then
!			allocate(profile(0:int(yrise_fall)))
!			do j = 0,int(yrise_fall)
!				ypos = real(j)
!				ytemp = (yrise_fall-ypos) / yrise_fall
!				profile(j) = 10.*(ytemp**3) - 15.*(ytemp**4) + 6.*(ytemp**5)
!			enddo
!		endif
!					
!		tfac = 1.
!		if (t < timerise) then
!			tfac = t / timerise
!		else if (t < timerise + timeflat) then
!			tfac = 1.
!		else if (t < timerise + timeflat + timefall) then
!			tfac = 1. - (t - (timerise+timeflat))/timefall
!		endif
!		
!		if (t < timerise + timeflat + timefall) then
!		
!			ylow = real(ny/2-yrise_fall)
!			yhigh = real(ny/2+yrise_fall)
!			do j = 1,nypmx
!				
!				ypos = real(proc_pos(idproc+1) + j - 2)
!				if (ypos >= ylow .AND. ypos <= yhigh) then
!					yarr = abs( ny/2 - int(ypos))
!		
!					tempamp = amp * tfac * profile(yarr)
!					yfac = tempamp*(ypos-real(ny/2))*kW2inv
!		
!				! Array position 2 corresponds to x=0 cause of gaurdcells
!					do i = 2, nx
!						xpos = real(i)-2.
!						
!						fxye(1,i,j,1) = fxye(1,i,j,1) + tempamp*sin(wavek*xpos-wavew*t)
!						fxye(2,i,j,1) = fxye(2,i,j,1) + yfac*cos(wavek*xpos-wavew*t)
!					enddo
!				endif
!			enddo
!		endif
	end subroutine gauss_tran_per_wavelen_pois
	
	subroutine rect_tran_finite_wavelen(fxye,t,nx,nxe,ny,nypmx,nvp,idproc)
		integer :: nx,nxe,ny,nypmx,nvp,idproc
		real :: t
		real,dimension(:,:,:,:) :: fxye
		integer :: i,j,stoppos, nyproc,ystart,yend
		integer, dimension(nvp) :: proc_pos
		real :: xstart,xpos,wavek,riseinv,tfac,xfac
		real :: ylow, yhigh, ypos
		real :: tempamp, tempamp2,ytemp
		real :: mfallinv,falloffset,r_f,r_f_f
		
!		wavek = real(nx)/wavemode
!		wavek = 6.283185307/wavek
!		r_f = rise + flat
!		r_f_f = rise + flat + fall
!		riseinv = 1. / rise
!		mfallinv = -1. / fall
!		falloffset = rise/fall + flat/fall + 1.
!		tfac = 1.
!		
!		nyproc = ny / nvp
!		do i = 0, nvp-1
!			proc_pos(i+1) = ny / nvp * (i)
!		enddo
!		
!		tfac = 1.
!		if (t < timerise) then
!			tfac = t / timerise
!		else if (t < timerise + timeflat) then
!			tfac = 1.
!		else if (t < timerise + timeflat + timefall) then
!			tfac = 1. - (t - (timerise+timeflat))/timefall
!		endif
!		
!		if (t < timerise + timeflat + timefall) then
!			
!			ylow = real(ny/2-yrise_fall)
!			yhigh = real(ny/2+yrise_fall)
!			do j = 1,nypmx
!		
!				tempamp = amp * tfac
!		
!				ypos = real(proc_pos(idproc+1) + j - 2)
!				if (ypos >= ylow .AND. ypos <= yhigh) then
!					ytemp = 1.
!				else
!					ytemp = 0.
!				endif
!		
!				tempamp = tempamp * ytemp
!		
!				! Array position 2 corresponds to x=0 cause of gaurdcells
!				stoppos = floor(r_f_f)+3
!				do i = 2, stoppos-2
!					xpos = real(i)-2.
!					if (xpos < rise) then
!						xfac = xpos * riseinv
!					else if (xpos < r_f) then
!						xfac = 1.
!					else if (xpos < r_f_f) then
!						xfac = mfallinv*xpos + falloffset
!					endif
!					
!					tempamp2 = tempamp*xfac  
!		
!					fxye(1,i,j,1) = fxye(1,i,j,1) + &
!						& tempamp2*sin(wavek*xpos-wavew*t)
!				enddo
!			enddo
!		endif
	end subroutine rect_tran_finite_wavelen
	
	subroutine rect_tran_per_wavelen(fxye,t,nx,nxe,ny,nypmx,nvp,idproc)
	integer :: nx,nxe,ny,nypmx,nvp,idproc
	real :: t
	real,dimension(:,:,:,:) :: fxye
	integer :: i,j,stoppos, nyproc,ystart,yend
	integer, dimension(nvp) :: proc_pos
	real :: xstart,xpos,wavek,riseinv,tfac,xfac
	real :: ylow, yhigh, ypos
	real :: tempamp, tempamp2,ytemp
	real :: mfallinv,falloffset,r_f,r_f_f
	
!	wavek = real(nx)/wavemode
!	wavek = 6.283185307/wavek
!	r_f = rise + flat
!	r_f_f = rise + flat + fall
!	riseinv = 1. / rise
!	mfallinv = -1. / fall
!	falloffset = rise/fall + flat/fall + 1.
!	tfac = 1.
!	
!	nyproc = ny / nvp
!	do i = 0, nvp-1
!		proc_pos(i+1) = ny / nvp * (i)
!	enddo
!	
!	if (t < timerise) then
!		tfac = t / timerise
!	else if (t < timerise + timeflat) then
!		tfac = 1.
!	else if (t < timerise + timeflat + timefall) then
!		tfac = 1. - (t - (timerise+timeflat))/timefall
!	endif
!	
!	if (t < timerise + timeflat + timefall) then
!		
!		ylow = real(ny/2-yrise_fall)
!		yhigh = real(ny/2+yrise_fall)
!		do j = 1,nypmx
!	
!			tempamp = amp * tfac
!	
!			ypos = real(proc_pos(idproc+1) + j - 2)
!			if (ypos >= ylow .AND. ypos <= yhigh) then
!				ytemp = 1.
!			else
!				ytemp = 0.
!			endif
!	
!			tempamp = tempamp * ytemp
!	
!			! Array position 2 corresponds to x=0 cause of gaurdcells
!			stoppos = floor(r_f_f)+3
!			do i = 1, nxe
!				xpos = real(i)-2.						
!				fxye(1,i,j,1) = fxye(1,i,j,1) + tempamp*sin(wavek*xpos-wavew*t)
!			enddo
!		enddo
!	endif
	end subroutine rect_tran_per_wavelen

		
	subroutine laguerre_gaussian_ponderomotive_force(fxyze,t,nx,nxe,ny,nypmx,nz,nzpmx,nvpy,nvpz,idproc, l_number)
		integer :: nx,nxe,ny,nypmx,nz,nzpmx,nvpy,nvpz,idproc, l_number
		real :: t
		real,dimension(:,:,:,:,:) :: fxyze
		integer :: i,j,k, yproc, zproc
		integer, dimension(:,:,:),allocatable,save :: proc_pos
		real :: xstart,xpos,ypos, zpos, r2, spot_size2, z_R2, wavek,tempamp,fac, tfac

		wavek = real(nx)/wavemode
		wavek = 6.283185307/wavek

		if (.not. allocated(proc_pos)) then
			allocate(proc_pos(nvpy,nvpz,2))
			do i = 0, nvpy-1
				do j = 0, nvpz-1
					proc_pos(i+1,j+1,1) = ny / nvpy
					proc_pos(i+1,j+1,1) = proc_pos(i+1,j+1,1) * i
					proc_pos(i+1,j+1,2) = nz / nvpz
					proc_pos(i+1,j+1,2) = proc_pos(i+1,j+1,2) * j
				enddo
			enddo
		endif
		
		if (t < timerise) then
			tfac = t / timerise
		else if (t < timerise + timeflat) then
			tfac = 1.
		else if (t < timerise + timeflat + timefall) then
			tfac = 1. - (t - (timerise+timeflat))/timefall
		endif
		
		tfac = tfac * amp

		if (t < timerise + timeflat + timefall) then
		
			! Array position 2 corresponds to x=0 cause of guard cells
			do j = 1, nypmx
				do k = 1, nzpmx
					yproc = mod(idproc,nvpy) + 1
					zproc = idproc / nvpy + 1
					ypos = real(proc_pos(yproc,zproc,1) + j)
					zpos = real(proc_pos(yproc,zproc,2) + k)
					r2 = ypos ** 2. + zpos ** 2.
					spot_size2 = spot_size ** 2. 
                                        do i=1, nxe
                                                xpos = real(i) - 2.
						fxyze(1,i,j,k,1) = fxyze(1,i,j,k,1) + &
                                                &tfac * (2.*xpos/(xpos ** 2. + z_R2) *(amp &
                                                &/ (1. + xpos ** 2. / z_R2) * (2. * r2 / &
                                                &spot_size2) ** l_number * exp(-2.*r2/ &
                                                &spot_size2)))
						fxyze(2,i,j,k,1) = fxyze(2,i,j,k,1) + tfac * &
                                                &(-ypos *(2.*l_number/r2 - 4./spot_size2)*(amp &
                                                &/ (1. + xpos ** 2. / z_R2) * (2. * r2 / &
                                                &spot_size2) ** l_number * exp(-2.*r2/ &
                                                &spot_size2)))
						fxyze(3,i,j,k,1) = fxyze(3,i,j,k,1) + tfac * &
                                                &(-zpos * (2.*l_number/r2 - 4./spot_size2) *&
                                                &(amp / (1. + xpos ** 2. / z_R2) * (2. * r2 /&
                                                & spot_size2) ** l_number * exp(-2.*r2/ &
                                                &spot_size2)))
					enddo
				enddo
			enddo
		endif

	end subroutine laguerre_gaussian_ponderomotive_force

	subroutine two_laguerre_gaussian_beams(fxyze,t,nx,nxe,ny,nypmx,nz,nzpmx,nvpy,nvpz,idproc, l_number1, l_number2)
		integer :: nx,nxe,ny,nypmx,nz,nzpmx,nvpy,nvpz,idproc, l_number1, l_number2
		real :: t
		real,dimension(:,:,:,:,:) :: fxyze
		integer :: i,j,k, yproc, zproc
		integer, dimension(:,:,:),allocatable,save :: proc_pos
		real :: xstart,xpos,ypos, zpos, r2, z_R2, l_sum, l_diff, wavek,tempamp,fac, tfac, w_p, phi_l, spot_size2

		wavek = real(nx)/wavemode
		wavek = 6.283185307/wavek

		if (.not. allocated(proc_pos)) then
			allocate(proc_pos(nvpy,nvpz,2))
			do i = 0, nvpy-1
				do j = 0, nvpz-1
					proc_pos(i+1,j+1,1) = ny / nvpy
					proc_pos(i+1,j+1,1) = proc_pos(i+1,j+1,1) * i
					proc_pos(i+1,j+1,2) = nz / nvpz
					proc_pos(i+1,j+1,2) = proc_pos(i+1,j+1,2) * j
				enddo
			enddo
		endif
	

!               Exponential in time
!               Just to be extra careful, we have a cutoff
                if ((t < timerise + timeflat + timefall) .and. (t > 0)) then
                        tfac = exp(-((t-(timerise+0.5*timeflat))/timerise)**4.0)
                endif

!		if (t < timerise) then
!			tfac = t / timerise
!		else if (t < timerise + timeflat) then
!			tfac = 1.
!		else if (t < timerise + timeflat + timefall) then
!			tfac = 1. - (t - (timerise+timeflat))/timefall
!		endif
	
                !Here we use amp from the input deck and tfac depending on rise fall and flat.        
		tfac = tfac * amp
               
		if (t < timerise + timeflat + timefall) then
		        yproc = mod(idproc,nvpy) + 1
        	        zproc = idproc / nvpy + 1
		        do j = 1, nypmx 
			        do k = 1, nzpmx
				        ypos = real(proc_pos(yproc,zproc,1) + j - ny / 2)
				        zpos = real(proc_pos(yproc,zproc,2) + k - nz / 2)
                                        r2 = ypos ** 2. + zpos ** 2.
				        spot_size2 = spot_size ** 2. 
				        l_sum = real(l_number1 + l_number2)
				        l_diff = abs(real(l_number1 - l_number2))
                                        phi_l = l_diff * atan2(zpos, ypos)
                                        if (r2 /= 0.) then
                                                tempamp = tfac * (2. * r2 / spot_size2)&
                                                & ** (l_sum / 2.) * exp(-2.*r2 / spot_size2)
                                        else
                                                tempamp = tfac
                                        endif

                                        do i=1, nxe
                                                xpos = real(i) - 2.
                                                fxyze(1,i,j,k,1) = fxyze(1,i,j,k,1) - &
                                                &wavek * tempamp * sin(phi_l -wavek * &
                                                &xpos + wavew * t)
                                                if (r2 /= 0.) then
                                                       fxyze(2,i,j,k,1) = fxyze(2,i,j,k,1)&
                                                       & + tempamp * (-1 * zpos * l_diff / &
                                                       &r2 * sin(phi_l - wavek * xpos + wavew&
                                                       & * t) + ypos * (-l_sum / r2 + 4. / &
                                                       &spot_size2) * cos(phi_l-wavek * xpos &
                                                       &+ wavew * t))
				                       fxyze(3,i,j,k,1) = fxyze(3,i,j,k,1) +&
                                                       & tempamp * (ypos * l_diff / r2 * &
                                                       &sin(phi_l -wavek * xpos + wavew * t) +&
                                                       & zpos * (-l_sum / r2 + 4. / spot_size2)&
                                                       & * cos(phi_l - wavek * xpos + wavew * t))
                                               endif
                                         enddo
				enddo
			enddo
		endif

	end subroutine two_laguerre_gaussian_beams

end module ext_driver32_jf
