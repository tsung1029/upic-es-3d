!-----------------------------------------------------------------------
!
      module diag32_ie
      
! Has diagnostic subroutines written by Jay Fahlen for use with 
! Viktor Decyk's 2d BEPS code
! written Jan. 2008
			
      use pinit32d, only: idrun, indx, indy, indz, ntp, ntd, ntr, psolve, tend, dt,npx,npy,npz
      use p32d_jf
      use p0d
      use pinit32d_jf
	  use hdf5_write32_ie
      
      implicit none
			include 'mpif.h'
      private
      public :: idrun, indx, indy, ntp!, modesx, modesy
      public :: tend, dt
      public :: ntfield,nfieldxrec,nfieldyrec,nphvxxrec,nphvyxrec,nphvxyrec,nphvyyrec
      public :: write_jf_diag_file,mkdir_structure
      !public :: phase_slices_x,phase_slices_y,phase_space_vxy_vs_x_and_y
      !public :: phase_space_vxy_vs_x
      public :: phase_space_vxyz_vs_x
      public :: phase_space_vxyz_vs_y
      public :: phase_space_vxyz_vs_z
      !public :: write_Eoft_aty
!      public :: write_raw
		! name convention is xx mean vx vs x, xy means vx vs y
      public :: phsl_xx_nrec,phsl_yx_nrec,phsl_xx_funit,phsl_yx_funit
      public :: phsl_yy_nrec,phsl_xy_nrec,phsl_yy_funit,phsl_xy_funit
     	public :: phasevxxfunit,phasevyxfunit,phasevxyfunit,phasevyyfunit,nlinerec,linesfunit,nlines
      public :: ntdenfunit,ntpenefunit,ntjxfunit,ntjyfunit,ntdenrec,ntpenerec,ntjxrec,ntjyrec
      public :: j_diag_deposit,vdotE,ntfieldfunitx,ntfieldfunity,ntfieldfunitz
      public :: nvdotExrec,nvdotEyrec,nvdotExfunit,nvdotEyfunit,addtag
      public :: nESPoynt_xrec,ESPoynt_xfunit,nESPoynt_yrec,ESPoynt_yfunit
      public :: nESPoynt_sum_xrec,ESPoynt_sum_xfunit,nESPoynt_sum_yrec,ESPoynt_sum_yfunit
      public :: store_init_vel
      !public :: calc_ene_dist,calc_vx_dist,calc_vy_dist
      public :: kivx_loc,kivy_loc, vEx_loc,vEy_loc
      public :: vdotEx_partrec,vdotEy_partrec,vdotEx_partfunit,vdotEy_partfunit
      public :: vdotEx_sumrec,vdotEy_sumrec,vdotEx_sumfunit,vdotEy_sumfunit
      public :: fvx_xy_rec,fvy_xy_rec,fvx_xy_funit,fvy_xy_funit
      public :: set_part_array_vdotE,n_jE_u_poynt
      public :: vdotEx_follow_partrec,vdotEy_follow_partrec
      public :: vdotEx_follow_partfunit,vdotEy_follow_partfunit
      public :: ntfield_ene_int_xrec,ntfield_ene_int_yrec
      public :: ntfield_ene_int_funitx,ntfield_ene_int_funity
      public :: div_ESPoynt_int_rec,div_ESPoynt_int_funit
      public :: div_ESPoynt_rec,div_ESPoynt_funit
      public :: nt_phase_vx_vy_rec,nt_phase_vx_vy_funit
      !public :: phase_vx_vy
      !public :: phase_space_vxy_vs_y
      public :: get_ene_in_driver_region,get_Py_sumover_x,get_1D_sumover_x
!      public :: calc_ene_three_regions
      public :: get_2D_sumover_x
      public :: nt_write_div_P_vs_y_funit,nt_write_jE_sumover_x_funit,nt_write_Py_sumover_x_funit
      public :: nt_write_U_sumover_x_funit, bin_Ex_rec, bin_Ey_rec, bin_Ex_funit,bin_Ey_funit
      !public :: bin_Ex_and_Ey
      public :: ntpotfunit, npotrec
      public :: get_3D_sumover_x
      
      save
      integer :: nfieldxrec = 0, nfieldyrec = 0
      integer :: nphvxxrec = 0, nphvyxrec = 0,nphvxyrec=0,nphvyyrec=0,nlines
      integer,dimension(:),allocatable :: phsl_xx_nrec,phsl_yx_nrec
      integer,dimension(:),allocatable :: phsl_yy_nrec,phsl_xy_nrec
      integer,dimension(5) :: phsl_xx_funit,phsl_yx_funit
      integer,dimension(5) :: phsl_yy_funit,phsl_xy_funit
      integer,dimension(:),allocatable :: nlinerec
      integer,dimension(5) :: linesfunit
      integer :: phasevxxfunit,phasevyxfunit,phasevxyfunit,phasevyyfunit
      integer :: ntdenfunit,ntpenefunit,ntjxfunit,ntjyfunit,ntdenrec=0,ntpenerec=0
      integer :: ntjxrec=0,ntjyrec=0
      integer :: nvdotExrec=0,nvdotEyrec=0,nvdotExfunit,nvdotEyfunit
      integer :: addtag = 0, ntfieldfunitx,ntfieldfunity,ntfieldfunitz
      integer :: nESPoynt_xrec=0,ESPoynt_xfunit,nESPoynt_yrec=0,ESPoynt_yfunit
      integer :: nESPoynt_sum_xrec=0,ESPoynt_sum_xfunit,nESPoynt_sum_yrec=0,ESPoynt_sum_yfunit
      integer :: calc_ene_dist_rec=0,calc_ene_dist_funit,calc_vx_dist_rec=0,calc_vx_dist_funit
      integer :: calc_vy_dist_rec=0,calc_vy_dist_funit
      integer :: vdotEx_partrec=0,vdotEy_partrec=0,vdotEx_partfunit,vdotEy_partfunit
      integer :: vdotEx_sumrec=0,vdotEy_sumrec=0,vdotEx_sumfunit,vdotEy_sumfunit
      integer :: kivx_loc=0,kivy_loc=0, vEx_loc=0,vEy_loc=0
		integer :: fvx_xy_rec=0,fvy_xy_rec=0,fvx_xy_funit,fvy_xy_funit   
		integer :: set_part_array_vdotE = 0, n_jE_u_poynt = 0
      integer :: vdotEx_follow_partrec=0,vdotEy_follow_partrec=0
      integer :: vdotEx_follow_partfunit,vdotEy_follow_partfunit
      integer :: ntfield_ene_int_xrec=0,ntfield_ene_int_yrec=0
      integer :: ntfield_ene_int_funitx,ntfield_ene_int_funity
      integer :: div_ESPoynt_int_rec=0,div_ESPoynt_int_funit
      integer :: div_ESPoynt_rec=0,div_ESPoynt_funit
      integer :: nt_phase_vx_vy_rec=0,nt_phase_vx_vy_funit
      integer :: nt_write_div_P_vs_y_funit,nt_write_jE_sumover_x_funit,nt_write_Py_sumover_x_funit
      integer :: nt_write_U_sumover_x_funit
      integer :: bin_Ex_rec = 0, bin_Ey_rec = 0, bin_Ex_funit,bin_Ey_funit
      integer :: ntpotfunit, npotrec = 0


     
      contains
      	subroutine mkdir_structure(idproc)
      		implicit none
      		integer :: idproc
      		character(len=4) :: temp
      		
				integer :: ierr, i
      		
      		if (idproc == 0) then
						if (ntd .ne. 0) then
							call mkdir_f('./DIAG/Qe'//char(0),ierr)
                     call mkdir_f('./DIAG/Qi'//char(0),ierr)
						endif
						if (ntfield .ne. 0) then
							call mkdir_f('./DIAG/Ex'//char(0),ierr)
							call mkdir_f('./DIAG/Ey'//char(0),ierr)
							call mkdir_f('./DIAG/Ez'//char(0),ierr)
						endif
						if (ntESPoynt .ne. 0) then
							call mkdir_f('./DIAG/ESPoynt_x/'//char(0),ierr)
							call mkdir_f('./DIAG/ESPoynt_y/'//char(0),ierr)
						endif
						if (ntESPoynt_int .ne. 0) then
							call mkdir_f('./DIAG/ESPoynt_int_x/'//char(0),ierr)
							call mkdir_f('./DIAG/ESPoynt_int_y/'//char(0),ierr)
							call mkdir_f('./DIAG/init_cond_ESPoynt_int/'//char(0),ierr)							
						endif
						if (ntj .ne. 0) then
							call mkdir_f('./DIAG/jx'//char(0),ierr)
							call mkdir_f('./DIAG/jy'//char(0),ierr)
						endif
						if (nttrack .ne. 0) then
							call mkdir_f('./DIAG/TRACKS'//char(0),ierr)
						endif
						if (keep_init_vel .ne. 0) then
							call mkdir_f('./DIAG/EneBin'//char(0),ierr)
							call mkdir_f('./DIAG/VxBin'//char(0),ierr)
							call mkdir_f('./DIAG/VyBin'//char(0),ierr)
						endif
						if (nvdotE_int .ne. 0) then
							call mkdir_f('./DIAG/vdotEx_int'//char(0),ierr)
							call mkdir_f('./DIAG/vdotEy_int'//char(0),ierr)
							call mkdir_f('./DIAG/init_cond_vdotE_int'//char(0),ierr)
						endif
						if (nvdotE_part .ne. 0) then
							call mkdir_f('./DIAG/vdotEx_part'//char(0),ierr)
							call mkdir_f('./DIAG/vdotEy_part'//char(0),ierr)
						endif
						if (ntvdotE .ne. 0) then
							call mkdir_f('./DIAG/vdotEx'//char(0),ierr)
							call mkdir_f('./DIAG/vdotEy'//char(0),ierr)
						endif
						if (ntphxy .ne. 0) then
							call mkdir_f('./DIAG/fvx_xy'//char(0),ierr)
							call mkdir_f('./DIAG/fvy_xy'//char(0),ierr)
						endif
						if (nvdotE_follow_part .ne. 0) then
							call mkdir_f('./DIAG/vdotEx_follow_part'//char(0),ierr)
							call mkdir_f('./DIAG/vdotEy_follow_part'//char(0),ierr)
						endif
						if (ntfield_ene_int .ne. 0) then
							call mkdir_f('./DIAG/u_x_int'//char(0),ierr)
							call mkdir_f('./DIAG/u_y_int'//char(0),ierr)
							call mkdir_f('./DIAG/init_cond_u_int'//char(0),ierr)
						endif
						if (ntraw .ne. 0) then
							call mkdir_f('./DIAG/RAW'//char(0),ierr)
						endif
						if (nt_div_ESPoynt_int .ne. 0) then
							call mkdir_f('./DIAG/div_ESPoynt_int'//char(0),ierr)
						endif
						if (nt_div_ESPoynt .ne. 0) then
							call mkdir_f('./DIAG/div_ESPoynt'//char(0),ierr)
						endif
						if (nt_phase_vx_vy .ne. 0) then
							call mkdir_f('./DIAG/VxVy'//char(0),ierr)
						endif
						if (nt_bin_E .ne. 0) then
							call mkdir_f('./DIAG/bin_Ex'//char(0),ierr)
							call mkdir_f('./DIAG/bin_Ey'//char(0),ierr)
						endif

						if (nphxx .ne. 0) then
							call mkdir_f('./DIAG/Vx_x'//char(0),ierr)
						endif
						if (nphyx .ne. 0) then
							call mkdir_f('./DIAG/Vy_x'//char(0),ierr)
						endif
						if (nphzx .ne. 0) then
							call mkdir_f('./DIAG/Vz_x'//char(0),ierr)
						endif
						if (nphxy .ne. 0) then
							call mkdir_f('./DIAG/Vx_y'//char(0),ierr)
						endif
						if (nphyy .ne. 0) then
							call mkdir_f('./DIAG/Vy_y'//char(0),ierr)
						endif
						if (nphzy .ne. 0) then
							call mkdir_f('./DIAG/Vz_y'//char(0),ierr)
						endif
						if (nphxz .ne. 0) then
							call mkdir_f('./DIAG/Vx_z'//char(0),ierr)
						endif
						if (nphyz .ne. 0) then
							call mkdir_f('./DIAG/Vy_z'//char(0),ierr)
						endif
						if (nphzz .ne. 0) then
							call mkdir_f('./DIAG/Vz_z'//char(0),ierr)
						endif
						if (nt_write_Py_vs_y .ne. 0) then
							call mkdir_f('./DIAG/Py_sumover_x'//char(0),ierr)
						endif
						if (nt_write_div_P_vs_y .ne. 0) then
							call mkdir_f('./DIAG/div_P_sumover_x'//char(0),ierr)
						endif
						if (ntp .ne. 0) then
							call mkdir_f('./DIAG/Pot'//char(0),ierr)
						endif
						
						if (ntr .ne. 0) then
							call mkdir_f('./RST'//char(0),ierr)
						endif

						
						do i = 1, num_phsl_x
							write(temp,'(i4)') phsl_x_pos(i)
							temp = adjustl(temp)
							call mkdir_f('./DIAG/vx_x/'//trim(temp)//char(0),ierr)
							call mkdir_f('./DIAG/vy_x/'//trim(temp)//char(0),ierr)
						enddo
						do i = 1, num_phsl_y
							write(temp,'(i4)') phsl_y_pos(i)
							temp = adjustl(temp)
							call mkdir_f('./DIAG/vx_y/'//trim(temp)//char(0),ierr)
							call mkdir_f('./DIAG/vy_y/'//trim(temp)//char(0),ierr)
						enddo
				endif
      		
      	end subroutine mkdir_structure
      	
      	subroutine write_jf_diag_file(id0,cdrun)
      		implicit none
      		integer::id0
      		character(len=*) :: cdrun
      		character(len=32)::filename
      		integer :: funitnum
      		!integer::modesx=11,modesy=11

      		if (id0 == 0) then
      			filename='diagparams_jf'
      			funitnum = get_funit(20)
      			open(unit=funitnum,file=trim(filename),form='unformatted',status=&
     &'replace')
					write (funitnum) idrun,indx,indy,indz,ntfield,psolve,tend,dt,npx,npy,npz,&
						&fvxmax,fvymax,nphbx,nphby,nphxx,nphyx,ntlines,nlines,linepos,phsl_x_pos,&
						&phsl_x_thick,ntphsl_x,num_phsl_x,ntden,ntpene,ntj,ntvdotE
				endif
			end subroutine write_jf_diag_file
			
			subroutine get_1D_sumover_x(div_ESPoynt,div_P_sumover_x,nx,nxe,ny,nypmx,nvp,idproc,inorder)
				implicit none
				real,dimension(:,:,:) :: div_ESPoynt
				real,dimension(:) :: div_P_sumover_x
				real :: num_par_cell
				integer :: nx, nxe, ny, nypmx, nvp, idproc, inorder

				integer :: i, j, k, nyproc,nylow,nyhigh, cury
      	integer, dimension(nvp) :: proc_pos
      	
      	nyproc = ny / nvp
      	do i = 0, nvp-1
      		proc_pos(i+1) = ny / nvp * (i)
      	enddo
				do j = inorder, nyproc+inorder-1
					cury = proc_pos(idproc+1) + j
					div_P_sumover_x(cury-inorder+1) = sum(div_ESPoynt(inorder:inorder+nx,j,1))
				enddo
				call plsum(div_P_sumover_x)
			end subroutine get_1D_sumover_x

			subroutine get_2D_sumover_x(jE,jE_sumover_x,nx,nxe,ny,nypmx,nvp,idproc,inorder)
				implicit none
				real,dimension(:,:,:,:) :: jE
				real,dimension(:) :: jE_sumover_x
				real :: num_par_cell
				integer :: nx, nxe, ny, nypmx, nvp, idproc, inorder

				integer :: i, j, k, nyproc,nylow,nyhigh, cury
      	integer, dimension(nvp) :: proc_pos
      	
      	nyproc = ny / nvp
      	do i = 0, nvp-1
      		proc_pos(i+1) = ny / nvp * (i)
      	enddo
				do j = inorder, nyproc+inorder-1
					cury = proc_pos(idproc+1) + j
					jE_sumover_x(cury-inorder+1) = sum(jE(1,inorder:inorder+nx,j,1)) + &
																			&	 sum(jE(2,inorder:inorder+nx,j,1))
				enddo
				call plsum(jE_sumover_x)
			end subroutine get_2D_sumover_x

			subroutine get_3D_sumover_x(f,out,nx,nxe,ny,nypmx,nz,nzpmx,nvpy,&
				& nvpz,idproc,inorder)
				implicit none
				real,dimension(:,:,:,:,:) :: f
				real,dimension(:,:) :: out
				real :: num_par_cell
				integer :: nx, nxe, ny, nypmx, nz, nzpmx, nvpy, nvpz, idproc, inorder

				integer :: i, j, k, nyproc,nzproc,nylow,nyhigh, cury, curz, yoff, zoff
      	integer, dimension(:,:,:),allocatable,save :: proc_pos
      	
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
				
				nyproc = ny / nvpy
				nzproc = nz / nvpz
				yoff = proc_pos(mod(idproc,nvpy) + 1,idproc / nvpy + 1,1)
				zoff = proc_pos(mod(idproc,nvpy) + 1,idproc / nvpy + 1,2)

				do i = inorder, nyproc+inorder-1
					do j = inorder, nzproc+inorder-1
						cury = i+yoff
						curz = j+zoff

						out(cury-inorder+1,curz-inorder+1) = &
							& sum(f(1,inorder:inorder+nx,i,j,1)) + &
							&	sum(f(2,inorder:inorder+nx,i,j,1)) + &
							& sum(f(3,inorder:inorder+nx,i,j,1))
					enddo
				enddo
				call plsum(out)
			end subroutine get_3D_sumover_x

			
			subroutine get_Py_sumover_x(ESPoynt,Py_sumover_x,nx,nxe,ny,nypmx,nvp,idproc,inorder)
				implicit none
				real,dimension(:,:,:,:) :: ESPoynt
				real,dimension(:) :: Py_sumover_x
				real :: num_par_cell
				integer :: nx, nxe, ny, nypmx, nvp, idproc, inorder

				integer :: i, j, k, nyproc,nylow,nyhigh, cury
      	integer, dimension(nvp) :: proc_pos
      	
      	nyproc = ny / nvp
      	do i = 0, nvp-1
      		proc_pos(i+1) = ny / nvp * (i)
      	enddo
				do j = inorder, nyproc+inorder-1
					cury = proc_pos(idproc+1) + j
					Py_sumover_x(cury-inorder+1) = sum(ESPoynt(2,inorder:inorder+nx,j,1))
				enddo
				call plsum(Py_sumover_x)
			end subroutine get_Py_sumover_x
			
			!This returns the field energy in the rectagle -yrise_fall < y - ny/2 < yrise_fall and all x
			subroutine get_ene_in_driver_region(efield,field_ene,nx,nxe,ny,nypmx,nvp,idproc,inorder)
				implicit none
				integer :: nx,nxe,ny,nypmx,nvp,idproc,inorder
				real :: field_ene
				real,dimension(:,:,:,:) :: efield
				
				integer :: i, j, nyproc,nylow,nyhigh
      	integer, dimension(nvp) :: proc_pos
      	real,dimension(1) :: temp
				
      	nyproc = ny / nvp
      	do i = 0, nvp-1
      		proc_pos(i+1) = ny / nvp * (i)
      	enddo
      	field_ene = 0.
      	
!      	print*,idproc,"a"
      	do j = inorder, nyproc+inorder
      		if ( (proc_pos(idproc+1) + j-1 >= ny/2-yrise_fall) .AND. &
      			&  (proc_pos(idproc+1) + j-1 <= ny/2+yrise_fall) ) then
!						if (idproc ==0) then
!						print*,j
!						endif
						do i = inorder, inorder+nx
							field_ene = field_ene + 0.5*(efield(1,i,j,1)**2 + efield(2,i,j,1)**2)
						enddo
					endif      		
      	enddo
!      	print*,idproc,"b"
      	
      	temp(1) = field_ene
      	call plsum(temp)
      	field_ene = temp(1)

			end subroutine get_ene_in_driver_region
			
			subroutine store_init_vel(part,npp)
				implicit none
				real,dimension(:,:,:) :: part
				integer,dimension(:),pointer :: npp
				
				integer :: i
				
				do i = 1, npp(1)
					part(kivx_loc,i,1) = part(3,i,1)
					part(kivy_loc,i,1) = part(4,i,1)
				enddo
			end subroutine store_init_vel
			
			! This will calculate the current from the particle data 
			subroutine j_diag_deposit(part,arr,npp,xory,ny,nvp,idproc)
				implicit none
				real,dimension(:,:,:), pointer :: part, arr
				integer,dimension(:) :: npp
				integer :: xory,ny,nvp,idproc
				
				integer :: i,xpos,ypos,dir,procpos
				
				procpos = ny /nvp * idproc
				procpos = procpos - 2 !to get correct array pos with guard cells

				dir = 2+xory
				do i = 1, npp(1)
					xpos = floor(part(1,i,1)) + 2
					ypos = floor(part(2,i,1)) - procpos
					arr(xpos,ypos,1) = arr(xpos,ypos,1) + part(dir,i,1)
				enddo
			end subroutine j_diag_deposit
			
			!This will calculate v_i E_i, where i is specified by xory
			subroutine vdotE(part,fxye,arr,npp,xory,ny,nvp,idproc)
				implicit none
				real,dimension(:,:,:), pointer :: part,arr
				real,dimension(:,:,:,:) :: fxye
				integer,dimension(:) :: npp
				integer :: xory,ny,nvp,idproc
				
				integer :: i,xpos,ypos,dir,procpos

				procpos = ny /nvp * idproc
				procpos = procpos - 2 !to get correct array pos with guard cells

				dir = 2+xory
				do i = 1, npp(1)
					xpos = floor(part(1,i,1)) + 2
					ypos = floor(part(2,i,1)) - procpos
					arr(xpos,ypos,1) = arr(xpos,ypos,1) + part(dir,i,1) * fxye(xory,xpos,ypos,1)
				enddo
			end subroutine vdotE

			! This calculates the phase space at a given time for 1 direction
			! specified by xoryorz(x=1,y=2,z=3),nxyz is nx or ny or nz depending on xoryorz
			! fvxyz = phase space for x or y or z
			subroutine phase_space_vxyz_vs_x(part,fvxyz,xoryorz,nxyz,npp,dt,it,nblok,idproc)
				implicit none
				integer :: nblok,xoryorz,nxyz,idproc
	      	real,dimension(:,:,:) :: part
				integer :: it
				real :: dt
				real,dimension(:,:,:),pointer :: fvxyz
				real,dimension(:,:,:),pointer :: sfieldtemp
				integer,dimension(nblok) :: npp
				integer :: i,varrpos,xarrpos,pos
				real :: dv,low,high,v,vrange,dvinv
				character(len=32) :: fname
				integer,dimension(2) :: int_tag
				real :: real_tag
				equivalence (real_tag,int_tag)
				
				if (xoryorz == 1) then
					vrange = 2.*fvxmax
					dv = vrange / real(nphbx)
					dvinv = 1. / dv
					low = -1.*fvxmax
					high = fvxmax
	
!					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
						
							pos = floor(part(1,i,1))+1
							v = part(4,i,1)
							varrpos = floor( (v + fvxmax)*dvinv + 0.5) + 1
							
							if (varrpos .le. nphbx) then 
								fvxyz(pos,varrpos,1) = fvxyz(pos,varrpos,1) + 1.
							endif
						enddo
!					else
!						do i = 1, npp(1)
!							real_tag = part(addtag_loc,i,1)
!							if (int_tag(1) == -1) then
!								pos = floor(part(1,i,1))+1
!								v = part(3,i,1)
!								varrpos = floor( (v + fvxmax)*dvinv + 0.5) + 1
!								
!								if (varrpos .le. nphbx) then 
!									fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
!								endif
!							endif
!						enddo
!					endif
						
					call plsum(fvxyz(:,:,1))
!					if (idproc == 0) then
					call writef(fvxyz,nxyz,nphbx,it,dt,&
								&PH_VX_X,1.,dv,0.,low)		
!					endif
!					if (idproc == 0) then
!					
!						if (nphvxxrec == 0) then
!							nphvxxrec = -1
!							fname = './DIAG/Vx_x/vx_x'
!							phasevxxfunit = get_funit(20)
!							call write_jf(fvxy,nxy,nphbx,phasevxxfunit,nphvxxrec,&
!								&trim(fname),PH_VX_X,time,it,1.,dv,0.,low)					
!		!					call writef(fvxvy,nphbx,nphby,nt_phase_vx_vy_funit,-1,itime,itime*dt,VXVY,trim(fname),inorder)
!						else 
!							fname = './DIAG/Vx_x/vx_x'
!							call write_jf(fvxy,nxy,nphbx,phasevxxfunit,nphvxxrec,&
!								&trim(fname),PH_VX_X,time,it,1.,dv,0.,low)					
!						endif
!					endif
				else if (xoryorz == 2) then
					vrange = 2.*fvymax
					dv = vrange / real(nphby)
					dvinv = 1. / dv
					low = -1.*fvymax
					high = fvymax
	
!					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
							pos = floor(part(1,i,1))+1
							v = part(5,i,1)
							varrpos = floor( (v + fvymax)*dvinv + 0.5) + 1
							if (varrpos .le. nphby) then 
								fvxyz(pos,varrpos,1) = fvxyz(pos,varrpos,1) + 1.
							endif
						enddo
!					else
!						do i = 1, npp(1)
!							real_tag = part(addtag_loc,i,1)
!							if (int_tag(1) == -1) then
!								pos = floor(part(1,i,1))+1
!								v = part(4,i,1)
!								varrpos = floor( (v + fvymax)*dvinv + 0.5) + 1
!								
!								if (varrpos .le. nphby) then 
!									fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
!								endif
!							endif
!						enddo
!					endif
					call plsum(fvxyz(:,:,1))
!					if (idproc == 0) then
!						call phasespace_writef(fvxy,nxy,nphby,it,dt,&
!								&PH_VY_X,1.,dv,0.,low)					
!					endif
					call writef(fvxyz,nxyz,nphby,it,dt,&
								&PH_VY_X,1.,dv,0.,low)		
!					if (idproc == 0) then
!
!						if (nphvyxrec == 0) then
!							nphvyxrec = -1
!							fname = './DIAG/Vy_x/vy_x'
!							phasevyxfunit = get_funit(20)
!							call write_jf(fvxy,nxy,nphby,phasevyxfunit,nphvyxrec,&
!								&trim(fname),PH_VY_X,time,it,1.,dv,0.,low)					
!						else 
!							fname = './DIAG/Vy_x/vy_x'
!							call write_jf(fvxy,nxy,nphby,phasevxyfunit,nphvxyrec,&
!								&trim(fname),PH_VY_X,time,it,1.,dv,0.,low)					
!						endif
!					endif
				else if (xoryorz == 3) then
					vrange = 2.*fvzmax
					dv = vrange / real(nphbz)
					dvinv = 1. / dv
					low = -1.*fvzmax
					high = fvzmax
	
!					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
						
							pos = floor(part(1,i,1))+1
							v = part(6,i,1)
							varrpos = floor( (v + fvzmax)*dvinv + 0.5) + 1
							
							if (varrpos .le. nphbz) then 
								fvxyz(pos,varrpos,1) = fvxyz(pos,varrpos,1) + 1.
							endif
						enddo
!					else
!						do i = 1, npp(1)
!							real_tag = part(addtag_loc,i,1)
!							if (int_tag(1) == -1) then
!								pos = floor(part(1,i,1))+1
!								v = part(3,i,1)
!								varrpos = floor( (v + fvxmax)*dvinv + 0.5) + 1
!								
!								if (varrpos .le. nphbx) then 
!									fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
!								endif
!							endif
!						enddo
!					endif
						
					call plsum(fvxyz(:,:,1))
!					if (idproc == 0) then
					call writef(fvxyz,nxyz,nphbz,it,dt,&
								&PH_VZ_X,1.,dv,0.,low)		
!					endif
!					if (idproc == 0) then
!					
!						if (nphvxxrec == 0) then
!							nphvxxrec = -1
!							fname = './DIAG/Vx_x/vx_x'
!							phasevxxfunit = get_funit(20)
!							call write_jf(fvxy,nxy,nphbx,phasevxxfunit,nphvxxrec,&
!								&trim(fname),PH_VX_X,time,it,1.,dv,0.,low)					
!		!					call writef(fvxvy,nphbx,nphby,nt_phase_vx_vy_funit,-1,itime,itime*dt,VXVY,trim(fname),inorder)
!						else 
!							fname = './DIAG/Vx_x/vx_x'
!							call write_jf(fvxy,nxy,nphbx,phasevxxfunit,nphvxxrec,&
!								&trim(fname),PH_VX_X,time,it,1.,dv,0.,low)					
!						endif
!					endif

				endif

			end subroutine phase_space_vxyz_vs_x

			! This calculates the phase space at a given time for 1 direction
			! specified by xoryorz(x=1,y=2,z=3),nxyz is nx or ny or nz depending on xoryorz
			! fvxyz = phase space for x or y or z
			subroutine phase_space_vxyz_vs_y(part,fvxyz,xoryorz,nxyz,npp,dt,it,nblok,idproc)
				implicit none
				integer :: nblok,xoryorz,nxyz,idproc
	      	real,dimension(:,:,:) :: part
				integer :: it
				real :: dt
				real,dimension(:,:,:),pointer :: fvxyz
				real,dimension(:,:,:),pointer :: sfieldtemp
				integer,dimension(nblok) :: npp
				integer :: i,varrpos,xarrpos,pos
				real :: dv,low,high,v,vrange,dvinv
				character(len=32) :: fname

				integer,dimension(2) :: int_tag
				real :: real_tag
				equivalence (real_tag,int_tag)
								
				if (xoryorz == 1) then
					vrange = 2.*fvxmax
					dv = vrange / real(nphbx)
					dvinv = 1. / dv
					low = -1.*fvxmax
					high = fvxmax
					
!					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
							pos = floor(part(2,i,1))+1
							v = part(4,i,1)
							varrpos = floor( (v + fvxmax)*dvinv + 0.5) + 1
							
							if ((varrpos .le. nphbx) .and. (varrpos .gt. 0)) then 
								fvxyz(pos,varrpos,1) = fvxyz(pos,varrpos,1) + 1.
							endif
						enddo
!					else
!						do i = 1, npp(1)
!							real_tag = part(addtag_loc,i,1)
!							if (int_tag(1) == -1) then
!								pos = floor(part(2,i,1))+1
!								v = part(3,i,1)
!								varrpos = floor( (v + fvxmax)*dvinv + 0.5) + 1
!								
!								if ((varrpos .le. nphbx) .and. (varrpos .gt. 0)) then 
!									fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
!								endif
!							endif
!						enddo
!					endif
						
					call plsum(fvxyz(:,:,1))
!					if (idproc == 0) then
!						call phasespace_writef(fvxy,nxy,nphby,it,dt,&
!								&PH_VX_Y,1.,dv,0.,low)					
!					endif
					call writef(fvxyz,nxyz,nphbx,it,dt,&
								&PH_VX_Y,1.,dv,0.,low)		
!					if (idproc == 0) then
!					
!						if (nphvxyrec == 0) then
!							nphvxyrec = -1
!							fname = './DIAG/Vx_y/vx_y'
!							phasevxyfunit = get_funit(20)
!							call write_jf(fvxy,nxy,nphbx,phasevxyfunit,nphvxyrec,&
!								&trim(fname),PH_VX_Y,time,it,1.,dv,0.,low)					
!						else 
!							fname = './DIAG/Vx_y/vx_y'
!							call write_jf(fvxy,nxy,nphbx,phasevxyfunit,nphvxyrec,&
!								&trim(fname),PH_VX_Y,time,it,1.,dv,0.,low)					
!						endif
!					endif
				else if (xoryorz == 2) then
					vrange = 2.*fvymax
					dv = vrange / real(nphby)
					dvinv = 1. / dv
					low = -1.*fvymax
					high = fvymax
	
!					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
							pos = floor(part(2,i,1))+1
							v = part(5,i,1)
							varrpos = floor( (v + fvymax)*dvinv + 0.5) + 1
							
							if (varrpos .le. nphby) then 
								fvxyz(pos,varrpos,1) = fvxyz(pos,varrpos,1) + 1.
							endif
						enddo
!					else
!						do i = 1, npp(1)
!							real_tag = part(addtag_loc,i,1)
!							if (int_tag(1) == -1) then
!								pos = floor(part(2,i,1))+1
!								v = part(4,i,1)
!								varrpos = floor( (v + fvymax)*dvinv + 0.5) + 1
!								
!								if (varrpos .le. nphby) then 
!									fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
!								endif
!							endif
!						enddo
!					endif
						
					call plsum(fvxyz(:,:,1))
!					if (idproc == 0) then
!						call phasespace_writef(fvxy,nxy,nphby,it,dt,&
!								&PH_VY_Y,1.,dv,0.,low)					
!					endif
					call writef(fvxyz,nxyz,nphby,it,dt,&
								&PH_VY_Y,1.,dv,0.,low)		
!					if (idproc == 0) then
!						if (nphvyyrec == 0) then
!							nphvyyrec = -1
!							fname = './DIAG/Vy_y/vy_y'
!							phasevyyfunit = get_funit(20)
!							call write_jf(fvxy,nxy,nphby,phasevyyfunit,nphvyyrec,&
!								&trim(fname),PH_VY_Y,time,it,1.,dv,0.,low)					
!						else 
!							fname = './DIAG/Vy_y/vy_y'
!							call write_jf(fvxy,nxy,nphby,phasevyyfunit,nphvyyrec,&
!								&trim(fname),PH_VY_Y,time,it,1.,dv,0.,low)					
!						endif
!					endif
				else if (xoryorz == 3) then
					vrange = 2.*fvzmax
					dv = vrange / real(nphbz)
					dvinv = 1. / dv
					low = -1.*fvzmax
					high = fvzmax
	
!					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
							pos = floor(part(2,i,1))+1
							v = part(6,i,1)
							varrpos = floor( (v + fvzmax)*dvinv + 0.5) + 1
							
							if (varrpos .le. nphbz) then 
								fvxyz(pos,varrpos,1) = fvxyz(pos,varrpos,1) + 1.
							endif
						enddo
!					else
!						do i = 1, npp(1)
!							real_tag = part(addtag_loc,i,1)
!							if (int_tag(1) == -1) then
!								pos = floor(part(2,i,1))+1
!								v = part(4,i,1)
!								varrpos = floor( (v + fvymax)*dvinv + 0.5) + 1
!								
!								if (varrpos .le. nphby) then 
!									fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
!								endif
!							endif
!						enddo
!					endif
						
					call plsum(fvxyz(:,:,1))
!					if (idproc == 0) then
!						call phasespace_writef(fvxy,nxy,nphby,it,dt,&
!								&PH_VY_Y,1.,dv,0.,low)					
!					endif
					call writef(fvxyz,nxyz,nphbz,it,dt,&
								&PH_VZ_Y,1.,dv,0.,low)		
!					if (idproc == 0) then
!						if (nphvyyrec == 0) then
!							nphvyyrec = -1
!							fname = './DIAG/Vy_y/vy_y'
!							phasevyyfunit = get_funit(20)
!							call write_jf(fvxy,nxy,nphby,phasevyyfunit,nphvyyrec,&
!								&trim(fname),PH_VY_Y,time,it,1.,dv,0.,low)					
!						else 
!							fname = './DIAG/Vy_y/vy_y'
!							call write_jf(fvxy,nxy,nphby,phasevyyfunit,nphvyyrec,&
!								&trim(fname),PH_VY_Y,time,it,1.,dv,0.,low)					
!						endif
!					endif
				endif

			end subroutine phase_space_vxyz_vs_y
			
			! This calculates the phase space at a given time for 1 direction
			! specified by xoryorz(x=1,y=2,z=3),nxyz is nx or ny or nz depending on xoryorz
			! fvxyz = phase space for x or y or z
			subroutine phase_space_vxyz_vs_z(part,fvxyz,xoryorz,nxyz,npp,dt,it,nblok,idproc)
				implicit none
				integer :: nblok,xoryorz,nxyz,idproc
	      	real,dimension(:,:,:) :: part
				integer :: it
				real :: dt
				real,dimension(:,:,:),pointer :: fvxyz
				real,dimension(:,:,:),pointer :: sfieldtemp
				integer,dimension(nblok) :: npp
				integer :: i,varrpos,xarrpos,pos
				real :: dv,low,high,v,vrange,dvinv
				character(len=32) :: fname

				integer,dimension(2) :: int_tag
				real :: real_tag
				equivalence (real_tag,int_tag)
								
				if (xoryorz == 1) then
					vrange = 2.*fvxmax
					dv = vrange / real(nphbx)
					dvinv = 1. / dv
					low = -1.*fvxmax
					high = fvxmax
					
!					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
							pos = floor(part(3,i,1))+1
							v = part(4,i,1)
							varrpos = floor( (v + fvxmax)*dvinv + 0.5) + 1
							
							if ((varrpos .le. nphbx) .and. (varrpos .gt. 0)) then 
								fvxyz(pos,varrpos,1) = fvxyz(pos,varrpos,1) + 1.
							endif
						enddo
!					else
!						do i = 1, npp(1)
!							real_tag = part(addtag_loc,i,1)
!							if (int_tag(1) == -1) then
!								pos = floor(part(2,i,1))+1
!								v = part(3,i,1)
!								varrpos = floor( (v + fvxmax)*dvinv + 0.5) + 1
!								
!								if ((varrpos .le. nphbx) .and. (varrpos .gt. 0)) then 
!									fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
!								endif
!							endif
!						enddo
!					endif
						
					call plsum(fvxyz(:,:,1))
!					if (idproc == 0) then
!						call phasespace_writef(fvxy,nxy,nphby,it,dt,&
!								&PH_VX_Y,1.,dv,0.,low)					
!					endif
					call writef(fvxyz,nxyz,nphbx,it,dt,&
								&PH_VX_Z,1.,dv,0.,low)		
!					if (idproc == 0) then
!					
!						if (nphvxyrec == 0) then
!							nphvxyrec = -1
!							fname = './DIAG/Vx_y/vx_y'
!							phasevxyfunit = get_funit(20)
!							call write_jf(fvxy,nxy,nphbx,phasevxyfunit,nphvxyrec,&
!								&trim(fname),PH_VX_Y,time,it,1.,dv,0.,low)					
!						else 
!							fname = './DIAG/Vx_y/vx_y'
!							call write_jf(fvxy,nxy,nphbx,phasevxyfunit,nphvxyrec,&
!								&trim(fname),PH_VX_Y,time,it,1.,dv,0.,low)					
!						endif
!					endif
				else if (xoryorz == 2) then
					vrange = 2.*fvymax
					dv = vrange / real(nphby)
					dvinv = 1. / dv
					low = -1.*fvymax
					high = fvymax
	
!					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
							pos = floor(part(3,i,1))+1
							v = part(5,i,1)
							varrpos = floor( (v + fvymax)*dvinv + 0.5) + 1
							
							if (varrpos .le. nphby) then 
								fvxyz(pos,varrpos,1) = fvxyz(pos,varrpos,1) + 1.
							endif
						enddo
!					else
!						do i = 1, npp(1)
!							real_tag = part(addtag_loc,i,1)
!							if (int_tag(1) == -1) then
!								pos = floor(part(2,i,1))+1
!								v = part(4,i,1)
!								varrpos = floor( (v + fvymax)*dvinv + 0.5) + 1
!								
!								if (varrpos .le. nphby) then 
!									fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
!								endif
!							endif
!						enddo
!					endif
						
					call plsum(fvxyz(:,:,1))
!					if (idproc == 0) then
!						call phasespace_writef(fvxy,nxy,nphby,it,dt,&
!								&PH_VY_Y,1.,dv,0.,low)					
!					endif
					call writef(fvxyz,nxyz,nphby,it,dt,&
								&PH_VY_Z,1.,dv,0.,low)		
!					if (idproc == 0) then
!						if (nphvyyrec == 0) then
!							nphvyyrec = -1
!							fname = './DIAG/Vy_y/vy_y'
!							phasevyyfunit = get_funit(20)
!							call write_jf(fvxy,nxy,nphby,phasevyyfunit,nphvyyrec,&
!								&trim(fname),PH_VY_Y,time,it,1.,dv,0.,low)					
!						else 
!							fname = './DIAG/Vy_y/vy_y'
!							call write_jf(fvxy,nxy,nphby,phasevyyfunit,nphvyyrec,&
!								&trim(fname),PH_VY_Y,time,it,1.,dv,0.,low)					
!						endif
!					endif
				else if (xoryorz == 3) then
					vrange = 2.*fvzmax
					dv = vrange / real(nphbz)
					dvinv = 1. / dv
					low = -1.*fvzmax
					high = fvzmax
	
!					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
							pos = floor(part(3,i,1))+1
							v = part(6,i,1)
							varrpos = floor( (v + fvzmax)*dvinv + 0.5) + 1
							
							if (varrpos .le. nphbz) then 
								fvxyz(pos,varrpos,1) = fvxyz(pos,varrpos,1) + 1.
							endif
						enddo
!					else
!						do i = 1, npp(1)
!							real_tag = part(addtag_loc,i,1)
!							if (int_tag(1) == -1) then
!								pos = floor(part(2,i,1))+1
!								v = part(4,i,1)
!								varrpos = floor( (v + fvymax)*dvinv + 0.5) + 1
!								
!								if (varrpos .le. nphby) then 
!									fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
!								endif
!							endif
!						enddo
!					endif
						
					call plsum(fvxyz(:,:,1))
!					if (idproc == 0) then
!						call phasespace_writef(fvxy,nxy,nphby,it,dt,&
!								&PH_VY_Y,1.,dv,0.,low)					
!					endif
					call writef(fvxyz,nxyz,nphbz,it,dt,&
								&PH_VZ_Z,1.,dv,0.,low)		
!					if (idproc == 0) then
!						if (nphvyyrec == 0) then
!							nphvyyrec = -1
!							fname = './DIAG/Vy_y/vy_y'
!							phasevyyfunit = get_funit(20)
!							call write_jf(fvxy,nxy,nphby,phasevyyfunit,nphvyyrec,&
!								&trim(fname),PH_VY_Y,time,it,1.,dv,0.,low)					
!						else 
!							fname = './DIAG/Vy_y/vy_y'
!							call write_jf(fvxy,nxy,nphby,phasevyyfunit,nphvyyrec,&
!								&trim(fname),PH_VY_Y,time,it,1.,dv,0.,low)					
!						endif
!					endif
				endif

			end subroutine phase_space_vxyz_vs_z
			
			end module diag32_ie
