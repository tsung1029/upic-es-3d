!-----------------------------------------------------------------------
!
      module diag32_jf
      
! Has diagnostic subroutines written by Jay Fahlen for use with 
! Viktor Decyk's 2d BEPS code
! written Jan. 2008
			
      use pinit32d, only: idrun, indx, indy, indz, ntp, psolve, tend, dt,npx&
     &,npy,npz
      use p32d_jf
      use p0d
      use pinit32d_jf
      use hdf_write32_jf
      implicit none
			include 'mpif.h'
      private
      public :: idrun, indx, indy, ntp!, modesx, modesy
      public :: tend, dt
      public :: ntfield,nfieldxrec,nfieldyrec,nphvxxrec,nphvyxrec,nphvxyrec,nphvyyrec
      public :: write_jf_diag_file,phase_space_vxy_vs_x,write_Eoft_aty,mkdir_structure
      public :: phase_slices_x,phase_slices_y,phase_space_vxy_vs_x_and_y
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
      public :: store_init_vel,calc_ene_dist,calc_vx_dist,calc_vy_dist
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
      public :: phase_vx_vy,nt_phase_vx_vy_rec,nt_phase_vx_vy_funit
      public :: phase_space_vxy_vs_y
      public :: get_ene_in_driver_region,get_Py_sumover_x,get_1D_sumover_x
!      public :: calc_ene_three_regions
      public :: get_2D_sumover_x
      public :: nt_write_div_P_vs_y_funit,nt_write_jE_sumover_x_funit,nt_write_Py_sumover_x_funit
      public :: nt_write_U_sumover_x_funit, bin_Ex_rec, bin_Ey_rec, bin_Ex_funit,bin_Ey_funit
      public :: bin_Ex_and_Ey
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


!     HDF Function declaration.
			integer,external :: sfstart, sfcreate, sfwdata, sfsdtstr, sfsdmstr
			integer,external :: sfdimid, sfsdmname, sfsdscale, sfsblsz 
			integer,external :: sfendacc, sfend, sfsnatt, sfscompress

!     HDF Constant declaration
			integer, parameter :: DFACC_CREATE = 4
			integer, parameter :: DFACC_WRITE = 2
			integer, parameter :: DFNT_CHAR = 4
			integer, parameter :: DFNT_FLOAT32 = 5
			integer, parameter :: DFNT_FLOAT64 = 6
			integer, parameter :: DFNT_INT32 = 24
			integer, parameter :: COMP_CODE_DEFLATE = 4
			integer, parameter :: SD_UNLIMITED = 0
			integer, parameter :: SD_FILL = 0

            
!      character(len=*)::idrun
!      integer::modesx=11,modesy=11
!      real::t0,tend,dt,ceng
      
      contains
      	subroutine mkdir_structure(idproc)
      		implicit none
      		integer :: idproc
      		character(len=4) :: temp
      		
				integer :: ierr, i
      		
      		if (idproc == 0) then
						if (ntden .ne. 0) then
							call mkdir_f('./DIAG/den'//char(0),ierr)
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
						if (nphxy .ne. 0) then
							call mkdir_f('./DIAG/Vx_y'//char(0),ierr)
						endif
						if (nphyy .ne. 0) then
							call mkdir_f('./DIAG/Vy_y'//char(0),ierr)
						endif
						if (nt_write_Py_vs_y .ne. 0) then
							call mkdir_f('./DIAG/Py_sumover_x'//char(0),ierr)
						endif
						if (nt_write_div_P_vs_y .ne. 0) then
							call mkdir_f('./DIAG/div_P_sumover_x'//char(0),ierr)
						endif
						if (ntp .ne. 0) then
							call mkdir_f('./DIAG/pot'//char(0),ierr)
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
			
!			subroutine calc_ene_three_regions(grad_phi_int,ESPoynt_int,jE_sum,div_ESPoynt_sum,totU,Pflow,&
!				&totjE,div_tot,num_par_cell,dt,nx,nxe,ny,nypmx,nvp,idproc,inorder)
!				implicit none
!				real,dimension(:,:,:,:) :: grad_phi_int,ESPoynt_int,jE_sum
!				real,dimension(:,:,:) :: div_ESPoynt_sum
!				real,dimension(:,:) :: Pflow
!				real,dimension(:) :: totjE, totU, div_tot
!				real :: num_par_cell,dt
!				integer :: nx, nxe, ny, nypmx, nvp, idproc, inorder
!				
!				integer :: i, j, k, nyproc,nylow,nyhigh, cury
!      	integer, dimension(nvp) :: proc_pos
!      	
!      	nyproc = ny / nvp
!      	do i = 0, nvp-1
!      		proc_pos(i+1) = ny / nvp * (i)
!      	enddo
!!				print*,nyproc,proc_pos(idproc+1),idproc,inorder
!				do j = inorder, nyproc+inorder-1
!!					print*,idproc, proc_pos(idproc+1)+j
!					cury = proc_pos(idproc+1) + j
!					!Do bottom box
!					if (cury .ge. int_pos(1) .AND. cury .lt. int_pos(2)) then
!						do i = inorder, inorder + nx
!							totU(1) = totU(1) + grad_phi_int(1,i,j,1) + grad_phi_int(2,i,j,1)
!							totjE(1) = totjE(1) + je_sum(1,i,j,1) + je_sum(2,i,j,1)
!							div_tot(1) = div_tot(1) + div_ESPoynt_sum(i,j,1)
!						enddo
!						if (cury .eq. int_pos(1)) then
!							do i = inorder, inorder + nx
!								Pflow(1,1) = Pflow(1,1) - 0.5*( ESPoynt_int(2,i,j-1,1) + ESPoynt_int(2,i,j,1) )
!!								Pflow(1,1) = Pflow(1,1) - 0.5*( ESPoynt_int(1,i,j-1,1) + ESPoynt_int(1,i,j,1) &
!!																						&	+ ESPoynt_int(2,i,j-1,1) + ESPoynt_int(2,i,j,1) )
!							enddo
!						endif
!						if (cury .eq. int_pos(2)-1) then
!							do i = inorder, inorder + nx
!								Pflow(1,3) = Pflow(1,3) + 0.5*( ESPoynt_int(2,i,j+1,1) + ESPoynt_int(2,i,j,1) )
!!								Pflow(1,3) = Pflow(1,3) + 0.5*( ESPoynt_int(1,i,j+1,1) + ESPoynt_int(1,i,j,1) &
!!																						&	+ ESPoynt_int(2,i,j+1,1) + ESPoynt_int(2,i,j,1) )
!							enddo
!						endif
!					endif
!					!Do middle box
!					if (cury .ge. int_pos(2) .AND. cury .lt. int_pos(3)) then
!						do i = inorder, inorder + nx
!							totU(2) = totU(2) + grad_phi_int(1,i,j,1) + grad_phi_int(2,i,j,1)
!							totjE(2) = totjE(2) + je_sum(1,i,j,1) + je_sum(2,i,j,1)
!							div_tot(2) = div_tot(2) + div_ESPoynt_sum(i,j,1)
!						enddo
!						if (cury .eq. int_pos(2)) then
!							do i = inorder, inorder + nx
!								Pflow(2,1) = Pflow(2,1) - 0.5*( ESPoynt_int(2,i,j-1,1) + ESPoynt_int(2,i,j,1) )
!!								Pflow(2,1) = Pflow(2,1) - 0.5*( ESPoynt_int(1,i,j-1,1) + ESPoynt_int(1,i,j,1) &
!!																						&	+ ESPoynt_int(2,i,j-1,1) + ESPoynt_int(2,i,j,1) )
!							enddo
!						endif
!						if (cury .eq. int_pos(3)-1) then
!							do i = inorder, inorder + nx
!								Pflow(2,3) = Pflow(2,3) + 0.5*( ESPoynt_int(2,i,j+1,1) + ESPoynt_int(2,i,j,1) )
!!								Pflow(2,3) = Pflow(2,3) + 0.5*( ESPoynt_int(1,i,j+1,1) + ESPoynt_int(1,i,j,1) &
!!																						&	+ ESPoynt_int(2,i,j+1,1) + ESPoynt_int(2,i,j,1) )
!							enddo
!						endif
!					endif
!					!Do top box
!					if (cury .ge. int_pos(3) .AND. cury .lt. int_pos(4)) then
!						do i = inorder, inorder + nx
!							totU(3) = totU(3) + grad_phi_int(1,i,j,1) + grad_phi_int(2,i,j,1)
!							totjE(3) = totjE(3) + je_sum(1,i,j,1) + je_sum(2,i,j,1)
!							div_tot(3) = div_tot(3) + div_ESPoynt_sum(i,j,1)
!						enddo
!						if (cury .eq. int_pos(3)) then
!							do i = inorder, inorder + nx
!								Pflow(3,1) = Pflow(3,1) - 0.5*( ESPoynt_int(2,i,j-1,1) + ESPoynt_int(2,i,j,1) )
!!								Pflow(3,1) = Pflow(3,1) - 0.5*( ESPoynt_int(1,i,j-1,1) + ESPoynt_int(1,i,j,1) &
!!																						&	+ ESPoynt_int(2,i,j-1,1) + ESPoynt_int(2,i,j,1) )
!							enddo
!						endif
!						if (cury .eq. int_pos(4)-1) then
!							do i = inorder, inorder + nx
!								Pflow(3,3) = Pflow(3,3) + 0.5*( ESPoynt_int(2,i,j+1,1) + ESPoynt_int(2,i,j,1) )
!!								Pflow(3,3) = Pflow(3,3) + 0.5*( ESPoynt_int(1,i,j+1,1) + ESPoynt_int(1,i,j,1) &
!!																						&	+ ESPoynt_int(2,i,j+1,1) + ESPoynt_int(2,i,j,1) )
!							enddo
!						endif
!					endif
!				enddo
!				
!				call plsum(totU)
!				call plsum(totjE)
!				call plsum(Pflow)
!				call plsum(div_tot)
!				totU = totU*num_par_cell
!				totjE = totje
!				Pflow = Pflow
!				div_tot = div_tot
!			end subroutine calc_ene_three_regions
			
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
			
			!This function puts all particles in a grid of their initial energy
			! by the final energy and writes it to a file
			subroutine calc_ene_dist(part,npp,time,it,idproc)
				implicit none
				real,dimension(:,:,:) :: part
				integer,dimension(:),pointer :: npp
				integer :: it,idproc
				real :: time

				real,dimension(nEneBin_init,nEneBin_final) :: ene_bins
				real,dimension(nEneBin_init,nEneBin_final,1) :: ene_bins3d
				real :: dini,dfin,inv_dini,inv_dfin
				real :: i_ene,f_ene
				integer :: i, i_arr,f_arr,funitnum
				character(len=50) :: name,fname
				character(len=5) :: nlabel

				dini = initEneRange / real(nEneBin_init)
				dfin = finalEneRange / real(nEneBin_final)
				inv_dini = 1./ dini
				inv_dfin = 1./ dfin
				
				ene_bins = 0.
				
				
				do i = 1, npp(1)
					i_ene = 0.5*(part(kivx_loc,i,1)*part(kivx_loc,i,1) + part(kivy_loc,i,1)*part(kivy_loc,i,1))
					f_ene = 0.5*(part(3,i,1)*part(3,i,1) + part(4,i,1)*part(4,i,1))
					
					i_arr = floor( i_ene * inv_dini ) + 1
					f_arr = floor( f_ene * inv_dfin ) + 1
					
					if (i_arr > nEneBin_init) then 
						i_arr = nEneBin_init
					endif
					if (f_arr > nEneBin_final) then 
						f_arr = nEneBin_final
					endif
					
					ene_bins(i_arr,f_arr) = ene_bins(i_arr,f_arr) + 1.
				enddo

				call plsum(ene_bins)
				
!				write_jf3d(f,nx,kyp,iunit,nrec,name,order)
				
				if (idproc == 0) then
					if(calc_ene_dist_rec == 0) then
						calc_ene_dist_rec = -1
						write (nlabel,'(i5)') it + 10000
						name = './DIAG/EneBin/eneBin'
						fname = trim(name)//'_'//trim(adjustl(nlabel))//'.hdf'
						calc_ene_dist_funit = get_funit(20)
						ene_bins3d(:,:,1) = ene_bins
						call write_jf(ene_bins3d,nEneBin_init,nEneBin_final,calc_ene_dist_funit,&
							&calc_ene_dist_rec,name,FIN_VS_INIT_ENE,time,it,dini,dfin,0.,0.)
					else 
						write (nlabel,'(i5)') it + 10000
						name = './DIAG/EneBin/eneBin'
						fname = trim(name)//'_'//trim(adjustl(nlabel))//'.hdf'
						calc_ene_dist_funit = get_funit(20)
						ene_bins3d(:,:,1) = ene_bins
						call write_jf(ene_bins3d,nEneBin_init,nEneBin_final,calc_ene_dist_funit,&
							&calc_ene_dist_rec,name,FIN_VS_INIT_ENE,time,it,dini,dfin,0.,0.)
					endif
				endif

			end subroutine calc_ene_dist

			!This function puts all particles in a grid of their initial vx
			! by the final vx and writes it to a file
			subroutine calc_vx_dist(part,npp,time,it,idproc)
				implicit none
				real,dimension(:,:,:) :: part
				integer,dimension(:),pointer :: npp
				integer :: it,idproc
				real :: time

				real,dimension(nVxBin_init,nVxBin_final) :: vx_bins
				real,dimension(nVxBin_init,nVxBin_final,1) :: vx_bins3d
				real :: dini,dfin,inv_dini,inv_dfin
				real :: i_vx,f_vx
				integer :: i, i_arr,f_arr,funitnum
				character(len=50) :: name,fname
				character(len=5) :: nlabel

				dini = (initVxHigh-initVxLow) / real(nVxBin_init)
				dfin = (finalVxHigh-finalVxLow) / real(nVxBin_final)
				inv_dini = 1./ dini
				inv_dfin = 1./ dfin
				
				vx_bins = 0.
				
				
				do i = 1, npp(1)
					i_vx = part(kivx_loc,i,1)
					f_vx = part(3,i,1)
					
					i_arr = floor( (i_vx - initVxLow) * inv_dini ) + 1
					f_arr = floor( (f_vx - finalVxLow) * inv_dfin ) + 1
					
					if (i_arr > nVxBin_init) then 
						i_arr = nVxBin_init
					endif
					if (i_arr < 1) then 
						i_arr = 1
					endif
					if (f_arr > nVxBin_final) then 
						f_arr = nVxBin_final
					endif
					if (f_arr < 1) then 
						f_arr = 1
					endif
					
					vx_bins(i_arr,f_arr) = vx_bins(i_arr,f_arr) + 1.
				enddo

				call plsum(vx_bins)
				
!				write_jf3d(f,nx,kyp,iunit,nrec,name,order)
				if (idproc == 0) then
					if(calc_vx_dist_rec == 0) then
						calc_vx_dist_rec = -1
						write (nlabel,'(i5)') it + 10000
						name = './DIAG/VxBin/vxbin'
						fname = trim(name)//'_'//trim(adjustl(nlabel))//'.hdf'
						calc_vx_dist_funit = get_funit(20)
						vx_bins3d(:,:,1) = vx_bins
						call write_jf(vx_bins3d,nVxBin_init,nVxBin_final,calc_vx_dist_funit,&
							&calc_vx_dist_rec,name,FIN_VS_INIT_VX,time,it,dini,dfin,initVxLow,finalVxLow)
					else 
						write (nlabel,'(i5)') it + 10000
						name = './DIAG/VxBin/vxbin'
						fname = trim(name)//'_'//trim(adjustl(nlabel))//'.hdf'
						calc_vx_dist_funit = get_funit(20)
						vx_bins3d(:,:,1) = vx_bins
						call write_jf(vx_bins3d,nVxBin_init,nVxBin_final,calc_vx_dist_funit,&
							&calc_vx_dist_rec,name,FIN_VS_INIT_VX,time,it,dini,dfin,initVxLow,finalVxLow)
					endif
				endif

			end subroutine calc_vx_dist

			!This function puts all particles in a grid of their initial vy
			! by the final vy and writes it to a file
			subroutine calc_vy_dist(part,npp,time,it,idproc)
				implicit none
				real,dimension(:,:,:) :: part
				integer,dimension(:),pointer :: npp
				integer :: it,idproc
				real :: time

				real,dimension(nVyBin_init,nVyBin_final) :: vy_bins
				real,dimension(nVyBin_init,nVyBin_final,1) :: vy_bins3d
				real :: dini,dfin,inv_dini,inv_dfin
				real :: i_vy,f_vy
				integer :: i, i_arr,f_arr,funitnum
				character(len=50) :: name,fname
				character(len=5) :: nlabel

				dini = (initVyHigh-initVyLow) / real(nVyBin_init)
				dfin = (finalVyHigh-finalVyLow) / real(nVyBin_final)
				inv_dini = 1./ dini
				inv_dfin = 1./ dfin
				
				vy_bins = 0.
								
				do i = 1, npp(1)
					i_vy = part(kivy_loc,i,1)
					f_vy = part(4,i,1)
					
					i_arr = floor( (i_vy - initVyLow) * inv_dini ) + 1
					f_arr = floor( (f_vy - finalVyLow) * inv_dfin ) + 1
					
					if (i_arr > nVyBin_init) then 
						i_arr = nVyBin_init
					endif
					if (i_arr < 1) then 
						i_arr = 1
					endif
					if (f_arr > nVyBin_final) then 
						f_arr = nVyBin_final
					endif
					if (f_arr < 1) then 
						f_arr = 1
					endif
					
					vy_bins(i_arr,f_arr) = vy_bins(i_arr,f_arr) + 1.
				enddo

				call plsum(vy_bins)

!				write_jf3d(f,nx,kyp,iunit,nrec,name,order)
				if (idproc == 0) then
					if(calc_vy_dist_rec == 0) then
						calc_vy_dist_rec = -1
						write (nlabel,'(i5)') it + 10000
						name = './DIAG/VyBin/vybin'
						fname = trim(name)//'_'//trim(adjustl(nlabel))//'.hdf'
						calc_vy_dist_funit = get_funit(20)
						vy_bins3d(:,:,1) = vy_bins
						call write_jf(vy_bins3d,nVyBin_init,nVyBin_final,calc_vy_dist_funit,&
							&calc_vy_dist_rec,name,FIN_VS_INIT_VY,time,it,dini,dfin,initVyLow,finalVyLow)
					else 
						write (nlabel,'(i5)') it + 10000
						name = './DIAG/VyBin/vybin'
						fname = trim(name)//'_'//trim(adjustl(nlabel))//'.hdf'
						calc_vy_dist_funit = get_funit(20)
						vy_bins3d(:,:,1) = vy_bins
						call write_jf(vy_bins3d,nVyBin_init,nVyBin_final,calc_vy_dist_funit,&
							&calc_vy_dist_rec,name,FIN_VS_INIT_VY,time,it,dini,dfin,initVyLow,finalVyLow)
					endif
				endif
			end subroutine calc_vy_dist
			
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
			! specified by xory(x=1,y=2),nxy is nx or ny depending on xory
			! fvxy = phase space for x or y
			subroutine phase_space_vxy_vs_x(part,fvxy,xory,nxy,npp,time,it,nblok,idproc)
				implicit none
				integer :: nblok,xory,nxy,idproc
	      	real,dimension(:,:,:) :: part
				integer :: it
				real :: time
				real,dimension(:,:,:),pointer :: fvxy
				real,dimension(:,:,:),pointer :: sfieldtemp
				integer,dimension(nblok) :: npp
				integer :: i,varrpos,xarrpos,pos
				real :: dv,low,high,v,vrange,dvinv
				character(len=32) :: fname
				integer,dimension(2) :: int_tag
				real :: real_tag
				equivalence (real_tag,int_tag)
				
				if (xory == 1) then
					vrange = 2.*fvxmax
					dv = vrange / real(nphbx)
					dvinv = 1. / dv
					low = -1.*fvxmax
					high = fvxmax
	
!					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
						
							pos = floor(part(1,i,1))+1
							v = part(3,i,1)
							varrpos = floor( (v + fvxmax)*dvinv + 0.5) + 1
							
							if (varrpos .le. nphbx) then 
								fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
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
						
					call plsum(fvxy(:,:,1))
					if (idproc == 0) then
					
						if (nphvxxrec == 0) then
							nphvxxrec = -1
							fname = './DIAG/Vx_x/vx_x'
							phasevxxfunit = get_funit(20)
							call write_jf(fvxy,nxy,nphbx,phasevxxfunit,nphvxxrec,&
								&trim(fname),PH_VX_X,time,it,1.,dv,0.,low)					
		!					call writef(fvxvy,nphbx,nphby,nt_phase_vx_vy_funit,-1,itime,itime*dt,VXVY,trim(fname),inorder)
						else 
							fname = './DIAG/Vx_x/vx_x'
							call write_jf(fvxy,nxy,nphbx,phasevxxfunit,nphvxxrec,&
								&trim(fname),PH_VX_X,time,it,1.,dv,0.,low)					
						endif
					endif
				else
					vrange = 2.*fvymax
					dv = vrange / real(nphby)
					dvinv = 1. / dv
					low = -1.*fvymax
					high = fvymax
	
!					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
							pos = floor(part(1,i,1))+1
							v = part(4,i,1)
							varrpos = floor( (v + fvymax)*dvinv + 0.5) + 1
							if (varrpos .le. nphby) then 
								fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
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
					call plsum(fvxy(:,:,1))
					if (idproc == 0) then

						if (nphvyxrec == 0) then
							nphvyxrec = -1
							fname = './DIAG/Vy_x/vy_x'
							phasevyxfunit = get_funit(20)
							call write_jf(fvxy,nxy,nphby,phasevyxfunit,nphvyxrec,&
								&trim(fname),PH_VY_X,time,it,1.,dv,0.,low)					
						else 
							fname = './DIAG/Vy_x/vy_x'
							call write_jf(fvxy,nxy,nphby,phasevxyfunit,nphvxyrec,&
								&trim(fname),PH_VY_X,time,it,1.,dv,0.,low)					
						endif
					endif
				endif

			end subroutine phase_space_vxy_vs_x

			! This calculates the phase space at a given time for 1 direction
			! specified by xory(x=1,y=2),nxy is nx or ny depending on xory
			! fvxy = phase space for x or y
			subroutine phase_space_vxy_vs_y(part,fvxy,xory,nxy,npp,time,it,nblok,idproc)
				implicit none
				integer :: nblok,xory,nxy,idproc
	      	real,dimension(:,:,:) :: part
				integer :: it
				real :: time
				real,dimension(:,:,:),pointer :: fvxy
				real,dimension(:,:,:),pointer :: sfieldtemp
				integer,dimension(nblok) :: npp
				integer :: i,varrpos,xarrpos,pos
				real :: dv,low,high,v,vrange,dvinv
				character(len=32) :: fname

				integer,dimension(2) :: int_tag
				real :: real_tag
				equivalence (real_tag,int_tag)
								
				if (xory == 1) then
					vrange = 2.*fvxmax
					dv = vrange / real(nphbx)
					dvinv = 1. / dv
					low = -1.*fvxmax
					high = fvxmax
					
!					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
							pos = floor(part(2,i,1))+1
							v = part(3,i,1)
							varrpos = floor( (v + fvxmax)*dvinv + 0.5) + 1
							
							if ((varrpos .le. nphbx) .and. (varrpos .gt. 0)) then 
								fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
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
						
					call plsum(fvxy(:,:,1))
					if (idproc == 0) then
					
						if (nphvxyrec == 0) then
							nphvxyrec = -1
							fname = './DIAG/Vx_y/vx_y'
							phasevxyfunit = get_funit(20)
							call write_jf(fvxy,nxy,nphbx,phasevxyfunit,nphvxyrec,&
								&trim(fname),PH_VX_Y,time,it,1.,dv,0.,low)					
						else 
							fname = './DIAG/Vx_y/vx_y'
							call write_jf(fvxy,nxy,nphbx,phasevxyfunit,nphvxyrec,&
								&trim(fname),PH_VX_Y,time,it,1.,dv,0.,low)					
						endif
					endif
				else
					vrange = 2.*fvymax
					dv = vrange / real(nphby)
					dvinv = 1. / dv
					low = -1.*fvymax
					high = fvymax
	
!					if (phase_keep_only_tracked .eq. 0) then
						do i = 1, npp(1)
							pos = floor(part(2,i,1))+1
							v = part(4,i,1)
							varrpos = floor( (v + fvymax)*dvinv + 0.5) + 1
							
							if (varrpos .le. nphby) then 
								fvxy(pos,varrpos,1) = fvxy(pos,varrpos,1) + 1.
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
						
					call plsum(fvxy(:,:,1))
					if (idproc == 0) then
						if (nphvyyrec == 0) then
							nphvyyrec = -1
							fname = './DIAG/Vy_y/vy_y'
							phasevyyfunit = get_funit(20)
							call write_jf(fvxy,nxy,nphby,phasevyyfunit,nphvyyrec,&
								&trim(fname),PH_VY_Y,time,it,1.,dv,0.,low)					
						else 
							fname = './DIAG/Vy_y/vy_y'
							call write_jf(fvxy,nxy,nphby,phasevyyfunit,nphvyyrec,&
								&trim(fname),PH_VY_Y,time,it,1.,dv,0.,low)					
						endif
					endif
				endif

			end subroutine phase_space_vxy_vs_y
			
			subroutine bin_Ex_and_Ey(fxye,Ex_bin,Ey_bin,nx,ny,kyp,it,time,nvp,idproc,inorder)
				implicit none
				real,dimension(:,:,:,:) :: fxye
				real,dimension(:,:,:) :: Ex_bin, Ey_bin
				integer :: nx, ny, it, nvp, idproc, inorder,kyp
				real :: time
				
				integer :: i, j, nyproc, arr, ypos,xrange, yrange
				real :: Exrange, Eyrange, dEx, dEy,invdEx,invdEy
				integer, dimension(nvp) :: proc_pos
				character(len=32) :: fname
				
				xrange = bin_E_xrange(2) - bin_E_xrange(1)
				yrange = bin_E_yrange(2) - bin_E_yrange(1)
				Exrange = bin_E_Exrange(2) - bin_E_Exrange(1)
				Eyrange = bin_E_Eyrange(2) - bin_E_Eyrange(1)
				dEx = Exrange / real(bin_E)
				dEy = Exrange / real(bin_E)
				invdEx = 1./dEx
				invdEy = 1./dEy
				
				nyproc = ny / nvp
				do i = 0, nvp-1
					proc_pos(i+1) = ny / nvp * (i)
				enddo
				do i = 0, kyp-1
					do j = 0, nx-1
	
						ypos = proc_pos(idproc+1) + i
						
						if (ypos >= bin_E_yrange(1) .AND. ypos < bin_E_yrange(2)) then
							if (j >= bin_E_xrange(1) .AND. j < bin_E_xrange(2)) then
								arr = floor( (fxye(1,j+inorder+1,i+inorder+1,1) - bin_E_Exrange(1))*invdEx + 0.5) + 1
								if (arr > 0 .AND. arr < bin_E) then
									Ex_bin(j-bin_E_xrange(1)+1,ypos-bin_E_yrange(1)+1,arr) = 1.
								endif

								arr = floor( (fxye(2,j+inorder+1,i+inorder+1,1) - bin_E_Eyrange(1))*invdEy + 0.5) + 1
								if (arr > 0 .AND. arr < bin_E) then
									Ey_bin(j-bin_E_xrange(1)+1,ypos-bin_E_yrange(1)+1,arr) = 1.
								endif
							endif
						endif
					enddo
				enddo
				
				call plsum(Ex_bin)
				call plsum(Ey_bin)

				if (idproc == 0) then
					if(bin_Ex_rec == 0) then
						bin_Ex_rec = -1
						fname = './DIAG/bin_Ex/bin_Ex'
						bin_Ex_funit = get_funit(20)
						!Use LINEAR because no gaurd cells in fvy
						call write_jf(Ex_bin,xrange,yrange,bin_E,bin_Ex_funit,&
							&bin_Ex_rec,fname,BIN_EX_LABEL,time,it,1.,1.,dEx,real(bin_E_xrange(1)),&
							&real(bin_E_yrange(1)),bin_E_Exrange(1))
					else
						fname = './DIAG/bin_Ex/bin_Ex'
						call write_jf(Ex_bin,xrange,yrange,bin_E,bin_Ex_funit,&
							&bin_Ex_rec,fname,BIN_EX_LABEL,time,it,1.,1.,dEx,real(bin_E_xrange(1)),&
							&real(bin_E_yrange(1)),bin_E_Exrange(1))
					endif

					if(bin_Ey_rec == 0) then
						bin_Ey_rec = -1
						fname = './DIAG/bin_Ey/bin_Ey'
						bin_Ey_funit = get_funit(20)
						!Use LINEAR because no gaurd cells in fvy
						call write_jf(Ey_bin,xrange,yrange,bin_E,bin_Ey_funit,&
							&bin_Ey_rec,fname,BIN_EY_LABEL,time,it,1.,1.,dEy,real(bin_E_xrange(1)),&
							&real(bin_E_yrange(1)),bin_E_Eyrange(1))
					else
						fname = './DIAG/bin_Ey/bin_Ey'
						call write_jf(Ey_bin,xrange,yrange,bin_E,bin_Ey_funit,&
							&bin_Ey_rec,fname,BIN_EY_LABEL,time,it,1.,1.,dEy,real(bin_E_xrange(1)),&
							&real(bin_E_yrange(1)),bin_E_Eyrange(1))
					endif
				endif

			end subroutine bin_Ex_and_Ey

			! This calculates the phase space for vx or vy vs. x and y
			! specified by xory(x=1,y=2),nxy is nx or ny depending on xory
			! fvxy = phase space for x or y
			subroutine phase_space_vxy_vs_x_and_y(part,fvxy,xory,nx,ny,np,npp,it,time,idproc)
				implicit none
				integer :: np,xory,nx,ny,idproc,it
				real :: time
      	real,dimension(:,:,:) :: part
				real,dimension(:,:,:),pointer :: fvxy
				integer,dimension(:),pointer :: npp
				integer :: i,varrpos,xarrpos
!				real :: dv,low,high,v,vrange,dvinv
				real :: v
				character(len=32) :: fname
				real,dimension(2) :: fmax,vrange,low,dv,dvinv,r
				real :: dx,dy,invdx,invdy
				integer,dimension(2) :: pos
				
				fmax(1) = fvx_xy_max
				fmax(2) = fvy_xy_max
				vrange(1) = 2.*fvx_xy_max
				vrange(2) = 2.*fvy_xy_max
				low(1) = -1.*fmax(1)
				low(2) = -1.*fmax(2)
				dv(1) = vrange(1) / real(nphb_xy)
				dv(2) = vrange(2) / real(nphb_xy)
				dvinv = 1./dv
				
				dx = (fvxy_xy_xrange(2) - fvxy_xy_xrange(1)) / real(nfvxy_xy_xrange)
				dy = (fvxy_xy_yrange(2) - fvxy_xy_yrange(1)) / real(nfvxy_xy_yrange)
				invdx = 1./dx
				invdy = 1./dy
				
				do i = 1, npp(1)
					r = part(1:2,i,1)
					if (r(1) > fvxy_xy_xrange(1) .AND. r(1) < fvxy_xy_xrange(2)) then
						if (r(2) > fvxy_xy_yrange(1) .AND. r(2) < fvxy_xy_yrange(2)) then
							
							pos(1) = floor( (r(1) - fvxy_xy_xrange(1))*invdx + 0.5 ) + 1
							pos(2) = floor( (r(2) - fvxy_xy_yrange(1))*invdy + 0.5 ) + 1
							
							if (pos(1) > 0 .AND. pos(1) < nfvxy_xy_xrange+1) then
								if (pos(2) > 0 .AND. pos(2) < nfvxy_xy_yrange+1) then
									
									v = part(2+xory,i,1)
									varrpos = floor( (v + fmax(xory))*dvinv(xory) + 0.5) + 1
									if (varrpos .le. nphb_xy) then 
										fvxy(pos(1),pos(2),varrpos) = fvxy(pos(1),pos(2),varrpos) + 1.
									endif
								endif
							endif
						endif
					endif
				enddo
				
				call plsum(fvxy(:,:,:))				
				
				if (xory == 1) then
					if (idproc == 0) then
						if(fvx_xy_rec == 0) then
							fvx_xy_rec = -1
							fname = './DIAG/fvx_xy/fvx_xy'
							fvx_xy_funit = get_funit(20)
							!Use LINEAR because no gaurd cells in fvy
							call write_jf(fvxy,nfvxy_xy_xrange,nfvxy_xy_yrange,nphb_xy,fvx_xy_funit,&
								&fvx_xy_rec,fname,FVX_XY_LABEL,time,it,dx,dy,dv(xory),fvxy_xy_xrange(1),&
								&fvxy_xy_yrange(1),low(xory))
						else
							fname = './DIAG/fvx_xy/fvx_xy'
							call write_jf(fvxy,nfvxy_xy_xrange,nfvxy_xy_yrange,nphb_xy,fvx_xy_funit,&
								&fvx_xy_rec,fname,FVX_XY_LABEL,time,it,dx,dy,dv(xory),fvxy_xy_xrange(1),&
								&fvxy_xy_yrange(1),low(xory))
						endif
					endif
				else
					if (idproc == 0) then
						if(fvy_xy_rec == 0) then
							fvy_xy_rec = -1
							fname = './DIAG/fvy_xy/fvy_xy'
							fvy_xy_funit = get_funit(20)
							!Use LINEAR because no gaurd cells in fvxy
							call write_jf(fvxy,nfvxy_xy_xrange,nfvxy_xy_yrange,nphb_xy,fvy_xy_funit,&
								&fvy_xy_rec,fname,FVY_XY_LABEL,time,it,dx,dy,dv(xory),fvxy_xy_xrange(1),&
								&fvxy_xy_yrange(1),low(xory))
						else
							fname = './DIAG/fvy_xy/fvy_xy'
							call write_jf(fvxy,nfvxy_xy_xrange,nfvxy_xy_yrange,nphb_xy,fvy_xy_funit,&
								&fvy_xy_rec,fname,FVY_XY_LABEL,time,it,dx,dy,dv(xory),fvxy_xy_xrange(1),&
								&fvxy_xy_yrange(1),low(xory))
						endif
					endif
				endif
			end subroutine phase_space_vxy_vs_x_and_y

			!make vx vs. vy phase space plot
			subroutine phase_vx_vy(part,fvxvy,np,npp,time,it,idimp,nblok,idproc)
				implicit none
				integer :: np, idimp, nblok,idproc,it
      	real,dimension(:,:,:) :: part
				real,dimension(:,:,:),pointer :: fvxvy
				integer,dimension(nblok) :: npp
				real :: time
				integer :: i,vxarrpos,vyarrpos
				real :: dvx,lowx,highx,vx,vxrange,dvxinv
				real :: dvy,lowy,highy,vy,vyrange,dvyinv
				character(len=32) :: fname
				
				vxrange = 2.*fvxmax
				dvx = vxrange / real(nphbx)
				dvxinv = 1. / dvx
				lowx = -1.*fvxmax
				highx = fvxmax
				vyrange = 2.*fvymax
				dvy = vyrange / real(nphby)
				dvyinv = 1. / dvy
				lowy = -1.*fvymax
				highy = fvymax

				do i = 1, npp(1)
					vx = part(3,i,1)
					vxarrpos = floor( (vx + fvxmax)*dvxinv + 0.5) + 1
					vy = part(4,i,1)
					vyarrpos = floor( (vy + fvymax)*dvyinv + 0.5) + 1
					
					if ((vxarrpos .le. nphbx) .AND. (vyarrpos .le. nphby) .AND. &
						&	(vxarrpos .gt. 0) .AND. (vyarrpos .gt. 0)) then 
						fvxvy(vxarrpos,vyarrpos,1) = fvxvy(vxarrpos,vyarrpos,1) + 1.
					endif
				enddo
				call plsum(fvxvy(:,:,1))
				if (nt_phase_vx_vy_rec == 0) then
					nt_phase_vx_vy_rec = -1
					fname = './DIAG/VxVy/vx_vy'
					nt_phase_vx_vy_funit = get_funit(20)
					call write_jf(fvxvy,nphby,nphbx,nt_phase_vx_vy_funit,nt_phase_vx_vy_rec,&
						&trim(fname),VXVY,time,it,dvx,dvy,lowx,lowy)					
!					call writef(fvxvy,nphbx,nphby,nt_phase_vx_vy_funit,-1,itime,itime*dt,VXVY,trim(fname),inorder)
				else 
					fname = './DIAG/VxVy/vx_vy'
					call write_jf(fvxvy,nphby,nphbx,nt_phase_vx_vy_funit,nt_phase_vx_vy_rec,&
						&trim(fname),VXVY,time,it,dvx,dvy,lowx,lowy)
				endif

				
			end subroutine phase_vx_vy
					
			subroutine write_Eoft_aty(lines,fxyze,nx,ny,nz,nvpy,nvpz,idproc)
				implicit none
				integer :: nx,ny,nz,nvpy,nvpz,idproc
				real,dimension(:,:,:,:,:) :: fxyze
				real,dimension(:,:) :: lines
				integer :: nyproc,nzproc,i, j, zproc, yproc, yprocpos, zprocpos
				integer, dimension(:,:,:), allocatable :: proc_pos
				real,dimension(:),allocatable,save :: sfieldtemp1
				character(len=10) :: num
				character(len=30) :: fname
				
				if (.not. allocated(sfieldtemp1)) then
					allocate(sfieldtemp1(nx))
				endif
				
				nyproc = ny / nvpy
				nzproc = nz / nvpz
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
				
				yproc = mod(idproc,nvpy) + 1; zproc = idproc / nvpy + 1
				yprocpos = proc_pos(yproc,zproc,1)
				zprocpos = proc_pos(yproc,zproc,2)
				
				lines = 0.
				do i = 1, nlines
!					linepos_3d = ny / nlines * (i-1)
					if (linepos_3d(i,1) < yprocpos + nyproc .AND. linepos_3d(i,1) >= yprocpos ) then
						if (linepos_3d(i,2) < zprocpos + nzproc .AND. linepos_3d(i,2) >= zprocpos ) then
							lines(i,:) = fxyze(1,2:nx+1,linepos_3d(i,1)-yproc+2,linepos_3d(i,2)-zproc+2,1)
						endif
					endif
				enddo
				call plsum(lines)
				if (idproc == 0) then
					do i = 1, nlines
!						linepos_3d = ny / nlines * (i-1)
						if (nlinerec(i) == 0) then
							nlinerec(i) = -1
							write (num,'(i10)') linepos_3d(i,1)
							num = adjustl(num)
							fname = 'line.'//num
							write (num,'(i10)') linepos_3d(i,2)
							num = adjustl(num)
							fname = trim(fname)//'_'//num
							linesfunit(i) = get_funit(20)
							sfieldtemp1 = lines(i,:)
							!Use LINEAR because no gaurd cells in fvxy
							call write_jf(sfieldtemp1,nx,linesfunit(i),nlinerec(i),fname,LINEAR)
						else 
							write (num,'(i10)') linepos_3d(i,1)
							num = adjustl(num)
							fname = 'line.'//num
							write (num,'(i10)') linepos_3d(i,2)
							num = adjustl(num)
							fname = trim(fname)//'_'//num
							sfieldtemp1 = lines(i,:)
							!Use LINEAR because no gaurd cells in fvxy
							call write_jf(sfieldtemp1,nx,linesfunit(i),nlinerec(i),fname,LINEAR)
						endif
					enddo
				endif
			end subroutine	write_Eoft_aty

			subroutine phase_slices_y(part,fv,xory,nx,ny,np,npp,time,it,idimp,nblok,nvp,idproc)
				implicit none
				integer :: np, idimp, nblok,xory,nx,ny,nvp,idproc,it
				real :: time
	      	real,dimension(:,:,:) :: part
!	      	real,dimension(idimp,np,nblok) :: part
				real,dimension(:,:,:),pointer :: fv
				real,dimension(:,:,:),allocatable,save :: sfieldtempx,sfieldtempy
				integer,dimension(nblok) :: npp
				integer, dimension(nvp) :: proc_pos
				integer :: i,j,varrpos,xarrpos,pos,nyproc
				real :: dv,low,high,posx,v,vrange,dvinv
				character(len=100) :: fname,num
				real :: pposhigh,pposlow
				
				if (.not. allocated(sfieldtempx)) then
					allocate(sfieldtempx(ny,nphbx,1))
				endif
				if (.not. allocated(sfieldtempy)) then
					allocate(sfieldtempy(ny,nphby,1))
				endif
				nyproc = ny / nvp
				do i = 0, nvp-1
					proc_pos(i+1) = ny / nvp * (i)
				enddo

				if (xory == 1) then
					vrange = 2.*fvxmax
					dv = vrange / real(nphbx)
					dvinv = 1. / dv
					low = -1.*fvxmax
					high = fvxmax
					do i = 1, npp(1)
						pos = floor(part(2,i,1))+1
						v = part(3,i,1)
						do j = 1, num_phsl_y
							pposlow = real(phsl_y_pos(j))
							pposhigh = real(phsl_y_pos(j)+phsl_y_thick)
							if ((part(1,i,1) .ge. pposlow) .AND. &
								& (part(1,i,1) .lt. pposhigh) ) then
								
								varrpos = floor( (v + fvxmax)*dvinv + 0.5) + 1
								if (varrpos .le. nphbx) then 
									fv(pos,varrpos,j) = fv(pos,varrpos,j) + 1.
								endif
							endif
						enddo
					enddo
					! Since each slice of this is zero on all other processors, just add them all
					do i = 1, num_phsl_y
						call plsum(fv(:,:,i))
					enddo
					if (idproc == 0) then
						do i = 1, num_phsl_y
							if (phsl_xy_nrec(i) == 0) then
								phsl_xy_nrec(i) = -1
								write (num,'(i10)') phsl_y_pos(i)
								num = adjustl(num)
								fname = './DIAG/vx_y/'//trim(num)//'/ph_xy_slice_'//trim(num)
								sfieldtempx(:,:,1) = fv(:,:,i)
								phsl_xy_funit(i) = get_funit(20)
								!Use LINEAR because no gaurd cells in fvxy
! HDF write
								call write_jf(sfieldtempx,ny,nphbx,phsl_xy_funit(i),phsl_xy_nrec(i),&
									&trim(fname),PHSL_XY,time,it,1.,1.,0.,0.)
!old write
!								call write_jf(sfieldtempx,ny,nphbx,phsl_xy_funit(i),phsl_xy_nrec(i),fname,LINEAR)
							else 
								write (num,'(i10)') phsl_y_pos(i)
								num = adjustl(num)
								fname = './DIAG/vx_y/'//trim(num)//'/ph_xy_slice_'//trim(num)
								sfieldtempx(:,:,1) = fv(:,:,i)
								!Use LINEAR because no gaurd cells in fvxy
								call write_jf(sfieldtempx,ny,nphbx,phsl_xy_funit(i),phsl_xy_nrec(i),&
									&trim(fname),PHSL_XY,time,it,1.,1.,0.,0.)
!								call write_jf(sfieldtempx,ny,nphbx,phsl_xy_funit(i),phsl_xy_nrec(i),fname,LINEAR)
							endif
						enddo
					endif
				else
					vrange = 2.*fvymax
					dv = vrange / real(nphby)
					dvinv = 1. / dv
					low = -1.*fvymax
					high = fvymax
					
					do i = 1, npp(1)
						pos = floor(part(2,i,1))+1
						v = part(4,i,1)
						do j = 1, num_phsl_y
							pposlow = real(phsl_y_pos(j))
							pposhigh = real(phsl_y_pos(j)+phsl_y_thick)
							if ((part(1,i,1) .ge. pposlow) .AND. &
								& (part(1,i,1) .lt. pposhigh) ) then
								
								varrpos = floor( (v + fvymax)*dvinv + 0.5) + 1
								if (varrpos .le. nphby) then 
									fv(pos,varrpos,j) = fv(pos,varrpos,j) + 1.
								endif
							endif
						enddo
					enddo
					! Since each slice of this is zero on all other processors, just add them all
					do i = 1, num_phsl_y
						call plsum(fv(:,:,i))
					enddo
					if (idproc == 0) then
						do i = 1, num_phsl_y
							if (phsl_yy_nrec(i) == 0) then
								phsl_yy_nrec(i) = -1
								write (num,'(i10)') phsl_y_pos(i)
								num = adjustl(num)
								fname = './DIAG/vy_y/'//trim(num)//'/ph_yy_slice_'//trim(num)
!								fname = './DIAG/vy_y/ph_yy_slice.'//num
								sfieldtempy(:,:,1) = fv(:,:,i)
								phsl_yy_funit(i) = get_funit(20)
								!Use LINEAR because no gaurd cells in fvxy
!old write
!								call write_jf(sfieldtempy,ny,nphby,phsl_yy_funit(i),phsl_yy_nrec(i),fname,LINEAR)
!HDF write
								call write_jf(sfieldtempy,ny,nphby,phsl_yy_funit(i),phsl_yy_nrec(i),&
									&trim(fname),PHSL_YY,time,it,1.,1.,0.,0.)
							else 
								write (num,'(i10)') phsl_y_pos(i)
								num = adjustl(num)
								fname = './DIAG/vy_y/'//trim(num)//'/ph_yy_slice_'//trim(num)
!								fname = './DIAG/vy_y/ph_yy_slice.'//num
								sfieldtempy(:,:,1) = fv(:,:,i)
								!Use LINEAR because no gaurd cells in fvxy
!								call write_jf(sfieldtempy,ny,nphby,phsl_yy_funit(i),phsl_yy_nrec(i),fname,LINEAR)
								call write_jf(sfieldtempy,ny,nphby,phsl_yy_funit(i),phsl_yy_nrec(i),&
									&trim(fname),PHSL_YY,time,it,1.,1.,0.,0.)
							endif
						enddo
					endif
				endif
			end subroutine phase_slices_y
			
			subroutine phase_slices_x(part,fv,xory,nx,ny,np,npp,time,it,idimp,nblok,nvp,idproc)
				implicit none
				integer :: np, idimp, nblok,xory,nx,ny,nvp,idproc,it
				real :: time
	      	real,dimension(:,:,:) :: part
!	      	real,dimension(idimp,np,nblok) :: part
				real,dimension(:,:,:),pointer :: fv
				real,dimension(:,:,:),allocatable,save :: sfieldtempx,sfieldtempy
				integer,dimension(nblok) :: npp
				integer, dimension(nvp) :: proc_pos
				integer :: i,j,varrpos,xarrpos,pos,nyproc
				real :: dv,low,high,posx,v,vrange,dvinv
				character(len=100) :: fname,num
				real :: pposhigh,pposlow
				
				if (.not. allocated(sfieldtempx)) then
					allocate(sfieldtempx(nx,nphbx,1))
				endif
				if (.not. allocated(sfieldtempy)) then
					allocate(sfieldtempy(nx,nphby,1))
				endif
				nyproc = ny / nvp
				do i = 0, nvp-1
					proc_pos(i+1) = ny / nvp * (i)
				enddo
				if (xory == 1) then
					vrange = 2.*fvxmax
					dv = vrange / real(nphbx)
					dvinv = 1. / dv
					low = -1.*fvxmax
					high = fvxmax
	
					do j = 1, num_phsl_x
						! If the jth slice is on this processor then do it
						if (phsl_x_pos(j) >= proc_pos(idproc+1) .AND. &
							&	phsl_x_pos(j) < proc_pos(idproc+1) + nyproc) then
							
							pposlow = real(phsl_x_pos(j))
							pposhigh = real(phsl_x_pos(j)+phsl_x_thick)
							do i = 1, npp(1)
								! If this particle falls in slice thickness
								if (part(2,i,1) < pposhigh .AND. part(2,i,1) >= pposlow) then
									pos = floor(part(1,i,1))+1
									v = part(3,i,1)
									varrpos = floor( (v + fvxmax)*dvinv + 0.5) + 1
									if (varrpos .le. nphbx) then 
										fv(pos,varrpos,j) = fv(pos,varrpos,j) + 1.
									endif
								endif
							enddo
						endif
					enddo
					
					! Since each slice of this is zero on all other processors, just add them all
					do i = 1, num_phsl_x
						call plsum(fv(:,:,i))
					enddo

					if (idproc == 0) then
						do i = 1, num_phsl_x
							if (phsl_xx_nrec(i) == 0) then
								phsl_xx_nrec(i) = -1
								write (num,'(i10)') phsl_x_pos(i)
								num = adjustl(num)
!								fname = './DIAG/vx_x/ph_xx_slice.'//num
								fname = './DIAG/vx_x/'//trim(num)//'/ph_xx_slice_'//trim(num)
								sfieldtempx(:,:,1) = fv(:,:,i)
								phsl_xx_funit(i) = get_funit(20)
								!Use LINEAR because no gaurd cells in fvxy
!								call write_jf(sfieldtempx,nx,nphbx,phsl_xx_funit(i),phsl_xx_nrec(i),fname,LINEAR)
								call write_jf(sfieldtempx,nx,nphbx,phsl_xx_funit(i),phsl_xx_nrec(i),&
									&trim(fname),PHSL_XX,time,it,1.,1.,0.,0.)
							else 
								write (num,'(i10)') phsl_x_pos(i)
								num = adjustl(num)
!								fname = './DIAG/vx_x/ph_xx_slice.'//num
								fname = './DIAG/vx_x/'//trim(num)//'/ph_xx_slice_'//trim(num)
								sfieldtempx(:,:,1) = fv(:,:,i)
								!Use LINEAR because no gaurd cells in fvxy
								call write_jf(sfieldtempx,nx,nphbx,phsl_xx_funit(i),phsl_xx_nrec(i),&
									&trim(fname),PHSL_XX,time,it,1.,1.,0.,0.)
!								call write_jf(sfieldtempx,nx,nphbx,phsl_xx_funit(i),phsl_xx_nrec(i),fname,LINEAR)
							endif
						enddo
					endif

				else
					vrange = 2.*fvymax
					dv = vrange / real(nphby)
					dvinv = 1. / dv
					low = -1.*fvymax
					high = fvymax
					do j = 1, num_phsl_x
						! If the jth slice is on this processor then do it
						if (phsl_x_pos(j) >= proc_pos(idproc+1) .AND. &
							&	phsl_x_pos(j) < proc_pos(idproc+1) + nyproc) then							
							pposlow = real(phsl_x_pos(j))
							pposhigh = real(phsl_x_pos(j)+phsl_x_thick)

							do i = 1, npp(1)
								! If this particle falls in slice thickness
								if (part(2,i,1) < pposhigh .AND. part(2,i,1) >= pposlow) then
									pos = floor(part(1,i,1))+1
									v = part(3,i,1)
									varrpos = floor( (v + fvymax)*dvinv + 0.5) + 1
									if (varrpos .le. nphby) then 
										fv(pos,varrpos,j) = fv(pos,varrpos,j) + 1.
									endif
								endif
							enddo
						endif
					enddo

					! Since each slice of this is zero on all other processors, just add them all
					do i = 1, num_phsl_x
						call plsum(fv(:,:,i))
					enddo
	
					if (idproc == 0) then
						do i = 1, num_phsl_x
							if (phsl_yx_nrec(i) == 0) then
								phsl_yx_nrec(i) = -1
								write (num,'(i10)') phsl_x_pos(i)
								num = adjustl(num)
!								fname = './DIAG/vy_x/ph_yx_slice.'//num
								fname = './DIAG/vy_x/'//trim(num)//'/ph_yx_slice_'//trim(num)
								sfieldtempy(:,:,1) = fv(:,:,i)
								phsl_yx_funit(i) = get_funit(20)
								!Use LINEAR because no gaurd cells in fvxy
!								call write_jf(sfieldtempy,nx,nphby,phsl_yx_funit(i),phsl_yx_nrec(i),fname,LINEAR)
								call write_jf(sfieldtempy,nx,nphby,phsl_yx_funit(i),phsl_yx_nrec(i),&
									&trim(fname),PHSL_YX,time,it,1.,1.,0.,0.)
							else 
								write (num,'(i10)') phsl_x_pos(i)
								num = adjustl(num)
!								fname = './DIAG/vy_x/ph_yx_slice.'//num
								fname = './DIAG/vy_x/'//trim(num)//'/ph_yx_slice_'//trim(num)
								sfieldtempy(:,:,1) = fv(:,:,i)
								!Use LINEAR because no gaurd cells in fvxy
								call write_jf(sfieldtempy,nx,nphby,phsl_yx_funit(i),phsl_yx_nrec(i),&
									&trim(fname),PHSL_YX,time,it,1.,1.,0.,0.)
!								call write_jf(sfieldtempy,nx,nphby,phsl_yx_funit(i),phsl_yx_nrec(i),fname,LINEAR)
							endif
						enddo
					endif
				endif

			end subroutine phase_slices_x

!Doesn't work for 3D yet!
! Taken from OSIRIS particle tracking
!subroutine write_raw( spec, no_co, g_space, n, t, n_aux, coordinates )
!subroutine write_raw( part, nx, ny, npp, nvp, t, iter, ndump, idproc, trackset )
!!-------------------------------------------------------------------------------
!  implicit none
!  
!  ! dummy variables
!	real,dimension(:,:,:) :: part
!	integer :: nx, ny, nvp, iter, ndump, idproc  !ndump is the dump number for the raw files itime/nt_raw
!	integer,dimension(:) :: npp
!	real :: t
!	type(t_track_set) :: trackset
!	
!
!!       local variables
!  integer :: sdfileID               ! HDF SD file ID number  
!!       local variables
!
!  character(80)   :: path
!  character(80)   :: full_name, file_name
!
!	real :: raw_eval
!	real, dimension(size(part,1)) :: raw_var
!  integer :: i, j, l, lbuf,iCnt
!
!  ! Buffer to hold particle data
!  real, allocatable, dimension(:,:) :: buffer
!  real, allocatable, dimension(:) :: write_buffer
!
!  integer, allocatable, dimension(:,:) :: buffer_tag
!  
!  integer :: num_par_buffer, buffer_dim
!
!  real, dimension(2) :: box_limits
!  integer, dimension(2) :: lperiodic, lmove_c
!
!  integer::npE,npar_xid,nbase_x,nbase_p,nbase_q,nbase_pE,nbase_par_xid,nbuffer_size
!!        particles written to the buffer on all nodes         
!  integer, allocatable, dimension(:) :: num_par_buffer_all
!
!!        local variables for HDF          
!  character(20) :: sdname          
!  integer       :: sdsTagID, hdf_type
!
!	!This should be 5 for the space (2) and velocity (3 for idl script), plus addtag, plus 1 for q
!  integer, dimension(5+addtag+1) :: sdsDataID !+1 because adding extra velocity dim (all 0)
!  																								! to make idl diagnostic happy
!!  integer, dimension(size(part,1)+1) :: sdsDataID !+1 because adding extra velocity dim (all 0)
!  																								! to make idl diagnostic happy
!
!  integer, dimension(1) :: hdfstart, stride, edge, dimsize
!  integer, dimension(2) :: hdfstartTag, strideTag, edgeTag, dimsizeTag
!  
!  ! local variables for MPI
!  integer :: handle, mpi_type
!  integer, dimension(mpi_status_size):: stat
!  integer :: count
!  integer :: dest, source, tag
!  integer :: dim, bnd
!	
!  integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
!  common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
!
!  integer,dimension(2) :: tags,node_conf
!  real :: temp_real
!  equivalence (temp_real,tags)
!  character(10) :: num
!	!For this case, make vdim 3 so that the third vel dim is zero to make idl script happy
!  integer :: xdim = 2, vdim = 3
!  
!  real :: left, right
!  real, dimension(2) :: lbnd, rbnd
!
!  real, external :: randum
!
!  integer :: ierr
!   
!	mpi_type = MPI_DOUBLE_PRECISION
!	hdf_type = DFNT_FLOAT64
!  
!
!  ! initialize the buffers
!  buffer_dim = npp(1)
!  num_par_buffer = 0
!  nbuffer_size = xdim + vdim
!  nbase_x=0;
!  nbase_p=nbase_x+xdim;
!  nbase_q=nbase_p+vdim;
!  nbase_pE=1+nbase_q;
!  nbase_par_xid=nbase_pE+npE
!
!  ! buffer must be organized like this so that each quantity is contiguous in memory
!  allocate( buffer(nbuffer_size, buffer_dim), stat = ierr)
!  if (ierr /= 0) then 
!		print*,"Unable to allocate memory for raw dump diagnostic : 1"
!  endif
!
!  if ( addtag /= 0) then
!    allocate( buffer_tag(2, buffer_dim), stat = ierr)
!    if (ierr /= 0) then 
!			print*,"Unable to allocate memory for raw dump diagnostic : 2"
!    endif
!  endif
!  
!  ! loop over all particles
!  l = 1
!  lbuf = 0
!
!  ! Copy all the particles to the temp buffer
!  if (raw_range(1,1) < 0) then
!  	
!  	if (raw_dump_tagged .eq. 0) then
!			do                                
!				if (l > npp(1))  exit             
!					 
!				raw_eval = 1.
!				
!				if (raw_eval > 0.) then
!				
!					if ( randum() <= raw_fraction ) then
!					 
!						lbuf = lbuf + 1
!			
!						do i=1, xdim
!						buffer(i+nbase_x, lbuf) = part(i,l,1)
!						enddo
!						do i=1, vdim-1
!						 buffer(i+nbase_p, lbuf) = part(i+nbase_p,l,1)
!						enddo
!						buffer(vdim+nbase_p, lbuf) = 0.
!						
!						if ( addtag /= 0) then
!						 temp_real = part(5,l,1)
!						 buffer_tag(1,lbuf) = tags(1)
!						 buffer_tag(2,lbuf) = tags(2)
!						endif
!				
!					endif  
!				
!				endif              
!				l = l+1
!			enddo
!		else !this if if raw_dump_tagged .ne. 0 but raw_range(1,1) < 0
!			do                                
!				if (l > npp(1)) exit
!					temp_real = part(addtag_loc,l,1)
!				if (tags(1) .eq. -1) then
!					raw_eval = 1.
!					
!					if (raw_eval > 0.) then
!					
!						if ( randum() <= raw_fraction ) then
!						 
!							lbuf = lbuf + 1
!											
!							do i=1, xdim
!							buffer(i+nbase_x, lbuf) = part(i,l,1)
!							enddo
!							do i=1, vdim-1
!							 buffer(i+nbase_p, lbuf) = part(i+nbase_p,l,1)
!							enddo
!							buffer(vdim+nbase_p, lbuf) = 0.
!							
!							if ( addtag /= 0) then
!							 temp_real = part(addtag_loc,l,1)
!							 buffer_tag(1,lbuf) = tags(1)
!							 buffer_tag(2,lbuf) = tags(2)
!							endif
!					
!						endif  
!					endif
!				endif
!				l = l+1
!			enddo
!		endif
!	else	!This is if raw_range(1,1) >= 0
!		if (raw_cycle_x_move > 0.) then
!			left = raw_range(1,1) + t*raw_cycle_x_move
!			right = raw_range(1,2) + t*raw_cycle_x_move
!			left = modulo(left,real(nx))
!			right = modulo(right,real(nx))
!			
!			!This block accounts for when the window wraps around the periodic box and there are
!			! two sections of the window.  This happens when right < left
!			if (right > left) then
!				lbnd(1) = left
!				lbnd(2) = right
!				rbnd(1) = -1.
!				rbnd(2) = -1.
!			else
!				lbnd(1) = 0.
!				lbnd(2) = right
!				rbnd(1) = left
!				rbnd(2) = nx+2.
!			endif
!		else	!raw_cycle_x_move < 0, so just use the regular raw_range
!			lbnd(1) = raw_range(1,1)
!			lbnd(2) = raw_range(1,2)
!			rbnd(1) = -1
!			rbnd(2) = -1
!		endif
!
!		do                                
!			if (l > npp(1)) exit
!			
!			! only keep if particle fits inside this box
!			if ( (part(2,l,1) .ge. raw_range(2,1)) .AND. (part(2,l,1) .le. raw_range(2,2)) ) then
!			if ( (part(3,l,1) .ge. raw_range(3,1)) .AND. (part(3,l,1) .le. raw_range(3,2)) ) then
!			if ( (part(4,l,1) .ge. raw_range(4,1)) .AND. (part(4,l,1) .le. raw_range(4,2)) ) then
!			if ( ((part(1,l,1) .ge. lbnd(1)) .AND. (part(1,l,1) .le. lbnd(2))) .OR. &
!				&  ((part(1,l,1) .ge. rbnd(1)) .AND. (part(1,l,1) .le. rbnd(2)))) then
!
!				if ( randum() <= raw_fraction ) then
!					lbuf = lbuf + 1
!		
!					do i=1, xdim
!					buffer(i+nbase_x, lbuf) = part(i,l,1)
!					enddo
!					do i=1, vdim-1
!					 buffer(i+nbase_p, lbuf) = part(i+nbase_p,l,1)
!					enddo
!					buffer(vdim+nbase_p, lbuf) = 0.
!					
!					if ( addtag /= 0) then
!					 temp_real = part(addtag_loc,l,1)
!					 buffer_tag(1,lbuf) = tags(1)
!					 buffer_tag(2,lbuf) = tags(2)
!					endif
!				endif
!			endif
!			endif
!			endif
!			endif
!			l = l+1
!		enddo
!	endif
!	
! 
!  num_par_buffer = lbuf
!  
!  ! get total number of particles
!  if (nvp > 1) then            
!                
!      ! get number of particles on other nodes
!      allocate(num_par_buffer_all(nvp))
!      
!      call MPI_GATHER(num_par_buffer, 1, mint, &
!                 num_par_buffer_all, 1, mint, &
!                 0, lworld, ierr) 
!            
!  endif
!
!  dim = 0
!  bnd = 0
!
!  if (idproc == 0) then ! root node
!    
!    ! Create the file   
!
!    ! Prepare path and file name
!    path  = './DIAG/RAW/' 
!  
!		write (num,'(i5)') ndump + 10000
!		num = adjustl(num)
!		file_name = 'RAW_'//trim(num)//'.hdf'
!
!    full_name = trim(path) // trim(file_name)
!         
!    ! Create directory if necessary
!    call mkdir_f( trim(path)//char(0), ierr )
!   
!    ! Create the HDF file
!    sdfileID = sfstart(full_name, DFACC_CREATE)
!    ! Add simulation time to the file
!    ierr = sfsnatt(sdfileID, 'TIME', DFNT_FLOAT64, 1, t) 
!    ! Add simulation timestep to the file 
!    ierr = sfsnatt(sdfileID, 'ITER', DFNT_INT32, 1, iter)
!    
!    ! Simulation box information
!    box_limits(:) = 0
!    ierr = sfsnatt(sdfileID, 'XMIN', DFNT_FLOAT64, xdim, box_limits) 
!    box_limits(1) = nx
!    box_limits(2) = ny
!    ierr = sfsnatt(sdfileID, 'XMAX', DFNT_FLOAT64, xdim, box_limits) 
!    
!		! write periodic boundaries information
!		lperiodic(:) = 1
!		ierr = sfsnatt(sdfileID, 'PERIODIC', DFNT_INT32, xdim, lperiodic) 
!	
!		! write moving window information
!		! No moving window in this code!
!		lmove_c(:)=0
!		ierr = sfsnatt(sdfileID, 'MOVE C', DFNT_INT32, xdim, lmove_c) 
!    
!    ! Add node number information to the file
!    ierr = sfsnatt(sdfileID, 'NODE NUMBER', DFNT_INT32, 1, 0)
!    ! Add node configuration information to the file
!    node_conf(1) = 1
!    node_conf(2) = nvp
!    dim = xdim
!    ierr = sfsnatt(sdfileID, 'NODE CONFIGURATION', DFNT_INT32, dim, node_conf )
!    dim = 0
!    ! Add gamma limit information to the file
!    ierr = sfsnatt(sdfileID, 'GAMMA LIMIT', DFNT_FLOAT64, 1, 100. )
!    ! Add particle fraction information to the file
!    ierr = sfsnatt(sdfileID, 'PARTICLE FRACTION', DFNT_FLOAT64, 1, &
!                   raw_fraction )
!    ! Add species ID to the file 
!    ierr = sfsnatt(sdfileID, 'SP_ID', DFNT_INT32, 1, 1)
!
!    ! Add species Name to the file
!    ierr = sfsnatt(sdfileID, 'NAME', DFNT_CHAR, 4, &
!                             trim('elec') )
!
!    ! Add rqm to the file
!    ierr = sfsnatt(sdfileID, 'RQM', DFNT_FLOAT64, 1, -1.)
!
!    ! get total number of particles to save
!    dimsize = num_par_buffer
!    do i = 2, nvp
!      dimsize = dimsize + num_par_buffer_all(i)
!    enddo
!
!    ! Create species datasets
!    do i=1, xdim
!      sdname = 'x'//char(ichar('0')+i)
!      sdsDataID(i) = sfcreate(sdfileID, sdname, hdf_type,1, dimsize)
!      ierr = sfsdtstr(sdsDataID(i), 'x_'//char(ichar('0')+i), 'k_D', 'F5.2', 'cartesian')      
!    enddo
!
!    do i=1, vdim
!      sdname = 'p'//char(ichar('0')+i)
!      sdsDataID(xdim+i) = sfcreate(sdfileID, sdname, hdf_type,1, dimsize) 
!      ierr = sfsdtstr(sdsDataID(xdim+i), 'p_'//char(ichar('0')+i), &
!                      'm_e v_th', 'F5.2', 'cartesian')
!    enddo
!
!    sdsDataID(xdim+vdim+1) = sfcreate(sdfileID, 'q', hdf_type,1, dimsize) 
!    ierr = sfsdtstr( sdsDataID(xdim+vdim+1), 'q', 'e', 'F5.2', 'cartesian')
!
!    ! write node 0 data to file
!    stride = 1
!    hdfstart = 0
!    hdfstartTag = 0
!  
!    if (num_par_buffer > 0) then
!       allocate( write_buffer(num_par_buffer) )
!       ! write data
!       edge = num_par_buffer
!       do j = 1, xdim+vdim
!				write_buffer = buffer(j,1:num_par_buffer)
!				ierr = sfwdata(sdsDataID(j), hdfstart, stride, edge, write_buffer)
!				if (ierr /= 0) then 
!				write(*,*) "Failed writing file" 
!	!			call abort_program( p_err_diagfile )
!				endif
!       enddo
!
!       deallocate( write_buffer )
!    endif
!    ! write pE information if needed
!!    if ( if_pE ) then
!!      do i=1, vdim
!!        sdname = 'pE'//char(ichar('0')+i)
!!        sdsDataID(nbase_pE+i) = sfcreate(sdfileID, sdname, hdf_type,1, dimsize) 
!!        ierr = sfsdtstr(sdsDataID(nbase_pE+i), 'pE_'//char(ichar('0')+i), &
!!                        'm_e c^2 !Mw!D0!N', 'F5.2', 'cartesian')
!!                        'm_e c^2 \omega_0', 'F5.2', 'cartesian')
!!	    if (ierr /= 0) then 
!!	      write(*,*) "Failed setting attributes" 
!!	      call abort_program( p_err_diagfile )
!!	   endif
! !     enddo
!      ! write node 0 data to file
!  
!!      if (num_par_buffer > 0) then
!!        allocate( write_buffer(num_par_buffer) )
!        ! write data
!!        edge = num_par_buffer
!!        do j = 1, vdim
!!		  write_buffer = buffer(nbase_pE+j,1:num_par_buffer)
!!		  ierr = sfwdata(sdsDataID(nbase_pE+j), hdfstart, stride, edge, write_buffer)
!!		  if (ierr /= 0) then 
!!			write(*,*) "Failed writing file" 
!!			call abort_program( p_err_diagfile )
!!		  endif
! !       enddo
!!        deallocate( write_buffer )
!!      endif
!!    endif
!    ! write par_xid information if needed
!!    if ( spec%if_par_xid ) then
!!      do i=1, xdim
!!        sdname = 'par_xid'//char(ichar('0')+i)
!!        sdsDataID(nbase_par_xid+i) = sfcreate(sdfileID, sdname, hdf_type,1, dimsize) 
!!        ierr = sfsdtstr(sdsDataID(nbase_par_xid+i), 'par_xid_'//char(ichar('0')+i), &
!!                        'c / !Mw!D0!N', 'F5.2', 'cartesian')
!!                        'c / \omega_0', 'F5.2', 'cartesian')
!!	    if (ierr /= 0) then 
!!	      write(*,*) "Failed setting attributes" 
!!	      call abort_program( p_err_diagfile )
!!	   endif
! !     enddo
!      ! write node 0 data to file
!  
!!      if (num_par_buffer > 0) then
!!        allocate( write_buffer(num_par_buffer) )
!        ! write data
!!        edge = num_par_buffer
!!        do j = 1, xdim
!!		  write_buffer = buffer(nbase_par_xid+j,1:num_par_buffer)
!!		  ierr = sfwdata(sdsdataID(nbase_par_xid+j), hdfstart, stride, edge, write_buffer)
!!		  if (ierr /= 0) then 
!!			write(*,*) "Failed writing file" 
!!			call abort_program( p_err_diagfile )
!!		  endif
! !       enddo
!!        deallocate( write_buffer )
!!      endif
!
! !   endif
!	   ! Update start position                   
!    hdfstart = hdfstart + num_par_buffer
!    
!    ! clear buffer
!    deallocate(buffer)
!    
!    ! write tag information if needed
!    if ( addtag /= 0) then
!	   dimsizeTag(1) = 2
!	   dimsizeTag(2) = dimsize(1)
!   
!	   sdname = 'tag'
!	   
!	   sdsTagID = sfcreate(sdfileID, sdname, DFNT_INT32 , 2, dimsizeTag) 
!						   
!	   ! write node 0 data to file
!	   strideTag = 1
!	   edgeTag(1) = dimsizeTag(1)
!
!	   if ( num_par_buffer > 0 ) then
!		  edgeTag(2) = num_par_buffer
!
!		  ierr = sfwdata(sdsTagID, hdfstartTag, strideTag, edgeTag, buffer_tag)
!		  if (ierr /= 0) then 
!			write(*,*) "Failed writing file" 
!!			call abort_program( p_err_diagfile )
!		  endif
!	   
!	   endif
!	   
!	   deallocate( buffer_tag )
!	   hdfstartTag(2) = hdfstartTag(2) + num_par_buffer
!    endif
!    
!    if (nvp > 1) then            
!    
!      ! loop through remaining nodes
!    
!      do i=2, nvp             
!
!		 if (num_par_buffer_all(i) > 0) then
!			
!			! prepare buffer to receive data
!			allocate( buffer(nbuffer_size, num_par_buffer_all(i)), stat = ierr)
!			if (ierr /= 0) then 
!			  write(*,*) "Unable to allocate memory for diagnostic" 
!			endif
!		  
!			source = i-1
!			tag = 1000+i
!			count = num_par_buffer_all(i) * nbuffer_size 
!		  
!			call MPI_RECV( buffer, count , mpi_type, source, &
!							  tag, lworld, stat, ierr )
!   			
!			! write data
!			edge = num_par_buffer_all(i)
!			
!			allocate( write_buffer(num_par_buffer_all(i)) )
!			do j = 1, xdim+vdim
!				 write_buffer = buffer(j,1:num_par_buffer_all(i))
!				 ierr = sfwdata(sdsDataID(j), hdfstart, stride, edge, write_buffer)
!			   if (ierr /= 0) then 
!				 write(*,*) "Failed writing file" 
!			   endif
!			enddo
!!			if ( if_pE ) then
!!			  do j = 1, npE
! !                           write_buffer = buffer(nbase_pE+j,1:num_par_buffer_all(i))
!  !                          ierr = sfwdata(sdsDataID(nbase_pE+j), hdfstart, stride, edge, write_buffer)
!!			    if (ierr /= 0) then 
!!			      write(*,*) "Failed writing file" 
!!			      call abort_program( p_err_diagfile )
!!			    endif
!!			  enddo
!!			endif
!!			if ( spec%if_par_xid ) then
!!			  do j = 1, npar_xid
! !                           write_buffer = buffer(nbase_par_xid+j,1:num_par_buffer_all(i))
!  !                          ierr = sfwdata(sdsDataID(nbase_par_xid+j), hdfstart, stride, edge, write_buffer)
!!			    if (ierr /= 0) then 
!!			      write(*,*) "Failed writing file" 
!!			      call abort_program( p_err_diagfile )
!!			    endif
!!			  enddo
!!			endif
!
!			deallocate( write_buffer )
!			
!			! clear buffer
!			deallocate(buffer)
!			
!			! Update start position        
!			hdfstart = hdfstart + num_par_buffer_all(i)
!			
!			if ( addtag /=0) then
!			   ! prepare buffer to receive data
!			   allocate( buffer_tag(2,num_par_buffer_all(i)))
!			 
!			   ! get data
!			   source = i-1
!				 tag = 1000+i
!			   count = 2 * num_par_buffer_all(i)
!			 
!			   call MPI_RECV( buffer_tag, count , MPI_INTEGER, source, &
!								 tag, lworld, stat, ierr )
!	  
!			   ! save date from node i
!			   edgeTag(2) = num_par_buffer_all(i)
!			   ierr = sfwdata(sdsTagID, hdfstartTag, strideTag, edgeTag, buffer_tag)
!			   if (ierr /= 0) then 
!				  write(*,*) "Failed writing file" 
!			   endif
!			   
!			   ! clear buffer
!			   deallocate(buffer_tag)
!			   
!			   hdfstartTag(2) = hdfstartTag(2) + num_par_buffer_all(i)
!   
!			endif
!  
!		 endif 
!       
!      enddo 
!       
!      deallocate(num_par_buffer_all) 
!    
!    endif
!    
!    ! Close the SDS
!    do j = 1, nbuffer_size
!      ierr = sfendacc(sdsDataID(j))
!    enddo
!    
!    if ( addtag /=0) then
!      ierr = sfendacc(sdsTagID)
!    endif
!    
!    ! Close HDF file
!    ierr = sfend(sdfileID)
! 
!  else  ! other nodes        
!    
!    deallocate(num_par_buffer_all) 
!    
!    if (num_par_buffer > 0) then             
!      ! send particles on local node to node 1                          
!			dest = 0
!			tag = 1000 + idproc + 1
!      count =  num_par_buffer * nbuffer_size
!      
!      call MPI_ISEND( buffer, count, mpi_type, dest, tag, &
!&                        lworld, handle, ierr )
!      call MPI_WAIT( handle, stat, ierr )
!      
!      if (addtag /= 0) then
!         count = 2 * num_par_buffer
!         
!         call MPI_ISEND( buffer_tag, count, MPI_INTEGER, dest, tag, &
!   &                        lworld, handle, ierr )
!         call MPI_WAIT( handle, stat, ierr )
!      endif
!      
!    endif
!    
!    ! clear buffer
!    deallocate(buffer)
!  endif
!
!
!end subroutine write_raw
!-------------------------------------------------------------------------------

			end module diag32_jf