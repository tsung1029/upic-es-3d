!-----------------------------------------------------------------------
! This module does multi-domain I/O, writing output in LLNL's PDB format.
! The PDB format is chosen because it is easy to compile on Blue Gene computers.
! PACT wouldn't compile on Blue Gene with MPI support, so each I/O node writes its own file.


      module hdf5_write32_ie
      
      use globals, only: LINEAR,QUADRATIC
      use pinit32d_jf
      use m_pdata_utils

      implicit none
      include 'mpif.h'

      private
      public :: writef
      public :: EX,EY,JX,JY,DEN,VDOTEX,VDOTEY,PHSL_XX,PHSL_XY,PHSL_YX,PHSL_YY
      public :: ESPOYNT_X,ESPOYNT_Y,FIN_VS_INIT_ENE,FIN_VS_INIT_VX,FIN_VS_INIT_VY
      public :: VDOTEX_PART,VDOTEY_PART,VDOTEX_INT,VDOTEY_INT,FVX_XY_LABEL,FVY_XY_LABEL
      public :: VDOTEX_FOLLOW_PART,VDOTEY_FOLLOW_PART,ESPOYNT_INT_X,ESPOYNT_INT_Y
      public :: EX_ENE_INT,EY_ENE_INT,DIV_ESPOYNTINT,VXVY
      public :: PH_VX_X,PH_VY_X,PH_VZ_X,PH_VX_Y,PH_VY_Y,PH_VZ_Y,PH_VX_Z,PH_VY_Z,PH_VZ_Z
      public :: DIVESPOYNT,BIN_EX_LABEL,BIN_EY_LABEL,POT,EZ
      public :: E_CHARGE, I_CHARGE, SQN_CHARGE
		
		save
      integer,parameter :: EX = 1,EY = 2,JX =3,JY = 4,DEN = 5,VDOTEX = 6,VDOTEY = 7
      integer,parameter :: PHSL_XX = 8,PHSL_XY=9,PHSL_YX=10,PHSL_YY=11
      integer,parameter :: ESPOYNT_X=12,ESPOYNT_Y=13,FIN_VS_INIT_ENE=14,FIN_VS_INIT_VX=15
      integer,parameter :: FIN_VS_INIT_VY=16,VDOTEX_PART=17,VDOTEY_PART=18
      integer,parameter :: VDOTEX_INT=19,VDOTEY_INT=20,FVX_XY_LABEL=21,FVY_XY_LABEL=22
      integer,parameter :: VDOTEX_FOLLOW_PART = 23,VDOTEY_FOLLOW_PART=24
      integer,parameter :: ESPOYNT_INT_X=25,ESPOYNT_INT_Y=26,EX_ENE_INT=27,EY_ENE_INT=28
      integer,parameter :: DIV_ESPOYNTINT=29,VXVY=30,PH_VX_X=31,PH_VY_X=32,PH_VX_Y=33,PH_VY_Y=34
      integer,parameter :: DIVESPOYNT=35,BIN_EX_LABEL=36,BIN_EY_LABEL=37,POT=38,EZ=39
      integer,parameter :: E_CHARGE=40, I_CHARGE=41, SQN_CHARGE=42
      integer,parameter :: PH_VZ_X=43, PH_VZ_Y=44, PH_VX_Z=45, PH_VY_Z=46, PH_VZ_Z=47

		interface writef
			 module procedure ipwrite3_hdf5
          	 module procedure ipwrite2_hdf5
		end interface
      
      contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! this is the wrapper function for pwrite_hdf3d
			subroutine ipwrite3_hdf5(f,nxv,nypmx,nzpmx,nvpy,nvpz,it,dt,label_code,inorder,stride)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! this subroutine collects distributed real 3d data f 
! and it writes to a file.  Only some of the data can be written.
! this = ufield2d descriptor of data
! f = input data to be written
! name = file name (used only if new file is to be opened)
! it = iteration
! inorder = interpolation order, 1 = LINEAR, 2 = QUADRATIC
			implicit none
			real, dimension(:,:,:,:) :: f
			integer :: nxv, nypmx, nzpmx,nvpy, nvpz, nx, ny, nz
			integer :: it, label_code
			integer, optional :: inorder
			real :: dt
			integer, dimension(3), optional :: stride
         integer, dimension(3) :: in_stride

         

			! local data
			integer :: mx, lrec, order,nblok,i,temp
         character(len=40) :: label, units
         character(len=40), dimension(3) :: xname, xlabel, xunits


			order = QUADRATIC
			if (present(inorder)) then
				if ((inorder >= 1) .and. (inorder <= 2)) then
					 order = inorder
				endif
			endif

         if ( .NOT. present(stride) ) then
            in_stride(1:3) = 1
         else
            in_stride(1:3) = stride(1:3)
         endif
			
			if (order==QUADRATIC) then
            nx = nxv-4
            ny = nypmx-3
            nz = nzpmx-3
			else
            nx = nxv-1
            ny = nypmx-3
            nz = nzpmx-3
			endif
         
         call do_labels_3D(label_code, label, units, xname, xlabel, xunits)
         call pwrite_data(f(3:nx+2, 3:ny+2, 3:nz+2,1), label, it, dt, &
               & nvpy, nvpz, in_stride, label, units, &
               & xlabel, xunits, (/0.0, 0.0, 0.0/), path="DIAG/"//trim(label)//"/")

		end subroutine ipwrite3_hdf5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! this is the wrapper function for pwrite_hdf2d
			subroutine ipwrite2_hdf5(f,nxv,nypmx,it,dt,label_code,&
							&ben_delta1,ben_delta2,ben_min1,ben_min2,inorder,stride)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! this subroutine collects distributed real 3d data f 
! and it writes to a file.  Only some of the data can be written.
! this = ufield2d descriptor of data
! f = input data to be written
! name = file name (used only if new file is to be opened)
! it = iteration
! inorder = interpolation order, 1 = LINEAR, 2 = QUADRATIC
			implicit none
			real, dimension(:,:,:) :: f
			integer :: nxv, nypmx, nzpmx,nvpy, nvpz, nx, ny, nz
			integer :: it, label_code
			integer, optional :: inorder
			real :: dt
			integer, dimension(3), optional :: stride
         integer, dimension(3) :: in_stride
         real :: ben_delta1,ben_delta2,ben_min1,ben_min2
         real, dimension(2) :: benmin, bendelta

			! local data
			integer :: mx, lrec, order,nblok,i,temp
         character(len=40) :: label, units
         character(len=40), dimension(3) :: xname, xlabel, xunits

		benmin(1) = ben_min1
		benmin(2) = ben_min2
		bendelta(1) = ben_delta1
		bendelta(2) = ben_delta2
		
			order = QUADRATIC
			if (present(inorder)) then
				if ((inorder >= 1) .and. (inorder <= 2)) then
					 order = inorder
				endif
			endif

         if ( .NOT. present(stride) ) then
            in_stride(1:3) = 1
         else
            in_stride(1:3) = stride(1:3)
         endif
			
			if (order==QUADRATIC) then
            nx = nxv-4
            ny = nypmx-3
!            nz = nzpmx-3
			else
            nx = nxv-1
            ny = nypmx-3
 !           nz = nzpmx-3
			endif

!!!! I don't know why the split between 3d and 2d labels, but hopefully this works         
         call do_labels_3D(label_code, label, units, xname, xlabel, xunits)
!
         if ((label_code .eq. PH_VY_X) .OR. (label_code .eq. PH_VX_X) &
         &.OR. (label_code .eq. PH_VZ_X)) then
         	nx = nxv-4
			ny = nypmx
		 endif
         if ((label_code .eq. PH_VY_Y) .OR. (label_code .eq. PH_VX_Y) &
         &.OR. (label_code .eq. PH_VZ_Y)) then
         	nx = nxv-3
			ny = nypmx
         endif
         if ((label_code .eq. PH_VY_Z) .OR. (label_code .eq. PH_VX_Z) &
         &.OR. (label_code .eq. PH_VZ_Z)) then
         	nx = nxv-3
			ny = nypmx
         endif

      call pwrite_data_direct(f(3:nx+2,:,1), 0, label, it,&
                          &dt,.false., 1, 0, label,units&
                          &,xlabel,xunits,benmin,bendelta,                      &
                          &.false.,path="DIAG/"//trim(label)//"/")
!      subroutine pwrite_hdf2d_direct(datain, idsend, dataname, timestep,&
!                          &dt,doappend, nappend, append_pos, label,units&
!                          &,xlabel,xunits,xmin,dx,                      &
!                          &docompress,path)

		end subroutine ipwrite2_hdf5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!takes the label_code and generates titles  and units for all the different diagnostics
      subroutine do_labels_3D(label_code,label,units,xname,xlabel,xunits)
      	implicit none
      	integer :: label_code
      	character(len=*) :: label,units
         character(len=*), dimension(3) :: xname, xlabel, xunits
         integer :: ierror, ierror2
      	
      	select case(label_code)
      	
      	case(EX)
      		label = 'Ex'
            units = 'm!Mw!N!D0!N!U2!N!9d/e'
            xname(1) = 'x1'
            xname(2) = 'x2'
            xname(3) = 'x3'
      		xlabel(1) = 'X'
            xlabel(2) = 'Y'
            xlabel(3) = 'Z'
      		xunits = '!9d'
      	case(EY)
      		label = 'Ey'
            units = 'm!Mw!N!D0!N!U2!N!9d/e'
            xname(1) = 'x1'
            xname(2) = 'x2'
            xname(3) = 'x3'
      		xlabel(1) = 'X'
            xlabel(2) = 'Y'
            xlabel(3) = 'Z'
      		xunits = '!9d'
      	case(EZ)
      		label = 'Ez'
            units = 'm!Mw!N!D0!N!U2!N!9d/e'
            xname(1) = 'x1'
            xname(2) = 'x2'
            xname(3) = 'x3'
      		xlabel(1) = 'X'
            xlabel(2) = 'Y'
            xlabel(3) = 'Z'
      		xunits = '!9d'
      	!case(JX)
      	!	title = 'jx'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'j!Dx!N'
      	!	zunits = 'env!Dth!N'
      	!case(JY)
      	!	title = 'jy'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'j!Dy!N'
      	!	zunits = 'env!Dth!N'
      	!case(VDOTEX)
      	!	title = 'vx Ex'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'v!Dx!N E!Dx!N'
      	!	zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	!case(VDOTEY)
      	!	title = 'vy Ey'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'v!Dy!N E!Dy!N'
      	!	zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	!case(DEN)
      	!	label = 'Density'
         !   units = 'n!D0!N'
         !   xname(1) = 'x1'
         !   xname(2) = 'x2'
         !   xname(3) = 'x3'
      	!	xlabel(1) = 'X'
         !   xlabel(2) = 'Y'
         !   xlabel(3) = 'Z'
      	!	xunits = '!9d'
      	!case(PHSL_XX)
      	!	title = 'Phase Space Slice - vx vs. x'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'v!Dx!N'
      	!	yunits = '1/v!Dth!N'
      	!	ztitle = '# of particles'
      	!	zunits = ' '
      	!case(PHSL_XY)
      	!	title = 'Phase Space Slice - vx vs. y'
      	!	xtitle = 'y'
      	!	xunits = '!9d'
      	!	ytitle = 'v!Dx!N'
      	!	yunits = '1/v!Dth!N'
      	!	ztitle = '# of particles'
      	!	zunits = ' '
      	!case(PHSL_YX)
      	!	title = 'Phase Space Slice - vy vs. x'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'v!Dy!N'
      	!	yunits = '1/v!Dth!N'
      	!	ztitle = '# of particles'
      	!	zunits = ' '
      	!case(PHSL_YY)
      	!	title = 'Phase Space Slice - vy vs. y'
      	!	xtitle = 'y'
      	!	xunits = '!9d'
      	!	ytitle = 'v!Dy!N'
      	!	yunits = '1/v!Dth!N'
      	!	ztitle = '# of particles'
      	!	zunits = ' '
      	!case(ESPOYNT_X)
      	!	title = 'ES Comp. Poynting Vector Px'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'P!Dx!N'
      	!	zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	!case(ESPOYNT_Y)
      	!	title = 'ES Comp. Poynting Vector Py'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'P!Dy!N'
      	!	zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	!case(FIN_VS_INIT_ENE)
      	!	title = 'Initial vs. Final Energy'
      	!	xtitle = 'Initial Ene'
      	!	xunits = 'mv!Dth!N!U2!N'
      	!	ytitle = 'Final Ene'
      	!	yunits = 'mv!Dth!N!U2!N'
      	!	ztitle = 'Number of Particles'
      	!	zunits = ' '
      	!case(FIN_VS_INIT_VX)
      	!	title = 'Initial vs. Final vx'
      	!	xtitle = 'Initial v!Dx!N'
      	!	xunits = 'v!Dth!N'
      	!	ytitle = 'Final v!Dx!N'
      	!	yunits = 'v!Dth!N'
      	!	ztitle = 'Number of Particles'
      	!	zunits = ' '
      	!case(FIN_VS_INIT_VY)
      	!	title = 'Initial vs. Final vy'
      	!	xtitle = 'Initial v!Dy!N'
      	!	xunits = 'v!Dth!N'
      	!	ytitle = 'Final v!Dy!N'
      	!	yunits = 'v!Dth!N'
      	!	ztitle = 'Number of Particles'
      	!	zunits = ' '
      	!case(VDOTEX_PART)
      	!	title = 'vx Ex'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'v!Dx!N E!Dx!N'
      	!	zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	!case(VDOTEY_PART)
      	!	title = 'vy Ey'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'v!Dy!N E!Dy!N'
      	!	zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	!case(VDOTEX_INT)
      	!	title = 'Time Integrated vx Ex'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'v!Dx!N E!Dx!N'
      	!	zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	!case(VDOTEY_INT)
      	!	title = 'Time Integrated vy Ey'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'v!Dy!N E!Dy!N'
      	!	zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	!case(FVX_XY_LABEL)
      	!	title = 'Phase Space vx vs x and y'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'v!Dx!N'
      	!	zunits = '1/v!Dth!N'
      	!case(FVY_XY_LABEL)
      	!	title = 'Phase Space vy vs x and y'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'v!Dy!N'
      	!	zunits = '1/v!Dth!N'
      	!case(VXVY)
      	!	title = 'Phase Space vx vs vy'
      	!	xtitle = 'v!Dx!N'
      	!	xunits = '1/v!Dth!N'
      	!	ytitle = 'v!Dy!N'
      	!	yunits = '1/v!Dth!N'
      	!	ztitle = 'Number of Particles'
      	!	zunits = ' '
      	!case(VDOTEX_FOLLOW_PART)
      	!	title = 'Summed by following particle vx Ex'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'v!Dx!N E!Dx!N'
      	!	zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	!case(VDOTEY_FOLLOW_PART)
      	!	title = 'Summed by following particle vy Ey'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'v!Dy!N E!Dy!N'
      	!	zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	!case(ESPOYNT_INT_X)
      	!	title = 'Time Integrated ES Comp. Poynting Vector Px'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'P!Dx!N'
      	!	zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	!case(ESPOYNT_INT_Y)
      	!	title = 'Time Integrated ES Comp. Poynting Vector Py'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'P!Dy!N'
      	!	zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	!case(EX_ENE_INT)
      	!	title = 'x Component of ES Field Energy'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'Ex!U2!N'
      	!	zunits = '(!Ne/m!Mw!N!Dp!Nv!Dth!N)!U2!N'
      	!case(EY_ENE_INT)
      	!	title = 'y Component of ES Field Energy'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'Ey!U2!N'
      	!	zunits = '(!Ne/m!Mw!N!Dp!Nv!Dth!N)!U2!N'
      	!case(DIV_ESPOYNTINT)
      	!	title = 'Time Integrated div ES Poynting Vector'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'div P!D!N'
      	!	zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	!case(DIVESPOYNT)
      	!	title = 'div ES Poynting Vector'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'div P!D!N'
      	!	zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	case(PH_VX_X)
      		label = 'Vx_x'
            units = 'm!Mw!N!D0!N!U2!N!9d/e'
      		xname(1) = 'x1'
			xlabel(1) = 'X'
      		xunits(1) = '!9d'
      		xname(2) = 'v!Dx!N'
			xlabel(2) = 'V_X'
      		xunits(2) = '1/v!Dth!N'
      		xname(3) = '# of particles'
			xlabel(3) = '# of particles'
      		xunits(3) = ' '
      	case(PH_VY_X)
      		label = 'Vy_x'
            units = 'm!Mw!N!D0!N!U2!N!9d/e'
      		xname(1) = 'x1'
			xlabel(1) = 'X'
      		xunits(1) = '!9d'
      		xname(2) = 'v!Dy!N'
			xlabel(2) = 'V_Y'
      		xunits(2) = '1/v!Dth!N'
      		xname(3) = '# of particles'
			xlabel(3) = '# of particles'
      		xunits(3) = ' '
      	case(PH_VZ_X)
      		label = 'Vz_x'
            units = 'm!Mw!N!D0!N!U2!N!9d/e'
      		xname(1) = 'x1'
			xlabel(1) = 'X'
      		xunits(1) = '!9d'
      		xname(2) = 'v!Dz!N'
			xlabel(2) = 'V_Z'
      		xunits(2) = '1/v!Dth!N'
      		xname(3) = '# of particles'
			xlabel(3) = '# of particles'
      		xunits(3) = ' '
      	case(PH_VX_Y)
      		label = 'Vx_y'
            units = 'm!Mw!N!D0!N!U2!N!9d/e'
      		xname(1) = 'x2'
			xlabel(1) = 'Y'
      		xunits(1) = '!9d'
      		xname(2) = 'v!Dx!N'
			xlabel(2) = 'V_X'
      		xunits(2) = '1/v!Dth!N'
      		xname(3) = '# of particles'
			xlabel(3) = '# of particles'
      		xunits(3) = ' '
      	case(PH_VY_Y)
      		label = 'Vy_y'
            units = 'm!Mw!N!D0!N!U2!N!9d/e'
      		xname(1) = 'x2'
			xlabel(1) = 'Y'
      		xunits(1) = '!9d'
      		xname(2) = 'v!Dy!N'
			xlabel(2) = 'V_Y'
      		xunits(2) = '1/v!Dth!N'
      		xname(3) = '# of particles'
			xlabel(3) = '# of particles'
      		xunits(3) = ' '
      	case(PH_VZ_Y)
      		label = 'Vz_y'
            units = 'm!Mw!N!D0!N!U2!N!9d/e'
      		xname(1) = 'x2'
			xlabel(1) = 'Y'
      		xunits(1) = '!9d'
      		xname(2) = 'v!Dz!N'
			xlabel(2) = 'V_Z'
      		xunits(2) = '1/v!Dth!N'
      		xname(3) = '# of particles'
			xlabel(3) = '# of particles'
      		xunits(3) = ' '
      	case(PH_VX_Z)
      		label = 'Vx_z'
            units = 'm!Mw!N!D0!N!U2!N!9d/e'
      		xname(1) = 'x3'
			xlabel(1) = 'Z'
      		xunits(1) = '!9d'
      		xname(2) = 'v!Dx!N'
			xlabel(2) = 'V_X'
      		xunits(2) = '1/v!Dth!N'
      		xname(3) = '# of particles'
			xlabel(3) = '# of particles'
      		xunits(3) = ' '
      	case(PH_VY_Z)
      		label = 'Vy_z'
            units = 'm!Mw!N!D0!N!U2!N!9d/e'
      		xname(1) = 'x3'
			xlabel(1) = 'Z'
      		xunits(1) = '!9d'
      		xname(2) = 'v!Dy!N'
			xlabel(2) = 'V_Y'
      		xunits(2) = '1/v!Dth!N'
      		xname(3) = '# of particles'
			xlabel(3) = '# of particles'
      		xunits(3) = ' '
      	case(PH_VZ_Z)
      		label = 'Vz_z'
            units = 'm!Mw!N!D0!N!U2!N!9d/e'
      		xname(1) = 'x3'
			xlabel(1) = 'Z'
      		xunits(1) = '!9d'
      		xname(2) = 'v!Dz!N'
			xlabel(2) = 'V_Z'
      		xunits(2) = '1/v!Dth!N'
      		xname(3) = '# of particles'
			xlabel(3) = '# of particles'
      		xunits(3) = ' '
      	!case(BIN_EX_LABEL)
      	!	title = 'Ex'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'E!Dx!N'
      	!	zunits = '!Ne/m!Mw!N!Dp!Nv!Dth!N'
      	!case(BIN_EY_LABEL)
      	!	title = 'Ey'
      	!	xtitle = 'x'
      	!	xunits = '!9d'
      	!	ytitle = 'y'
      	!	yunits = '!9d'
      	!	ztitle = 'E!Dy!N'
      	!	zunits = '!Ne/m!Mw!N!Dp!Nv!Dth!N'
      	case(POT)
      		label = 'Pot'
      		units = 'm!Mw!N!D0!N!U2!N!9d!U2!N/e'
            xname(1) = 'x1'
            xname(2) = 'x2'
            xname(3) = 'x3'
      		xlabel(1) = 'X'
            xlabel(2) = 'Y'
            xlabel(3) = 'Z'
      		xunits = '!9d'
        case(E_CHARGE)
      		label = 'Qe'
      		units = 'en!D0!N'
            xname(1) = 'x1'
            xname(2) = 'x2'
            xname(3) = 'x3'
      		xlabel(1) = 'X'
            xlabel(2) = 'Y'
            xlabel(3) = 'Z'
      		xunits = '!9d'
        case(I_CHARGE)
      		label = 'Qi'
      		units = 'en!D0!N'
            xname(1) = 'x1'
            xname(2) = 'x2'
            xname(3) = 'x3'
      		xlabel(1) = 'X'
            xlabel(2) = 'Y'
            xlabel(3) = 'Z'
      		xunits = '!9d'
        case(SQN_CHARGE)
      		label = 'QiSQN'
      		units = 'en!D0!N'
            xname(1) = 'x1'
            xname(2) = 'x2'
            xname(3) = 'x3'
      		xlabel(1) = 'X'
            xlabel(2) = 'Y'
            xlabel(3) = 'Z'
      		xunits = '!9d'
        case default
            print*, 'You need to setup labels for this diagnostic'
            call MPI_ABORT(MPI_COMM_WORLD, ierror, ierror2)
        end select

      end subroutine do_labels_3D

      subroutine do_labels_2D(label_code,label,units,xlabel,xunits)
      	implicit none
      	integer :: label_code
      	character(len=*) :: label,units
         character(len=*), dimension(3) :: xlabel, xunits
         integer :: ierror, ierror2
      	
         print*, '2D labels not setup yet...you need to do that yourself.'
         call MPI_ABORT(MPI_COMM_WORLD, ierror, ierror2)

      end subroutine do_labels_2D

		end module hdf5_write32_ie
