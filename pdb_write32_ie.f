!-----------------------------------------------------------------------
! This module does multi-domain I/O, writing output in LLNL's PDB format.
! The PDB format is chosen because it is easy to compile on Blue Gene computers.
! PACT wouldn't compile on Blue Gene with MPI support, so each I/O node writes its own file.


      module pdb_write32_ie
      
      use globals, only: LINEAR,QUADRATIC
      use pinit32d_jf
      implicit none
      include 'mpif.h'

      private
      public :: writef
      public :: EX,EY,JX,JY,DEN,VDOTEX,VDOTEY,PHSL_XX,PHSL_XY,PHSL_YX,PHSL_YY
      public :: ESPOYNT_X,ESPOYNT_Y,FIN_VS_INIT_ENE,FIN_VS_INIT_VX,FIN_VS_INIT_VY
      public :: VDOTEX_PART,VDOTEY_PART,VDOTEX_INT,VDOTEY_INT,FVX_XY_LABEL,FVY_XY_LABEL
      public :: VDOTEX_FOLLOW_PART,VDOTEY_FOLLOW_PART,ESPOYNT_INT_X,ESPOYNT_INT_Y
      public :: EX_ENE_INT,EY_ENE_INT,DIV_ESPOYNTINT,VXVY,PH_VX_X,PH_VY_X,PH_VX_Y,PH_VY_Y
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

		interface writef
			 module procedure ipwrite3_pdb
		end interface
      
      contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
	subroutine PWRITE3_PDB(f,nx,nxv,kyp,kzp,nvpy,nvpz,name,label_code,&
		& iter,time,iorder,in_stride)
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	! this subroutine collects distributed real data f and writes it to disk
        ! The data is written to multiple files, which must be combined when doing analysis,
        ! The number of files depends on the number of nodes being used in the computation.
        ! Data is transfered in chunks to nodes doing the I/O to conserve memory,
        ! especially on machines with small amounts of it, such as BGL.

        ! As Jim and Dave did in ddcMD, this funciton will simply dump the data to disks.
        ! Each task's array piece will be stored separately.  A post-processor needs to
        ! be used to assemble the data into a single, unified array.  I'm simply dumping
        ! data this way for convenience.

	! f = input data to be written
	! nx/kyp = length of data f in x/y on each processor to write
	! nxv = first dimension of data array f, must be >= nx
	! name = file name (used only if nren < 0)
	! input: f, nx, kyp, nxv, fname
	implicit none
	integer nx,kyp,kzp, nxv,iter,nvpy,nvpz,in_stride
	real f
	character*(*) name
	dimension f(nxv,kyp,kzp,1)
	dimension in_stride(3)
	real time
	integer iorder, label_code
	! common block for parallel processing
	integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
        ! nio = number of tasks doing io
        integer nio, ioLeader, ntpg
        ! counters
        integer ii,jj,kk
	! lstat = length of status array
	parameter(lstat=10)
        integer, dimension(1) :: readyflag = 1


	! nproc = number of real or virtual processors obtained
	! lgrp = current communicator
	! lworld = MPI_COMM_WORLD communicator
	common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld
	! local data
	integer istatus, lrec, nvp, idproc, np, ioff, id, i, j, k,l,temp
	integer ierr,rc,ny,nz
	integer per
	dimension per(2)
	dimension istatus(lstat)

        character*40 title
        character*40 xtitle,xunits
        character*40 ytitle,yunits
        character*40 ztitle,zunits

	!**** PDB stuff *******************************************
        ! variables
        integer :: pfid, scratchid
        integer, parameter :: errlen=255
        character(len=errlen) :: errmsg
        character(len=70) :: fname
        character(len=10) :: fnum   ! the number of the file (num files=num I/O tasks)
        character(len=20) :: tempstr
        character(len=20) :: titlestr  ! for the variable name
        character(len=70) :: blockstr  ! for telling PDB the dimensions of the array to write out
        character(len=70) :: varstr  ! the string describing the variable to write

        ! arrays needed for PDB write
	integer,dimension(:), allocatable :: xaxis,yaxis, zaxis
	real,dimension(:), allocatable :: ftemp

        ! PDB functions
        integer :: pfopen, pfdefa, pfwrta, pfflsh, pfclos, pfgerr, pftrgt, pfcd, pfread

	integer start(3), nodestart(3)
	
        integer nyout, nzout
	real*8 time_double

	ny = kyp-3; nz = kzp-3

        ! just being lazy and recycling Jay's code
        call do_labels(label_code,title,xtitle,xunits,ytitle,yunits,ztitle,zunits)
        
	! this segment is used for mpi computers

	! determine the rank of the calling process in the communicator
	call MPI_COMM_RANK(lworld,idproc,ierr)
	! determine the size of the group associated with a communicator
	call MPI_COMM_SIZE(lworld,nvp,ierr)

        ! determine the number of tasks doing I/O
        if ( nvp <= 16 ) then
            nio = nvp
        else
            nio = max((nvp+64-1)/64,16)
            nio = min(nio,256)
        endif

        ! Now group nodes for I/O.  The member with the lowest rank will do I/O
        ! this system does not take into account network mapping
        ntpg = (nvp+nio-1)/nio      ! max number of tasks per group
        ioLeader = (idproc/ntpg)*ntpg

        ! open the files and write the leader tasks data
        ! we have to use a separate file for each iteration because the program can crash
        ! during writes, making the files useless, and I haven't found out how to make backups in Fortran
        if (idproc==ioLeader) then

            ! write out field data from this task to a scratch file
            write(fnum,'(i6.6)') idproc/ntpg
            tempstr = 'it'//trim(adjustl(tempstr))
            fname = trim('scratch-')//trim(adjustl(fnum))//'.pdb'
            fname = trim(adjustl(fname))

            scratchid = pfopen(len(trim(fname)),fname,'w')
            ! Abort if we couldn't open the file
            if (scratchid==0) then
                ierr = pfgerr(errlen, errmsg)
                print*, 'PWRITE3_PDB: could not open PDB file ', idproc, trim(errmsg)
                call MPI_ABORT(lworld,rc,ierr)
            endif

            ! write the data
            write(tempstr,'(i15)') 1
            varstr = trim(title)//'('//trim(adjustl(tempstr))//':'
            write(tempstr,'(i15)') nx
            varstr = trim(varstr)//trim(adjustl(tempstr))//','
            write(tempstr,'(i15)') 1
            varstr = trim(varstr)//trim(adjustl(tempstr))//':'
            write(tempstr,'(i15)') kyp
            varstr = trim(varstr)//trim(adjustl(tempstr))
            write(tempstr,'(i15)') 1
            varstr = trim(varstr)//','//trim(adjustl(tempstr))//':'
            write(tempstr,'(i15)') kzp
            varstr = trim(varstr)//trim(adjustl(tempstr))//')'

            if (mreal==MPI_REAL) then
                ierr = pfwrta(scratchid, len(trim(varstr)), varstr, 5, 'float', f(3:nx+2,:,:,1))
                if (ierr==0) then
                    ierr = pfgerr(errlen,errmsg)
                    print*,'PWRITE3_PDB: Could not write the current bunch of data ',idproc,errmsg
                    call MPI_ABORT(lworld,rc,ierr)
                endif
            elseif (mreal==MPI_DOUBLE_PRECISION) then
                ierr = pfwrta(scratchid, len(trim(varstr)), varstr, 6, 'double', f(3:nx+2,:,:,1))
                if (ierr==0) then
                    ierr = pfgerr(errlen,errmsg)
                    print*,'PWRITE3_PDB: Could not write the current bunch of data ',idproc,errmsg
                    call MPI_ABORT(lworld,rc,ierr)
                endif
            else
                print*, 'PWRITE3_PDB: could not determine size of real type'
                call MPI_ABORT(lworld,rc,ierr)
            endif

            ierr = pfflsh(scratchid)
            if (ierr == 0) then
                ierr = pfgerr(errlen, errmsg)
                print*, 'PWRITE3_PDB: Could not flush the file ', idproc, errmsg
                call MPI_ABORT(lworld, rc, ierr)
            endif

            write(fnum,'(i6.6)') idproc/ntpg
            write(tempstr,'(i8.8)') iter
            tempstr = 'it'//trim(adjustl(tempstr))
            fname = trim(name)//'#'//trim(adjustl(tempstr))//'#'//trim(adjustl(fnum))//'.pdb'
            fname = trim(adjustl(fname))

            pfid = pfopen(len(trim(fname)),fname,'w')

            ! Abort if we couldn't open the file
            if (pfid==0) then
                ierr = pfgerr(errlen, errmsg)
                print*, 'PWRITE3_PDB: could not open PDB file ', idproc, trim(errmsg)
                call MPI_ABORT(lworld,rc,ierr)
            endif

            if ( (nx / in_stride(1) * in_stride(1) + 1) == nx) then
                allocate(xaxis(nx/in_stride(1)+1))
            else
                allocate(xaxis(nx/in_stride(1)))
            endif

            ! we only need to write the x-axis data once
            if (idproc==0) then
                do ii=1,size(xaxis)
                    xaxis(ii)=(ii-1)*in_stride(1) + 1
                enddo

                write(tempstr,'(i15)') size(xaxis)
                varstr = 'xaxis(1:'//trim(adjustl(tempstr))//')'
                varstr = trim(varstr)
            
                ierr = pfwrta(pfid, len(trim(varstr)), varstr, 7, 'integer', xaxis)
                if (ierr==0) then
                    ierr = pfgerr(errlen,errmsg)
                    print*,'PWRITE3_PDB: Could not write x-axis data ',idproc,errmsg
                    call MPI_ABORT(lworld,rc,ierr)
                endif

                !write out the iteration
                tempstr = 'iter(1)'
                ierr = pfwrta(pfid, len(trim(tempstr)), tempstr, 7, 'integer', iter)
                if (ierr==0) then
                    ierr = pfgerr(errlen,errmsg)
                    print*,'PWRITE3_PDB: Could not write the iteration ',idproc,errmsg
                    call MPI_ABORT(lworld,rc,ierr)
                endif

                ! write out the time
                tempstr = 'time(1)'
                if (mreal==MPI_REAL) then
                    ierr = pfwrta(pfid, len(trim(tempstr)), tempstr, 5, 'float', time)
                    if (ierr==0) then
                        ierr = pfgerr(errlen,errmsg)
                        print*,'PWRITE3_PDB: Could not write the time ',idproc,errmsg
                        call MPI_ABORT(lworld,rc,ierr)
                    endif
                elseif (mreal==MPI_DOUBLE_PRECISION) then
                    ierr = pfwrta(pfid, len(trim(tempstr)), tempstr, 6, 'double', time)
                    if (ierr==0) then
                        ierr = pfgerr(errlen,errmsg)
                        print*,'PWRITE3_PDB: Could not write the time ',idproc,errmsg
                        call MPI_ABORT(lworld,rc,ierr)
                    endif
                else
                    print*, 'PWRITE3_PDB: could not determine size of real type'
                    call MPI_ABORT(lworld,rc,ierr)
                endif

            endif

            allocate(yaxis(ny/in_stride(2)+1),zaxis(nz/in_stride(3)+1))

            ! loop over the nodes in the group
            do ii=idproc,min(ntpg+idproc-1,nvp-1)

                ! write out the y-axis data for the node
                nyout = 0
                start(2) = ( ((ii - ii / nvpy * nvpy) * ny + in_stride(2) - 1) / in_stride(2) ) * in_stride(2) + 1
                do jj = start(2), (ii - ii / nvpy * nvpy + 1) * ny, in_stride(2)
                    nyout = nyout+1
                    yaxis(nyout) = jj
                enddo

                write(tempstr,'(i10)') ii
                varstr = 'yaxis'//trim(adjustl(tempstr))//'('
                write(tempstr,'(i15)') 1
                varstr = trim(varstr)//trim(adjustl(tempstr))//':'
                write(tempstr,'(i15)') nyout
                varstr = trim(varstr)//trim(adjustl(tempstr))//')'
                varstr = trim(varstr)

                ierr = pfwrta(pfid, len(trim(varstr)), varstr, 7, 'integer', yaxis)
                if (ierr==0) then
                    ierr = pfgerr(errlen,errmsg)
                    print*,'PWRITE3_PDB: Could not write the y-axis data ',idproc,errmsg
                    call MPI_ABORT(lworld,rc,ierr)
                endif

                ! write out the z-axis data for the node
                nzout = 0;
                start(3) = ( (ii / nvpy * nz + in_stride(3) - 1) / in_stride(3) ) * in_stride(3) + 1
                do jj = start(3), ( ii / nvpy + 1 ) * nz, in_stride(3)
                    nzout = nzout+1
                    zaxis(nzout) = jj
                enddo

                write(tempstr,'(i10)') ii
                varstr = 'zaxis'//trim(adjustl(tempstr))//'('
                write(tempstr,'(i15)') 1
                varstr = trim(varstr)//trim(adjustl(tempstr))//':'
                write(tempstr,'(i15)') nzout
                varstr = trim(varstr)//trim(adjustl(tempstr))//')'
                varstr = trim(varstr)

                ierr = pfwrta(pfid, len(trim(varstr)), varstr, 7, 'integer', zaxis)
                if (ierr==0) then
                    ierr = pfgerr(errlen,errmsg)
                    print*,'PWRITE3_PDB: Could not write the z-axis data ',idproc,errmsg
                    call MPI_ABORT(lworld,rc,ierr)
                endif

                ! calculate the offset for this node
                !if (ii == idproc) then
                    nodestart(2) = start(2) - (ii - ii / nvpy * nvpy) * ny
                    nodestart(3) = start(3) - ii / nvpy * nz
                !endif

                ! receive this node's data if necessary
                if (ii /= idproc) then ! some other tasks's data
                    ! send is necessary to stop buffer overflow on uBGL
                    call MPI_SEND(readyflag, size(readyflag), mint, ii, 99, lworld, ierr)
                    call MPI_RECV(f(3:nx+2,:,:,1),size(f(3:nx+2,:,:,1)),mreal,ii,99,lworld,istatus,ierr)
                endif

                ! write out the data for this node
                write(tempstr,'(i10)') ii
                titlestr = trim(title)//trim(adjustl(tempstr))

                ! We'll recycle this part below
                write(tempstr,'(i15)') 1
                blockstr = '('//trim(adjustl(tempstr))//':'
                write(tempstr,'(i15)') size(xaxis)
                blockstr = trim(blockstr)//trim(adjustl(tempstr))//','

                ! Allocate space in the file for this node's data
                write(tempstr,'(i15)') 1
                varstr = trim(adjustl(tempstr))//':'
                write(tempstr,'(i15)') nyout
                varstr = trim(varstr)//trim(adjustl(tempstr))
                write(tempstr,'(i15)') 1
                varstr = trim(varstr)//','//trim(adjustl(tempstr))//':'
                write(tempstr,'(i15)') nzout
                varstr = trim(varstr)//trim(adjustl(tempstr))//')'

                varstr = trim(titlestr)//trim(blockstr)//trim(varstr)

                ! write the data
                if (mreal==MPI_REAL) then
                    ierr = pfwrta(pfid, len(trim(varstr)), varstr, 5, 'float', &
                        & f(3:nx+2:in_stride(1), 2+nodestart(2):ny+2:in_stride(2), 2+nodestart(3):nz+2:in_stride(3), 1))
                    if (ierr==0) then
                        ierr = pfgerr(errlen,errmsg)
                        print*,'PWRITE3_PDB: Could not write the current bunch of data ',idproc,errmsg
                        call MPI_ABORT(lworld,rc,ierr)
                    endif
                elseif (mreal==MPI_DOUBLE_PRECISION) then
                    ierr = pfwrta(pfid, len(trim(varstr)), varstr, 6, 'double', &
                        & f(3:nx+2:in_stride(1), 2+nodestart(2):ny+2:in_stride(2), 2+nodestart(3):nz+2:in_stride(3), 1))
                    if (ierr==0) then
                        ierr = pfgerr(errlen,errmsg)
                        print*,'PWRITE3_PDB: Could not write the current bunch of data ',idproc,errmsg
                        call MPI_ABORT(lworld,rc,ierr)
                    endif
                else
                    print*, 'PWRITE3_PDB: could not determine size of real type'
                    call MPI_ABORT(lworld,rc,ierr)
                endif
            enddo
            
            ierr = pfclos(pfid)
            if (ierr==0) then
                ierr = pfgerr(errlen,errmsg)
                print*,'PWRITE3_PDB: Oh, woe is me, I can''t even close a file. Good-bye cruel world! ',idproc,errmsg
                call MPI_ABORT(lworld,rc,ierr)
            endif
            deallocate(xaxis,yaxis,zaxis)

            
            ! read back in the data from the scratch file
            write(tempstr,'(i15)') 1
            varstr = trim(title)//'('//trim(adjustl(tempstr))//':'
            write(tempstr,'(i15)') nx
            varstr = trim(varstr)//trim(adjustl(tempstr))//','
            write(tempstr,'(i15)') 1
            varstr = trim(varstr)//trim(adjustl(tempstr))//':'
            write(tempstr,'(i15)') kyp
            varstr = trim(varstr)//trim(adjustl(tempstr))
            write(tempstr,'(i15)') 1
            varstr = trim(varstr)//','//trim(adjustl(tempstr))//':'
            write(tempstr,'(i15)') kzp
            varstr = trim(varstr)//trim(adjustl(tempstr))//')'

            ierr = pfread(scratchid, len(trim(varstr)), varstr, f(3:nx+2,:,:,1))
            if (ierr == 0) then
                ierr = pfgerr(errlen, errmsg)
                print*, 'PWRITE3_PDB: Could not read field data from the scratch file ', idproc, errmsg
                call MPI_ABORT(lworld, rc, ierr)
            endif

            ierr = pfclos(scratchid)
            if (ierr == 0) then
                ierr = pfgerr(errlen, errmsg)
                print*, 'PWRITE3_PDB: Oh, woe is me, I can''t even close a scratch file. Good-bye cruel world! ', idproc, errmsg
                call MPI_ABORT(lworld, rc, ierr)
            endif

        else ! send the data to the leader node
            call MPI_RECV(readyflag, size(readyflag), mint, ioLeader, 99, lworld, istatus, ierr)
            call MPI_SEND(f(3:nx+2,:,:,1),size(f(3:nx+2,:,:,1)),mreal,ioLeader,99,lworld,ierr)
        endif
        
	return
	
    end subroutine pwrite3_pdb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! this is the wrapper function for pwrite3_pdb
			subroutine ipwrite3_pdb(f,nxv,nypmx,nzpmx,nvpy,nvpz,it,time,label_code,name,inorder,stride)
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
			integer :: nxv, nypmx, nzpmx,nvpy, nvpz
			character(len=*), optional :: name
			integer :: it, label_code
			integer, optional :: inorder
			real :: time
			integer, dimension(3), optional :: stride
                        integer, dimension(3) :: in_stride
			! local data
			integer :: mx, lrec, order,nblok,i,temp
			character(len=8) :: nlabel
			character(len=60) :: fname
			character(len=10), save :: sname = ':ipwrite3:'
			nblok = 1	!This is only appropriate for when there is no shared memory
			fname = ' '
			if (present(name)) then
				fname = name
			!write (nlabel,'(i8.8)') it
			!fname = trim(name)//'#it'//trim(adjustl(nlabel))
			endif
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
				call PWRITE3_PDB(f,nxv-4,nxv,nypmx,nzpmx,nvpy,nvpz,trim(fname), &
				& label_code,it,time,QUADRATIC,in_stride)
			else
				call PWRITE3_PDB(f,nxv-1,nxv,nypmx,nzpmx,nvpy,nvpz,trim(fname), &
				& label_code,it,time,LINEAR,in_stride)
			endif
		end subroutine ipwrite3_pdb

!takes the label_code and generates titles  and units for all the different diagnostics
      subroutine do_labels(label_code,title,xtitle,xunits,ytitle,yunits,ztitle,zunits)
      	implicit none
      	integer :: label_code
      	character(len=*) :: title,xtitle,xunits,ytitle,yunits,ztitle,zunits
      	
      	select case(label_code)
      	
      	case(EX)
      		title = 'Ex'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'z'
      		zunits = 'k!N!DD!N'
      	case(EY)
      		title = 'Ey'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'E!Dy!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!N'
      	case(EZ)
      		title = 'Ez'
      		xtitle = 'z'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'E!Dx!N'
      		zunits = '!Ne/m!Mw!N!Dp!Nv!Dth!N'
      	case(JX)
      		title = 'jx'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'j!Dx!N'
      		zunits = 'env!Dth!N'
      	case(JY)
      		title = 'jy'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'j!Dy!N'
      		zunits = 'env!Dth!N'
      	case(VDOTEX)
      		title = 'vx Ex'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dx!N E!Dx!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(VDOTEY)
      		title = 'vy Ey'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dy!N E!Dy!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(DEN)
      		title = 'den'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'Density'
      		zunits = 'n!D0!N'
      	case(PHSL_XX)
      		title = 'Phase Space Slice - vx vs. x'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dx!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(PHSL_XY)
      		title = 'Phase Space Slice - vx vs. y'
      		xtitle = 'y'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dx!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(PHSL_YX)
      		title = 'Phase Space Slice - vy vs. x'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dy!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(PHSL_YY)
      		title = 'Phase Space Slice - vy vs. y'
      		xtitle = 'y'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dy!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(ESPOYNT_X)
      		title = 'ES Comp. Poynting Vector Px'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'P!Dx!N'
      		zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	case(ESPOYNT_Y)
      		title = 'ES Comp. Poynting Vector Py'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'P!Dy!N'
      		zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	case(FIN_VS_INIT_ENE)
      		title = 'Initial vs. Final Energy'
      		xtitle = 'Initial Ene'
      		xunits = 'mv!Dth!N!U2!N'
      		ytitle = 'Final Ene'
      		yunits = 'mv!Dth!N!U2!N'
      		ztitle = 'Number of Particles'
      		zunits = ' '
      	case(FIN_VS_INIT_VX)
      		title = 'Initial vs. Final vx'
      		xtitle = 'Initial v!Dx!N'
      		xunits = 'v!Dth!N'
      		ytitle = 'Final v!Dx!N'
      		yunits = 'v!Dth!N'
      		ztitle = 'Number of Particles'
      		zunits = ' '
      	case(FIN_VS_INIT_VY)
      		title = 'Initial vs. Final vy'
      		xtitle = 'Initial v!Dy!N'
      		xunits = 'v!Dth!N'
      		ytitle = 'Final v!Dy!N'
      		yunits = 'v!Dth!N'
      		ztitle = 'Number of Particles'
      		zunits = ' '
      	case(VDOTEX_PART)
      		title = 'vx Ex'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dx!N E!Dx!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(VDOTEY_PART)
      		title = 'vy Ey'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dy!N E!Dy!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(VDOTEX_INT)
      		title = 'Time Integrated vx Ex'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dx!N E!Dx!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(VDOTEY_INT)
      		title = 'Time Integrated vy Ey'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dy!N E!Dy!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(FVX_XY_LABEL)
      		title = 'Phase Space vx vs x and y'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dx!N'
      		zunits = '1/v!Dth!N'
      	case(FVY_XY_LABEL)
      		title = 'Phase Space vy vs x and y'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dy!N'
      		zunits = '1/v!Dth!N'
      	case(VXVY)
      		title = 'Phase Space vx vs vy'
      		xtitle = 'v!Dx!N'
      		xunits = '1/v!Dth!N'
      		ytitle = 'v!Dy!N'
      		yunits = '1/v!Dth!N'
      		ztitle = 'Number of Particles'
      		zunits = ' '
      	case(VDOTEX_FOLLOW_PART)
      		title = 'Summed by following particle vx Ex'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dx!N E!Dx!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(VDOTEY_FOLLOW_PART)
      		title = 'Summed by following particle vy Ey'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'v!Dy!N E!Dy!N'
      		zunits = 'e/m!Mw!N!Dp!Nv!Dth!U2!N'
      	case(ESPOYNT_INT_X)
      		title = 'Time Integrated ES Comp. Poynting Vector Px'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'P!Dx!N'
      		zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	case(ESPOYNT_INT_Y)
      		title = 'Time Integrated ES Comp. Poynting Vector Py'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'P!Dy!N'
      		zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	case(EX_ENE_INT)
      		title = 'x Component of ES Field Energy'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'Ex!U2!N'
      		zunits = '(!Ne/m!Mw!N!Dp!Nv!Dth!N)!U2!N'
      	case(EY_ENE_INT)
      		title = 'y Component of ES Field Energy'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'Ey!U2!N'
      		zunits = '(!Ne/m!Mw!N!Dp!Nv!Dth!N)!U2!N'
      	case(DIV_ESPOYNTINT)
      		title = 'Time Integrated div ES Poynting Vector'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'div P!D!N'
      		zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	case(DIVESPOYNT)
      		title = 'div ES Poynting Vector'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'div P!D!N'
      		zunits = 'n!D0!Nmv!Dt!N!U2!N'
      	case(PH_VX_X)
      		title = 'Phase Space - vx vs. x'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dx!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(PH_VY_X)
      		title = 'Phase Space - vy vs. x'
      		xtitle = 'y'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dx!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(PH_VX_Y)
      		title = 'Phase Space - vx vs. y'
      		xtitle = 'y'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dx!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(PH_VY_Y)
      		title = 'Phase Space - vy vs. y'
      		xtitle = 'y'
      		xunits = 'k!N!DD!N'
      		ytitle = 'v!Dx!N'
      		yunits = '1/v!Dth!N'
      		ztitle = '# of particles'
      		zunits = ' '
      	case(BIN_EX_LABEL)
      		title = 'Ex'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'E!Dx!N'
      		zunits = '!Ne/m!Mw!N!Dp!Nv!Dth!N'
      	case(BIN_EY_LABEL)
      		title = 'Ey'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'E!Dy!N'
      		zunits = '!Ne/m!Mw!N!Dp!Nv!Dth!N'
      	case(POT)
      		title = 'pot'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'Potential'
      		zunits = '!Ne/m!Mw!N!Dp!Nv!Dth!N!U2!N'
        case(E_CHARGE)
      		title = 'qe'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'Charge Density'
      		zunits = 'en!D0!N'
        case(I_CHARGE)
      		title = 'qi'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'Charge Density'
      		zunits = 'en!D0!N'
        case(SQN_CHARGE)
      		title = 'qiSQN'
      		xtitle = 'x'
      		xunits = 'k!N!DD!N'
      		ytitle = 'y'
      		yunits = 'k!N!DD!N'
      		ztitle = 'Charge Density'
      		zunits = 'en!D0!N'


      	end select
      end subroutine do_labels


		end module pdb_write32_ie
