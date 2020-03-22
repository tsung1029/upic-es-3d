!-----------------------------------------------------------------------
!
       module psimul32d_ie
       ! Higher level subroutines for electrostatics
       ! written by viktor k. decyk, ucla
       ! copyright 1999, regents of the university of california
       ! update: may 2, 2008
       !
       use globals, only: LINEAR, QUADRATIC
       use pdiag32d, only: vdist, bfopen
       !     use espush2d, only: dpost, push, rpush, pushzf, rpushzf, pushglx, &
       !    &rpushglx, pushgl, rpushgl, fft-
       use pespush32d_ie, only: get_funit, plbcast, plsum, plmax, wrdata, &
       &rddata, writebf, paguard, pcguard, fft
       !     use field2d, only: sguard, aguard, cguard, spois, pois, gtmodes
       use pfield32d, only: aguard, cguard, pois, gtmodes

       use m_h5_diagnostic_utilities, only: add_h5_atribute, detect_precision, start_hdf5, stop_hdf5
       use m_system

       use par_track_ie

       use HDF5

       implicit none
       include 'mpif.h'

       private
       public :: restart_write, restart_exists
       public :: restart_read
       public :: ion_write
       !     public :: dpostg, pushg, pushzfg, pushglg, pushglxg
       public :: initmodediag, initveldiag
       public :: dendiag, veldiag, potdiag, esenergy

       character(len=70), save :: restart_name_old

        interface write_dataset
          module procedure write_dataset1d
          module procedure write_dataset1d_int
          module procedure write_dataset2d
          module procedure write_dataset2d_int
          module procedure write_dataset3d
        end interface

        interface read_dataset
          module procedure read_dataset1d
          module procedure read_dataset1d_int
          module procedure read_dataset2d
          module procedure read_dataset2d_int
          module procedure read_dataset3d
        end interface

   contains
       logical function restart_exists()
           ! this function checks for the existence of restart files and
           integer :: idproc, ierr, rc

           ! common block for parallel processing
           integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
           ! nproc = number of real or virtual processors obtained
           ! lgrp = current communicator
           ! lworld = MPI_COMM_WORLD communicator
           common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld

           ! variables
           integer :: fid
           character(len = 70) :: fname

           ! determine the rank of the calling process in the communicator
           call MPI_COMM_RANK(lworld, idproc, ierr)

           call start_hdf5(ierr)

           ! try opening the restart files to see if they exist
           if (idproc==0) then
               fname = 'RST/rstrtdat_f000000.h5'

               !first, check that the file is present
               inquire( file=fname, exist=restart_exists )

               if (restart_exists) then
                  !now, check that we can open it
                  call h5fopen_f(fname,H5F_ACC_RDONLY_F, fid, ierr)
                  if (ierr /= 0) then 
                     restart_exists = .false.
                  else       
                     call h5fclose_f(fid, ierr)
                     if (ierr /= 0) then
                        print*, 'restart_exists: could not close HDF5 file', idproc
                        call MPI_ABORT(lworld, rc, ierr)
                     endif         
                  endif
               endif
           endif

           call MPI_BCAST(restart_exists, 1, mint, 0, lworld, ierr)

           call stop_hdf5(ierr)

       end function restart_exists

       subroutine restart_write(itime, itime0, npp, part, movion, nppi, parti, sqn, nppiSQN, partiSQN, &
                     & nppeTest, parteTest, nppiTest, partiTest, itw, wt, tracks)
           ! write file for basic restart or continuation
           implicit none
           integer, intent(in) :: itime, itime0, movion
           real, dimension(:,:,:), pointer :: part, parti, partiSQN, parteTest, partiTest
           real, dimension(:,:) :: sqn
           integer, dimension(:), pointer :: npp, nppi, nppiSQN, nppeTest, nppiTest
           integer :: itw
           real, dimension(:,:), pointer :: wt
           type(t_track_set), optional :: tracks

           ! local data
           integer :: ntime
           integer :: ii, j, it
           integer nio, idproc, ioLeader, ntpg, nvp
           integer istatus, ierror, ierror2, rc
           integer, dimension(1) :: readyflag = 1
           integer, dimension(4) :: msid=0, istatuses=0

           integer, dimension(:), pointer :: tracks_npoints => null(), tracks_savedpoints => null(), tracks_part_idx => null()
           integer, dimension(:,:), pointer :: tracks_tags => null(), tracks_n => null()
           real, dimension(:,:,:), pointer :: tracks_data => null()
           integer :: max_npoints = 0

           ! common block for parallel processing
           integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
           ! lstat = length of status array
           parameter(lstat = 10)
           ! nproc = number of real or virtual processors obtained
           ! lgrp = current communicator
           ! lworld = MPI_COMM_WORLD communicator
           common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld

           !**** file stuff *******************************************
           ! variables
           integer(HID_T) :: file_id, rootID, tracks_gid
           
           integer(hsize_t), dimension(1) :: dims1
           integer(hsize_t), dimension(2) :: dims2, max_dims
           integer(hsize_t), dimension(3) :: dims3

           character(len = 70) :: fname, name, linkname
           character(len = 10) :: fnum ! the number of the file (num files=num I/O tasks)
           character(len = 70) :: tempstr
           character(len = 70) :: varstr ! the string describing the variable to write


           ! this segment is used for mpi computers without shared memory!

           ntime = itime + itime0

           ! determine the rank of the calling process in the communicator
           call MPI_COMM_RANK(lworld, idproc, ierror)
           ! determine the size of the group associated with a communicator
           call MPI_COMM_SIZE(lworld, nvp, ierror)

           call start_hdf5(ierror)

           ! determine the number of tasks doing I/O
           if (nvp <= 16) then
               nio = nvp
           else
               nio = max((nvp + 64 - 1)/64, 16)
               nio = min(nio,256)
           endif

           ! Now group nodes for I/O.  The member with the lowest rank will do I/O
           ! this system does not take into account network mapping
           ntpg = (nvp + nio - 1)/nio ! max number of tasks per group
           ioLeader = (idproc/ntpg) * ntpg

           ! get track data ready....which everybody needs to do
           if (present(tracks)) then
              allocate(tracks_npoints(tracks%ntracks),tracks_savedpoints(tracks%ntracks))
              allocate(tracks_part_idx(tracks%ntracks),tracks_tags(tracks%ntracks,2))
              allocate(tracks_n(tracks%ntracks,tracks%maxpoints))
              allocate(tracks_data(tracks%ntracks,size(tracks%tracks(1)%data,1),tracks%maxpoints))

              tracks_npoints = 0
              tracks_savedpoints = 0
              tracks_part_idx = 0
              tracks_tags = 0
              tracks_n = 0

              do ii = 1,tracks%ntracks
                 tracks_npoints(ii) = tracks%tracks(ii)%npoints
                 tracks_savedpoints(ii) = tracks%tracks(ii)%savedpoints
                 tracks_part_idx(ii) = tracks%tracks(ii)%part_idx
                 tracks_tags(ii,:) = tracks%tracks(ii)%tag(:)
                 tracks_n(ii,1:tracks_npoints(ii)) = tracks%tracks(ii)%n(1:tracks_npoints(ii))
                 tracks_data(ii,:,1:tracks_npoints(ii)) = tracks%tracks(ii)%data(:,1:tracks_npoints(ii))
              enddo
           endif

           ! open the files and write the leader tasks' data
           if (idproc == ioLeader) then
               write(fnum, '(i6.6)') idproc/ntpg
               write(tempstr, '(i6.6)') ntime
               name = 'rst-' // trim(adjustl(tempstr))
               name = trim(adjustl(name)) // '_f' // trim(adjustl(fnum))
               name = trim(adjustl(name))
               fname = 'RST/' // trim(adjustl(name)) // '.h5'
               fname = trim(adjustl(fname))

               call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, ierror)
               if (ierror /= 0) then
                   print*, 'restart_write: could not open HDF5 file for restart', idproc
                   !call MPI_ABORT(lworld, rc, ierror)
               endif

               call h5gopen_f( file_id, '/', rootID, ierror )
               if (ierror.ne.0) then 
                  print*, "restart_write: h5gopen_f encountered an error", idproc
                  !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif               

               call add_h5_atribute( rootID, 'name', name ) 
               call add_h5_atribute( rootID, 'itime', itime )
               call add_h5_atribute( rootID, 'itime0', itime0 ) 
               call add_h5_atribute( rootID, 'movion', movion ) 
               call add_h5_atribute( rootID, 'nvp', nvp ) 
               call add_h5_atribute( rootID, 'nfiles', nio ) 

               ! save file_write for tracks....do it here so we don't do it for every node
               if (present(tracks)) then
                  call add_h5_atribute( rootID, 'tracks_file', tracks%file_write ) 
               endif
     
               call h5gclose_f( rootID, ierror )
               if (ierror.ne.0) then 
                  print*, "restart_write: h5gclose_f encountered an error", idproc
                  !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif


               ! some stuff we only need to write out once
               if (idproc == 0) then
                   
                  call h5gopen_f( file_id, '/', rootID, ierror )
                  if (ierror.ne.0) then 
                     print*, "restart_write: Error in h5gopen_f", idproc
                     !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                  endif

                  !for wt
                  dims2(1) = itw
                  dims2(2) = 4
                  call write_dataset(rootID, dims2, wt(1:itw,:), "wt", 'itw', itw)
                  

                   ! SQN position, if we're using it
                   if (associated(nppiSQN)) then

                     !for sqn
                     dims2(1) = size(sqn,1)
                     dims2(2) = size(sqn,2)
                     call write_dataset(rootID, dims2, sqn, "sqn", 'no attribute', 0)

                   endif

                   call h5gclose_f(rootID, ierror)
                   if (ierror.ne.0) then 
                      print*, "restart_write: Error in h5gclose_f", idproc
                      !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                   endif
                   
               endif


               ! Now the fun begins.  We'll write out the ion data (if it exists) and the electron data
               ! for this node.  Then we'll receive the data from the other nodes in the group, overwriting
               ! the data for this node, and write that out.  When we're doine, we'll read in the data for
               ! this node.  It should save memory, but I don't know if it's better than chunking.

               ! deal with a problem I ran into with the Intel compiler
               if (npp(1)==0) part = 0
               if (nppeTest(1)==0) parteTest = 0

               if (movion /= 0) then
                  if (nppi(1)==0) part = 0
                  if (nppiTest(1)==0) partiTest = 0
               endif

               ! loop over the nodes in the group
               do ii = idproc, min(ntpg + idproc - 1, nvp - 1)

                  call h5gopen_f( file_id, '/', rootID, ierror )
                  if (ierror.ne.0) then 
                     print*, "restart_write: Error in h5gopen_f", idproc
                     !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                  endif


                   ! first, electron data
                   if (ii /= idproc) then
                       ! npp
                       call MPI_RECV(npp, size(npp), mint, ii, 99, lworld, istatus, ierror)
                       !part
                       ! send is necessary to stop buffer overflow on uBGL
                       !call MPI_SEND(readyflag, size(readyflag), mint, ii, 99, lworld, ierror)
                       call MPI_RECV(part, size(part, 1) * npp(1), mreal, ii, 99, lworld, istatus, ierror)
                   endif


                  !for part
                  write(tempstr, '(i10)') ii
                  varstr = 'part' // trim(adjustl(tempstr))

                  dims2(1) = size(part,1)
                  if (npp(1) == 0) then
                     dims2(2) = 1
                     part(1,1,1) = -1.0
                  else
                     dims2(2) = npp(1)
                  endif
                  call write_dataset(rootID, dims2, part(:,1:npp(1), 1), varstr, 'npp', npp(1))

                  ! test electrons
                  if (ii /= idproc) then
                       ! nppeTest
                       call MPI_RECV(nppeTest, size(nppeTest), mint, ii, 99, lworld, istatus, ierror)
                       !part
                       ! send is necessary to stop buffer overflow on uBGL
                       !call MPI_SEND(readyflag, size(readyflag), mint, ii, 99, lworld, ierror)
                       call MPI_RECV(parteTest, size(parteTest, 1) * nppeTest(1), mreal, ii, 99, lworld, istatus, ierror)
                   endif

                  !for parteTest
                  write(tempstr, '(i10)') ii
                  varstr = 'parteTest' // trim(adjustl(tempstr))

                  dims2(1) = size(parteTest,1)
                  if (nppeTest(1) == 0) then
                     ! Because HDF5 gets pissed if a dim is 0
                     dims2(2) = 1
                     parteTest(1,1,1) = -1.0
                  else
                     dims2(2) = nppeTest(1)
                  endif
                  call write_dataset(rootID, dims2, parteTest(:,1:nppeTest(1), 1), varstr, 'nppeTest', nppeTest(1))
                  

                   ! SQN data, if it exists
                   if (associated(nppiSQN)) then
                       ! nppiSQN
                       if (ii /= idproc) then
                           call MPI_RECV(nppiSQN, size(nppiSQN), mint, ii, 99, lworld, istatus, ierror)

                           !partiSQN
                           ! send is necessary to stop buffer overflow on uBGL
                           !call MPI_SEND(readyflag, size(readyflag), mint, ii, 99, lworld, ierror)
                           call MPI_RECV(partiSQN, size(partiSQN, 1) * nppiSQN(1), mreal, ii, 99, lworld, istatus, ierror)
                       endif

                      !for partiSQN
                      write(tempstr, '(i10)') ii
                      varstr = 'partiSQN' // trim(adjustl(tempstr))

                      dims2(1) = size(partiSQN,1)
                        if (nppiSQN(1) == 0) then
                           dims2(2) = 1
                           partiSQN(1,1,1) = -1.0
                        else
                           dims2(2) = nppiSQN(1)
                        endif
                      call write_dataset(rootID, dims2, partiSQN(:,1:nppiSQN(1), 1), varstr, 'nppiSQN', nppiSQN(1))
                      
                   endif

                   ! then, maybe ion data
                   if (movion /= 0) then
                       !nppi
                       if (ii /= idproc) then
                           call MPI_RECV(nppi, size(nppi), mint, ii, 99, lworld, istatus, ierror)

                           !parti
                           !call MPI_SEND(readyflag, size(readyflag), mint, ii, 99, lworld, ierror)
                           call MPI_RECV(parti, size(parti, 1) * nppi(1), mreal, ii, 99, lworld, istatus, ierror)
                       endif

                      !for parti
                      write(tempstr, '(i10)') ii
                      varstr = 'parti' // trim(adjustl(tempstr))

                      dims2(1) = size(parti,1)
                     if (nppi(1) == 0) then
                        dims2(2) = 1
                        parti(1,1,1) = -1.0
                     else
                        dims2(2) = nppi(1)
                     endif
                      call write_dataset(rootID, dims2, parti(:,1:nppi(1), 1), varstr, 'nppi', nppi(1))   


                       !nppiTest
                       if (ii /= idproc) then
                           call MPI_RECV(nppiTest, size(nppiTest), mint, ii, 99, lworld, istatus, ierror)

                           !parti
                           !call MPI_SEND(readyflag, size(readyflag), mint, ii, 99, lworld, ierror)
                           call MPI_RECV(partiTest, size(partiTest, 1) * nppiTest(1), mreal, ii, 99, lworld, istatus, ierror)
                       endif

                      !for partiTest
                      write(tempstr, '(i10)') ii
                      varstr = 'partiTest' // trim(adjustl(tempstr))

                      dims2(1) = size(partiTest,1)
                     if (nppiTest(1) == 0) then
                        dims2(2) = 1
                        partiTest(1,1,1) = -1.0
                     else
                        dims2(2) = nppiTest(1)
                     endif

                     call write_dataset(rootID, dims2, partiTest(:,1:nppiTest(1), 1), varstr, 'nppiTest', nppiTest(1))   
                      
                   endif

                   ! now tracks, if applicable
                   if (present(tracks)) then

                      ! create group for this node's tracks
                      write(tempstr, '(i10)') ii
                      varstr = 'tracks' // trim(adjustl(tempstr))

                      call h5gcreate_f(rootID, varstr, tracks_gid, ierror)
                      if (ierror.ne.0) then 
                         print*, "restart_write: Error in h5gcreate_f", idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif
                      
                      ! receive the crap
                      if (ii/=idproc) then 
                        call MPI_IRECV(tracks_npoints,tracks%ntracks,mint,ii,98,lworld,msid(1),ierror)
                        call MPI_IRECV(tracks_part_idx,tracks%ntracks,mint,ii,97,lworld,msid(2),ierror)
                        call MPI_IRECV(tracks_tags,size(tracks_tags),mint,ii,96,lworld,msid(3),ierror)

                        ! wait to receive the rest and start writing
                        call MPI_WAIT(msid(1),istatus,ierror)
                      endif

                     max_npoints = maxval(tracks_npoints) ! limit size of send

                     if (ii/=idproc) then 
                        call MPI_IRECV(tracks_n(:,1:max_npoints),size(tracks_n(:,1:max_npoints)),mint,ii,95,lworld,msid(1),ierror)
                        call MPI_IRECV(tracks_data(:,:,1:max_npoints),size(tracks_data(:,:,1:max_npoints)),mreal,ii,94,lworld,msid(4),ierror)
                     endif

                     ! write the crap

                     !npoints
                     dims1(1) = tracks%ntracks
                     varstr = 'npoints'
                     call write_dataset(tracks_gid, dims1, tracks_npoints, varstr, 'ntracks', tracks%ntracks)                      

                     !savedpoints if leader (rest don't need it)
                     if (ii==idproc) then
                        varstr = 'savedpoints'
                        call write_dataset(tracks_gid, dims1, tracks_savedpoints, varstr, 'ntracks', tracks%ntracks)
                     endif

                     ! part_idx
                     if (ii/=idproc) call MPI_WAIT(msid(2),istatus,ierror)
                     varstr = 'part_idx'
                     call write_dataset(tracks_gid, dims1, tracks_part_idx, varstr, 'ntracks', tracks%ntracks)

                     ! tags
                     if (ii/=idproc) call MPI_WAIT(msid(3),istatus,ierror)
                     varstr = 'tags'
                     dims2(1) = tracks%ntracks
                     dims2(2) = 2
                     call write_dataset(tracks_gid, dims2, tracks_tags, varstr, 'ntracks', tracks%ntracks)

                     ! write out rest of data only if anything relevant
                     if (max_npoints>0) then
                        !n
                        if (ii/=idproc) call MPI_WAIT(msid(1),istatus,ierror)
                        varstr = 'n'
                        dims2(1) = tracks%ntracks
                        dims2(2) = max_npoints
                        call write_dataset(tracks_gid, dims2, tracks_n(:,1:max_npoints), varstr, 'ntracks', tracks%ntracks)

                        !data
                        if (ii/=idproc) call MPI_WAIT(msid(4),istatus,ierror)
                        varstr = 'data'
                        dims3(1) = tracks%ntracks
                        dims3(2) = size(tracks_data,2)
                        dims3(3) = max_npoints
                        call write_dataset(tracks_gid, dims3, tracks_data(:,:,1:max_npoints), varstr, 'ntracks', tracks%ntracks)
                     elseif (ii/=idproc) then
                        call MPI_WAIT(msid(1),istatus,ierror)
                        call MPI_WAIT(msid(4),istatus,ierror)
                     endif

                      call h5gclose_f(tracks_gid, ierror)
                      if (ierror.ne.0) then 
                         print*, "restart_write: Error closing a tracks group", idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

                   endif

               enddo
   
               call h5gclose_f(rootID, ierror)
               if (ierror.ne.0) then 
                   print*, "restart_write: Error in h5gclose_f", idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif
               
               ! Read back in this task's data
               ! first, electron data
               ! open part
               write(tempstr, '(i10)') idproc
               varstr = '/part' // trim(adjustl(tempstr))

               call read_dataset(file_id, dims2, part(:,:,1), varstr)

               if (part(1,1,1)<0) then
                  npp(1) = 0
               else
                  npp(1) = dims2(2)
               endif


               ! test electron data
               ! open parteTest
               write(tempstr, '(i10)') idproc
               varstr = '/parteTest' // trim(adjustl(tempstr))

               call read_dataset(file_id, dims2, parteTest(:,:,1), varstr)

               if (parteTest(1,1,1)<0) then
                  nppeTest(1) = 0
               else
                  nppeTest(1) = dims2(2)
               endif

               ! SQN data, if present
               if (associated(nppiSQN)) then
                  write(tempstr, '(i10)') idproc
                  varstr = '/partiSQN' // trim(adjustl(tempstr))

                  call read_dataset(file_id, dims2, partiSQN(:,:,1), varstr)

                  if (partiSQN(1,1,1)<0) then
                     nppiSQN(1) = 0
                  else
                     nppiSQN(1) = dims2(2)
                  endif
               endif

               ! then, maybe ion data
               if (movion /= 0) then
                   write(tempstr, '(i10)') idproc
                   varstr = '/parti' // trim(adjustl(tempstr))
               
                   call read_dataset(file_id, dims2, parti(:,:,1), varstr)

                  if (parti(1,1,1)<0) then
                     nppi(1) = 0
                  else
                     nppi(1) = dims2(2)
                  endif

                   !test ions
                   write(tempstr, '(i10)') idproc
                   varstr = '/partiTest' // trim(adjustl(tempstr))
               
                   call read_dataset(file_id, dims2, partiTest(:,:,1), varstr)

                  if (partiTest(1,1,1)<0) then
                     nppiTest(1) = 0
                  else
                     nppiTest(1) = dims2(2)
                  endif

               endif

               ! close the file
               call h5fclose_f(file_id, ierror)
               if (ierror.ne.0) then 
                   print*, "restart_write: Error closing HDF file", idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif

               ! link to the file
               fname = trim(adjustl(name)) // '.h5'
               linkname = 'RST/' // 'rstrtdat' // '_f' // trim(adjustl(fnum)) // '.h5'
               call link( trim(fname), trim(linkname), ierror )
               if (ierror .ne. 0) then
	               print*, "error creating restart file link on node", idproc
	               call MPI_ABORT(lworld,ierror,ierror2)
               endif

           else ! send data to the leader node

               !print*, 'idproc=', idproc, 'tracks_part_idx=', tracks_part_idx

               call MPI_SEND(npp, size(npp), mint, ioLeader, 99, lworld, ierror)
               !call MPI_RECV(readyflag, size(readyflag), mint, ioLeader, 99, lworld, istatus, ierror)
               call MPI_SEND(part, size(part, 1) * npp(1), mreal, ioLeader, 99, lworld, ierror)

               call MPI_SEND(nppeTest, size(nppeTest), mint, ioLeader, 99, lworld, ierror)
               !call MPI_RECV(readyflag, size(readyflag), mint, ioLeader, 99, lworld, istatus, ierror)
               call MPI_SEND(parteTest, size(parteTest, 1) * nppeTest(1), mreal, ioLeader, 99, lworld, ierror)

               if (associated(nppiSQN)) then
                   call MPI_SEND(nppiSQN, size(nppiSQN), mint, ioLeader, 99, lworld, ierror)
                   !call MPI_RECV(readyflag, size(readyflag), mint, ioLeader, 99, lworld, istatus, ierror)
                   call MPI_SEND(partiSQN, size(partiSQN, 1) * nppiSQN(1), mreal, ioLeader, 99, lworld, ierror)
               endif

               if (movion /= 0) then
                   call MPI_SEND(nppi, size(nppi), mint, ioLeader, 99, lworld, ierror)
                   !call MPI_RECV(readyflag, size(readyflag), mint, ioLeader, 99, lworld, istatus, ierror)
                   call MPI_SEND(parti, size(parti, 1) * nppi(1), mreal, ioLeader, 99, lworld, ierror)

                   call MPI_SEND(nppiTest, size(nppiTest), mint, ioLeader, 99, lworld, ierror)
                   !call MPI_RECV(readyflag, size(readyflag), mint, ioLeader, 99, lworld, istatus, ierror)
                   call MPI_SEND(partiTest, size(partiTest, 1) * nppiTest(1), mreal, ioLeader, 99, lworld, ierror)
               endif

               ! tracks, if applicable
               if (present(tracks)) then

                  !send the stuff so far using non-blocking sends
                  
                  call MPI_ISEND(tracks_npoints,tracks%ntracks,mint,ioLeader,98,lworld,msid(1),ierror)
                  call MPI_ISEND(tracks_part_idx,tracks%ntracks,mint,ioLeader,97,lworld,msid(2),ierror)
                  call MPI_ISEND(tracks_tags,size(tracks_tags),mint,ioLeader,96,lworld,msid(3),ierror)
   
                  max_npoints = maxval(tracks_npoints) ! limit size of send

                  !if (max_npoints>0) then
                  !   print*, 'idproc=', idproc, 'npoints=', tracks%tracks(:)%npoints
                  !   print*, 'idproc=', idproc, 'track1,n=', tracks%tracks(1)%n(1:tracks%tracks(1)%npoints)
                  !   print*, 'idproc=', idproc, 'track2,n=', tracks%tracks(2)%n(1:tracks%tracks(2)%npoints)
                     
                  !   print*, 'idproc=', idproc, 'max_npoints=', max_npoints, 'itime=', itime
                  !   print*, 'idproc=', idproc, 'n=', tracks_n(:,1:max_npoints)
                  !endif

                  ! wait to send the rest...so the receiving node has info to be ready to receive (necessary to wait?)         
                  call MPI_WAIT(msid(1),istatus,ierror)

                  call MPI_ISEND(tracks_n(:,1:max_npoints),size(tracks_n(:,1:max_npoints)),mint,ioLeader,95,lworld,msid(1),ierror)
                  call MPI_ISEND(tracks_data(:,:,1:max_npoints),size(tracks_data(:,:,1:max_npoints)),mreal,ioLeader,94,lworld,msid(4),ierror)

                  ! wait to finish sending
                  !call MPI_WAITALL(msid,istatuses,ierror)
                  call MPI_WAIT(msid(2),istatus,ierror)
                  call MPI_WAIT(msid(3),istatus,ierror)
                  call MPI_WAIT(msid(1),istatus,ierror)
                  call MPI_WAIT(msid(4),istatus,ierror)

               endif

           endif

           call stop_hdf5(ierror)

           ! remove old restart
           call MPI_BARRIER(lworld, ierror)  !synchronise nodes (in case some of them fail to write/link)
           if (idproc == ioLeader) then
              fname = 'RST/' // trim(restart_name_old) // '.h5'
              call remove(trim(fname),ierror)
              restart_name_old = name
           endif

           ! deallocate track storage memory
           if (present(tracks)) then
              deallocate(tracks_npoints, tracks_savedpoints)
              deallocate(tracks_part_idx, tracks_tags)
              deallocate(tracks_n, tracks_data)
           endif

       end subroutine restart_write


       subroutine ion_write(idrun, nppi, parti, nppiTest, partiTest, itime, itime0)
           ! write out ions at the beginning of a run or when they're frozen
           ! ...even though I'm too lazy to deal with reading back in frozen ions
           implicit none
           integer, intent(in) :: idrun, itime, itime0
           real, dimension(:,:,:), pointer :: parti, partiTest
           integer, dimension(:), pointer :: nppi, nppiTest

           ! local data
           integer :: ii, j, it
           integer :: nio, idproc, ioLeader, ntpg, nvp
           integer :: istatus, ierror, ierror2, rc
           integer, dimension(1) :: readyflag = 1

           ! common block for parallel processing
           integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
           ! lstat = length of status array
           parameter(lstat = 10)
           ! nproc = number of real or virtual processors obtained
           ! lgrp = current communicator
           ! lworld = MPI_COMM_WORLD communicator
           common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld

           integer :: ntime

           !**** file stuff *******************************************
           ! variables
           integer(HID_T) :: file_id, rootID
           
           integer, parameter :: rank = 2 
           integer(hsize_t), dimension(rank) :: dims, max_dims

           character(len = 70) :: fname, name, linkname
           character(len = 10) :: fnum ! the number of the file (num files=num I/O tasks)
           character(len = 20) :: tempstr
           character(len = 70) :: varstr ! the string describing the variable to write

           ! this segment is used for mpi computers without shared memory!

           ntime = itime + itime0

           ! determine the rank of the calling process in the communicator
           call MPI_COMM_RANK(lworld, idproc, ierror)
           ! determine the size of the group associated with a communicator
           call MPI_COMM_SIZE(lworld, nvp, ierror)

           call start_hdf5(ierror)

           ! determine the number of tasks doing I/O
           if (nvp <= 16) then
               nio = nvp
           else
               nio = max((nvp + 64 - 1)/64, 16)
               nio = min(nio,256)
           endif

           ! Now group nodes for I/O.  The member with the lowest rank will do I/O
           ! this system does not take into account network mapping
           ntpg = (nvp + nio - 1)/nio ! max number of tasks per group
           ioLeader = (idproc/ntpg) * ntpg

           ! open the files and write the leader tasks data
           ! we have to use a separate file for each iteration because the program can crash
           ! during writes, making the files useless, and I haven't found out how to make backups in Fortran
           if (idproc == ioLeader) then
               write(fnum, '(i6.6)') idproc/ntpg
               write(tempstr, '(i6.6)') ntime
               name = 'ion-rst-' // trim(adjustl(tempstr))
               name = trim(adjustl(name)) // '_f' // trim(adjustl(fnum))
               name = trim(adjustl(name))
               fname = 'RST/' // trim(adjustl(name)) // '.h5'
               fname = trim(adjustl(fname))

               call h5fcreate_f(fname, H5F_ACC_TRUNC_F, file_id, ierror)
               if (ierror /= 0) then
                   print*, 'ion_write: could not open HDF5 file for restart', idproc
                   !call MPI_ABORT(lworld, rc, ierror)
               endif

               call h5gopen_f( file_id, '/', rootID, ierror )
               if (ierror.ne.0) then 
                  print*, "ion_write: h5gopen_f encountered an error", idproc
                  !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif               

               call add_h5_atribute( rootID, 'name', name ) 
               call add_h5_atribute( rootID, 'itime', itime )
               call add_h5_atribute( rootID, 'itime0', itime0 ) 
               call add_h5_atribute( rootID, 'nvp', nvp ) 
               call add_h5_atribute( rootID, 'nfiles', nio ) 
     
               call h5gclose_f( rootID, ierror )
               if (ierror.ne.0) then 
                  print*, "ion_write: h5gclose_f encountered an error", idproc
                  !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif


               ! deal with a problem I ran into with the Intel compiler
               if (nppi(1)==0) parti = 0
               if (nppiTest(1)==0) partiTest = 0
      

               ! Now the fun begins.  We'll write out the ion data.

               ! loop over the nodes in the group
               do ii = idproc, min(ntpg + idproc - 1, nvp - 1)
                  call h5gopen_f( file_id, '/', rootID, ierror )
                  if (ierror.ne.0) then 
                     print*, "ion_write: Error in h5gopen_f", idproc
                     !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                  endif

                   ! ion data
                   if (ii /= idproc) then
                       ! nppi
                       call MPI_RECV(nppi, size(nppi), mint, ii, 99, lworld, istatus, ierror)
                       !parti
                       ! send is necessary to stop buffer overflow on uBGL
                       !call MPI_SEND(readyflag, size(readyflag), mint, ii, 99, lworld, ierror)
                       call MPI_RECV(parti, size(parti, 1) * nppi(1), mreal, ii, 99, lworld, istatus, ierror)
                   endif


                  !for parti
                  write(tempstr, '(i10)') ii
                  varstr = 'parti' // trim(adjustl(tempstr))

                  dims(1) = size(parti,1)
                     if (nppi(1) == 0) then
                        dims(2) = 1
                        parti(1,1,1) = -1.0
                     else
                        dims(2) = nppi(1)
                     endif

                  call write_dataset(rootID, dims, parti(:,1:nppi(1), 1), varstr, 'nppi', nppi(1))


                   ! test ion data
                   if (ii /= idproc) then
                       ! nppi
                       call MPI_RECV(nppiTest, size(nppiTest), mint, ii, 99, lworld, istatus, ierror)
                       !parti
                       ! send is necessary to stop buffer overflow on uBGL
                       !call MPI_SEND(readyflag, size(readyflag), mint, ii, 99, lworld, ierror)
                       call MPI_RECV(partiTest, size(partiTest, 1) * nppiTest(1), mreal, ii, 99, lworld, istatus, ierror)
                   endif


                  !for parti
                  write(tempstr, '(i10)') ii
                  varstr = 'partiTest' // trim(adjustl(tempstr))

                  dims(1) = size(partiTest,1)
                  if (nppiTest(1) == 0) then
                     dims(2) = 1
                     partiTest(1,1,1) = -1.0
                  else
                     dims(2) = nppiTest(1)
                  endif

                  call write_dataset(rootID, dims, partiTest(:,1:nppiTest(1), 1), varstr, 'nppiTest', nppiTest(1))

               enddo

               call h5gclose_f(rootID, ierror)
               if (ierror.ne.0) then 
                   print*, "ion_write: Error in h5gclose_f", idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif

               ! Read back in this task's data (in case we're dumping initial electron positions)
               write(tempstr, '(i10)') idproc
               varstr = '/parti' // trim(adjustl(tempstr))

               call read_dataset(file_id, dims, parti(:,:,1), varstr)

                  if (parti(1,1,1)<0) then
                     nppi(1) = 0
                  else
                     nppi(1) = dims(2)
                  endif


               write(tempstr, '(i10)') idproc
               varstr = '/partiTest' // trim(adjustl(tempstr))

               call read_dataset(file_id, dims, partiTest(:,:,1), varstr)

               if (partiTest(1,1,1)<0) then
                  nppiTest(1) = 0
               else
                  nppiTest(1) = dims(2)
               endif


               ! close the file
               call h5fclose_f(file_id, ierror)
               if (ierror.ne.0) then 
                   print*, "ion_write: Error closing HDF file", idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif

               ! link to the file
               fname = trim(adjustl(name)) // '.h5'
               linkname = 'RST/' // 'ion-rstrtdat' // '_f' // trim(adjustl(fnum)) // '.h5'
               call link( trim(fname), trim(linkname), ierror )
               if (ierror .ne. 0) then
	               print*, "error creating restart file link on node", idproc
	               call MPI_ABORT(lworld,ierror,ierror2)
               endif

           else ! send data to the leader node
               call MPI_SEND(nppi, size(nppi), mint, ioLeader, 99, lworld, ierror)
               !call MPI_RECV(readyflag, size(readyflag), mint, ioLeader, 99, lworld, istatus, ierror)
               call MPI_SEND(parti, size(parti, 1) * nppi(1), mreal, ioLeader, 99, lworld, ierror)

               !test particles
               call MPI_SEND(nppiTest, size(nppiTest), mint, ioLeader, 99, lworld, ierror)
               !call MPI_RECV(readyflag, size(readyflag), mint, ioLeader, 99, lworld, istatus, ierror)
               call MPI_SEND(partiTest, size(partiTest, 1) * nppiTest(1), mreal, ioLeader, 99, lworld, ierror)
           endif

           call stop_hdf5(ierror)

       end subroutine ion_write

       subroutine restart_read(idrun, npp, part, movion, nppi, parti, sqn, nppiSQN, partiSQN, &
                     & nppeTest, parteTest, nppiTest, partiTest, smoothIon, itime, itime0, itw, &
                     & wt, tracks)

           ! read file for restart or continuation
           implicit none
           integer, intent(in) :: movion, idrun, smoothIon
           integer, intent(out), optional :: itime, itime0, itw
           real, dimension(:,:), pointer, optional :: wt
           real, dimension(:,:,:), pointer :: part, parti, partiSQN, parteTest, partiTest
           real, dimension(:,:) :: sqn
           integer, dimension(:), pointer :: npp, nppi, nppiSQN, nppeTest, nppiTest
           type(t_track_set), optional :: tracks

           ! local data
           integer :: ii
           integer :: restartNum
           integer :: nvp, nio, idproc, ioLeader, ntpg
           integer :: istatus, ierr, rc
           integer :: tempvar
           integer, dimension(1) :: readyflag = 1
           integer, dimension(4) :: msid=0, istatuses=0

           integer, dimension(:), pointer :: tracks_npoints => null(), tracks_savedpoints => null(), tracks_part_idx => null()
           integer, dimension(:,:), pointer :: tracks_tags => null(), tracks_n => null()
           real, dimension(:,:,:), pointer :: tracks_data => null()
           integer :: max_npoints = 0

           ! common block for parallel processing
           integer nproc, lgrp, lstat, mreal, mint, mcplx, mdouble, lworld
           ! lstat = length of status array
           parameter(lstat = 10)
           ! nproc = number of real or virtual processors obtained
           ! lgrp = current communicator
           ! lworld = MPI_COMM_WORLD communicator
           common /PPARMS/ nproc, lgrp, mreal, mint, mcplx, mdouble, lworld

           !**** HDF5 stuff *******************************************
           ! variables
           integer(HID_T) :: file_id, rootID, attr_id, typeID
           integer(size_t) :: strlen
           integer(hsize_t), dimension(1) :: attr_dims
           integer(hsize_t), dimension(1) :: dims1
           integer(hsize_t), dimension(2) :: dims2
           integer(hsize_t), dimension(3) :: dims3
           character(len = 70) :: fname, name
           character(len = 10) :: fnum ! the number of the file (num files=num I/O tasks)
           character(len = 20) :: tempstr
           character(len = 70) :: varstr ! the string describing the variable to write

           ! this segment is used for mpi computers without shared memory!

           ! determine the rank of the calling process in the communicator
           call MPI_COMM_RANK(lworld, idproc, ierr)
           ! determine the size of the group associated with a communicator
           call MPI_COMM_SIZE(lworld, nvp, ierr)

           ! allocate memory for track storage if necessary
           if (present(tracks)) then
              allocate(tracks_npoints(tracks%ntracks),tracks_savedpoints(tracks%ntracks))
              allocate(tracks_part_idx(tracks%ntracks),tracks_tags(tracks%ntracks,2))
              allocate(tracks_n(tracks%ntracks,tracks%maxpoints))
              allocate(tracks_data(tracks%ntracks,size(tracks%tracks(1)%data,1),tracks%maxpoints))

              tracks_npoints = 0
              tracks_savedpoints = 0
              tracks_part_idx = 0
              tracks_tags = 0
           endif

           call start_hdf5(ierr)

           ! determine the number of tasks doing I/O
           if (nvp <= 16) then
               nio = nvp
           else
               nio = max((nvp + 64 - 1)/64, 16)
               nio = min(nio,256)
           endif

           ! Now group nodes for I/O.  The member with the lowest rank will do I/O
           ! this system does not take into account network mapping
           ntpg = (nvp + nio - 1)/nio ! max number of tasks per group
           ioLeader = (idproc/ntpg) * ntpg

           ! determine most recent restart file
           if (idproc == 0) then
               print*, 'Begin reading restarts'

               write(fnum, '(i6.6)') idproc/ntpg
               name = 'rstrtdat_f' // trim(adjustl(fnum))
               fname = 'RST/' // trim(adjustl(name)) // '.h5'

               call h5fopen_f(fname,H5F_ACC_RDONLY_F, file_id, ierr)
               if (ierr .ne. 0) then
                   print*, 'restart_read: Error opening restart file', idproc
                   call MPI_ABORT(lworld, rc, ierr)
               endif

               ! read in itime from restart file
               attr_dims = 1
               if (present(itime)) then
                  call h5aopen_by_name_f(file_id, "/","itime",attr_id, ierr)
                  call h5aread_f(attr_id, H5T_NATIVE_INTEGER, itime, attr_dims, ierr) 
                  if (ierr .ne. 0) then
                      print*, 'restart_read: Error reading itime', idproc
                      call MPI_ABORT(lworld, rc, ierr)
                  endif 

                  call h5aclose_f(attr_id, ierr) 
                  if (ierr .ne. 0) then
                      print*, 'restart_read: Error closing attribute', idproc
                      call MPI_ABORT(lworld, rc, ierr)
                  endif 
               endif

               ! itime0
               if (present(itime0)) then
                  call h5aopen_by_name_f(file_id, "/","itime0",attr_id, ierr)
                  call h5aread_f(attr_id, H5T_NATIVE_INTEGER, itime0, attr_dims, ierr) 
                  if (ierr .ne. 0) then
                      print*, 'restart_read: Error reading itime0', idproc
                      call MPI_ABORT(lworld, rc, ierr)
                  endif 
                  
                  call h5aclose_f(attr_id, ierr) 
                  if (ierr .ne. 0) then
                      print*, 'restart_read: Error closing attribute', idproc
                      call MPI_ABORT(lworld, rc, ierr)
                  endif 
               endif

               ! nvp must agree or we have a problem
               call h5aopen_by_name_f(file_id, "/","nvp",attr_id, ierr)
               call h5aread_f(attr_id, H5T_NATIVE_INTEGER, tempvar, attr_dims, ierr) 
               if (ierr .ne. 0) then
                   print*, 'restart_read: Error reading nvp', idproc
                   call MPI_ABORT(lworld, rc, ierr)
               endif 
                
               call h5aclose_f(attr_id, ierr) 
               if (ierr .ne. 0) then
                   print*, 'restart_read: Error closing attribute', idproc
                   call MPI_ABORT(lworld, rc, ierr)
               endif 

               if (tempvar /= nvp) then
                   print*, 'restart_read: Changing the number of tasks between restarts is not supported.'
                   call MPI_ABORT(lworld, rc, ierr)
               endif

               ! movion must agree or we have a problem
               call h5aopen_by_name_f(file_id, "/","movion",attr_id, ierr)
               call h5aread_f(attr_id, H5T_NATIVE_INTEGER, tempvar, attr_dims, ierr) 
               if (ierr .ne. 0) then
                   print*, 'restart_read: Error reading movion', idproc
                   call MPI_ABORT(lworld, rc, ierr)
               endif 
                
               call h5aclose_f(attr_id, ierr) 
               if (ierr .ne. 0) then
                   print*, 'restart_read: Error closing attribute', idproc
                   call MPI_ABORT(lworld, rc, ierr)
               endif

               if (tempvar /= movion) then
                   print*, 'restart_read: Changin movion between restarts is not supoorted.'
                   print*, 'This error may be occuring because you froze the ions, which is not supported either.'
                   call MPI_ABORT(lworld, rc, ierr)
               endif

               ! energy history stuff
               if (present(itw)) then

                  call read_dataset(file_id, dims2, wt, "/wt")
                  itw = dims2(1)

               endif

               ! SQN position if applicable
               if (associated(nppiSQN)) then
                  call read_dataset(file_id, dims2, sqn, "/sqn")
               endif

               ! track file name if present
               if (present(tracks)) then
                  call h5aopen_by_name_f(file_id, "/","tracks_file",attr_id, ierr)
                  call h5tcopy_f(H5T_NATIVE_CHARACTER, typeID, ierr)
                  strlen = len(tracks%file_write)
                  call h5tset_size_f(typeID, strlen, ierr)
                  call h5aread_f(attr_id, typeID, tracks%file_write, attr_dims, ierr) 
                  if (ierr .ne. 0) then
                      print*, 'restart_read: Error reading tracks file name', idproc
                      call MPI_ABORT(lworld, rc, ierr)
                  endif
                      
                  call h5tclose_f( typeID, ierr )
                  if (ierr .ne. 0) then
                      print*, 'restart_read: Error closing attribute type', idproc
                      call MPI_ABORT(lworld, rc, ierr)
                  endif

                  call h5aclose_f(attr_id, ierr) 
                  if (ierr .ne. 0) then
                      print*, 'restart_read: Error closing attribute', idproc
                      call MPI_ABORT(lworld, rc, ierr)
                  endif

               endif


           endif

           ! broadcast info from task 0 to the other tasks
           ! so we can keep track in the new run if necessary
           if (present(itime0)) then
               call MPI_BCAST(itime0, 1, mint, 0, lworld, ierr)
           endif
           if (present(itime)) then
               call MPI_BCAST(itime, 1, mint, 0, lworld, ierr)
           endif
           if (present(itw)) then
               call MPI_BCAST(itw, 1, mint, 0, lworld, ierr)
               call MPI_BCAST(wt, itw*4, mreal, 0, lworld, ierr)
           endif
           if (associated(nppiSQN)) then
               call MPI_BCAST(sqn, size(sqn,1)*size(sqn,2), mreal, 0, lworld, ierr)
           endif
           if (present(tracks)) then
               call MPI_BCAST(tracks%file_write, len(tracks%file_write), MPI_CHARACTER, 0, lworld, ierr)
           endif

           ! leaders read in data and send it to other tasks
           if (idproc == ioLeader) then
               ! open files on rest of leader tasks
               if (idproc == 0) then
                   print*, 'reading electron data'
               else
                  write(fnum, '(i6.6)') idproc/ntpg
                  name = 'rstrtdat_f' // trim(adjustl(fnum))
                  fname = 'RST/' // trim(adjustl(name)) // '.h5'

                  call h5fopen_f(fname,H5F_ACC_RDONLY_F, file_id, ierr)
                  if (ierr .ne. 0) then
                      print*, 'restart_read: Error opening restart file', idproc
                      call MPI_ABORT(lworld, rc, ierr)
                  endif

               endif

               ! read in previous file name
               call h5aopen_by_name_f(file_id, "/","name",attr_id, ierr)
               call h5tcopy_f(H5T_NATIVE_CHARACTER, typeID, ierr)
               strlen = len(name)
               call h5tset_size_f(typeID, strlen, ierr)
               call h5aread_f(attr_id, typeID, name, attr_dims, ierr) 
               if (ierr .ne. 0) then
                   print*, 'restart_read: Error reading name', idproc
                   call MPI_ABORT(lworld, rc, ierr)
               endif
                   
               call h5tclose_f( typeID, ierr )
               if (ierr .ne. 0) then
                   print*, 'restart_read: Error closing attribute type', idproc
                   call MPI_ABORT(lworld, rc, ierr)
               endif

               call h5aclose_f(attr_id, ierr) 
               if (ierr .ne. 0) then
                   print*, 'restart_read: Error closing attribute', idproc
                   call MPI_ABORT(lworld, rc, ierr)
               endif

               restart_name_old = name

               ! loop over the nodes in the group backwards (so we have the data for this node at the end)
               do ii = min(ntpg + idproc - 1, nvp - 1), idproc, -1
                   ! first, electron data
                   write(tempstr, '(i10)') ii
                   varstr = '/part' // trim(adjustl(tempstr))

                   call read_dataset(file_id, dims2, part(:,:,1), trim(varstr))

                  if (part(1,1,1)<0) then
                     npp(1) = 0
                  else
                     npp(1) = dims2(2)
                  endif

                   if (ii /= idproc) then
                       call MPI_SEND(npp, size(npp), mint, ii, 99, lworld, ierr)
                       call MPI_SEND(part, size(part, 1) * npp(1), mreal, ii, 99, lworld, ierr)
                       ! necessary for some reason or things hang on bgl...
                       !call MPI_RECV(readyflag, size(readyflag), mint, ii, 99, lworld, istatus, ierr)
                   endif

                   ! test electron data
                   ! open parteTest
                   write(tempstr, '(i10)') ii
                   varstr = '/parteTest' // trim(adjustl(tempstr))

                   call read_dataset(file_id, dims2, parteTest(:,:,1), varstr)

                  if (parteTest(1,1,1)<0) then
                     nppeTest(1) = 0
                  else
                     nppeTest(1) = dims2(2)
                  endif

                   if (ii /= idproc) then
                       call MPI_SEND(nppeTest, size(nppeTest), mint, ii, 99, lworld, ierr)
                       call MPI_SEND(parteTest, size(parteTest, 1) * nppeTest(1), mreal, ii, 99, lworld, ierr)
                       ! necessary for some reason or things hang on bgl...
                       !call MPI_RECV(readyflag, size(readyflag), mint, ii, 99, lworld, istatus, ierr)
                   endif

                   ! nppiSQN, if necessary
                   if (associated(nppiSQN)) then
                      write(tempstr, '(i10)') ii
                      varstr = '/partiSQN' // trim(adjustl(tempstr))

                      call read_dataset(file_id, dims2, partiSQN(:,:,1), trim(varstr))
                     if (partiSQN(1,1,1)<0) then
                        nppiSQN(1) = 0
                     else
                        nppiSQN(1) = dims2(2)
                     endif

                      if (ii /= idproc) then
                          call MPI_SEND(nppiSQN, size(nppiSQN), mint, ii, 99, lworld, ierr)
                          call MPI_SEND(partiSQN, size(partiSQN, 1) * nppiSQN(1), mreal, ii, 99, lworld, ierr)
                          ! necessary for some reason or things hang...
                          !call MPI_RECV(readyflag, size(readyflag), mint, ii, 99, lworld, istatus, ierr)
                      endif
                   endif

                   ! now tracks, if applicable
                   if (present(tracks)) then

                      write(tempstr, '(i10)') ii
                      tempstr = '/tracks' // trim(adjustl(tempstr))

                      !npoints                      
                      varstr = trim(adjustl(tempstr)) // '/npoints'
                      call read_dataset(file_id, dims1, tracks_npoints, trim(varstr))
                      if (ii/=idproc) call MPI_ISEND(tracks_npoints,tracks%ntracks,mint,ii,98,lworld,msid(1),ierr)

                      !savedpoints
                      if (ii==idproc) then
                        varstr = trim(adjustl(tempstr)) // '/savedpoints'
                        call read_dataset(file_id, dims1, tracks_savedpoints, trim(varstr))
                      endif

                      !part_idx                     
                      varstr = trim(adjustl(tempstr)) // '/part_idx'
                      call read_dataset(file_id, dims1, tracks_part_idx, trim(varstr))
                      if (ii/=idproc) call MPI_ISEND(tracks_part_idx,tracks%ntracks,mint,ii,97,lworld,msid(2),ierr)

                      !tags
                      varstr = trim(adjustl(tempstr)) // '/tags'
                      call read_dataset(file_id, dims2, tracks_tags, trim(varstr))
                      if (ii/=idproc) call MPI_ISEND(tracks_tags,size(tracks_tags),mint,ii,96,lworld,msid(3),ierr)

                      ! wait to make sure receiver has gotten necessary data (is this necessary?)
                      if (ii/=idproc) call MPI_WAIT(msid(1),istatus,ierr)

                      max_npoints = maxval(tracks_npoints) ! limit size of send

                      
                     ! read rest of data only if anything relevant
                     if (max_npoints>0) then
                         !n
                         varstr = trim(adjustl(tempstr)) // '/n'
                         call read_dataset(file_id, dims2, tracks_n(:,1:max_npoints), trim(varstr))
                         if (ii/=idproc) call MPI_ISEND(tracks_n(:,1:max_npoints),size(tracks_n(:,1:max_npoints)),mint,ii,95,lworld,msid(1),ierr)

                         !data
                         varstr = trim(adjustl(tempstr)) // '/data'
                         call read_dataset(file_id, dims3, tracks_data(:,:,1:max_npoints), trim(varstr))
                         if (ii/=idproc) call MPI_ISEND(tracks_data(:,:,1:max_npoints),size(tracks_data(:,:,1:max_npoints)),mreal,ii,94,lworld,msid(4),ierr)
                         if (ii/=idproc) then
                           !call MPI_WAITALL(msid,istatuses,ierr)
                           call MPI_WAIT(msid(2),istatus,ierr)
                           call MPI_WAIT(msid(3),istatus,ierr)
                           call MPI_WAIT(msid(1),istatus,ierr)
                           call MPI_WAIT(msid(4),istatus,ierr)
                         endif
                     else
                        if (ii/=idproc) then
                            call MPI_WAIT(msid(2),istatus,ierr)
                            call MPI_WAIT(msid(3),istatus,ierr)
                        endif
                     endif
                      
                      ! put the data in the right places at the end...save code (it's OK since we loop backwards)

                   endif
               enddo

               ! read in ions if they're moving or we're not using a smooth neutralizing background
               if (movion==1 .or. (movion==0 .and. smoothIon==0)) then
                  ! open the ion restart file if necessary
                  if (movion==0 .and. smoothIon==0) then
                     ! close previous HDF file
                     call h5fclose_f(file_id, ierr)
                     if (ierr.ne.0) then 
                         print*, "restart_read: Error closing HDF file", idproc
                         call MPI_ABORT(lworld, rc, ierr)
                     endif

                     write(fnum, '(i6.6)') idproc/ntpg
                     name = 'ion-rstrtdat_f' // trim(adjustl(fnum))
                     fname = 'RST/' // trim(adjustl(name)) // '.h5'

                     call h5fopen_f(fname,H5F_ACC_RDONLY_F, file_id, ierr)
                     if (ierr .ne. 0) then
                         print*, 'restart_read: Error opening ion restart file', idproc
                         call MPI_ABORT(lworld, rc, ierr)
                     endif

                     ! nvp must agree or we have a problem
                     if (idproc==0) then
                        call h5aopen_by_name_f(file_id, "/","nvp",attr_id, ierr)
                        call h5aread_f(attr_id, H5T_NATIVE_INTEGER, tempvar, attr_dims, ierr) 
                        if (ierr .ne. 0) then
                            print*, 'restart_read: Error reading nvp', idproc
                            call MPI_ABORT(lworld, rc, ierr)
                        endif 
                         
                        call h5aclose_f(attr_id, ierr) 
                        if (ierr .ne. 0) then
                            print*, 'restart_read: Error closing attribute', idproc
                            call MPI_ABORT(lworld, rc, ierr)
                        endif 

                        if (tempvar /= nvp) then
                            print*, 'restart_read: Number of tasks in ion file does not match number in this run.'
                            call MPI_ABORT(lworld, rc, ierr)
                        endif
                     endif
                  endif

                   if (idproc==0) print*, 'reading ion data'

                   ! now, read in stuff
                   do ii = min(ntpg + idproc - 1, nvp - 1), idproc, -1
                      ! first, electron data
                      write(tempstr, '(i10)') ii
                      varstr = '/parti' // trim(adjustl(tempstr))

                     call read_dataset(file_id, dims2, parti(:,:,1), trim(varstr))
                     if (parti(1,1,1)<0) then
                        nppi(1) = 0
                     else
                        nppi(1) = dims2(2)
                     endif

                       if (ii /= idproc) then
                           call MPI_SEND(nppi, size(nppi), mint, ii, 99, lworld, ierr)
                           call MPI_SEND(parti, size(parti, 1) * nppi(1), mreal, ii, 99, lworld, ierr)
                           ! necessary for some reason or things hang...
                           !call MPI_RECV(readyflag, size(readyflag), mint, ii, 99, lworld, istatus, ierr)
                       endif

                       !test ions
                       write(tempstr, '(i10)') ii
                       varstr = '/partiTest' // trim(adjustl(tempstr))
                  
                       call read_dataset(file_id, dims2, partiTest(:,:,1), varstr)
                        if (partiTest(1,1,1)<0) then
                           nppiTest(1) = 0
                        else
                           nppiTest(1) = dims2(2)
                        endif

                       if (ii /= idproc) then
                           call MPI_SEND(nppiTest, size(nppiTest), mint, ii, 99, lworld, ierr)
                           call MPI_SEND(partiTest, size(partiTest, 1) * nppiTest(1), mreal, ii, 99, lworld, ierr)
                           ! necessary for some reason or things hang...
                           !call MPI_RECV(readyflag, size(readyflag), mint, ii, 99, lworld, istatus, ierr)
                       endif

                   enddo
               endif

               call h5fclose_f(file_id, ierr)
               if (ierr.ne.0) then 
                   print*, "restart_read: Error closing HDF file", idproc
                   call MPI_ABORT(lworld, rc, ierr)
               endif

           else ! receive the sends
               call MPI_RECV(npp, size(npp), mint, ioLeader, 99, lworld, istatus, ierr)
               call MPI_RECV(part, size(part,1) * npp(1), mreal, ioLeader, 99, lworld, istatus, ierr)
               ! necessary for some reason or things hang...
               !call MPI_SEND(readyflag, size(readyflag), mint, ioLeader, 99, lworld, ierr)

               call MPI_RECV(nppeTest, size(nppeTest), mint, ioLeader, 99, lworld, istatus, ierr)
               call MPI_RECV(parteTest, size(parteTest,1) * nppeTest(1), mreal, ioLeader, 99, lworld, istatus, ierr)
               ! necessary for some reason or things hang...
               !call MPI_SEND(readyflag, size(readyflag), mint, ioLeader, 99, lworld, ierr)

               if (associated(nppiSQN)) then
                   call MPI_RECV(nppiSQN, size(nppiSQN), mint, ioLeader, 99, lworld, istatus, ierr)
                   call MPI_RECV(partiSQN, size(partiSQN,1) * nppiSQN(1), mreal, ioLeader, 99, lworld, istatus, ierr)
                   ! necessary for some reason or things hang...
                   !call MPI_SEND(readyflag, size(readyflag), mint, ioLeader, 99, lworld, ierr)
               endif

               if (movion==1 .or. (movion==0 .and. smoothIon==0)) then
                   call MPI_RECV(nppi, size(nppi), mint, ioLeader, 99, lworld, istatus, ierr)
                   call MPI_RECV(parti, size(parti,1) * nppi(1), mreal, ioLeader, 99, lworld, istatus, ierr)
                   ! necessary for some reason or things hang...
                   !call MPI_SEND(readyflag, size(readyflag), mint, ioLeader, 99, lworld, ierr)

                   call MPI_RECV(nppiTest, size(nppiTest), mint, ioLeader, 99, lworld, istatus, ierr)
                   call MPI_RECV(partiTest, size(partiTest,1) * nppiTest(1), mreal, ioLeader, 99, lworld, istatus, ierr)
                   ! necessary for some reason or things hang...
                   !call MPI_SEND(readyflag, size(readyflag), mint, ioLeader, 99, lworld, ierr)
               endif

               if (present(tracks)) then
                  ! receive the crap
                  call MPI_IRECV(tracks_npoints,tracks%ntracks,mint,ioLeader,98,lworld,msid(1),ierr)
                  call MPI_IRECV(tracks_part_idx,tracks%ntracks,mint,ioLeader,97,lworld,msid(2),ierr)
                  call MPI_IRECV(tracks_tags,size(tracks_tags),mint,ioLeader,96,lworld,msid(3),ierr)

                  ! wait to receive the rest and start writing
                  call MPI_WAIT(msid(1),istatus,ierr)

                  max_npoints = maxval(tracks_npoints) ! limit size of send

                  if (max_npoints>0) then
                     !print*, 'idproc=', idproc, 'maxnpoints=', max_npoints
                     call MPI_IRECV(tracks_n(:,1:max_npoints),size(tracks_n(:,1:max_npoints)),mint,ioLeader,95,lworld,msid(1),ierr)
                     !print*, 'idproc=', idproc, 'after n', ' maxnpoints=', max_npoints
                     call MPI_IRECV(tracks_data(:,:,1:max_npoints),size(tracks_data(:,:,1:max_npoints)),mreal,ioLeader,94,lworld,msid(4),ierr)
                     !print*, 'idproc=', idproc, 'after data', ' maxnpoints=', max_npoints
                     !call MPI_WAITALL(msid,istatuses,ierr)
                     call MPI_WAIT(msid(2),istatus,ierr)
                     call MPI_WAIT(msid(3),istatus,ierr)
                     call MPI_WAIT(msid(1),istatus,ierr)
                     call MPI_WAIT(msid(4),istatus,ierr)

                     !print*, 'idproc=', idproc, 'max_npoints=', max_npoints, 'itime=', itime
                     !print*, 'idproc=', idproc, 'n=', tracks_n(:,1:max_npoints)

                  else
                     call MPI_WAIT(msid(2),istatus,ierr)
                     call MPI_WAIT(msid(3),istatus,ierr)
                  endif

               endif

           endif

           call stop_hdf5(ierr)

           ! put tracks data in structure, if applicable
           if (present(tracks)) then
              do ii = 1,tracks%ntracks
                 tracks%tracks(ii)%npoints = tracks_npoints(ii)
                 tracks%tracks(ii)%savedpoints = tracks_savedpoints(ii)  ! junk if not received
                 tracks%tracks(ii)%part_idx = tracks_part_idx(ii)
                 tracks%tracks(ii)%tag(:) = tracks_tags(ii,:)
                 tracks%tracks(ii)%n(1:tracks_npoints(ii)) = tracks_n(ii,1:tracks_npoints(ii))
                 tracks%tracks(ii)%data(:,1:tracks_npoints(ii)) = tracks_data(ii,:,1:tracks_npoints(ii))
              enddo

                 !if( max_npoints >0 ) then
                 !    print*, 'idproc=', idproc, 'npoints=', tracks%tracks(:)%npoints
                 !    print*, 'idproc=', idproc, 'track1,n=', tracks%tracks(1)%n(1:tracks%tracks(1)%npoints)
                 !    print*, 'idproc=', idproc, 'track2,n=', tracks%tracks(2)%n(1:tracks%tracks(2)%npoints)
                 !endif

           endif

       end subroutine restart_read

      subroutine write_dataset1d(rootID, dims, dataset, varstr, attrstr, attrval)

           character(len = *), intent(in) :: varstr, attrstr
           
           integer, parameter :: rank = 1
           integer(hsize_t), dimension(rank), intent(in) :: dims
           real, dimension(:) :: dataset
           integer, intent(in) :: attrval
           integer(HID_T), intent(in) :: rootID

           !local variables
            integer :: ierror, idproc
           integer(HID_T) :: dspace_id, dset_id, d_float

           ! determine the rank of the calling process in the communicator
           call MPI_COMM_RANK(MPI_COMM_WORLD, idproc, ierror)


            ! determine real type precision
            d_float = detect_precision()

                     call h5screate_simple_f(rank, dims, dspace_id, ierror) 
                      if (ierror.ne.0) then
                         print*, "write_dataset1d: Error creating dataspace for",trim(varstr), idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

                      call h5dcreate_f(rootID, trim(varstr), d_float, dspace_id, dset_id, ierror)
                      if (ierror.ne.0) then 
                         print*, "write_dataset1d: Error in h5dcreate_f", idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif
                      
                      ! write data and attributes
                      call h5dwrite_f(dset_id, d_float, dataset, dims, ierror)
                      if (ierror.ne.0) then 
                         print*, "write_dataset1d: Error writing",trim(varstr), idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

                      ! write nppi as an attribute of parti
                      call add_h5_atribute(dset_id,attrstr, attrval)

                      call h5dclose_f(dset_id, ierror)
                      if (ierror.ne.0) then 
                         print*, "write_dataset1d: Error closing dataset", idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

                      call h5sclose_f(dspace_id, ierror)
                      if (ierror.ne.0) then 
                         print*, "write_dataset1d: Error closing dataspace", idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

   end subroutine write_dataset1d

      subroutine write_dataset1d_int(rootID, dims, dataset, varstr, attrstr, attrval)

           character(len = *), intent(in) :: varstr, attrstr
           
           integer, parameter :: rank = 1
           integer(hsize_t), dimension(rank), intent(in) :: dims
           integer, dimension(:) :: dataset
           integer, intent(in) :: attrval
           integer(HID_T), intent(in) :: rootID

           !local variables
            integer :: ierror, idproc
           integer(HID_T) :: dspace_id, dset_id

           ! determine the rank of the calling process in the communicator
           call MPI_COMM_RANK(MPI_COMM_WORLD, idproc, ierror)

                     call h5screate_simple_f(rank, dims, dspace_id, ierror) 
                      if (ierror.ne.0) then
                         print*, "write_dataset1d: Error creating dataspace for",trim(varstr), idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

                      call h5dcreate_f(rootID, trim(varstr), H5T_NATIVE_INTEGER, dspace_id, dset_id, ierror)
                      if (ierror.ne.0) then 
                         print*, "write_dataset1d: Error in h5dcreate_f", idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif
                      
                      ! write data and attributes
                      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, dataset, dims, ierror)
                      if (ierror.ne.0) then 
                         print*, "write_dataset1d: Error writing",trim(varstr), idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

                      ! write nppi as an attribute of parti
                      call add_h5_atribute(dset_id,attrstr, attrval)

                      call h5dclose_f(dset_id, ierror)
                      if (ierror.ne.0) then 
                         print*, "write_dataset1d: Error closing dataset", idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

                      call h5sclose_f(dspace_id, ierror)
                      if (ierror.ne.0) then 
                         print*, "write_dataset1d: Error closing dataspace", idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

   end subroutine write_dataset1d_int


       subroutine write_dataset2d(rootID, dims, dataset, varstr, attrstr, attrval)

           character(len = *), intent(in) :: varstr, attrstr
           
           integer, parameter :: rank = 2 
           integer(hsize_t), dimension(rank), intent(in) :: dims
           real, dimension(:,:) :: dataset
           integer, intent(in) :: attrval
           integer(HID_T), intent(in) :: rootID

           !local variables
            integer :: ierror, idproc
           integer(HID_T) :: dspace_id, dset_id, d_float

           ! determine the rank of the calling process in the communicator
           call MPI_COMM_RANK(MPI_COMM_WORLD, idproc, ierror)


            ! determine real type precision
            d_float = detect_precision()

                     call h5screate_simple_f(rank, dims, dspace_id, ierror) 
                      if (ierror.ne.0) then
                         print*, "write_dataset2d: Error creating dataspace for",trim(varstr), idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

                      call h5dcreate_f(rootID, trim(varstr), d_float, dspace_id, dset_id, ierror)
                      if (ierror.ne.0) then 
                         print*, "write_dataset2d: Error in h5dcreate_f", idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif
                      
                      ! write data and attributes
                      call h5dwrite_f(dset_id, d_float, dataset, dims, ierror)
                      if (ierror.ne.0) then 
                         print*, "write_dataset2d: Error writing",trim(varstr), idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

                      ! write nppi as an attribute of parti
                      call add_h5_atribute(dset_id,attrstr, attrval)

                      call h5dclose_f(dset_id, ierror)
                      if (ierror.ne.0) then 
                         print*, "write_dataset2d: Error closing dataset", idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

                      call h5sclose_f(dspace_id, ierror)
                      if (ierror.ne.0) then 
                         print*, "write_dataset2d: Error closing dataspace", idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

   end subroutine write_dataset2d

       subroutine write_dataset2d_int(rootID, dims, dataset, varstr, attrstr, attrval)

           character(len = *), intent(in) :: varstr, attrstr
           
           integer, parameter :: rank = 2 
           integer(hsize_t), dimension(rank), intent(in) :: dims
           integer, dimension(:,:) :: dataset
           integer, intent(in) :: attrval
           integer(HID_T), intent(in) :: rootID

           !local variables
            integer :: ierror, idproc
           integer(HID_T) :: dspace_id, dset_id

           ! determine the rank of the calling process in the communicator
           call MPI_COMM_RANK(MPI_COMM_WORLD, idproc, ierror)

                     call h5screate_simple_f(rank, dims, dspace_id, ierror) 
                      if (ierror.ne.0) then
                         print*, "write_dataset2d: Error creating dataspace for",trim(varstr), idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

                      call h5dcreate_f(rootID, trim(varstr), H5T_NATIVE_INTEGER, dspace_id, dset_id, ierror)
                      if (ierror.ne.0) then 
                         print*, "write_dataset2d: Error in h5dcreate_f", idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif
                      
                      ! write data and attributes
                      call h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, dataset, dims, ierror)
                      if (ierror.ne.0) then 
                         print*, "write_dataset2d: Error writing",trim(varstr), idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

                      ! write nppi as an attribute of parti
                      call add_h5_atribute(dset_id,attrstr, attrval)

                      call h5dclose_f(dset_id, ierror)
                      if (ierror.ne.0) then 
                         print*, "write_dataset2d: Error closing dataset", idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

                      call h5sclose_f(dspace_id, ierror)
                      if (ierror.ne.0) then 
                         print*, "write_dataset2d: Error closing dataspace", idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

   end subroutine write_dataset2d_int

       subroutine write_dataset3d(rootID, dims, dataset, varstr, attrstr, attrval)

           character(len = *), intent(in) :: varstr, attrstr
           
           integer, parameter :: rank = 3
           integer(hsize_t), dimension(rank), intent(in) :: dims
           real, dimension(:,:,:) :: dataset
           integer, intent(in) :: attrval
           integer(HID_T), intent(in) :: rootID

           !local variables
            integer :: ierror, idproc
           integer(HID_T) :: dspace_id, dset_id, d_float

           ! determine the rank of the calling process in the communicator
           call MPI_COMM_RANK(MPI_COMM_WORLD, idproc, ierror)


            ! determine real type precision
            d_float = detect_precision()

                     call h5screate_simple_f(rank, dims, dspace_id, ierror) 
                      if (ierror.ne.0) then
                         print*, "write_dataset3d: Error creating dataspace for",trim(varstr), idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

                      call h5dcreate_f(rootID, trim(varstr), d_float, dspace_id, dset_id, ierror)
                      if (ierror.ne.0) then 
                         print*, "write_dataset3d: Error in h5dcreate_f", idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif
                      
                      ! write data and attributes
                      call h5dwrite_f(dset_id, d_float, dataset, dims, ierror)
                      if (ierror.ne.0) then 
                         print*, "write_dataset3d: Error writing",trim(varstr), idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

                      ! write nppi as an attribute of parti
                      call add_h5_atribute(dset_id,attrstr, attrval)

                      call h5dclose_f(dset_id, ierror)
                      if (ierror.ne.0) then 
                         print*, "write_dataset3d: Error closing dataset", idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

                      call h5sclose_f(dspace_id, ierror)
                      if (ierror.ne.0) then 
                         print*, "write_dataset3d: Error closing dataspace", idproc
                         !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
                      endif

   end subroutine write_dataset3d

   subroutine read_dataset1d(file_id, dims, dataset, varstr)
           character(len = *), intent(in) :: varstr
           
           integer, parameter :: rank = 1 
           integer(hsize_t), dimension(rank), intent(out) :: dims
           real, dimension(:) :: dataset
           integer(HID_T), intent(in) :: file_id


           !local variables
            integer :: ierror, idproc, rankr
           integer(hsize_t), dimension(rank) :: max_dims

           integer(HID_T) :: dspace_id, dset_id, d_float

           ! determine the rank of the calling process in the communicator
           call MPI_COMM_RANK(MPI_COMM_WORLD, idproc, ierror)

            ! determine real type precision
            d_float = detect_precision()

               call h5dopen_f(file_id, trim(varstr), dset_id, ierror, H5P_DEFAULT_F)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not open dataset to read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif
   
               ! get dataset properties       
               call H5Dget_space_f(dset_id, dspace_id, ierror)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not open dataspace to read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif    

               call H5Sget_simple_extent_ndims_f(dspace_id, rankr, ierror)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not get rank of dataset to read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif  

               ! the dimensions include the npp info
               call H5Sget_simple_extent_dims_f(dspace_id, dims, max_dims, ierror)
               if (ierror .ne. rankr) then
                   print*, 'read_dataset: Could not get dims of dataset to read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif  

               ! read in part
               call H5Dread_f(dset_id, d_float, dataset(1:dims(1)), dims, ierror)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif  

               call h5dclose_f(dset_id, ierror)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not close dataset',trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif  


               call h5sclose_f(dspace_id, ierror)
               if (ierror.ne.0) then 
                    print*, "read_dataset: Error closing dataspace",trim(varstr), idproc
                    !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif
   end subroutine read_dataset1d

   subroutine read_dataset1d_int(file_id, dims, dataset, varstr)
           character(len = *), intent(in) :: varstr
           
           integer, parameter :: rank = 1 
           integer(hsize_t), dimension(rank), intent(out) :: dims
           integer, dimension(:) :: dataset
           integer(HID_T), intent(in) :: file_id


           !local variables
            integer :: ierror, idproc, rankr
           integer(hsize_t), dimension(rank) :: max_dims

           integer(HID_T) :: dspace_id, dset_id

           ! determine the rank of the calling process in the communicator
           call MPI_COMM_RANK(MPI_COMM_WORLD, idproc, ierror)

               call h5dopen_f(file_id, trim(varstr), dset_id, ierror, H5P_DEFAULT_F)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not open dataset to read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif
   
               ! get dataset properties       
               call H5Dget_space_f(dset_id, dspace_id, ierror)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not open dataspace to read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif    

               call H5Sget_simple_extent_ndims_f(dspace_id, rankr, ierror)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not get rank of dataset to read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif  

               ! the dimensions include the npp info
               call H5Sget_simple_extent_dims_f(dspace_id, dims, max_dims, ierror)
               if (ierror .ne. rankr) then
                   print*, 'read_dataset: Could not get dims of dataset to read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif  

               ! read in part
               call H5Dread_f(dset_id, H5T_NATIVE_INTEGER, dataset(1:dims(1)), dims, ierror)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif  

               call h5dclose_f(dset_id, ierror)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not close dataset',trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif  


               call h5sclose_f(dspace_id, ierror)
               if (ierror.ne.0) then 
                    print*, "read_dataset: Error closing dataspace",trim(varstr), idproc
                    !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif
   end subroutine read_dataset1d_int

   subroutine read_dataset2d(file_id, dims, dataset, varstr)
           character(len = *), intent(in) :: varstr
           
           integer, parameter :: rank = 2 
           integer(hsize_t), dimension(rank), intent(out) :: dims
           real, dimension(:,:) :: dataset
           integer(HID_T), intent(in) :: file_id


           !local variables
            integer :: ierror, idproc, rankr
           integer(hsize_t), dimension(rank) :: max_dims

           integer(HID_T) :: dspace_id, dset_id, d_float

           ! determine the rank of the calling process in the communicator
           call MPI_COMM_RANK(MPI_COMM_WORLD, idproc, ierror)

            ! determine real type precision
            d_float = detect_precision()

               call h5dopen_f(file_id, trim(varstr), dset_id, ierror, H5P_DEFAULT_F)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not open dataset to read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif
   
               ! get dataset properties       
               call H5Dget_space_f(dset_id, dspace_id, ierror)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not open dataspace to read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif    

               call H5Sget_simple_extent_ndims_f(dspace_id, rankr, ierror)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not get rank of dataset to read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif  

               ! the dimensions include the npp info
               call H5Sget_simple_extent_dims_f(dspace_id, dims, max_dims, ierror)
               if (ierror .ne. rankr) then
                   print*, 'read_dataset: Could not get dims of dataset to read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif  

               ! read in part
               call H5Dread_f(dset_id, d_float, dataset(1:dims(1),1:dims(2)), dims, ierror)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif  

               call h5dclose_f(dset_id, ierror)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not close dataset',trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif  


               call h5sclose_f(dspace_id, ierror)
               if (ierror.ne.0) then 
                    print*, "read_dataset: Error closing dataspace",trim(varstr), idproc
                    !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif
   end subroutine read_dataset2d

   subroutine read_dataset2d_int(file_id, dims, dataset, varstr)
           character(len = *), intent(in) :: varstr
           
           integer, parameter :: rank = 2 
           integer(hsize_t), dimension(rank), intent(out) :: dims
           integer, dimension(:,:) :: dataset
           integer(HID_T), intent(in) :: file_id


           !local variables
            integer :: ierror, idproc, rankr
           integer(hsize_t), dimension(rank) :: max_dims

           integer(HID_T) :: dspace_id, dset_id

           ! determine the rank of the calling process in the communicator
           call MPI_COMM_RANK(MPI_COMM_WORLD, idproc, ierror)

               call h5dopen_f(file_id, trim(varstr), dset_id, ierror, H5P_DEFAULT_F)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not open dataset to read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif
   
               ! get dataset properties       
               call H5Dget_space_f(dset_id, dspace_id, ierror)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not open dataspace to read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif    

               call H5Sget_simple_extent_ndims_f(dspace_id, rankr, ierror)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not get rank of dataset to read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif  

               ! the dimensions include the npp info
               call H5Sget_simple_extent_dims_f(dspace_id, dims, max_dims, ierror)
               if (ierror .ne. rankr) then
                   print*, 'read_dataset: Could not get dims of dataset to read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif  

               ! read in part
               call H5Dread_f(dset_id, H5T_NATIVE_INTEGER, dataset(1:dims(1),1:dims(2)), dims, ierror)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif  

               call h5dclose_f(dset_id, ierror)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not close dataset',trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif  


               call h5sclose_f(dspace_id, ierror)
               if (ierror.ne.0) then 
                    print*, "read_dataset: Error closing dataspace",trim(varstr), idproc
                    !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif
   end subroutine read_dataset2d_int

   subroutine read_dataset3d(file_id, dims, dataset, varstr)
           character(len = *), intent(in) :: varstr
           
           integer, parameter :: rank = 3 
           integer(hsize_t), dimension(rank), intent(out) :: dims
           real, dimension(:,:,:) :: dataset
           integer(HID_T), intent(in) :: file_id


           !local variables
            integer :: ierror, idproc, rankr
           integer(hsize_t), dimension(rank) :: max_dims

           integer(HID_T) :: dspace_id, dset_id, d_float

           ! determine the rank of the calling process in the communicator
           call MPI_COMM_RANK(MPI_COMM_WORLD, idproc, ierror)

            ! determine real type precision
            d_float = detect_precision()

               call h5dopen_f(file_id, trim(varstr), dset_id, ierror, H5P_DEFAULT_F)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not open dataset to read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif
   
               ! get dataset properties       
               call H5Dget_space_f(dset_id, dspace_id, ierror)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not open dataspace to read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif    

               call H5Sget_simple_extent_ndims_f(dspace_id, rankr, ierror)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not get rank of dataset to read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif  

               ! the dimensions include the npp info
               call H5Sget_simple_extent_dims_f(dspace_id, dims, max_dims, ierror)
               if (ierror .ne. rankr) then
                   print*, 'read_dataset: Could not get dims of dataset to read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif  

               ! read in part
               call H5Dread_f(dset_id, d_float, dataset(1:dims(1),1:dims(2),1:dims(3)), dims, ierror)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not read back in', trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif  

               call h5dclose_f(dset_id, ierror)
               if (ierror .ne. 0) then
                   print*, 'read_dataset: Could not close dataset',trim(varstr), idproc
                   !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif  


               call h5sclose_f(dspace_id, ierror)
               if (ierror.ne.0) then 
                    print*, "read_dataset: Error closing dataspace",trim(varstr), idproc
                    !call MPI_ABORT(MPI_COMM_WORLD,ierror,ierror2)
               endif
   end subroutine read_dataset3d

   subroutine initmodediag(dent, ntd, id0, nxh, nyh, nzh, kxyp, kyzp, mode&
       &sxd, modesyd, modeszd, jblok, mblok, iud, ndrec, fdname)
       ! initialize mode diagnostic
       implicit none
       integer :: ntd, id0, nxh, nyh, nzh, kxyp, kyzp
       integer :: modesxd, modesyd, modeszd, jblok, mblok, iud, ndrec
       character(len = *) :: fdname
       complex, dimension(:,:,:,:), pointer :: dent
       ! local data
       integer :: modesz2d
       if (ntd <= 0) return
       if (modesxd > nxh) modesxd = nxh
       if (modesyd > nyh) modesyd = nyh
       if (modeszd > nzh) modeszd = nzh
       modesz2d = 2 * modeszd - 1
       allocate(dent(modesz2d, min(modesxd, kxyp), min(modesyd, kyzp), jblo&
       &k * mblok))
       ! open output file
       if (id0 == 0) then
           if (ndrec == 0) then
               iud = get_funit(iud); ndrec = -1
               call bfopen(dent, modesz2d, iud, ndrec, trim(fdname))
           endif
       else
           if (ndrec == 0) ndrec = 1
       endif
   end subroutine initmodediag
   !
   subroutine initveldiag(fv, fvm, vtx, vty, vtz, ntv, ndv, id0, nmv, nblok&
       &, iuv, fvname)
       ! initialize velocity diagnostic
       implicit none
       integer :: ntv, ndv, id0, nmv, nblok, iuv
       real :: vtx, vty, vtz
       real, dimension(:,:,:), pointer :: fv, fvm
       character(len = *) :: fvname
       if ((ntv <= 0) .and. (ndv <= 0)) return
       allocate(fv(2 * nmv + 2, 3, nblok), fvm(3, 3, nblok))
       ! fix velocity range
       fv(1,:,:) = 8. * max(vtx, vty, vtz)
       if (ntv > 0) then
           if (id0 == 0) then
               iuv = get_funit(iuv)
               open(unit = iuv, file = trim(fvname), form = 'formatted', status = '&
               &unknown')
               ! write captions
               write (iuv, *) 'it vdx vdy vdz vtx vty vtz sk'
           endif
       endif
   end subroutine initveldiag
   !
   subroutine dendiag(qt, qs, qi, sfield, dent, sfieldt, ffc, nyzp, mixup, &
       &sct, tfft, ntd, ndd, nx, ny, nz, modesxd, modesyd, modeszd, iud, ndrec, indx, i&
       &ndy, indz, ntime, nvpy, nvpz, kstrt, kxyp, kyp, kyzp, kzp, ngds, iblok, jblok, &
       &kblok, mblok, inorder)
       ! ion density diagnostic
       implicit none
       integer :: ntd, ndd, nx, ny, nz, modesxd, modesyd, modeszd
       integer :: iud, ndrec, indx, indy, indz, ntime, nvpy, nvpz
       integer :: kstrt, kxyp, kyp, kyzp, kzp, ngds
       integer :: iblok, jblok, kblok, mblok
       real, dimension(2) :: tfft
       integer, optional :: inorder
       real, dimension(:,:,:,:), pointer :: qi, sfield
       complex, dimension(:,:,:,:), pointer :: qt, qs, dent, sfieldt
       complex, dimension(:,:,:,:), pointer :: ffc
       integer, dimension(:,:), pointer :: nyzp
       integer, dimension(:), pointer :: mixup
       complex, dimension(:), pointer :: sct
       ! local data
       integer :: it, jt, modesz2d, isign
   !    if ((ntd > 0) .or. (ndd > 0)) then
           it = -1; if (ntd > 0) it = ntime - ntd * (ntime/ntd)
           jt = -1; if (ndd > 0) jt = ntime - ndd * (ntime/ndd)
   !        if ((it == 0) .or. (jt == 0)) then
               sfield = qi
               ! add guard cells for ion density in x
   !            call aguard(sfield, nyzp, nx, inorder)
   !            ! add guard cells for ion density in y and z
   !            call paguard(sfield, kstrt, nvpy, nvpz, nx, kyp, kzp, ngds, iblok&
   !            &, inorder)
               ! transform ion density to fourier space
               isign = -1
               call fft(sfield, qs, qt, isign, mixup, sct, tfft, indx, indy, indz&
               &, kstrt, kxyp, kyp, kyzp, kzp, kblok, mblok, inorder)
   !             calculate smoothing in fourier space
               call pois(qt, sfieldt, ffc, nx, ny, nz, kstrt, jblok)
               ! store selected fourier modes
   !            if (it == 0) then
   !                modesz2d = 2 * modeszd - 1
   !                call gtmodes(sfieldt, dent, nx, ny, nz, modesxd, modesyd, mod&
   !                &eszd, kstrt, jblok)
   !                ! write diagnostic output
   !                call writebf(dent, modesxd, modesyd, modesz2d, kxyp, kyzp, j&
   !                &blok, iud, ndrec)
   !            endif
               ! transform ion density to real space
   !            if (jt == 0) then
                   isign = 1
                   call fft(sfield, qs, sfieldt, isign, mixup, sct, tfft, indx, i&
                   &ndy, indz, kstrt, kxyp, kyp, kyzp, kzp, kblok, mblok, inorder)
                   call pcguard(sfield, kstrt, nvpy, nvpz, kyp, kzp, ngds, iblok&
                   &, inorder)
                   call cguard(sfield, nyzp, nx, inorder)
   !            endif
   !        endif
   !    endif
   end subroutine dendiag
   !
   subroutine veldiag(part, fv, fvm, npp, ntv, ndv, id0, nmv, iuv, ntime)
       ! velocity diagnostic
       implicit none
       integer :: ntv, ndv, id0, nmv, iuv, ntime
       real, dimension(:,:,:), pointer :: part
       real, dimension(:,:,:), pointer :: fv, fvm
       integer, dimension(:), pointer :: npp
       ! local data
       integer :: it, jt, nt
       if ((ntv > 0) .or. (ndv > 0)) then
           it = -1; if (ntv > 0) it = ntime - ntv * (ntime/ntv)
           jt = -1; if (ndv > 0) jt = ntime - ndv * (ntime/ndv)
           if ((it == 0) .or. (jt == 0)) then
               ! calculate particle distribution function and moments
               call vdist(part, fv, fvm, npp, nmv)
               call plsum(fv(:,:, 1))
               ! print out velocity moments
               if (it == 0) then
                   if (id0 == 0) then
                       nt = ntime/ntv
                       write (iuv, *) nt, fvm(1,:, 1), fvm(2,:, 1), sum(fvm(3&
                       &,:, 1))
                   endif
               endif
           endif
       endif
   end subroutine veldiag
   !
   subroutine potdiag(qt, qs, sfield, pott, sfieldt, ffc, nyzp, mixup, sct&
       &, tfft, ntp, ndp, nx, ny, nz, modesxp, modesyp, modeszp, iup, nprec, indx, indy&
       &, indz, ntime, nvpy, nvpz, kstrt, kxyp, kyp, kyzp, kzp, ngds, iblok, jblok, kbl&
       &ok, mblok, inorder)
       ! potential diagnostic
       implicit none
       integer :: ntp, ndp, nx, ny, nz, modesxp, modesyp, modeszp
       integer :: iup, nprec, indx, indy, indz, ntime, nvpy, nvpz
       integer :: kstrt, kxyp, kyp, kyzp, kzp, ngds
       integer :: iblok, jblok, kblok, mblok
       real, dimension(2) :: tfft
       integer, optional :: inorder
       real, dimension(:,:,:,:), pointer :: sfield
       complex, dimension(:,:,:,:), pointer :: qt, qs, pott, sfieldt
       complex, dimension(:,:,:,:), pointer :: ffc
       integer, dimension(:,:), pointer :: nyzp
       integer, dimension(:), pointer :: mixup
       complex, dimension(:), pointer :: sct
       ! local data
       integer :: it, jt, modesz2p, isign
       real :: we
       if ((ntp > 0) .or. (ndp > 0)) then
           it = -1; if (ntp > 0) it = ntime - ntp * (ntime/ntp)
           jt = -1; if (ndp > 0) jt = ntime - ndp * (ntime/ndp)
           if ((it == 0) .or. (jt == 0)) then
               ! calculate potential in fourier space
               call pois(qt, sfieldt, ffc, we, nx, ny, nz, kstrt, jblok)
               ! store selected fourier modes
               if (it == 0) then
                   modesz2p = 2 * modeszp - 1
                   call gtmodes(sfieldt, pott, nx, ny, nz, modesxp, modesyp, mod&
                   &eszp, kstrt, jblok)
                   ! write diagnostic output
                   !      keeps running into div by 0 errors
                   !                   call writebf(pott,modesxp,modesyp,modesz2p,kxyp,kyzp,j&
                   !      &blok,iup,nprec)
               endif
               ! transform potential to real space
               if (jt == 0) then
                   isign = 1
                   call fft(sfield, qs, sfieldt, isign, mixup, sct, tfft, indx, i&
                   &ndy, indz, kstrt, kxyp, kyp, kyzp, kzp, kblok, mblok, inorder)
                   call pcguard(sfield, kstrt, nvpy, nvpz, kyp, kzp, ngds, iblok&
                   &, inorder)
                   call cguard(sfield, nyzp, nx, inorder)
               endif
           endif
       endif
   end subroutine potdiag
   !
   subroutine esenergy(wt, wtot, msg, we, wke, wki, ntw, ndw, id0, itw, iuot&
       &, ntime)
       ! electrostatic energy diagnostic
       implicit none
       integer :: ntw, ndw, id0, itw, iuot, ntime
       real :: we, wke, wki
       real, dimension(:,:), pointer :: wt
       real, dimension(4) :: wtot
       double precision, dimension(:) :: msg
       ! local data
       integer :: it, jt
       992 format (' field, kinetic, total energies = ', 3e14.7)
       if ((ntw > 0) .or. (ndw > 0)) then
           it = -1; if (ntw > 0) it = ntime - ntw * ((ntime + 1)/ntw) + 1
           jt = -1; if (ndw > 0) jt = ntime - ndw * ((ntime + 1)/ndw) + 1
           if ((it == 0) .or. (jt == 0)) then
               wtot(1) = we
               wtot(2) = wke
               wtot(3) = wki
               wtot(4) = we + wke + wki
               call plsum(wtot)
               ! send energy values to diagnostic node
               msg(1:4) = wtot
               call HARTBEAT(msg, 4)
               if (it == 0) then
                   if (id0 == 0) write (iuot, 992) wtot(1), wtot(2), wtot(4)
               endif
               if (jt == 0) then
                   itw = itw + 1
                   wt(itw,:) = wtot
               endif
           endif
       endif
   end subroutine esenergy
   !
   end module psimul32d_ie
