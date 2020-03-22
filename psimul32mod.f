!-----------------------------------------------------------------------
!
      module psimul32d
! Higher level subroutines for electrostatics
! written by viktor k. decyk, ucla
! copyright 1999, regents of the university of california
! update: may 2, 2008
!
      use globals, only: LINEAR, QUADRATIC
      use pdiag32d, only: vdist, bfopen
!     use espush2d, only: dpost, push, rpush, pushzf, rpushzf, pushglx, &
!    &rpushglx, pushgl, rpushgl, fft-
      use pespush32d, only: get_funit, plbcast, plsum, plmax, wrdata,   &
     &rddata, writebf, paguard, pcguard, fft
!     use field2d, only: sguard, aguard, cguard, spois, pois, gtmodes
      use pfield32d, only: aguard, cguard, pois, gtmodes
      implicit none
      private
      public :: restart_open, restart_bwrite, restart_dwrite
      public :: restart_bread, restart_dread
!     public :: dpostg, pushg, pushzfg, pushglg, pushglxg
      public :: initmodediag, initveldiag
      public :: dendiag, veldiag, potdiag, esenergy
!
      contains
!
         subroutine restart_open(nustrt,ntr,idrun,iur1,iur2,iuer)
! open restart files
         implicit none
         integer, intent(in) :: nustrt, ntr, idrun
         integer :: iur1, iur2, iuer
! local data
         integer :: ierr
         character(len=10) :: cdrun
         character(len=32) :: fname
! create string from idrun
         write (cdrun,'(i10)') idrun
         cdrun = adjustl(cdrun)
         iur1 = 0; iur2 = 0
! open old restart files
         if ((nustrt==2).or.(nustrt==0)) then
            iur1 = get_funit(16)
            fname = 'rstrt1.'//cdrun
            open(unit=iur1,file=trim(fname),form='unformatted',status='o&
     &ld',iostat=ierr)
            if (ierr /= 0)  then
               iur1 = -1
               write (iuer,*) 'Cannot open restart file1=',fname
            endif
            iur2 = get_funit(17)
            fname = 'rstrt2.'//cdrun
            open(unit=iur2,file=trim(fname),form='unformatted',status='o&
     &ld',iostat=ierr)
            if (ierr /= 0)  then
               iur2 = -1
               write (iuer,*) 'Cannot open restart file2=',fname
            endif
! open new restart files
         else if (ntr > 0) then
            iur1 = get_funit(16)
            fname = 'rstrt1.'//cdrun
            open(unit=iur1,file=trim(fname),form='unformatted',status='u&
     &nknown')
            iur2 = get_funit(17)
            fname = 'rstrt2.'//cdrun
            open(unit=iur2,file=trim(fname),form='unformatted',status='u&
     &nknown')
         endif
         end subroutine restart_open
!
         subroutine restart_bwrite(kunit,id0,itime,itime0,nvp,npp,part,m&
     &ovion,nppi,parti,qi,ef,bf)
! write file for basic restart or continuation
         implicit none
         integer, intent(in) :: kunit, id0, itime, itime0, nvp, movion
         real, dimension(:,:,:), pointer :: part, parti
         integer, dimension(:), pointer :: npp, nppi
         real, dimension(:,:,:,:), pointer :: qi
         complex, dimension(:,:,:,:,:), pointer, optional :: ef, bf
! local data
         integer :: j, it
         character(len=1), dimension(8), save :: arch = ' '
         integer, external :: NDIAN, NDPREC, IDPREC
! determine architecture information
         if (arch(1)==' ') then
! determine if number format is big or little endian
            it = NDIAN()
            if (it==0) then
               arch(1) = 'L'
            else if (it==1) then
               arch(1) = 'B'
            endif
! determine if default reals are double precision
            it = NDPREC()
            if (it==0) then
               arch(2) = 'S'
            else if (it==1) then
               arch(2) = 'D'
            endif
! determine if default integers are double precision
            it = IDPREC()
            if (it==0) then
               arch(3) = 'S'
            else if (it==1) then
               arch(3) = 'D'
            endif
         endif
! write out architecture information
         if (id0==0) write (kunit) (arch(j),j=1,8)
! write out current and initial time
         if (id0==0) write (kunit) itime, itime0
! write out number of processors
         if (id0==0) write (kunit) nvp
! write out size of particle array
         if (id0==0) write (kunit) size(part,1)
! write out electrons, if non-zero
         call wrdata(part,npp,kunit)
! write out if ions are moving
         if (id0==0) write (kunit) movion
         if (movion > 0) then
! write out size of ion array
            if (id0==0) write (kunit) size(parti,1)
! write out number of ions, size of particle array
            call wrdata(parti,nppi,kunit)
         else if (movion==0) then
! write out ion density, if ions are not moving
            if (id0==0) write (kunit) size(qi)
            call wrdata(qi,nvp,kunit)
         endif
! write out electromagnetic fields, if present
         it = 0
         if (present(ef)) then
            if (id0==0) write (kunit) size(ef)
            call wrdata(ef,nvp,kunit)
         else
            if (id0==0) write (kunit) it
         endif
         if (present(bf)) then
            if (id0==0) write (kunit) size(bf)
            call wrdata(bf,nvp,kunit)
         else
            if (id0==0) write (kunit) it
         endif
! write out current time for later confirmation
         if (id0==0) write (kunit) itime
         end subroutine restart_bwrite
!
         subroutine restart_dwrite(kunit,id0,itime,itw,wt,ndrec,fdname,n&
     &prec,fpname,narec,faname,nerec,fename)
! write diagnostic file for restart or continuation
         implicit none
         integer, intent(in) :: kunit, id0, itime
         integer :: itw, ndrec, nprec
         real, dimension(:,:), pointer :: wt
         character(len=*), intent(in) :: fdname, fpname
         integer, optional :: narec, nerec
         character(len=*), intent(in), optional :: faname, fename
! local data
         integer :: it, na, ne, irc
         character(len=32) :: fname
         na = 0; ne = 0
         if ((present(narec)).and.(present(faname))) na = 1
         if ((present(nerec)).and.(present(fename))) ne = 1
! write out number of diagnostics in file
         it = 5
         if (id0==0) write (kunit) it
! write out current energy time step
         if (id0==0) write (kunit) itw
         if (itw > 0) then
! write out size of energy array
            if (id0==0) write (kunit) size(wt,2)
! write out energy values
            call wrdata(wt,1,kunit)
         endif
! write out ion density diagnostic write location
         if (id0==0) write (kunit) ndrec
! write out record length
         if (ndrec > 0) then
            if (id0==0) then
               inquire(file=fdname,recl=it,iostat=irc)
               if (irc /= 0) it = 0
               write (kunit) it
               if (it > 0) then
                  fname = fdname
                  write (kunit) fname
               endif
            endif
         endif
! write out potential diagnostic write location
         if (id0==0) write (kunit) nprec
! write out record length
         if (nprec > 0) then
            if (id0==0) then
               inquire(file=fpname,recl=it,iostat=irc)
               if (irc /= 0) it = 0
               write (kunit) it
               if (it > 0) then
                  fname = fpname
                  write (kunit) fname
               endif
            endif
         endif
! write out vector potential diagnostic write location
         if (na==1) then
            if (id0==0) write (kunit) narec
! write out record length
            if (narec > 0) then
               if (id0==0) then
                  inquire(file=faname,recl=it,iostat=irc)
                  if (irc /= 0) it = 0
                  write (kunit) it
                  if (it > 0) then
                     fname = faname
                     write (kunit) fname
                  endif
               endif
            endif
         else
            if (id0==0) write (kunit) na
         endif
! write out electromagnetic diagnostic write location
         if (ne==1) then
            if (id0==0) write (kunit) nerec
! write out record length
            if (nerec > 0) then
               if (id0==0) then
                  inquire(file=fename,recl=it,iostat=irc)
                  if (irc /= 0) it = 0
                  write (kunit) it
                  if (it > 0) then
                     fname = fename
                  write (kunit) fname
                  endif
               endif
            endif
         else
            if (id0==0) write (kunit) ne
         endif
! write current time for later confirmation
         if (id0==0) then
            write (kunit) itime
            end file kunit
            rewind kunit
         endif
         end subroutine restart_dwrite
!
         subroutine restart_bread(iur1,iur2,id0,kunit,itime,itime0,nvp,n&
     &pp,part,movion,nppi,parti,qi,irc,iuer,ef,bf)
! read file for basic restart or continuation
         implicit none
         integer, intent(in) :: iur1, iur2, id0, nvp, movion
         integer :: iuer, kunit, itime, itime0, irc
         real, dimension(:,:,:), pointer :: part, parti
         integer, dimension(:), pointer :: npp, nppi
         real, dimension(:,:,:,:), pointer :: qi
         complex, dimension(:,:,:,:,:), pointer, optional :: ef, bf
! local data
         integer :: j, it
         integer, dimension(2) :: ktime
         character(len=1), dimension(8), save :: arch = ' '
         integer, external :: NDIAN, NDPREC, IDPREC
         if (id0==0) then
! determine most recent restart file
            if (kunit <= 0) then
               ktime(1) = -1
               read (iur1,iostat=irc) (arch(j),j=1,8)
               if (irc==0) then
                  read (iur1,iostat=irc) ktime(1)
                  if (irc /= 0) ktime(1) = -1
                  rewind iur1
               endif
               ktime(2) = -1
               read (iur2,iostat=irc) (arch(j),j=1,8)
               if (irc==0) then
                  read (iur2,iostat=irc) ktime(2)
                  if (irc /= 0) ktime(2) = -1
                  rewind iur2
               endif
               if (ktime(1) > ktime(2)) then
                  ktime(2) = iur1
               else if (ktime(2) >= 0) then
                  ktime(1) = ktime(2)
                  ktime(2) = iur2
               endif
! switch restart files
            else
               if (kunit==iur1) then
                  ktime(1) = 0
                  ktime(2) = iur2
               else if (kunit==iur2) then
                  ktime(1) = 0
                  ktime(2) = iur1
               else
                  ktime(1) = 0
                  ktime(2) = -1
               endif
            endif
         endif
         call plbcast(ktime)
! check if unit number is valid
         irc = 99; if (ktime(2) < 0) return
         kunit = ktime(2)
! read in architecture information
         if (id0==0)  then
            read (kunit,err=10) (arch(j),j=1,8)
! convert restart architecture information to integer
            ktime(1) = -1
! determine if number format is big or little endian
            if (arch(1)=='L') then
              ktime(1) = 0
            else if (arch(1)=='B') then
              ktime(1) = 1
            endif
! determine if default reals are double precision
            if (arch(2)=='S') then
               ktime(1) = 2*ktime(1)
            else if (arch(2)=='D') then
               ktime(1) = 2*ktime(1) + 1
            endif
! determine if default integers are double precision
            if (arch(3)=='S') then
               ktime(1) = 2*ktime(1)
            else if (arch(3)=='D') then
               ktime(1) = 2*ktime(1) + 1
            endif
! convert current architecture information to integer
            it = 4*NDIAN() + 2*NDPREC() + IDPREC()
            ktime(1) = abs(ktime(1)-it)
         endif
         call plbcast(ktime)
         irc = 1; if (ktime(1) > 0) go to 10
! read in current and initial time
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1), ktime(2)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = 2; if (ktime(1) < 0) go to 10
         itime = ktime(1); itime0 = ktime(2)
! read in number of processors
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = 3; if (ktime(1) /= nvp) go to 10
! read in size of particle array
         if (id0==0) then
           read (kunit,iostat=irc) ktime(1)
           if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = 4; if (ktime(1) /= size(part,1)) go to 10
! read in electrons, if non-zero
         call rddata(part,npp,kunit,it)
         irc = 5; if (it /= 0) go to 10
! read in if ions are moving
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = 6; if (ktime(1) /= movion) go to 10
         if (movion > 0) then
! read in size of ion array
            if (id0==0) then
               read (kunit,iostat=irc) ktime(1)
               if (irc /= 0) ktime(1) = -1
            endif
            call plbcast(ktime)
            irc = 7; if (ktime(1) /= size(parti,1)) go to 10
! read in ions, if non-zero
            call rddata(parti,nppi,kunit,it)
            irc = 8; if (it /= 0) go to 10
! read in ion density, if ions are not moving
         else if (movion==0) then
            if (id0==0) then
               read (kunit,iostat=irc) ktime(1)
               if (irc /= 0) ktime(1) = -1
            endif
            call plbcast(ktime)
            irc = 9; if (ktime(1) > size(qi)) go to 10
            call rddata(qi,nvp,kunit,it)
            irc = 10; if (it /= 0) go to 10
         endif
! read in first electromagnetic field, if present
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = 11; if (ktime(1) < 0) go to 10
         if (ktime(1) > 0) then
            if (.not.present(ef)) go to 10
            if (ktime(1) > size(ef)) go to 10
            irc = 12; call rddata(ef,nvp,kunit,it)
            if (it /= 0) go to 10
         endif
! read in second electromagnetic field, if present
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = 13; if (ktime(1) < 0) go to 10
         if (ktime(1) > 0) then
            if (.not.present(bf)) go to 10
            if (ktime(1) > size(bf)) go to 10
            irc = 14; call rddata(bf,nvp,kunit,it)
            if (it /= 0) go to 10
         endif
! read in current time for confirmation
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = 15; if (ktime(1) /= itime) go to 10
         irc = 0
         return
! write out errors
   10    if (id0==0) then
            write (iuer,*) 'Basic Restart Error, irc = ', irc
            if (irc==1) then
               write (iuer,*) 'Architecture=', (arch(j),j=1,8)
            else if (irc==3) then
               write (iuer,*) 'nvp=', ktime(1)
            else if (irc==4) then
               write (iuer,*) 'size(part,1)=', ktime(1)
            else if (irc==6) then
               write (iuer,*) 'movion=', ktime(1)
            else if (irc==7) then
               write (iuer,*) 'size(parti,1)=', ktime(1)
            else if (irc==9) then
               write (iuer,*) 'size(qi)=', ktime(1)
            else if (irc==11) then
               write (iuer,*) 'size(ef)=', ktime(1)
            else if (irc==13) then
               write (iuer,*) 'size(bf)=', ktime(1)
            else if (irc==15) then
               write (iuer,*) 'confirmation itime,it=',itime,ktime(1)
            endif
         endif
         end subroutine restart_bread
!
         subroutine restart_dread(kunit,id0,itime,itw,wt,iud,ndrec,fdnam&
     &e,iup,nprec,fpname,irc,iuer,iua,narec,faname,iue,nerec,fename)
! read diagnostic file for restart or continuation
         implicit none
         integer, intent(in) :: kunit, id0, itime
         integer :: itw, iud, ndrec, iup, nprec, irc, iuer
         real, dimension(:,:), pointer :: wt
         character(len=*) :: fdname, fpname
         integer, optional :: iua, narec, iue, nerec
         character(len=*), optional :: faname, fename
! local data
         integer :: it, nt, na, ne, ierr
         integer, dimension(1) :: ktime
         character(len=32) :: fname
         na = 0; ne = 0
         if ((present(iua)).and.(present(narec)).and.(present(faname))) &
     &na = 1
         if ((present(iue)).and.(present(nerec)).and.(present(fename))) &
     &ne = 1
! check if unit number is valid
         irc = -99; if (kunit < 0) return
! read in number of diagnostics in file
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = -98; if (ktime(1) < 0) go to 20
         nt = ktime(1)
! read in current energy time step
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = -1; if (ktime(1) < 0) go to 20
         itw = ktime(1)
         if (itw > 0) then
! read in size of energy array
            if (id0==0) then
               read (kunit,iostat=irc) ktime(1)
               if (irc /= 0) ktime(1) = -1
            endif
            call plbcast(ktime)
            irc = -2
            if ((itw > size(wt,1)).or.(ktime(1) > size(wt,2))) go to 20
! read in energy values
            irc = -3; call rddata(wt,1,kunit,it)
            if (it /= 0) go to 20
         endif
         if (nt==1) go to 10
! read in ion density diagnostic write location
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = -4; if (ktime(1) < 0) go to 20
         ndrec = ktime(1)
! read in record length and open file
         if (ndrec > 0) then
            if (id0==0) then
               read (kunit,iostat=irc) ktime(1)
               if (irc /= 0) ktime(1) = -1
            endif
            call plbcast(ktime)
            irc = -5; if (ktime(1) < 1) go to 20
            it = ktime(1)
            if (id0==0) then
               ktime(1) = 0
               read (kunit,iostat=irc) fname
               if (irc /= 0) ktime(1) = -1
            endif
            call plbcast(ktime)
            irc = -6; if (ktime(1) < 0) go to 20
            if (id0==0) then
               fdname = fname
               iud = get_funit(iud)
               open(unit=iud,file=fdname,form='unformatted',access='dire&
     &ct',recl=it,status='old',iostat=ktime(1))
            endif
            call plbcast(ktime)
            irc = -7; if (ktime(1) /= 0) go to 20
         endif
         if (nt==2) go to 10
! read in potential diagnostic write location
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = -8; if (ktime(1) < 0) go to 20
         nprec = ktime(1)
! read in record length and open file
         if (nprec > 0) then
            if (id0==0) then
               read (kunit,iostat=irc) ktime(1)
               if (irc /= 0) ktime(1) = -1
            endif
            call plbcast(ktime)
            irc = -9; if (ktime(1) < 1) go to 20
            it = ktime(1)
            if (id0==0) then
               ktime(1) = 0
               read (kunit,iostat=irc) fname
               if (irc /= 0) ktime(1) = -1
            endif
            call plbcast(ktime)
            irc = -10; if (ktime(1) < 0) go to 20
            if (id0==0) then
               fpname = fname
               iup = get_funit(iup)
               open(unit=iup,file=fpname,form='unformatted',access='dire&
     &ct',recl=it,status='old',iostat=ktime(1))
            endif
            call plbcast(ktime)
            irc = -11; if (ktime(1) /= 0) go to 20
         endif
         if (nt==3) go to 10
! read in vector potential diagnostic write location
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = -12; if (ktime(1) < 0) go to 20
         it = ktime(1)
         if (na==1) then
            narec = it
! read in record length and open file
            if (narec > 0) then
               if (id0==0) then
                  read (kunit,iostat=irc) ktime(1)
                  if (irc /= 0) ktime(1) = -1
               endif
               call plbcast(ktime)
               irc = -13; if (ktime(1) < 1) go to 20
               it = ktime(1)
               if (id0==0) then
                  ktime(1) = 0
                  read (kunit,iostat=irc) fname
                  if (irc /= 0) ktime(1) = -1
               endif
               call plbcast(ktime)
               irc = -14; if (ktime(1) < 0) go to 20
               if (id0==0) then
                  faname = fname
                  iua = get_funit(iua)
                  open(unit=iua,file=faname,form='unformatted',access='d&
     &irect',recl=it,status='old',iostat=ktime(1))
               endif
               call plbcast(ktime)
               irc = -15; if (ktime(1) /= 0) go to 20
            endif
         else
            if (it > 0) go to 20
         endif
         if (nt==4) go to 10
! read in electromagnetic diagnostic write location
         if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = -16; if (ktime(1) < 0) go to 20
         it = ktime(1)
         if (ne==1) then
            nerec = it
! read in record length and open file
            if (nerec > 0) then
               if (id0==0) then
                  read (kunit,iostat=irc) ktime(1)
                  if (irc /= 0) ktime(1) = -1
               endif
               call plbcast(ktime)
               irc = -17; if (ktime(1) < 1) go to 20
               it = ktime(1)
               if (id0==0) then
                  ktime(1) = 0
                  read (kunit,iostat=irc) fname
                  if (irc /= 0) ktime(1) = -1
               endif
               call plbcast(ktime)
               irc = -18; if (ktime(1) < 0) go to 20
               if (id0==0) then
                  fename = fname
                  iue = get_funit(iue)
                  open(unit=iue,file=fename,form='unformatted',access='d&
     &irect',recl=it,status='old',iostat=ierr)
               endif
               call plbcast(ktime)
               irc = -19; if (ktime(1) /= 0) go to 20
            endif
         else
            if (it > 0) go to 20
         endif
! read in current time for confirmation
   10    if (id0==0) then
            read (kunit,iostat=irc) ktime(1)
            if (irc /= 0) ktime(1) = -1
         endif
         call plbcast(ktime)
         irc = -20; if (ktime(1) /= itime) go to 20
         irc = 0
         if (id0==0) rewind kunit
         return
! write out errors
   20    if (id0==0) then
            write (iuer,*) 'Diagnostic Restart Error, irc = ', irc
            if (irc==(-1)) then
               write (iuer,*) 'itw=', itw
            else if (irc==(-2)) then
               write (iuer,*) 'itw,size(wt,2)=', itw, it
            else if ((irc==(-5)).or.(irc==(-9)).or.(irc==(-13)).or.(irc=&
     &=(-17))) then
               write (iuer,*) 'recl=', it
            else if ((irc==(-7)).or.(irc==(-11)).or.(irc==(-15)).or.(irc&
     &==(-19))) then
               write (iuer,*) 'fname=', fname
            else if (irc==(-20)) then
               write (iuer,*) 'confirmation itime,it=', itime, it
            endif
         endif
         end subroutine restart_dread
!
         subroutine initmodediag(dent,ntd,id0,nxh,nyh,nzh,kxyp,kyzp,mode&
     &sxd,modesyd,modeszd,jblok,mblok,iud,ndrec,fdname)
! initialize mode diagnostic
         implicit none
         integer :: ntd, id0, nxh, nyh, nzh, kxyp, kyzp
         integer :: modesxd, modesyd, modeszd, jblok, mblok, iud, ndrec
         character(len=*) :: fdname
         complex, dimension(:,:,:,:), pointer :: dent
! local data
         integer :: modesz2d
         if (ntd <= 0) return
         if (modesxd > nxh) modesxd = nxh
         if (modesyd > nyh) modesyd = nyh
         if (modeszd > nzh) modeszd = nzh
         modesz2d = 2*modeszd - 1
         allocate(dent(modesz2d,min(modesxd,kxyp),min(modesyd,kyzp),jblo&
     &k*mblok))
! open output file
         if (id0==0) then
            if (ndrec==0) then
               iud = get_funit(iud); ndrec = -1
               call bfopen(dent,modesz2d,iud,ndrec,trim(fdname))
            endif
         else
            if (ndrec==0) ndrec = 1
         endif
         end subroutine initmodediag
!
         subroutine initveldiag(fv,fvm,vtx,vty,vtz,ntv,ndv,id0,nmv,nblok&
     &,iuv,fvname)
! initialize velocity diagnostic
         implicit none
         integer :: ntv, ndv, id0, nmv, nblok, iuv
         real :: vtx, vty, vtz
         real, dimension(:,:,:), pointer :: fv, fvm
         character(len=*) :: fvname
         if ((ntv <= 0) .and. (ndv <= 0)) return
         allocate(fv(2*nmv+2,3,nblok),fvm(3,3,nblok))
! fix velocity range
         fv(1,:,:) = 8.*max(vtx,vty,vtz)
         if (ntv > 0) then
            if (id0==0) then
               iuv = get_funit(iuv)
               open(unit=iuv,file=trim(fvname),form='formatted',status='&
     &unknown')
! write captions
               write (iuv,*) 'it vdx vdy vdz vtx vty vtz sk'
            endif
         endif
         end subroutine initveldiag
!
         subroutine dendiag(qt,qs,qi,sfield,dent,sfieldt,ffc,nyzp,mixup,&
     &sct,tfft,ntd,ndd,nx,ny,nz,modesxd,modesyd,modeszd,iud,ndrec,indx,i&
     &ndy,indz,ntime,nvpy,nvpz,kstrt,kxyp,kyp,kyzp,kzp,ngds,iblok,jblok,&
     &kblok,mblok,inorder)
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
         if ((ntd > 0) .or. (ndd > 0)) then
            it = -1; if (ntd > 0) it = ntime - ntd*(ntime/ntd)
            jt = -1; if (ndd > 0) jt = ntime - ndd*(ntime/ndd)
            if ((it==0) .or. (jt==0)) then
               sfield = qi
! add guard cells for ion density in x
               call aguard(sfield,nyzp,nx,inorder)
! add guard cells for ion density in y and z
               call paguard(sfield,kstrt,nvpy,nvpz,nx,kyp,kzp,ngds,iblok&
     &,inorder)
! transform ion density to fourier space
               isign = -1
               call fft(sfield,qs,qt,isign,mixup,sct,tfft,indx,indy,indz&
     &,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
! calculate smoothing in fourier space
               call pois(qt,sfieldt,ffc,nx,ny,nz,kstrt,jblok)
! store selected fourier modes
               if (it==0) then
                  modesz2d = 2*modeszd - 1
                  call gtmodes(sfieldt,dent,nx,ny,nz,modesxd,modesyd,mod&
     &eszd,kstrt,jblok)
! write diagnostic output
                  call writebf(dent,modesxd,modesyd,modesz2d,kxyp,kyzp,j&
     &blok,iud,ndrec)
               endif
! transform ion density to real space
               if (jt==0) then
                  isign = 1
                  call fft(sfield,qs,sfieldt,isign,mixup,sct,tfft,indx,i&
     &ndy,indz,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
                  call pcguard(sfield,kstrt,nvpy,nvpz,kyp,kzp,ngds,iblok&
     &,inorder)
                  call cguard(sfield,nyzp,nx,inorder)
               endif
            endif
         endif
         end subroutine dendiag
!
         subroutine veldiag(part,fv,fvm,npp,ntv,ndv,id0,nmv,iuv,ntime)
! velocity diagnostic
         implicit none
         integer :: ntv, ndv, id0, nmv, iuv, ntime
         real, dimension(:,:,:), pointer :: part
         real, dimension(:,:,:), pointer :: fv, fvm
         integer, dimension(:), pointer :: npp
! local data
         integer :: it, jt, nt
         if ((ntv > 0) .or. (ndv > 0)) then
            it = -1; if (ntv > 0) it = ntime - ntv*(ntime/ntv)
            jt = -1; if (ndv > 0) jt = ntime - ndv*(ntime/ndv)
            if ((it==0) .or. (jt==0)) then
! calculate particle distribution function and moments
               call vdist(part,fv,fvm,npp,nmv)
               call plsum(fv(:,:,1))
! print out velocity moments
               if (it==0) then
                  if (id0==0) then
                     nt = ntime/ntv
                     write (iuv,*) nt, fvm(1,:,1), fvm(2,:,1), sum(fvm(3&
     &,:,1))
                  endif
               endif
            endif
         endif
         end subroutine veldiag
!
         subroutine potdiag(qt,qs,sfield,pott,sfieldt,ffc,nyzp,mixup,sct&
     &,tfft,ntp,ndp,nx,ny,nz,modesxp,modesyp,modeszp,iup,nprec,indx,indy&
     &,indz,ntime,nvpy,nvpz,kstrt,kxyp,kyp,kyzp,kzp,ngds,iblok,jblok,kbl&
     &ok,mblok,inorder)
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
            it = -1; if (ntp > 0) it = ntime - ntp*(ntime/ntp)
            jt = -1; if (ndp > 0) jt = ntime - ndp*(ntime/ndp)
            if ((it==0) .or. (jt==0)) then
! calculate potential in fourier space
               call pois(qt,sfieldt,ffc,we,nx,ny,nz,kstrt,jblok)
! store selected fourier modes
               if (it==0) then
                  modesz2p = 2*modeszp - 1
                  call gtmodes(sfieldt,pott,nx,ny,nz,modesxp,modesyp,mod&
     &eszp,kstrt,jblok)
! write diagnostic output
                  call writebf(pott,modesxp,modesyp,modesz2p,kxyp,kyzp,j&
     &blok,iup,nprec)
               endif
! transform potential to real space
               if (jt==0) then
                  isign = 1
                  call fft(sfield,qs,sfieldt,isign,mixup,sct,tfft,indx,i&
     &ndy,indz,kstrt,kxyp,kyp,kyzp,kzp,kblok,mblok,inorder)
                  call pcguard(sfield,kstrt,nvpy,nvpz,kyp,kzp,ngds,iblok&
     &,inorder)
                  call cguard(sfield,nyzp,nx,inorder)
               endif
            endif
         endif
         end subroutine potdiag
!
         subroutine esenergy(wt,wtot,msg,we,wke,wki,ntw,ndw,id0,itw,iuot&
     &,ntime)
! electrostatic energy diagnostic
         implicit none
         integer :: ntw, ndw, id0, itw, iuot, ntime
         real :: we, wke, wki
         real, dimension(:,:), pointer :: wt
         real, dimension(4) :: wtot
         double precision, dimension(:) :: msg
! local data
         integer :: it, jt
  992    format (' field, kinetic, total energies = ',3e14.7)
         if ((ntw > 0) .or. (ndw > 0)) then
            it = -1; if (ntw > 0) it = ntime - ntw*((ntime+1)/ntw) + 1
            jt = -1; if (ndw > 0) jt = ntime - ndw*((ntime+1)/ndw) + 1
            if ((it==0) .or. (jt==0)) then
               wtot(1) = we
               wtot(2) = wke
               wtot(3) = wki
               wtot(4) = we + wke + wki
               call plsum(wtot)
! send energy values to diagnostic node
               msg(1:4) = wtot
               call HARTBEAT(msg,4)
               if (it==0) then
                  if (id0==0) write (iuot,992) wtot(1), wtot(2), wtot(4)
               endif
               if (jt==0) then
                  itw = itw + 1
                  wt(itw,:) = wtot
               endif
            endif
         endif
         end subroutine esenergy
!
      end module psimul32d
