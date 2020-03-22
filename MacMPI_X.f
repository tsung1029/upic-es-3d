c-----------------------------------------------------------------------
c Partial MPI library based on Open Transport in the Macintosh OS,
c using TCP/IP protocol.
c No local buffering of messages is implemented, so that all messages
c must be received in the order sent, and receives with wildcard
c sources are not supported. the following subroutines are implemented:
c MPI_INIT, MPI_FINALIZE, MPI_SEND, MPI_RECV, MPI_ISEND, MPI_IRECV
c MPI_TEST, MPI_WAIT, MPI_SENDRECV, MPI_SSEND, MPI_ISSEND, MPI_WAITALL
c MPI_WAITANY, MPI_GET_COUNT, MPI_INITIALIZED, MPI_COMM_SIZE
c MPI_COMM_RANK, MPI_COMM_DUP, MPI_COMM_SPLIT, MPI_COMM_FREE
c MPI_CART_CREATE, MPI_CART_COORDS, MPI_CART_GET, MPI_CART_SHIFT
c MPI_CART_RANK, MPI_CART_SUB, MPI_DIMS_CREATE
c MPI_BCAST, MPI_BARRIER, MPI_REDUCE, MPI_SCAN
c MPI_ALLREDUCE, MPI_GATHER, MPI_ALLGATHER, MPI_SCATTER, MPI_ALLTOALL
c MPI_GATHERV, MPI_ALLGATHERV, MPI_SCATTERV, MPI_ALLTOALLV
c MPI_REDUCE_SCATTER, MPI_ABORT, MPI_WTIME, MPI_WTICK, MPI_TYPE_EXTENT
c MPI_REQUEST_FREE, MPI_GET_PROCESSOR_NAME, MPI_ERRHANDLER_SET
c Open Transport is described in Inside Macintosh: Networking with Open
c Transport, verson 1.3 [Apple Computer, Cupertino, CA, 1997], and at:
c http://developer.apple.com/techpubs/mac/NetworkingOT/NetworkingWOT-2.
c html
c The Message Passing Interface (MPI) is described in the reference,
c M. Snir, S. Otto, S. Huss-Lederman, D. Walker, and J. Dongarra,
c MPI: The Complete Reference [MIT Press, Cambridge, MA,1996].
c Fortran unit 2 is used throughout for error messages, and
c Fortran units 3 and 4 are used in MPI_INIT.
c written by viktor k. decyk, ucla
c copyright 1999, regents of the university of california.
c all rights reserved.
c no warranty for proper operation of this software is given or implied.
c software or information may be copied, distributed, and used at own
c risk; it may not be distributed without this notice included verbatim
c with each file. 
c update: april 1, 2004
      block data
      implicit none
c declare internal mpi common block
      integer nproc, idproc, cfig0, stime, mapcomm, notifierUPP
      integer MAXS, MAXC, MAXD, MAXQ, epref, ioc, nevents
      parameter(MAXS=32,MAXC=10,MAXD=6,MAXQ=MAXC*(MAXS+MAXD+3))
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c cfig0 = OTConfiguration structure template
c MAXS = maximum number of nodes connected
c MAXC = maximum number of communicators
c MAXD = maximum number of topology dimensions
c LDM = (0,1) = (NO,YES) register AppleTalk port for Launch Den Mother
c epref = array of endpoint references for each participating node
c ioc = array of context pointers for notifier function
c ioc(1,i) = endpoint reference for notifier for endpoint i
c ioc(2,i) = processor id for listener for endpoint i
c ioc(3,i) = handle for current receive from endpoint i
c ioc(4,i) = handle for current send to endpoint i
c nevents = log of unknown notifier events
c mapcomm = communicator map
c mapcomm(1:nproc,i) = actual proc id for given rank in communicator i
c mapcomm(nproc+1:MAXS,i) = MPI_UNDEFINED
c mapcomm(MAXS+1,i) = number of processes in comm i
c mapcomm(MAXS+2,i) = rank for this node in comm i
c mapcomm(MAXS+3,i) = ndims = number of dimensions in topology in comm i
c mapcomm(MAXS+4:MAXS+3+ndims,i) = size of dimension in comm i,
c negative if non-periodic
c mapcomm(MAXS+4+ndims:MAXS+3+MAXD,i) = 0
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm, notifierUPP
c declare common block for non-blocking messages
      integer MAXM, MAXMS, MAXP
      integer curreq, header, rwrec, monitor, trash, mqueue
      parameter(MAXM=2*MAXS,MAXMS=5*MAXM,MAXP=2*(MAXS+1))
      dimension curreq(5,MAXM), header(10,MAXM), rwrec(14,MAXM)
      dimension trash(256), mqueue(2,MAXS+1)
c MAXM = maximum number of outstanding messages on a node
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
c curreq = request record for transmission parameters, see rwstat
c header = message envelope
c header(1,i) = communicator, header(2,i) = tag, header(3,i) = datatype
c header(4,i) = length of data (in bytes) for message handle i
c rwrec = read/write record for asynchronous messages, see rwstat
c mqueue = message request queue
c mqueue(1,i) = end of message queue for receives from endpoint i
c mqueue(2,i) = end of message queue for sends to endpoint i
      common /mpisendrec/ monitor, curreq, header, rwrec, trash, mqueue
c common block for message window
c cpptr = pointer to window structure
c nsp = amount of space between boxes
c nbx = size of box
c nds = number of message sizes monitored
c mbs = maximum maximum bandwidth in MB/sec expected
      integer cpptr
      integer*2 crect(4), nsp, nbx, nds, mbs
      common /winmess/ cpptr, crect, nsp, nbx, nds, mbs
c declare adsp common block
      integer nameid
c nameid = adsp portname ID
      common /adsp/ nameid
c declare error handler common block
      integer errh
      dimension errh(MAXC)
c errh = error handler
      common /mpierrh/ errh
      save /mpiparms/, /mpisendrec/, /winmess/, /adsp/, /mpierrh/
      data nproc, cfig0, notifierUPP /-1,0,0/
      data epref /MAXS*0,0/
      data nevents /MAXS*0,0/
      data mapcomm /MAXQ*0/
      data curreq /MAXMS*0/
      data monitor /1/
      data mqueue /MAXP*0/
      data cpptr, nsp, nbx, nds, mbs /0,8,16,24,10/
      data nameid /0/
      data errh /MAXC*1/
      end
c-----------------------------------------------------------------------
      subroutine MPI_INIT(ierror)
c initialize the MPI execution environment
c ierror = error indicator
c input: none, output: ierror
      implicit none
      integer ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, LDM, epref, ioc, nevents, stime, mapcomm
      integer notifierUPP
      parameter(MAXS=32,MAXC=10,MAXD=6,LDM=0)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c idproc = processor id
c cfig0 = OTConfiguration structure template
c epref = array of endpoint references for each participating node
c ioc = array of context pointers for notifier function
c stime = first time stamp if MPI_Init successful
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm, notifierUPP
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c function declarations
      integer InitOpenTransportInContext, OTCreateConfiguration
      integer OTCloneConfiguration, OTOpenEndpointInContext
      integer OTInetGetInterfaceInfo, OTCloseProvider
      integer OTInetHostToString, OTInetStringToHost, NewOTNotifyUPP
      integer OTConnect, OTElapsedMilliseconds, OTGetEndpointState
      integer otpinit, ioresult, lencstr, lenstr
      logical checkesc
      external InitOpenTransportInContext, OTCreateConfiguration
      external OTCloneConfiguration, OTOpenEndpointInContext
      external OTInetGetInterfaceInfo, OTCloseProvider
      external OTInetHostToString, OTInetStringToHost, NewOTNotifyUPP
      external OTConnect, OTElapsedMilliseconds, OTGetEndpointState
      external otpinit, ioresult, lencstr, lenstr, checkesc, notifier
c OT constants
      integer kInitOTForApplicationMask
      integer kOTNotFoundErr, AF_INET, kOTNoAddressErr
      integer kOTAnyInetAddress, kOTOutStateErr, kOTBadAddressErr
      integer kDefaultInetInterface
      parameter(kInitOTForApplicationMask=1)
      parameter(kOTNotFoundErr=-3201,AF_INET=2,kOTNoAddressErr=-3154)
      parameter(kOTAnyInetAddress=0,kOTOutStateErr=-3155)
      parameter(kOTBadAddressErr=-3150,kDefaultInetInterface=-1)
c MPI constants
      integer MPI_UNDEFINED, MPI_STATUS_SIZE, MPI_INTEGER
      parameter(MPI_UNDEFINED=-1,MPI_STATUS_SIZE=5,MPI_INTEGER=18)
c local data
      integer*2 intw(2), info(276), addrb(8), addrl(8), addrn(8)
      integer longw, address, nerr, nofile, nv, i, oss, response, oldadd
      integer cfig, portnum, pshift, epref0
      integer pnumid(3), reqad(4), repad(4), scall(10), enum(8)
      integer ptime(2), delta(2), epinfo(8)
      integer stat(MPI_STATUS_SIZE)
      character*36 ename
      character*16 myself, location, fname
c longw and intw are used to convert types between integer*2 and integer
      equivalence(longw,intw)
      data cfig, portnum, pshift, epref0 /0,0,0,0/
c set MPI_COMM_WORLD and MPI_COMM_SELF mapping
      do 5 i = 1, MAXS
      mapcomm(i,1) = i - 1
      mapcomm(i,2) = MPI_UNDEFINED
    5 continue
      if (monitor.eq.2) write (2,*) 'MPI_INIT started'
c
c Open Open Transport
c
c Get information about the operating environment
      call Gestalt(VAL4('sysv'),response)
c Check if MacOS 8.6 or higher is installed
      if (response.lt.2144) then
         write (2,*) 'MacMPI_X requires MacOS 8.6 or higher'
         nv = response/256
         response = response - 256*nv
         nerr = response/16
         response = response - 16*nerr
         write (ename,'(i1,".",2i1)') nv, nerr, response
         write (2,*) 'MacOS detected: ', ename
         ierror = -7
         return
      endif
c Initialize Open Transport for use by application
      oss = InitOpenTransportInContext(val4(kInitOTForApplicationMask),v
     1al4(0))
      if (oss.ne.0) then
         write (2,*) 'Cannot Open Transport, oss = ', oss
         ierror = oss
         return
      endif
c Obtain the current time stamp
      call OTGetTimeStamp(stime)
c Create a structure defining a provider's configuration
      cfig0 = OTCreateConfiguration('tcp'//char(0))
      if (cfig0.eq.0) then
         write (2,*) 'Insufficient Memory for tcp Endpoint'
         ierror = -4
         call CloseOpenTransportInContext(val4(0))
         return
      elseif (cfig0.eq.(-1)) then
         write (2,*) 'Invalid Configuration for tcp Endpoint'
         ierror = -4
         call CloseOpenTransportInContext(val4(0))
         return
      endif
c Copy an OTConfiguration structure
      cfig = OTCloneConfiguration(val4(cfig0))
c Create a dummy endpoint provider just to make sure TCP/IP is active
      epref0 = OTOpenEndpointInContext(val4(cfig),val4(0),epinfo,oss,val
     14(0))
c Obtain information about the Internet environment
      oss = OTInetGetInterfaceInfo(info,val4(kDefaultInetInterface))
      if (oss.ne.0) then
         write (2,*) 'INS GetInterfaceInfo Error, oss = ', oss
         if (oss.eq.kOTNotFoundErr) then
            write (2,*) 'Requested information does not exist'
            write (2,*) 'TCP/IP may not be active'
         endif
         ierror = oss
         oss = OTCloseProvider(val4(epref0))
         call CloseOpenTransportInContext(val4(0))
         return
      endif
c Close the dummy endpoint provider
      oss = OTCloseProvider(val4(epref0))
      intw(1) = info(1)
      intw(2) = info(2)
      address = longw
      oldadd = address
c
c Everyone opens a port
c
c Open file containing portnum (and possibly participating nodes)
c first line in nodelist file on all nodes contains common portnum
c if the file is missing or empty, a default number of 5013 is used
      open(unit=3,file='nodelist_ip',form='formatted',status='old',iosta
     1t=nofile)
      if (nofile.ne.0) then
         open(unit=3,file='nodelist',form='formatted',status='old',iosta
     1t=nofile)
      endif
      if (nofile.eq.0) then
         read (3,'(i10)',iostat=nerr) portnum
         if ((nerr.ne.0).or.(portnum.lt.5000).or.(portnum.gt.49152)) the
     1n
            portnum = 5013
         endif
      else
         portnum = 5013
      endif
c Construct InetAddress of listener
      addrn(1) = AF_INET
      addrn(2) = portnum
      addrn(3) = info(1)
      addrn(4) = info(2)
c Convert an address into a character string
      intw(1) = addrn(3)
      intw(2) = addrn(4)
      call OTInetHostToString(val4(longw),myself)
      nv = lencstr(myself)
      myself(nv+1:16) = ' '
      ename = myself
c Establish port connection end
c Set TNetbuf for specifying address
      reqad(2) = 16
      reqad(3) = loc(addrn)
      reqad(4) = 1
c Set TNetbuf for address which is returned
      repad(1) = 16
      repad(3) = loc(addrb)
c Set TNetbuf for specifying address for local computer
      scall(2) = 16
      scall(3) = loc(addrn)
      nproc = 0
c Create Univeral Proc Ptr for notifier function
      notifierUPP = NewOTNotifyUPP(notifier)
c Initialize synchronous listener endpoint provider
      oss = otpinit(MAXS+1,reqad,repad)
      if (oss.ne.0) then
         if (oss.eq.kOTNoAddressErr) then
            write (2,*) 'portnum likely already exists = ', portnum
         endif
         ierror = oss
         call MPI_FINALIZE(nerr)
         return
      endif
c Save copy of portnum
      pnumid(1) = addrn(2)
c Register adsp endpoint for puppy
      if (LDM.ne.0) then
         call regport(portnum,ierror)
         if (ierror.ne.0) then
            write (2,*) 'Unable to register adsp port for puppy'
c           call MPI_FINALIZE(nerr)
         endif
      endif
c Pass processor id to listener notifier
      ioc(2,MAXS+1) = nproc + 1
c debug
      if (monitor.eq.2) then
         write (2,*) 'local host=', myself(1:nv), ',port=', addrb(2)
      endif
c
c Determine if node is master (idproc=0) or slave (idproc>0).
c on the master node, the second and subsequent lines of nodelist file
c contain IP addresses of the nodes participating, in dotted-decimal
c format (for example, "12.13.14.15").
c if this list of nodes is missing, then the node is a slave.
c every node also makes a connection to itself.
c
      location = ' '
      if (nofile.eq.0) read (3,'(a16)',iostat=nerr) location
c must be slave
      if ((nofile.ne.0).or.(nerr.ne.0).or.(location.eq.' ')) then
         idproc = 1
c may be master
      else
         if (((location.eq.'self').or.(location.eq.myself)).and.(pnumid(
     11).eq.portnum)) then
            idproc = 0
            close(unit=2)
            fname = 'MPIERRS00'
            open(unit=2,file=fname,form='formatted',status='unknown')
         else
c must be slave 
            idproc = 1
         endif
      endif
c Update port number 
      portnum = pnumid(1)
c 
c * * * begin main iteration loop * * *
c
c Prepare to accept connection
c
c Establish connection end
c Connections to oneself requires a new socket
   10 addrb(2) = kOTAnyInetAddress
c Set TNetbuf for specifying address
      reqad(2) = 16
      reqad(3) = loc(addrb)
      reqad(4) = 0
c Set TNetbuf for address which is returned
      repad(1) = 16
      repad(3) = loc(addrl)
c Initialize synchronous endpoint provider
      oss = otpinit(nproc+1,reqad,repad)
      if (oss.ne.0) then
         ierror = oss
         call MPI_FINALIZE(nerr)
         return
      endif
c For connections to oneself, jump to OTConnect
      if (idproc.eq.nproc) then
         ioc(2,MAXS+1) = MAXS+1
         go to 70
      endif
c Obtain the current time stamp
      call OTGetTimeStamp(ptime)
c Wait for connection
   20 if (ioresult(ioc(1,MAXS+1)).gt.0) then
         if (checkesc(1)) then
            ierror = -9
            call writerrs('MPI_INIT: ',ierror)
            return
c Wait up to one minute for next connection
         else
            if (OTElapsedMilliseconds(ptime).lt.60000) then
               go to 20
            else
               write (2,*) 'OTListen Wait Time Exceeded'
               ierror = -5
               call MPI_FINALIZE(nerr)
               return
            endif
         endif
      else
         oss = ioresult(ioc(1,MAXS+1))
         if (oss.ne.0) then
            write (2,*) 'OTListen failed, oss = ', oss
            ierror = oss
            call MPI_FINALIZE(nerr)
            return
         endif
      endif
c Receive portnum for verification and processor id from remote node
      nproc = nproc + 1
      mapcomm(MAXS+1,1) = nproc
c Reset iocompletion flag to notifier
      ioc(2,MAXS+1) = nproc + 1
      call MPI_RECV(enum,3,MPI_INTEGER,nproc-1,3,0,stat,ierror)
c Extract processor id and base portnum and shift on first connection
      if (nproc.eq.1) then
         idproc = enum(2)
         pshift = portnum - enum(3)
         portnum = enum(3)
         mapcomm(MAXS+2,1) = idproc
         close(unit=2)
         if (idproc.lt.10) then
            write (fname,'(8hMPIERRS0,i1)') idproc
         else
            write (fname,'(7hMPIERRS,i2)') idproc
         endif
         open(unit=2,file=fname,form='formatted',status='unknown')
      endif
c Check if remote portnum agrees with local portnum
      nerr = pnumid(1) - enum(1)
c Send reject flag to remote node
      call MPI_SEND(nerr,1,MPI_INTEGER,nproc-1,4,0,ierror)
c Reject if portnums disagree
      if (nerr.ne.0) then
         write (2,*) 'Session rejected, idproc = ', idproc
         write (2,*) 'Portnums do not agree'
         ierror = 5
         write (2,*) 'remote portnum = ', enum(1)
         call MPI_FINALIZE(nerr)
         return
      endif
c Check for processor number overflow
      if (nproc.gt.MAXS) then
         write (2,*) 'processor number overflow, nproc = ', nproc
         ierror = 6
         call MPI_FINALIZE(nerr)
         return
      endif
c debug
      if (monitor.eq.2) then
         write (2,*) 'connection accepted with idproc = ', nproc-1
      endif
c Accept more connections
      if (idproc.ge.nproc) go to 10
c
c Master prepares to start connection
c
c Find internet address of requested node
   40 nv = lenstr(location)
c Convert a character string into an address
      oss = OTInetStringToHost(location(1:nv)//char(0),address)
      if (oss.ne.0) then
         write (2,*) 'OTInetStringToHost failed, oss = ', oss
         if (oss.eq.kOTBadAddressErr) then
            write (2,*) 'Invalid IP address = ', location(1:nv)
            write (2,*) 'Invalid nodelist file possibly being used'
         endif
         ierror = oss
         call MPI_FINALIZE(nerr)
         return
      endif
c second cpu on a node uses incremented port number
      if (address.eq.oldadd) then
         pshift = pshift + 1
      else
         pshift = 0
         oldadd = address
      endif
      pnumid(1) = portnum + pshift
      addrn(1) = AF_INET
      addrn(2) = pnumid(1)
      longw = address
      addrn(3) = intw(1)
      addrn(4) = intw(2)
c Convert an address into a character string
      call OTInetHostToString(val4(longw),ename)
      nv = lencstr(ename)
c Establish connection end
c Create a new socket for this endpoint
      addrb(2) = kOTAnyInetAddress
c Set TNetbuf for specifying address
      reqad(2) = 16
      reqad(3) = loc(addrb)
      reqad(4) = 0
c Set TNetbuf for address which is returned
      repad(1) = 16
      repad(3) = loc(addrl)
c Initialize synchronous endpoint provider
      oss = otpinit(nproc+1,reqad,repad)
      if (oss.ne.0) then
         ierror = oss
         call MPI_FINALIZE(nerr)
         return
      endif
c Set TNetbuf for specifying address for remote computer
      scall(2) = 16
      scall(3) = loc(addrn)
c Set TNetbuf for specifying options
   70 scall(5) = 0
      scall(6) = 0
c Set TNetbuf for specifying data
      scall(8) = 0
      scall(9) = 0
c debug
      if (monitor.eq.2) then
         write (2,*) 'requesting connection with=',ename(1:nv),pnumid(1)
      endif
c Pause to let OS get some time
      if (checkesc(1)) then
         ierror = -9
         call writerrs('MPI_INIT: ',ierror)
         return
      endif
c Obtain the current time stamp
      call OTGetTimeStamp(ptime)
      call OTGetTimeStamp(delta)
c Request a connection to a remote peer
      oss = OTConnect(val4(epref(nproc+1)),scall,val4(0))
c Wait for connection
   80 if (ioresult(ioc(1,nproc+1)).gt.0) then
         if (checkesc(1)) then
            ierror = -9
            call writerrs('MPI_INIT: ',ierror)
            return
c Wait up to one minute for next connection
         else
            if (OTElapsedMilliseconds(ptime).lt.60000) then
               response = OTGetEndpointState(val4(epref(nproc+1)))
c Try again every second if no response 
               if ((response.lt.5).and.(OTElapsedMilliseconds(delta).gt.
     11000)) then
                  oss = OTConnect(val4(epref(nproc+1)),scall,val4(0))
c Obtain the current time stamp
                  call OTGetTimeStamp(delta)
               endif
               go to 80
            else
               write (2,*) 'OTConnect Wait Time Exceeded'
               write (2,*) 'Trying to start location = ', ename(1:nv),  
     1pnumid(1)
               ierror = -6
               call MPI_FINALIZE(nerr)
               return
            endif
         endif
      else
         oss = ioresult(ioc(1,nproc+1))
         if (oss.ne.0) then
            write (2,*) 'OTConnect Error, oss = ', oss
            write (2,*) 'Trying to start location = ', ename(1:nv), pnum
     1id(1)
            if (oss.eq.kOTOutStateErr) then
               write (2,*) 'Endpoint not in an appropriate state'
               write (2,*) 'Remote node may be unknown or inaccessible'
            endif
            ierror = oss
            call MPI_FINALIZE(nerr)
            return
         endif
      endif
c debug
      if (monitor.eq.2) then
         write (2,*) 'tentative connection started with=', ename(1:nv),
     1pnumid(1)
      endif
c Send portnum for verification and processor id to remote node
      nerr = nproc
      nproc = nproc + 1
      mapcomm(MAXS+1,1) = nproc
      if (nerr.gt.idproc) then
c Set processor id
         pnumid(2) = nerr
         pnumid(3) = portnum
         call MPI_SEND(pnumid,3,MPI_INTEGER,nproc-1,3,0,ierror)
c Read and check reject flag
         call MPI_RECV(nerr,1,MPI_INTEGER,nproc-1,4,0,stat,ierror)
         if (nerr.ne.0) then
            write (2,*) 'Connection rejected, reject info = ', nerr
            write (2,*) 'Portnums do not agree, idproc = ', nproc-1
            ierror = 12
            call MPI_FINALIZE(nerr)
            return
         endif
      endif
c Check for processor number overflow
      if (nproc.gt.MAXS) then
         write (2,*) 'processor number overflow, nproc = ', nproc
         ierror = 6
         call MPI_FINALIZE(nerr)
         return
      endif
c debug
      if (monitor.eq.2) then
         write (2,*) 'connection confirmed with idproc = ', nproc-1
      endif
c Pass current location to next node
      if (nproc.gt.(idproc+2)) then
         call MPI_SEND(location,4,MPI_INTEGER,idproc+1,1,0,nerr)
      endif
c Read location of next node from file
      if (idproc.eq.0) then
         if ((nproc.ge.2).or.(location.eq.'self').or.(location.eq.myself
     1)) then
            if (nofile.ne.0) go to 90
            read (3,'(a16)',iostat=nerr) location
c End of file
            if ((nerr.ne.0).or.(location.eq.' ')) go to 90
         endif
c Receive location of next node from another processor
      else
         call MPI_RECV(location,4,MPI_INTEGER,idproc-1,1,0,stat,nerr)
c End of file marker received
         if (stat(4).eq.0) go to 90
      endif
c Start another connection
      go to 40
c 
c * * * end main iteration loop * * *
c
c All expected nodes activated
   90 nv = nproc - 1
c debug
      if (monitor.eq.2) then
         write (2,*) 'all nodes activated: idproc, nproc=', idproc,nproc
      endif
c Send null record to next processor
      if (idproc.lt.nv) then
         call MPI_SEND(location,0,MPI_INTEGER,idproc+1,1,0,nerr)
      endif
      if (nofile.eq.0) close(unit=3)
c Check number of processors
      if (idproc.eq.nv) then
         do 100 i = 1, nv
         call MPI_SEND(nproc,1,MPI_INTEGER,nv-i,2,0,nerr)
  100    continue
      else
         call MPI_RECV(response,1,MPI_INTEGER,nv,2,0,stat,nerr)
c Local processor does not agree with last processor on total number
         if (response.ne.nproc) then
            write (2,*) 'processor number error, local/remote nproc = '
     1, nproc, response
            ierror = 7
            call MPI_FINALIZE(nerr)
            return
         endif
      endif
c Clear unused MPI_COMM_WORLD mapping
      do 95 i = nproc+1, MAXS
      mapcomm(i,1) = MPI_UNDEFINED
   95 continue
      mapcomm(MAXS+3,1) = 0
c Set MPI_COMM_SELF
      mapcomm(1,2) = idproc
      mapcomm(MAXS+1,2) = 1
      mapcomm(MAXS+2,2) = 0
      mapcomm(MAXS+3,2) = 0
c Create window for showing MPI message status
      if (monitor.gt.0) then
         call messwin(nproc)
         call checkesc(1)
         if (monitor.eq.2) then
            write (2,*) 'MPI_INIT complete'
            write (2,*)
         endif
      endif
c Set error code to success
      ierror = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine regport(portnum,ierror)
c Finds name of current computer by first finding network address of
c current computer, then finding a registered entity with type name
c 'PPCToolBox' which resides at that address.  The object name
c corresponding to that entity is the current computer name.  This only
c works if the user has set Program Linking on in the File Sharing
c control panel.  If computer name is found, then a portname is created
c whose entity name consists of object = computer name, type = portnum
c If it already exists, it is first deleted, then reregistered.
      implicit none
      integer portnum, ierror
c declare adsp common block
      integer nameid
c nameid = adsp portname ID
      common /adsp/ nameid
c function declarations
      integer OTOpenAppleTalkServicesInContext, OTATalkGetInfo
      integer OTCloseProvider
      integer OTCreateConfiguration, OTOpenMapperInContext, OTIoctl
      integer OTDestroyConfiguration, OTSetBlocking, OTLookupName
      integer OTGetNBPEntityLengthAsAddress
      integer OTDeleteName, OTRegisterName
      integer lenstr, instr
      external OTOpenAppleTalkServicesInContext, OTATalkGetInfo
      external OTCloseProvider
      external OTCreateConfiguration, OTOpenMapperInContext, OTIoctl
      external OTDestroyConfiguration, OTSetBlocking, OTLookupName
      external OTGetNBPEntityLengthAsAddress
      external OTDeleteName, OTRegisterName
      external lenstr, instr
c OT constants
      integer kDefaultAppleTalkServicesPath
      integer kATalkFullSelfSend, kOTNoDataErr, kOTNoAddressErr
      parameter(kDefaultAppleTalkServicesPath=-3)
      parameter(kATalkFullSelfSend=21551,kOTNoDataErr=-3162)
      parameter(kOTNoAddressErr=-3154)
c local data
      integer*2 intw(2), atinfo(11), addrb(4), addrn(54)
      integer longw, aref, oss, address, nodeid, cfig, mref, nv, i, club
      integer ucptr
      integer info(3), req(9), rep(4), aname(25), treq(7), trep(4)
      character selfsend, rbuff(1888)
      character*36 ename
c longw and intw are used to convert types between integer*2 and integer
      equivalence(longw,intw)
      data ucptr, cfig, mref, oss /0,0,0,0/
      ierror = 0
      selfsend = char(1)
      aref = 0
c Open a synchronous AppleTalk service provider
c     aref = OTOpenAppleTalkServicesInContext(val4(kDefaultAppleTalkServ
c    1icesPath),val4(0),oss,val4(0))
      if (aref.eq.0) then
         write (2,*) 'Invalid ATS, oss = ', oss
         ierror = oss
         return
      endif
c Set TNetbuf to receive AppleTalk Info
      info(1) = 22
      info(2) = 0
      info(3) = loc(atinfo)
c Obtain information about the AppleTalk environment
      oss = OTATalkGetInfo(val4(aref),info)
      if (oss.ne.0) then
         write (2,*) 'ATS GetInfo Error, oss = ', oss
         ierror = oss
         call OTCloseProvider(val4(aref))
         return
      endif
      address = atinfo(2)
      intw(1) = 0
      intw(2) = atinfo(3)
      nodeid = longw/256
c Close a provider of any type
      oss = OTCloseProvider(val4(aref))
      if (oss.ne.0) write (2,*) 'Close ATS Provider Error, oss = ', oss
c Create a structure defining a provider's configuration
      cfig = OTCreateConfiguration('nbp'//char(0))
      if (cfig.eq.0) then
         write (2,*) 'Insufficient Memory for nbp Mapper'
         ierror = -4
         return
      elseif (cfig.eq.(-1)) then
         write (2,*) 'Invalid Configuration for nbp Mapper'
         ierror = -4
         return
      endif
c Create a synchronous mapper provider
      mref = OTOpenMapperInContext(val4(cfig),val4(0),oss,val4(0))
      if (mref.eq.0) then
         write (2,*) 'Invalid NBP Mapper, oss = ', oss
         ierror = oss
         call OTDestroyConfiguration(val4(cfig))
         return
      endif
c Enable ADSP SelfSend
c Send a module-specific command to Open Transport protocol module
      oss = OTIoctl(val4(mref),val4(kATalkFullSelfSend),val1(ichar(selfs
     1end)))
      if (oss.lt.0) then
         write (2,*) 'OTIoctl SelfSend Error, oss = ', oss
c Save old value
      else
         selfsend = char(oss)
      endif
c Set a provider to wait or block until function can complete
      oss = OTSetBlocking(val4(mref))
      if (oss.ne.0) then
         write (2,*) 'NBP Mapper Set Blocking Error, oss = ', oss
         ierror = oss
         call OTCloseProvider(val4(mref))
         return
      endif
c Look for names with type 'PPCToolBox'
      ename = '=:PPCToolBox@*'//char(0)
c Set TNetbuf for name to be looked up
      req(2) = 14
      req(3) = loc(ename)
c Set TNetbuf to use defaults for address search
      req(5) = 0
      req(6) = 0
c Set number of names you expect to be returned
      req(7) = 16
c Set timeout to default
      req(8) = 0
      req(9) = 0
c Set TNetbuf for names which are found
      rep(1) = 1888
      rep(3) = loc(rbuff)
c Find all addresses that correspond to a particular name or pattern
      oss = OTLookupName(val4(mref),req,rep)
      if (oss.ne.0) then
         write (2,*) 'OTLookup Error, rspcount, oss = ', rep(4), oss
         if (oss.eq.kOTNoDataErr) then
            write (2,*) 'Cannot find PPCToolBox on any computer'
         endif
         ierror = oss
      else
         nv = 1
         i = 1
    5    club = loc(rbuff(nv))
c Get network address of found entity
         addrb(2) = word(club+6)
         intw(1) = 0
         intw(2) = word(club+8)
c Check if network and nodeid agrees
         if ((addrb(2).ne.address).or.((longw/256).ne.nodeid)) then
            i = i + 1
            if (i.gt.rep(4)) then
               write (2,*) 'Cannot find PPCToolBox on current computer'
               ierror = 8
c Look up more addresses
            else
               nv = nv + ((word(club)-1)/4 + (word(club+2)-1)/4 + 3)*4
               go to 5
            endif
c Found correct address
         else
            ierror = 0
         endif
      endif
      if (ierror.ne.0) then
         write (2,*) 'Program Linking may not be enabled'
         call OTCloseProvider(val4(mref))
         return
      endif
c Extract object name
      nv = word(club)
      ucptr = club + 4 + nv
      nv = word(club+2)
      call OTSetNBPEntityFromAddress(aname,val4(ucptr),val4(nv))
c Check if port name is already being used
      write (ename,'(i5)') portnum
      nv = lenstr(ename) + 1
      ename(nv:nv) = char(0)
      nv = instr(ename)
c Set the type part of an NBP entity structure
      call OTSetNBPType(aname,ename(nv:))
c Set the zone part of an NBP entity structure
      call OTSetNBPZone(aname,'*'//char(0))
c Store an NBP entity structure as an NBP address string
      call OTSetAddressFromNBPEntity(addrn(2),aname)
      addrn(1) = 258
c Obtain the size of an NBP entity structure formatted as string
      nv = OTGetNBPEntityLengthAsAddress(aname)
      i = nv/2
      if (nv.gt.(2*i)) addrn(i+2) = 256*(addrn(i+2)/256)
c Set TNetbuf for name to be deleted
      info(2) = nv
      info(3) = loc(addrn(2))
c Remove a previously registered entity name
      oss = OTDeleteName(val4(mref),info)
      if (oss.eq.0) then
         write (2,'(1x,a22,16a2)') 'closed old adsp port: '
     1, (addrn(i+1),i=1,(nv-1)/2+1)
      endif
c Set TNetbufs for name to be registered
      treq(2) = nv
      treq(3) = loc(addrn(2))
      treq(5) = 0
      treq(6) = 0
      treq(7) = 0
      trep(1) = 8
      trep(3) = loc(addrb)
      trep(4) = 0
c Register a name on the network
      oss = OTRegisterName(val4(mref),treq,trep)
      nameid = trep(4)
      if (oss.ne.0) then
         write (2,*) 'OTRegisterName Error, oss = ', oss
      endif
c Close a provider of any type
      oss = OTCloseProvider(val4(mref))
      if (oss.ne.0) write (2,*) 'Close NBP Provider Error, oss=', oss
      return
      end
c-----------------------------------------------------------------------
      subroutine delport()
c this subroutine deletes adsp portname
      implicit none
c declare adsp common block
      integer nameid
c nameid = adsp portname ID
      common /adsp/ nameid
c function declarations
      integer OTCreateConfiguration, OTOpenMapperInContext
      integer OTSetBlocking, OTDeleteNameByID, OTCloseProvider
      external OTCreateConfiguration, OTOpenMapperInContext
      external OTSetBlocking, OTDeleteNameByID, OTCloseProvider
c local data
      integer cfig, oss, mref
      data cfig, mref /0,0/
c Create a structure defining a provider's configuration
      cfig = OTCreateConfiguration('nbp'//char(0))
      if (cfig.eq.0) then
         write (2,*) 'delport: Insufficient Memory for nbp Mapper'
         return
      elseif (cfig.eq.(-1)) then
         write (2,*) 'delport: Invalid Configuration for nbp Mapper'
         return
      endif
c Create a synchronous mapper provider
      mref = OTOpenMapperInContext(val4(cfig),val4(0),oss,val4(0))
      if (mref.eq.0) then
         write (2,*) 'delport: Invalid NBP Mapper, oss = ', oss
         call OTDestroyConfiguration(val4(cfig))
         return
      endif
c Set a provider to wait or block until function can complete
      oss = OTSetBlocking(val4(mref))
      if (oss.ne.0) then
         write (2,*) 'delport: NBP Mapper Set Blocking Error, oss=', oss
         call OTCloseProvider(val4(mref))
         return
      endif
c Remove a previously registered name as specified by its name ID
      oss = OTDeleteNameByID(val4(mref),val4(nameid))
      if (oss.ne.0) then
         write (2,*) 'OTDeleteNameByID Error, oss = ', oss
      endif
      oss = OTCloseProvider(val4(mref))
      if (oss.ne.0) then
         write (2,*) 'delport: Close NBP Provider Error, oss=', oss
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine notifier(context,event,result,cookie)
c notifier function for asynchronous and completion events
      implicit none
      integer context, event, result, cookie
      value context, event, result, cookie
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c epref = array of endpoint references for each participating node
c nevents = log of unknown notifier events
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer MAXM, curreq, header, rwrec, monitor, trash, mqueue
      parameter(MAXM=2*MAXS)
      dimension curreq(5,MAXM), header(10,MAXM), rwrec(14,MAXM)
      dimension trash(256), mqueue(2,MAXS+1)
c rwrec = read/write record for asynchronous messages
c trash = trash bin for unwanted data
c mqueue = message request queue
      common /mpisendrec/ monitor, curreq, header, rwrec, trash, mqueue
c function declarations
      integer OTListen, OTAccept, OTRcvConnect, OTRcvDisconnect
      integer OTRcvOrderlyDisconnect, OTGetEndpointState, OTIoctl
      integer OTSndOrderlyDisconnect, OTRcv, OTLook, OTSnd
      external OTListen, OTAccept, OTRcvConnect, OTRcvDisconnect
      external OTRcvOrderlyDisconnect, OTGetEndpointState, OTIoctl
      external OTSndOrderlyDisconnect, OTRcv, OTLook, OTSnd
c OT constants
      integer T_LISTEN, T_CONNECT, T_DISCONNECT, T_DISCONNECTCOMPLETE
      integer T_ACCEPTCOMPLETE, T_PASSCON, T_ORDREL, kStreamIoctlEvent
      integer T_DATA, T_GODATA, T_INREL, I_FLUSH, FLUSHRW, kOTNoDataErr
      integer kOTLookErr, kOTFlowErr, kENOMEMErr
      parameter(T_LISTEN=1,T_CONNECT=2,T_DISCONNECT=16)
      parameter(T_DISCONNECTCOMPLETE=536870917)
      parameter(T_ACCEPTCOMPLETE=536870915,T_PASSCON=4096,T_ORDREL=128)
      parameter(kStreamIoctlEvent=553648133,T_DATA=4,T_GODATA=256)
      parameter(T_INREL=7,I_FLUSH=16645,FLUSHRW=3,kOTNoDataErr=-3162)
      parameter(kOTLookErr=-3158,kOTFlowErr=-3161,kENOMEMErr=-3211)
c local data
      integer n, m, oss, tryagn, addrb(4), call(10)
      tryagn = 0
c Get reference to endpoint
      n = long(context)
c Check event
   10 select case (event)
         case (T_LISTEN)
c Set TNetbuf for address which is returned
            call(1) = 16
            call(3) = loc(addrb)
c Set TNetbuf for options which are returned
            call(4) = 0
            call(6) = 0
c Set TNetbuf for data which is returned
            call(7) = 0
            call(9) = 0
c Listen for an incoming connection request
            oss = OTListen(val4(epref(n)),call)
            if (oss.ne.0) then
c Set completion flag to error
               long(context+4) = oss
            else 
               m = long(context+4) + 1
c Accept an incoming connection request
               oss = OTAccept(val4(epref(n)),val4(epref(m-1)),call)
            endif
            return
         case (T_CONNECT)
c Read the status of an asynchronous call to OTConnect
            oss = OTRcvConnect(val4(epref(n)),val4(0))
c Set completion flag
            long(context+4) = oss
            return
         case (T_DISCONNECT)
c Identify cause of connection rejection and clear event
            oss = OTRcvDisconnect(val4(epref(n)),val4(0))
c Set completion flag
c           long(context+4) = oss
c Read pointer to current read record
            m = long(context+8)
c Set iocompletion flag, if there is a pending read
            if ((m.gt.0).and.(m.le.MAXM)) rwrec(2,m) = -2
c Read pointer to current send record
            m = long(context+12)
c Set iocompletion flag, if there is a pending send
            if ((m.gt.0).and.(m.le.MAXM)) rwrec(2,m) = -2
            return
         case (T_DISCONNECTCOMPLETE)
c Set completion flag
c           long(context+4) = oss
            return
         case (T_ACCEPTCOMPLETE)
            long(context+4) = result
            return
         case (T_PASSCON)
c Set completion flag
            long(context+4) = 0
            return
         case (T_ORDREL)
c Reset event to complete OTLookErr handling
            if (tryagn.ne.0) event = tryagn
c Clear an incoming orderly disconnect event
            oss = OTRcvOrderlyDisconnect(val4(epref(n)))
c Set completion flag
c           long(context+4) = oss
c Read pointer to current read record
            m = long(context+8)
c Set iocompletion flag, if there is a pending read
            if ((m.gt.0).and.(m.le.MAXM)) rwrec(2,m) = -1
c Read pointer to current send record
            m = long(context+12)
c Set iocompletion flag, if there is a pending send
            if ((m.gt.0).and.(m.le.MAXM)) rwrec(2,m) = -1
c Obtain the current state of an endpoint
            oss = OTGetEndpointState(val4(epref(n)))
            if (oss.eq.T_INREL) then
c Send a module-specific command to Open Transport protocol module
               oss = OTIoctl(val4(epref(n)),val4(I_FLUSH),val4(FLUSHRW))
            endif
            if (tryagn.eq.0) return
        case (kStreamIoctlEvent)
c Initiate an orderly disconnect
            oss = OTSndOrderlyDisconnect(val4(epref(n)))
c Set completion flag
            long(context+4) = oss
            return
         case (T_DATA)
            tryagn = 0
c Read pointer to current read record
            m = long(context+8)
c Read data sent from a remote peer
            if ((m.gt.0).and.(m.le.MAXM)) then
c Obtain the current time stamp
               call OTGetTimeStamp(rwrec(11,m))
   20          oss = OTRcv(val4(rwrec(1,m)),val4(rwrec(3,m)),val4(rwrec(
     14,m)),rwrec(5,m))
c Unexpected flag returned
               if (rwrec(5,m).gt.1) then
                  rwrec(2,m) = -3
c Clear pointer to current read record
                  long(context+8) = 0
c Process data which arrived
               elseif (oss.ge.0) then
c Clear more flag
                  rwrec(5,m) = 0
c Clear non-fatal error code
                  rwrec(14,m) = 0
c Set actual length received
                  rwrec(8,m) = rwrec(8,m) + oss
c Check if all the data has arrived
                  if (rwrec(8,m).lt.header(10,m)) then
c Incomplete message
                    if (rwrec(4,m).gt.oss) then
c Readjust buffer pointer
                        rwrec(3,m) = rwrec(3,m) + oss
                        rwrec(4,m) = rwrec(4,m) - oss
                        if (rwrec(8,m).ge.0) rwrec(2,m) = rwrec(2,m) + 1
c Header is received, readjust parameters to receive data
                     elseif (rwrec(2,m).eq.1) then
                        rwrec(3,m) = rwrec(6,m)
                        rwrec(4,m) = min(header(10,m),rwrec(7,m))
                        rwrec(2,m) = 2
c Data is received and buffer is full
                     else
                        rwrec(3,m) = loc(trash)
                        rwrec(4,m) = 1024
                     endif
                     go to 20
c Message complete
                  else
c Obtain the current time stamp
                     call OTGetTimeStamp(rwrec(11,m))
c Set iocompletion flag
                     rwrec(2,m) = 0
c Get next message if messages are queued
                     if (rwrec(13,m).gt.0) then
                        m = rwrec(13,m)
                        if (m.eq.mqueue(1,n)) mqueue(1,n) = 0
                        long(context+8) = m
                        go to 20
c Clear pointer to current read record
                     else
                        long(context+8) = 0
                     endif
                  endif
c Check for errors
               else
c Quit if no data is available
                  if (oss.eq.kOTNoDataErr) then
                     return
c Determine cause of a kOTLookErr
                  elseif (oss.eq.kOTLookErr) then
                     event = OTLook(val4(rwrec(1,m)))
                     if ((event.eq.T_GODATA).or.(event.eq.T_ORDREL)) the
     1n
                        tryagn = T_DATA
c Store unprocessed event returned by OTLookErr
                     else
                        nevents(n) = event
                     endif
c Store non-fatal error code
                  elseif (oss.eq.kENOMEMErr) then
                     rwrec(14,m) = oss
                     return
                  endif
c Set iocompletion flag to error
                  if (tryagn.eq.0) then
                     rwrec(2,m) = oss
                     long(context+8) = 0
                  endif
               endif
            endif
            if (tryagn.eq.0) return
         case (T_GODATA)
c Read pointer to current send record
            m = long(context+12)
c Reset event to complete OTLookErr handling
            if (tryagn.ne.0) event = tryagn
c Send data to a remote peer
            if ((m.gt.0).and.(m.le.MAXM)) then
c Obtain the current time stamp
               call OTGetTimeStamp(rwrec(11,m))
   30          oss = OTSnd(val4(rwrec(1,m)),val4(rwrec(3,m)),val4(rwrec(
     14,m)),val4(rwrec(5,m)))
c Process data which has been sent
               if (oss.ge.0) then
c Clear non-fatal error code
                  rwrec(14,m) = 0
c Set actual length sent
                  rwrec(8,m) = rwrec(8,m) + oss
c Check for incomplete header
                  if (rwrec(8,m).lt.0) then
                     header(2,m) = header(2,m) + oss
                     header(3,m) = header(3,m) - oss
                     go to 30
c Check for incomplete data
                  elseif (rwrec(7,m).gt.rwrec(8,m)) then
c Header is sent, readjust parameters to send data
                     if (rwrec(2,m).eq.1) then
                        rwrec(3,m) = rwrec(6,m)
                        rwrec(4,m) = rwrec(7,m)
                        oss = oss - header(3,m)
                     endif
c Readjust buffer pointer
                     rwrec(3,m) = rwrec(3,m) + oss
                     rwrec(4,m) = rwrec(4,m) - oss
                     rwrec(2,m) = rwrec(2,m) + 1
                     go to 30
c Data is sent
                  else
c Obtain the current time stamp
                     call OTGetTimeStamp(rwrec(11,m))
c Set iocompletion flag
                     rwrec(2,m) = 0
c Get next message if messages are queued
                     if (rwrec(13,m).gt.0) then
                        m = rwrec(13,m)
                        if (m.eq.mqueue(2,n)) mqueue(2,n) = 0
                        long(context+12) = m
                        go to 30
c Clear pointer to current send record
                     else
                        long(context+12) = 0
                     endif
                  endif
c Check for errors
               else
c Determine cause of a kOTLookErr
                  if (oss.eq.kOTLookErr) then
                     event = OTLook(val4(rwrec(1,m)))
c Store unprocessed event returned by OTLookErr
                     nevents(n) = event
c Store non-fatal error code
                  elseif (oss.eq.kENOMEMErr) then
                     rwrec(14,m) = oss
                  endif
c Set iocompletion flag to error, unless no data can be sent
                  if ((oss.ne.kOTFlowErr).and.(oss.ne.kENOMEMErr)) then
                     rwrec(2,m) = oss
                     long(context+12) = 0
                  endif
               endif
            endif
         case default
c unknown event
            nevents(n) = event
            return
      end select
      if (tryagn.ne.0) go to 10
      return
      end
c-----------------------------------------------------------------------
      function lencstr(chr)
c this subroutine calculates the length of a C style string,
c chr = input characters
      implicit none
      integer lencstr, lc
      character*(*) chr
      lencstr = len(chr)
      if (lencstr.lt.1) return
      lc = 1
c find NULL
   10 if (chr(lc:lc).ne.char(0))  then
         lc = lc + 1
         if (lc.le.lencstr) go to 10
         lc = lc - 1
      endif
      lencstr = lc - 1
      return
      end
c-----------------------------------------------------------------------
      function lenstr(chr)
c this subroutine calculates the length of a Fortran style string,
c not counting trailing blanks
c chr = input characters
      implicit none
      integer lenstr
      character*(*) chr
      lenstr = len(chr)
      if (lenstr.lt.1) return
c omit trailing blanks
   10 if (chr(lenstr:lenstr).eq.' ')  then
         lenstr = lenstr - 1
         if (lenstr.gt.0) go to 10
      endif
      return
      end
c-----------------------------------------------------------------------
      function instr(chr)
c this subroutine calculates the start of a Fortran style string,
c not counting leading blanks
c chr = input characters
      implicit none
      integer instr, lenstr
      character*(*) chr
      instr = 1
      lenstr = len(chr)
      if (lenstr.lt.1) return
c omit leading blanks
   10 if (chr(instr:instr).eq.' ')  then
         instr = instr + 1
         if (instr.le.lenstr) go to 10
      endif
c special case of all blanks
      if (instr.gt.lenstr) instr = 1
      return
      end
c-----------------------------------------------------------------------
      function otpinit(np,reqad,repad)
c this subroutine initializes an Open Transport provider for
c index np, using address specified in reqad.  The cfig0
c OTConfiguration template is assumed to be already created
c provider is left in asynchronous, blocking mode
c np = index to endpoint reference array
c reqad = address to which endpoint is to be bound
c repad = address to which endpoint is actually bound
c returns OSStatus indicator
c input: np, reqad
      implicit none
      integer otpinit
      integer np, reqad(4), repad(4)
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      integer notifierUPP
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c cfig0 = OTConfiguration structure template
c epref = array of endpoint references for each participating node
c loc = array of context pointers for notifier function 
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm, notifierUPP
c function declarations
      integer OTCloneConfiguration, OTOpenEndpointInContext
      integer OTSetBlocking, OTInstallNotifier, OTBind
      integer OTSetAsynchronous, setopts
      external OTCloneConfiguration, OTOpenEndpointInContext
      external OTSetBlocking, OTInstallNotifier, OTBind
      external OTSetAsynchronous, setopts
c OT constants
      integer INET_TCP, TCP_NODELAY, kOTNoAddressErr
      parameter(INET_TCP=6,TCP_NODELAY=1,kOTNoAddressErr=-3154)
c local data
      integer i, oss, value, cfig, epinfo(8), laddrn
      otpinit = 0
c Copy an OTConfiguration structure
      cfig = OTCloneConfiguration(val4(cfig0))
c Create a synchronous endpoint provider
      epref(np) = OTOpenEndpointInContext(val4(cfig),val4(0),epinfo,oss,
     1val4(0))
      if (epref(np).eq.0) then
         write (2,*) 'Invalid Endpoint, oss = ', oss
         call OTDestroyConfiguration(val4(cfig))
         otpinit = oss
         return
      endif
c Set a provider to wait or block until function can complete
      oss = OTSetBlocking(val4(epref(np)))
      if (oss.ne.0) then
         write (2,*) 'Endpoint Set Blocking Error, oss = ', oss
         otpinit = oss
         return
      endif
c set TCP_NODELAY option
      value = 1
      oss = setopts(epref(np),INET_TCP,TCP_NODELAY,value)
      if (oss.ne.0) then
         write (2,*) 'TCP_NODELAY OTOptionManagement Error, oss = ', oss
      elseif (value.eq.0) then
         write (2,*) 'Info: TCP_NODELAY was not set'
      endif
c Install a notifier function
      oss = OTInstallNotifier(val4(epref(np)),val4(notifierUPP),ioc(1,np
     1))
      if (oss.ne.0) then
         write (2,*) 'OTInstall Notifier Error, oss = ', oss
         otpinit = oss
         return
      endif
      i = 1
      laddrn = reqad(3)
c Assign an address to an endpoint
   10 oss = OTBind(val4(epref(np)),reqad,repad)
      if (oss.ne.0) then
         if ((word(laddrn+2).ne.0).and.(i.le.16)) then
            word(laddrn+2) = word(laddrn+2) + 1
            i = i + 1
            go to 10
         endif
         write (2,*) 'OTBind Error, oss = ', oss
         if (oss.eq.kOTNoAddressErr) then
            write (2,*) 'unable to allocate address'
         endif
         otpinit = oss
         return
      endif
c Pass processor epref index and iocompletion to notifier
      ioc(1,np) = np
      ioc(2,np) = 1
c Clear read and send record pointers
      ioc(3,np) = 0
      ioc(4,np) = 0
c Set a provider to asynchronous mode
      oss = OTSetAsynchronous(val4(epref(np)))
      otpinit = oss
      return
      end
c-----------------------------------------------------------------------
      function setopts(eref,level,name,value)
c This function sets value of option
      implicit none
      integer setopts
      integer eref, level, name, value
c function declarations
      integer OTOptionManagement
      external OTOptionManagement
c OT constants
      integer T_NEGOTIATE, T_SUCCESS, T_NOTSUPPORT
      integer kOTBufferOverflowErr
      parameter(T_NEGOTIATE=4,T_SUCCESS=32,T_NOTSUPPORT=1024)
      parameter(kOTBufferOverflowErr=-3160)
c local data
      integer oss
      integer option(5)
      integer req(4), ret(4)
c Set option structure
      option(1) = 20
      option(2) = level
      option(3) = name
      option(4) = 0
      option(5) = value
c Set OptMgmt request structure
      req(2) = 20
      req(3) = loc(option)
      req(4) = T_NEGOTIATE
c Set OptMgmt reply structure
      ret(1) = 20
      ret(3) = loc(option)
c Determine an endpoint's option values or change the values
      oss = OTOptionManagement(val4(eref),req,ret)
      if (oss.eq.0) then
         if (option(4).eq.T_SUCCESS) then
            value = option(5)
         else
            oss = option(4)
            if (oss.eq.T_NOTSUPPORT) then
               write (2,*) 'Option Not Supported'
            endif
         endif
      else
         if (oss.eq.kOTBufferOverflowErr) then
            write (2,*) 'Option reply buffer is too small'
         endif
      endif
      setopts = oss
      return
      end
c-----------------------------------------------------------------------
      function ioresult(pblock)
c this function returns ioResult for asynchronous procedures
c input: pblock
      implicit none
      integer ioresult
      integer pblock(*)
      ioresult = pblock(2)
      return
      end
c-----------------------------------------------------------------------
      logical function checkesc(stk)
c this procedure allows user to abort a procedure by checking for
c escape, Cmd-. or Ctrl-C keystrokes.  Calling EventAvail also
c permits an idle procedure to time-share and checks for Quit Events
c returns true if an escape event occurred.
c recent keyboard events not processed are not removed from event queue
c stk = maximum number of sleepTicks (sixtieths of a second) that
c application agrees to relinquish the processor if no events are 
c pending for it.
c input: stk
      implicit none
      integer stk
c function declarations
      integer*2 EventAvail, GetNextEvent, WaitNextEvent, FindWindow
      integer GetMainDevice, TickCount, OTElapsedMilliseconds
      external EventAvail, GetNextEvent, WaitNextEvent, FindWindow
      external GetMainDevice, TickCount, OTElapsedMilliseconds
c common block for message window
c cpptr = pointer to window structure
c crect = current drag region
      integer cpptr
      integer*2 crect(4), nsp, nbx, nds, mbs
      common /winmess/ cpptr, crect, nsp, nbx, nds, mbs
c MPI constants
      integer MPI_COMM_WORLD
      parameter(MPI_COMM_WORLD=0)
c local data
      integer*2 mouseDown, keyDown, autoKey, updateEvt, kHighLevelEvent
      integer*2 osEvt, resus
c myEventMask looks for mouse, keyboard, update, and quit events
      parameter(mouseDown=1,keyDown=3,autoKey=5)
      parameter(updateEvt=6,kHighLevelEvent=23)
      parameter(osEvt=15,resus=256)
      integer*2 myEventMask
      integer*2 event(8), intw(2)
      integer key, longw, nvp, times(2)
c longw and intw are used to convert types between integer*2 and integer
      equivalence(longw,intw)
      save times
      data times /2*0/
      checkesc = .false.
c myEventMask looks for mouse, keyboard, and quit events
      myEventMask = 1086
c if monitor window is open, look for update events also
      if (cpptr.ne.0) myEventMask = myEventMask + 64
c get an event without removing it from the queue
      if (EventAvail(val2(myEventMask),event)) then
         if ((event(1).eq.keyDown).or.(event(1).eq.autoKey)) then
c check for escape key
            key = event(3) - 256*(event(3)/256)
            if (key.eq.27) then
               checkesc = .true.
c check for Cmd-.
            elseif (key.eq.46) then
               if ((event(8)/256).ne.(2*(event(8)/512))) then
                  checkesc = .true.
               endif
c check for Ctrl-C
            elseif (key.eq.3) then
               checkesc = .true.
            endif
c remove processed or keyboard event more than 3 seconds old
            intw(1) = event(4)
            intw(2) = event(5)
            if ((checkesc).or.((TickCount()-longw).gt.180)) then
c remove next available event
               call GetNextEvent(val2(myEventMask),event)
            endif
c check for 'QuitApplication' Apple Event
         elseif (event(1).eq.kHighLevelEvent) then
c remove high level events more than 5 seconds old
            intw(1) = event(4)
            intw(2) = event(5)
            if ((TickCount()-longw).gt.300) then
c remove next available event
               call GetNextEvent(val2(myEventMask),event)
            endif
c check if event(6) = 'qu' and event(7) = 'it'
            if ((event(6).eq.29045).and.(event(7).eq.26996)) then
               write (2,*) 'Quit Application Apple Event received'
               checkesc = .true.
            endif
c check for update events
         elseif (event(1).eq.updateEvt) then
c get window pointer
            intw(1) = event(2)
            intw(2) = event(3)
            if (cpptr.eq.longw) then
c remove next available event
               call GetNextEvent(val2(myEventMask),event)
c signal start of window update
               call BeginUpdate(val4(cpptr))
               call MPI_COMM_SIZE(MPI_COMM_WORLD,nvp,key)
               call messwin(nvp)
c signal end of update after BeginUpdate
               call EndUpdate(val4(cpptr))
            endif
c check for drag window event
         elseif (event(1).eq.mouseDown) then
c longw = position in global coordinates where mouse event took place
            intw(1) = event(6)
            intw(2) = event(7)
c see which window part, including menu bar, is at a point
            if (FindWindow(val4(longw),key).eq.4) then
               if (cpptr.eq.key) then
c remove next available event
                  call GetNextEvent(val2(myEventMask),event)
c track the mouse and move a window
                  call DragWindow(val4(cpptr),val4(longw),crect)
               endif
            endif
         endif
      endif
c yield time when time slice since last yield > 1 second
      if (stk.eq.0) then
c calculate time elapsed in milliseconds
         if (OTElapsedMilliseconds(times).gt.1000) then
c receive null event from event manager to relinquish the processor
            call WaitNextEvent(val2(0),event,val4(1),val4(0))
c Obtain the current time stamp
            call OTGetTimeStamp(times)
         endif
c yield time if requested
      else
c receive null event from event manager to relinquish the processor
         call WaitNextEvent(val2(0),event,val4(stk),val4(0))
c Obtain the current time stamp
         call OTGetTimeStamp(times)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_FINALIZE(ierror)
c terminate MPI execution environment
c ierror = error indicator
c output: ierror
      implicit none
      integer ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, LDM, epref, ioc, nevents, stime, mapcomm
      integer notifierUPP
      parameter(MAXS=32,MAXC=10,MAXD=6,LDM=0)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c cfig0 = OTConfiguration structure template
c epref = array of endpoint references for each participating node
c nevents = log of unknown notifier events
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm, notifierUPP
c declare common block for non-blocking messages
      integer MAXM, curreq, header, rwrec, monitor, trash, mqueue
      parameter(MAXM=2*MAXS)
      dimension curreq(5,MAXM), header(10,MAXM), rwrec(14,MAXM)
      dimension trash(256), mqueue(2,MAXS+1)
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
c curreq = request record for transmission parameters
      common /mpisendrec/ monitor, curreq, header, rwrec, trash, mqueue
c function declarations
      integer OTGetEndpointState, OTIoctl, OTRcvOrderlyDisconnect
      integer OTSetSynchronous, OTSndDisconnect, OTUnbind
      integer OTCloseProvider, OTElapsedMilliseconds
      logical checkesc
      external OTGetEndpointState, OTIoctl, OTRcvOrderlyDisconnect
      external OTSetSynchronous, OTSndDisconnect, OTUnbind
      external OTCloseProvider, OTElapsedMilliseconds
      external checkesc
c OT constants
      integer T_IDLE, T_DATAXFER, I_FLUSH, FLUSHRW, kOTOutStateErr
      parameter(T_IDLE=2,T_DATAXFER=5,I_FLUSH=16645,FLUSHRW=3)
      parameter(kOTOutStateErr=-3155)
c local data
      integer i, state, oss, ptime(2)
      ierror = 0
c MPI already finalized
      if (nproc.lt.0) then
         ierror = 1
         return
      endif
      if (monitor.eq.2) write (2,*) 'MPI_FINALIZE started'
c Dispose of an OTConfiguration structure
      call OTDestroyConfiguration(val4(cfig0))
c Flush endpoint providers
      do 10 i = 1, MAXS+1
      if (epref(i).ne.0) then
c Obtain the current state of an endpoint
         state = OTGetEndpointState(val4(epref(i)))
         if (state.eq.T_DATAXFER) then
c Send a module-specific command to Open Transport protocol module
            oss = OTIoctl(val4(epref(i)),val4(I_FLUSH),val4(FLUSHRW))
         endif
      endif
   10 continue
c Obtain the current time stamp
      call OTGetTimeStamp(ptime)
c Disconnect endpoints with orderly disconnect
      do 30 i = 1, MAXS+1
         if (epref(i).ne.0) then
c Wait 5 seconds for endpoints to acknowledge disconnect
   20       state = OTGetEndpointState(val4(epref(i)))
            if (state.gt.T_IDLE) then
c Calculate time elapsed in milliseconds
               if (OTElapsedMilliseconds(ptime).lt.5000) then
                  go to 20
               else
                  write (2,*) 'OrderlyDisconnect Timeout: i, state = '  
     1, i, state
c Tear down an open connection (abortive disconnect)
                  oss = OTSndDisconnect(val4(epref(i)),val4(0))
                  if (oss.ne.0) then
                     write (2,*) 'OTSndDisconnect Error, i, oss=', i,oss
                     if (oss.eq.kOTOutStateErr) then
                        write (2,*) 'Endpoint not in appropriate state='
     1, OTGetEndpointState(val4(epref(i)))
                     endif
                  endif
               endif
            endif
         endif
   30 continue
c Obtain the current time stamp
      call OTGetTimeStamp(ptime)
c Unbind and close providers
      do 50 i = 1, MAXS+1
         if (epref(i).ne.0) then
c Wait another 5 seconds for endpoints to acknowledge disconnect
   40       state = OTGetEndpointState(val4(epref(i)))
            if (state.gt.T_IDLE) then
c Calculate time elapsed in milliseconds
               if (OTElapsedMilliseconds(ptime).lt.5000) then
                  go to 40
               else
                  write (2,*) 'Disconnect Timeout: i, state=', i, state
               endif
            endif
c Set a provider to synchronous mode
            oss = OTSetSynchronous(val4(epref(i)))
c Dissociate an endpoint from its address
            oss = OTUnbind(val4(epref(i)))
            if (oss.ne.0) then
               write (2,*) 'Unbind Endpoint Error, i, oss = ', i, oss
               if (oss.eq.kOTOutStateErr) then
                     write (2,*) 'Endpoint not in an appropriate state='
     1, OTGetEndpointState(val4(epref(i)))
               endif
            endif
c Close a provider of any type
            oss = OTCloseProvider(val4(epref(i)))
            if (oss.ne.0) then
               write (2,*) 'Close Endpoint Error, i, oss = ', i, oss
               ierror = oss
            endif
         endif
   50 continue
c Dispose of Univeral Proc Ptr
      call DisposeOTNotifyUPP(val4(notifierUPP))
      notifierUPP = 0
      if (LDM.ne.0) call delport()
c Unregister your application connection to Open Transport
      call CloseOpenTransportInContext(val4(0))
c Write out any unknown notifier events
      do 60 i = 1, MAXS+1
      if (nevents(i).ne.0) then
         write (2,*) 'Unknown notifier event, endpoint= ', nevents(i), i
      endif
   60 continue
c Close window for showing MPI message status
      if (monitor.gt.0) then
         call logmess(0,0,-1,0,0)
         call delmess()
      endif
c Nullify nproc
      nproc = -1
c Nullify communicator mappings
      do 70 i = 1, MAXC
      mapcomm(MAXS+1,i) = 0
   70 continue
c Nullify endpoint references and unknown notifier events
      do 80 i = 1, MAXS+1
      epref(i) = 0
      nevents(i) = 0
   80 continue
c Check if any messages remain outstanding
      state = 0
      do 90 i = 1, MAXM
      if (curreq(1,i).ne.0) state = state + 1
   90 continue
      if (state.gt.0) then
         write (2,*) state, ' message(s) never cleared'
         do 100 i = 1, MAXM
         if (curreq(1,i).ne.0) then
            if (curreq(1,i).eq.(-1)) then
               write (2,*) 'transmission mode = send'
            elseif (curreq(1,i).eq.1) then
               write (2,*) 'transmission mode = receive'
            endif
            write (2,*) 'destination/source = ', curreq(2,i)
            write (2,*) 'communicator = ', curreq(3,i)
            write (2,*) 'tag = ', curreq(4,i)
            write (2,*) 'datatype = ', curreq(5,i)
            write (2,*)
         endif
  100    continue
      endif
      if (monitor.eq.2) write (2,*) 'MPI_FINALIZE complete'
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_SEND(buf,count,datatype,dest,tag,comm,ierror)
c blocking standard mode send
c buf = initial address of send buffer
c count = number of entries to send
c datatype = datatype of each entry
c dest = rank of destination
c tag = message tag
c comm = communicator
c ierror = error indicator
c input: buf, count, datatype, dest, tag, comm
c output: ierror
      implicit none
      integer buf(*)
      integer count, datatype, dest, tag, comm, ierror
c MPI constants
      integer MPI_STATUS_SIZE
      parameter(MPI_STATUS_SIZE=5)
c local data
      integer request, status
      dimension status(MPI_STATUS_SIZE)
      call MPI_ISEND(buf,count,datatype,dest,tag,comm,request,ierror)
      call MPI_WAIT(request,status,ierror)
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_RECV(buf,count,datatype,source,tag,comm,status,ierr
     1or)
c blocking receive
c buf = initial address of receive buffer
c count = maximum number of entries to receive
c datatype = datatype of each entry
c source = rank of source
c tag = message tag
c comm = communicator
c status = return status
c ierror = error indicator
c input: count, datatype, source, tag, comm
c output: buf, status, ierror
      implicit none
      integer buf(*), status(*)
      integer count, datatype, source, tag, comm, ierror
c local data
      integer request
      call MPI_IRECV(buf,count,datatype,source,tag,comm,request,ierror)
      call MPI_WAIT(request,status,ierror)
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_ISEND(buf,count,datatype,dest,tag,comm,request,ierr
     1or)
c start a non-blocking send
c buf = initial address of send buffer
c count = number of entries to send
c datatype = datatype of each entry
c dest = rank of destination
c tag = message tag
c comm = communicator
c request = request handle
c ierror = error indicator
c input: buf, count, datatype, dest, tag, comm
c output: request, ierror
      implicit none
      integer buf(*)
      integer count, datatype, dest, tag, comm, request, ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c idproc = processor id
c epref = array of endpoint references for each participating node
c ioc = array of context pointers for notifier function
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer MAXM, curreq, header, rwrec, monitor, trash, mqueue
      parameter(MAXM=2*MAXS)
      dimension curreq(5,MAXM), header(10,MAXM), rwrec(14,MAXM)
      dimension trash(256), mqueue(2,MAXS+1)
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
c curreq = request record for transmission parameters
c header = message envelope
c rwrec = read/write record for asynchronous messages
c mqueue = message request queue
      common /mpisendrec/ monitor, curreq, header, rwrec, trash, mqueue
c function declarations
      integer ioresult
      external ioresult
c OT constants
      integer kNetbufDataIsOTData
      parameter(kNetbufDataIsOTData=-2)
c MPI constants
      integer MPI_REQUEST_NULL, MPI_PROC_NULL
      parameter(MPI_REQUEST_NULL=-1,MPI_PROC_NULL=-3)
      integer MPI_INTEGER, MPI_REAL, MPI_DOUBLE_PRECISION
      parameter(MPI_INTEGER=18,MPI_REAL=19,MPI_DOUBLE_PRECISION=20)
      integer MPI_COMPLEX, MPI_DOUBLE_COMPLEX, MPI_BYTE
      parameter(MPI_COMPLEX=22,MPI_DOUBLE_COMPLEX=23,MPI_BYTE=2)
      integer MPI_2REAL, MPI_2DOUBLE_PRECISION, MPI_2INTEGER
      parameter(MPI_2REAL=35,MPI_2DOUBLE_PRECISION=36,MPI_2INTEGER=37)
c local data
      integer np, longw, i, response
      ierror = 0
      if (dest.eq.MPI_PROC_NULL)  then
         request = MPI_REQUEST_NULL
         return
      endif
c find space for record
      i = 0
   10 i = i + 1
      if (i.gt.MAXM) then
         write (2,*) 'too many sends waiting, dest, tag = ', dest, tag
         request = MPI_REQUEST_NULL
         ierror = 14
         call writerrs('MPI_ISEND: ',ierror)
         return
      elseif (curreq(1,i).ne.0) then
         go to 10
      endif
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c invalid count
      elseif (count.lt.0) then
         ierror = 3
c invalid destination
      elseif ((dest.lt.0).or.(dest.ge.nproc)) then
         write (2,*) 'destination = ', dest
         ierror = 4
c invalid tag
      elseif (tag.lt.(-1)) then
         ierror = 6
c communicator errors
      else
         longw = mapcomm(MAXS+1,comm+1)
         np = mapcomm(dest+1,comm+1)
c communicator not mapped
         if ((longw.le.0).or.(longw.gt.nproc)) then
            ierror = 2
c invalid destination
         elseif ((dest.lt.0).or.(dest.ge.longw)) then
            write (2,*) 'destination = ', dest
            ierror = 4
c invalid mapping
         elseif ((np.lt.0).or.(np.ge.nproc)) then
            write (2,*) 'Invalid mapping, destination, node = ', dest,np
            ierror = 2
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_ISEND: ',ierror)
         return
      endif
c Create header
c OTData structure for header
      header(1,i) = loc(header(4,i))
      header(2,i) = loc(header(7,i))
      header(3,i) = 16
c Save communicator
      header(7,i) = comm
c Save tag
      header(8,i) = tag
c Save datatype
      header(9,i) = datatype
c Set destination id for selfsends
      if (idproc.eq.np) then
         np = MAXS + 1
c Set normal destination id
      else
         np = np + 1
      endif
c Set endpoint reference pointer
      rwrec(1,i) = epref(np)
c Reset iocompletion flag
      rwrec(2,i) = 1
c Set pointer to header
      rwrec(3,i) = loc(header(1,i))
c Set buffer length for header
      rwrec(4,i) = kNetbufDataIsOTData
c Find buffer length for data
      if ((datatype.eq.MPI_INTEGER).or.(datatype.eq.MPI_REAL)) then
         longw = 4*count
      elseif ((datatype.eq.MPI_DOUBLE_PRECISION).or.(datatype.eq.MPI_COM
     1PLEX)) then
         longw = 8*count
      elseif (datatype.eq.MPI_DOUBLE_COMPLEX) then
         longw = 16*count
      elseif (datatype.eq.MPI_BYTE) then
         longw = count
      elseif ((datatype.eq.MPI_2INTEGER).or.(datatype.eq.MPI_2REAL)) the
     1n
         longw = 8*count
      elseif (datatype.eq.MPI_2DOUBLE_PRECISION) then
         longw = 16*count
c Invalid datatype
      else
         ierror = 7
         call writerrs('MPI_ISEND: ',ierror)
         return
      endif
c Set pointer to data buffer
      rwrec(6,i) = loc(buf)
c Set buffer lengths for data
      rwrec(7,i) = longw
      rwrec(8,i) = -header(3,i)
c Clear more flag
      rwrec(5,i) = 0
c Clear next message flag
      rwrec(13,i) = 0
c Clear non-fatal error code
      rwrec(14,i) = 0
c OTData structure for data
      header(4,i) = 0
      header(5,i) = rwrec(6,i)
      header(6,i) = longw
c Save length
      header(10,i) = longw
c Obtain the current time stamp
      call OTGetTimeStamp(rwrec(9,i))
c Limit notifications that can be sent to notifier
      call OTEnterNotifier(val4(epref(np)))
c Append this message to send queue if necessary
      if (ioc(4,np).gt.0) then
         if (mqueue(2,np).gt.0) then
            rwrec(13,mqueue(2,np)) = i
         else
            rwrec(13,ioc(4,np)) = i
         endif
         mqueue(2,np) = i
c Obtain the current time stamp
         call OTGetTimeStamp(rwrec(11,i))
         go to 30
      endif
c First send 4 word header, then data
      call sndmsgf(i,np)
      response = ioresult(rwrec(1,i))
c Set pointer to current send record
      if (response.gt.0) then
         ioc(4,np) = i
c Set iocompletion flag to error
      elseif (response.lt.0) then
         ierror = response
      endif
c Allow Open Transport to resume sending events
   30 call OTLeaveNotifier(val4(epref(np)))
c Handle read and write errors
      if (ierror.ne.0) then
c Write out readwrite record
         call rwstat(i,2)
         call wqueue(2)
         do 40 i = 1, MAXM
         if (curreq(1,i).ne.0) call rwstat(i,2)
   40    continue
         call writerrs('MPI_ISEND: ',ierror)
         return
      endif
c Find actual destination
      if (np.eq.(MAXS+1)) then
         np = idproc
      else
         np = np - 1
      endif
c log MPI message state change and display status
      if (monitor.gt.0) call logmess(np,1,rwrec(7,i),0,tag)
c save transmission mode as send
      curreq(1,i) = -1
c save destination/source id
      curreq(2,i) = np
c save communicator
      curreq(3,i) = comm
c save tag
      curreq(4,i) = tag
c save datatype
      curreq(5,i) = datatype
c assign request handle
      request = i
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_IRECV(buf,count,datatype,source,tag,comm,request,ie
     1rror)
c begin a non-blocking receive
c buf = initial address of receive buffer
c count = maximum number of entries to receive
c datatype = datatype of each entry
c source = rank of source
c tag = message tag
c comm = communicator
c request = request handle
c ierror = error indicator
c input: count, datatype, source, tag, comm
c output: buf, request, ierror
      implicit none
      integer buf(*)
      integer count, datatype, source, tag, comm, request, ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c epref = array of endpoint references for each participating node
c ioc = array of context pointers for notifier function
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer MAXM, curreq, header, rwrec, monitor, trash, mqueue
      parameter(MAXM=2*MAXS)
      dimension curreq(5,MAXM), header(10,MAXM), rwrec(14,MAXM)
      dimension trash(256), mqueue(2,MAXS+1)
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
c curreq = request record for transmission parameters
c header = message envelope
c rwrec = read/write record for asynchronous messages
c mqueue = message request queue
      common /mpisendrec/ monitor, curreq, header, rwrec, trash, mqueue
c function declarations
      integer ioresult
      external ioresult
c MPI constants
      integer MPI_ANY_SOURCE, MPI_REQUEST_NULL, MPI_PROC_NULL
      parameter(MPI_ANY_SOURCE=-1,MPI_REQUEST_NULL=-1,MPI_PROC_NULL=-3)
      integer MPI_INTEGER, MPI_REAL, MPI_DOUBLE_PRECISION
      parameter(MPI_INTEGER=18,MPI_REAL=19,MPI_DOUBLE_PRECISION=20)
      integer MPI_COMPLEX, MPI_DOUBLE_COMPLEX, MPI_BYTE
      parameter(MPI_COMPLEX=22,MPI_DOUBLE_COMPLEX=23,MPI_BYTE=2)
      integer MPI_2REAL, MPI_2DOUBLE_PRECISION, MPI_2INTEGER
      parameter(MPI_2REAL=35,MPI_2DOUBLE_PRECISION=36,MPI_2INTEGER=37)
c local data
      integer np, longw, i, response
      ierror = 0
      if (source.eq.MPI_PROC_NULL) then
         request = MPI_REQUEST_NULL
         return
      endif
c find space for record
      i = 0
   10 i = i + 1
      if (i.gt.MAXM) then
         write (2,*) 'too many receives waiting, source, tag = ', source
     1, tag
         request = MPI_REQUEST_NULL
         ierror = 15
         call writerrs('MPI_IRECV: ',ierror)
         return
      elseif (curreq(1,i).ne.0) then
         go to 10
      endif
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c invalid count
      elseif (count.lt.0) then
         ierror = 3
c invalid source
      elseif ((source.lt.0).or.(source.ge.nproc)) then
         if (source.eq.MPI_ANY_SOURCE) then
            write (2,*) 'MPI_ANY_SOURCE not supported'
         else
            write (2,*) 'source = ', source
         endif
         ierror = 5
c invalid tag
      elseif (tag.lt.(-1)) then
         ierror = 6
c communicator errors
      else
         longw = mapcomm(MAXS+1,comm+1)
         np = mapcomm(source+1,comm+1)
c communicator not mapped
         if ((longw.le.0).or.(longw.gt.nproc)) then
            ierror = 2
c invalid source
         elseif ((source.lt.0).or.(source.ge.longw)) then
            write (2,*) 'source = ', source
            ierror = 5
c invalid mapping
         elseif ((np.lt.0).or.(np.ge.nproc)) then
            write (2,*) 'Invalid mapping, source, node = ', source, np
            ierror = 2
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_IRECV: ',ierror)
         return
      endif
c Set source id
      np = np + 1
c Set endpoint reference pointer
      rwrec(1,i) = epref(np)
c Reset iocompletion flag for receive
      rwrec(2,i) = 1
c Set pointer to header
      rwrec(3,i) = loc(header(7,i))
c Set buffer length for header
      rwrec(4,i) = 16
c Set buffer length for data
      if ((datatype.eq.MPI_INTEGER).or.(datatype.eq.MPI_REAL)) then
         longw = 4*count
      elseif ((datatype.eq.MPI_DOUBLE_PRECISION).or.(datatype.eq.MPI_COM
     1PLEX)) then
         longw = 8*count
      elseif (datatype.eq.MPI_DOUBLE_COMPLEX) then
         longw = 16*count
      elseif (datatype.eq.MPI_BYTE) then
         longw = count
      elseif ((datatype.eq.MPI_2INTEGER).or.(datatype.eq.MPI_2REAL)) the
     1n
         longw = 8*count
      elseif (datatype.eq.MPI_2DOUBLE_PRECISION) then
         longw = 16*count
c invalid datatype
      else
         ierror = 7
         call writerrs('MPI_IRECV: ',ierror)
         return
      endif
c Set pointer to data buffer
      rwrec(6,i) = loc(buf)
c Set buffer lengths
      rwrec(7,i) = longw
      rwrec(8,i) = -rwrec(4,i)
c Clear more flag
      rwrec(5,i) = 0
c Clear next message flag
      rwrec(13,i) = 0
c Clear non-fatal error code
      rwrec(14,i) = 0
c Clear length
      header(10,i) = 0
c Obtain the current time stamp
      call OTGetTimeStamp(rwrec(9,i))
c Limit notifications that can be sent to notifier
      call OTEnterNotifier(val4(epref(np)))
c Append this message to receive queue if necessary
      if (ioc(3,np).gt.0) then
         if (mqueue(1,np).gt.0) then
            rwrec(13,mqueue(1,np)) = i
         else
            rwrec(13,ioc(3,np)) = i
         endif
         mqueue(1,np) = i
c Obtain the current time stamp
         call OTGetTimeStamp(rwrec(11,i))
         go to 40
      endif
c First receive 4 word header, then data
      call rcvmsgf(i,np)
      response = ioresult(rwrec(1,i))
c Set pointer to current receive record
      if (response.gt.0) then
         ioc(3,np) = i
c Set iocompletion flag to error
      elseif (response.lt.0) then
         ierror = response
      endif
c Allow Open Transport to resume sending events
   40 call OTLeaveNotifier(val4(epref(np)))
c Handle read and write errors
      if (ierror.ne.0) then
c Write out readwrite record
         call rwstat(i,2)
         call wqueue(2)
         do 50 i = 1, MAXM
         if (curreq(1,i).ne.0) call rwstat(i,2)
   50    continue
         call writerrs('MPI_IRECV: ',ierror)
         return
      endif
c Find actual source
      np = np - 1
c log MPI message state change and display status
      if (monitor.gt.0) call logmess(np,2,rwrec(7,i),0,tag)
c save transmission mode as receive
      curreq(1,i) = 1
c save destination/source id
      curreq(2,i) = np
c save communicator
      curreq(3,i) = comm
c save tag
      curreq(4,i) = tag
c save datatype
      curreq(5,i) = datatype
c assign request handle
      request = i
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_TEST(request,flag,status,ierror)
c check to see if a nonblocking send or receive operation has completed
c request = request handle
c flag = true if operation completed
c status = status object
c ierror = error indicator
c input: request
c output: request, flag, status, ierror
      implicit none
      integer status(*)
      integer request, ierror
      logical flag
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c ioc = array of context pointers for notifier function
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer MAXM, curreq, header, rwrec, monitor, trash, mqueue
      parameter(MAXM=2*MAXS)
      dimension curreq(5,MAXM), header(10,MAXM), rwrec(14,MAXM)
      dimension trash(256), mqueue(2,MAXS+1)
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
c curreq = request record for transmission parameters
c header = message envelope
c rwrec = read/write record for asynchronous messages
      common /mpisendrec/ monitor, curreq, header, rwrec, trash, mqueue
c function declarations
      integer ioresult, OTElapsedMilliseconds, OTTimeStampInMicroseconds
      logical checkesc
      external ioresult, checkesc
      external OTElapsedMilliseconds, OTTimeStampInMicroseconds
c OT constants
      integer kOTOutStateErr, kENOMEMErr, kOTFlowErr, kOTNoDataErr
      parameter(kOTOutStateErr=-3155,kENOMEMErr=-3211)
      parameter(kOTFlowErr=-3161,kOTNoDataErr=-3162)
c MPI constants
      integer MPI_COMM_WORLD, MPI_REQUEST_NULL
      parameter(MPI_COMM_WORLD=0,MPI_REQUEST_NULL=-1)
c local data
      integer i, dest, source, slen, comm, tag, datatype, j
      integer rlen, rcomm, rtag, rdatat, nerr, mstime, mptime
      integer delta
      real ts
      dimension delta(2)
c mstime = maximum time (msec) to wait for message to start arriving
c mptime = maximum time (msec) to wait for next part of message
      data mstime, mptime /300000,10000/
      save mstime, mptime
      ierror = 0
c check for error conditions
      i = request
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c null request
      elseif (request.lt.0) then
         flag = .true.
         return
c invalid request handle
      elseif ((i.lt.1).or.(i.gt.MAXM)) then
         ierror = 16
      elseif (curreq(1,i).eq.0) then
         ierror = 16
      endif
c handle errors
      if (ierror.ne.0) then
         status(3) = ierror
         call writerrs('MPI_TEST: ',ierror)
         return
      endif
c set status to empty
      status(1) = -1
      status(2) = -1
      status(3) = 0
      status(4) = 0
      status(5) = 0
      flag = .false.
c Check if data read or write has completed
      if (ioresult(rwrec(1,i)).gt.0) then
         if (checkesc(0)) then
            if (curreq(1,i).lt.0) then
               write (2,*) 'Send killed, dest,tag = ', curreq(2,i),
     1curreq(4,i)
            else
               write (2,*) 'Receive killed, source,tag = ', curreq(2,i),
     1curreq(4,i)
            endif
            ierror = -9
            call writerrs('MPI_TEST: ',ierror)
            return
c Check Timeout
         else
c Limit notifications that can be sent to notifier
            call OTEnterNotifier(val4(rwrec(1,i)))
c Calculate time elapsed in milliseconds
            nerr = OTElapsedMilliseconds(rwrec(11,i))
c Check if send or receive should be retried
            if ((nerr.lt.mptime).and.(rwrec(14,i).ne.kENOMEMErr)) then
c Allow Open Transport to resume sending events
               call OTLeaveNotifier(val4(rwrec(1,i)))
               return
            endif
c Check if message arrived during checkesc
            nerr = ioresult(rwrec(1,i))
            if (nerr.le.0) then
               write (2,*) 'Info: message arrived during checkesc'
               flag = .true.
               go to 20
            endif
c Check if any other messages need service
            do 10 j = 1, MAXM
            if ((j.ne.i).and.(curreq(1,j).lt.0)) then
c Limit notifications that can be sent to notifier
               if (rwrec(1,i).ne.rwrec(1,j)) then
                  call OTEnterNotifier(val4(rwrec(1,j)))
               endif
c Check if out of memory flag is set and message is incomplete
               if (rwrec(14,j).eq.kENOMEMErr) then
                  if (ioresult(rwrec(1,j)).gt.0) then
                     dest = curreq(2,j)
                     write (2,*) 'send endpt ', dest, 'out of memory'
                     if (idproc.eq.dest) dest = MAXS
                     dest = dest + 1
c Attempt to send next block of data
                     call sndmsgf(j,dest)
                     write (2,*) 'Info: request, response=', j,
     1ioresult(rwrec(1,j))
                  endif
               endif
c Allow Open Transport to resume sending events
               if (rwrec(1,i).ne.rwrec(1,j)) then
                  call OTLeaveNotifier(val4(rwrec(1,j)))
               endif
            elseif ((j.ne.i).and.(curreq(1,j).gt.0)) then
c Limit notifications that can be sent to notifier
               if (rwrec(1,i).ne.rwrec(1,j)) then
                  call OTEnterNotifier(val4(rwrec(1,j)))
               endif
c Check if out of memory flag is set and message is incomplete
               if (rwrec(14,j).eq.kENOMEMErr) then
                  if (ioresult(rwrec(1,j)).gt.0) then
                     source = curreq(2,j)
                     write (2,*) 'recv endpt ', source, 'out of memory'
                     source = source + 1
c Attempt to receive next block of data
                     call rcvmsgf(j,source)
                     write (2,*) 'Info: request, response=', j,
     1ioresult(rwrec(1,j))
                  endif
               endif
c Allow Open Transport to resume sending events
               if (rwrec(1,i).ne.rwrec(1,j)) then
                  call OTLeaveNotifier(val4(rwrec(1,j)))
               endif
            endif
   10       continue
c Try to determine cause of Timeout in Send
            if (curreq(1,i).lt.0) then
               dest = curreq(2,i)
               if (rwrec(14,i).eq.kENOMEMErr) then
                  write (2,*) 'OT Temporarily out of Memory'
               else
                  ts = .001*real(OTElapsedMilliseconds(rwrec(9,i)))
                  write (2,*) 'Send Timeout, ', ts, ' sec, Retrying...'
               endif
c Debug information
               write (2,*) 'destination=', dest, ' size=' , rwrec(7,i),
     1'tag=', curreq(4,i)
               if (idproc.eq.dest) dest = MAXS
               dest = dest + 1
c Attempt to send next block of data
               call sndmsgf(i,dest)
               nerr = ioresult(rwrec(1,i))
c Non-fatal errors
               if (rwrec(14,i).ne.0) then
                  if (rwrec(14,i).eq.kOTFlowErr) then
                     write (2,*) 'Flow Control prevents sending data'
                  endif
c Do not wait more than 5 minutes for message to start sending
                  if (OTElapsedMilliseconds(rwrec(9,i)).ge.mstime) then
                     nerr = rwrec(14,i)
                     if (nerr.eq.kENOMEMErr) then
                        write (2,*) 'Open Transport is out of memory'
                     endif
                     write (2,*) 'OTSnd Retry failed'
                     ierror = nerr
                     ioc(4,dest) = 0
                     flag = .true.
                  endif
c Fatal errors
               elseif (nerr.lt.0) then
                  ierror = nerr
                  flag = .true.
c Data successfully sent
               else
c Debug information
                  write (2,*) 'data block sent, current total = ',
     1rwrec(8,i)
                  if (nerr.eq.0) flag = .true.
               endif
c Try to determine cause of Timeout in Receive
            else
               source = curreq(2,i)
               if (rwrec(14,i).eq.kENOMEMErr) then
                  write (2,*) 'OT Temporarily out of Memory'
               else
                  ts = .001*real(OTElapsedMilliseconds(rwrec(9,i)))
                  write (2,*) 'Receive Timeout, ',ts,' sec, Retrying...'
               endif
c Debug information
               write (2,*) 'source=', source, ' size=', rwrec(7,i),
     1'tag=', curreq(4,i)
               source = source + 1
c Attempt to receive next block of data
               call rcvmsgf(i,source)
               nerr = ioresult(rwrec(1,i))
c Non-fatal errors
               if (rwrec(14,i).ne.0) then
                  if (rwrec(14,i).eq.kOTFlowErr) then
                     write (2,*) 'Flow Control prevents accepting data'
                  elseif (rwrec(14,i).eq.kOTNoDataErr) then
                     write (2,*) 'No data is currently available'
                  endif
c Do not wait more than 5 minutes for message to start sending
                  if (OTElapsedMilliseconds(rwrec(9,i)).ge.mstime) then
                     nerr = rwrec(14,i)
                     if (nerr.eq.kENOMEMErr) then
                        write (2,*) 'Open Transport is out of memory'
                     endif
                     write (2,*) 'OTRcv Retry failed'
                     ierror = nerr
                     ioc(3,source) = 0
                     flag = .true.
                  endif
c Fatal errors
               elseif (nerr.lt.0) then
                  ierror = nerr
                  flag = .true.
c Data successfully received
               else
c Debug information
                  write (2,*) 'data block received, current total=',
     1rwrec(8,i)
                  if (nerr.eq.0) flag = .true.
               endif
            endif
c Allow Open Transport to resume sending events
   20       call OTLeaveNotifier(val4(rwrec(1,i)))
            if (.not.flag) return
         endif
c Data read or write has completed
      else
         nerr = ioresult(rwrec(1,i))
         flag = .true.
      endif
c Get requested length
      slen = rwrec(7,i)
c Get actual length
      rlen = rwrec(8,i)
c Read current request record
      tag = curreq(4,i)
c Define length for MPI_GET_COUNT
      status(4) = rlen
c Check for send errors
      if (curreq(1,i).lt.0) then
         dest = curreq(2,i)
c Define type for MPI_GET_COUNT
         status(5) = curreq(5,i)
c Check for write errors
         if (nerr.lt.0) then
            if (ierror.eq.0) then
            write (2,*) 'OTSnd Error, nerr, dest, tag=', nerr, dest, tag
               if (nerr.eq.kOTOutStateErr) then
                  write (2,*) 'Endpoint not in an appropriate state'
               elseif (nerr.eq.kENOMEMErr) then
                  write (2,*) 'Open Transport has run out of memory'
               endif
               ierror = nerr
            endif
         elseif (rlen.ne.slen) then
            write (2,*) 'Send Length Error, dest, tag, requested/actual 
     1length = ', dest, tag, slen, rlen
            ierror = 8
         endif
c log MPI message state change and display status
         if (monitor.gt.0) then
            if (checkesc(0)) then
               ierror = -9
               call writerrs('MPI_TEST: ',ierror)
               return
            elseif (ierror.eq.0) then
c Subtract one timestamp value from another
               call OTSubtractTimeStamps(delta,rwrec(9,i),rwrec(11,i))
c Convert difference between time steps into microseconds
               nerr = OTTimeStampInMicroseconds(delta)
               call logmess(dest,-1,rlen,nerr,tag)
            endif
         endif
         go to 30
      endif
c Read current request record
      source = curreq(2,i)
      comm = curreq(3,i)
      datatype = curreq(5,i)
c Read header
c Get received tag from header
      rtag = header(8,i)
c Get received comm from header
      rcomm = header(7,i)
c Get received datatype from header
      rdatat = header(9,i)
c Define source, tag and type for MPI_GET_COUNT
      status(1) = source
      status(2) = rtag
      status(5) = rdatat
c Check for receive errors
      if (nerr.lt.0) then
         if (ierror.eq.0) then
         write (2,*) 'OTRcv Error, nerr, source, tag = ',nerr,source,tag
            if (nerr.eq.kOTOutStateErr) then
               write (2,*) 'Endpoint not in an appropriate state'
            elseif (nerr.eq.kENOMEMErr) then
               write (2,*) 'Open Transport is out of memory'
            endif
            ierror = nerr
         endif
c Length error
      elseif (rlen.gt.slen) then
         write (2,*) 'Read Length Error, source, tag, requested/actual='
     1, source, tag, slen, rlen
         write (2,*) 'Possibly attempting to receive data out of order'
         ierror = 13
c Check for incomplete message, this should never be able to happen
      elseif (rlen.lt.header(10,i)) then
         write (2,*) 'Incomplete Read, source, tag, expected/actual = '  
     1, source, tag, header(10,i), rlen
         ierror = 12
c Comm error
      elseif (rcomm.ne.comm) then
         write (2,*) 'Read Comm Error, source, tag, expected/received co
     1mm = ', source, tag, comm, rcomm
         ierror = 9
c Tag error
      elseif ((tag.ge.0).and.(rtag.ne.tag)) then
         write (2,*) 'Read Tag Error, source, expected/received tag = '
     1, source, tag, rtag
         write (2,*) 'Possibly attempting to receive data out of order'
         ierror = 10
c Type error
      elseif (rdatat.ne.datatype) then
         write (2,*) 'Read Type Error, source, tag, expected/received ty
     1pe = ', source, tag, datatype, rdatat
         ierror = 11
      endif
c log MPI message state change and display status
      if (monitor.gt.0) then
         if (checkesc(0)) then
            ierror = -9
            call writerrs('MPI_TEST: ',ierror)
            return
         elseif (ierror.eq.0) then
c Subtract one timestamp value from another
            call OTSubtractTimeStamps(delta,rwrec(9,i),rwrec(11,i))
c Convert difference between time steps into microseconds
            nerr = OTTimeStampInMicroseconds(delta)
            call logmess(source,-2,rlen,nerr,tag)
         endif
      endif
c Store error code
   30 status(3) = ierror
c Nullify transmission mode
      curreq(1,i) = 0
c Nullify request handle
      request = MPI_REQUEST_NULL
c Handle read and write errors
      if (ierror.ne.0) then
c Write out readwrite record
         call rwstat(i,2)
         call wqueue(2)
         do 40 i = 1, MAXM
         if (curreq(1,i).ne.0) call rwstat(i,2)
   40    continue
         call writerrs('MPI_TEST: ',ierror)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine sndmsgf(request,dest)
c send a message fragment
c request = request handle
c dest = rank of destination + 1
      implicit none
      integer request, dest
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c ioc = array of context pointers for notifier function
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer MAXM, curreq, header, rwrec, monitor, trash, mqueue
      parameter(MAXM=2*MAXS)
      dimension curreq(5,MAXM), header(10,MAXM), rwrec(14,MAXM)
      dimension trash(256), mqueue(2,MAXS+1)
c curreq = request record for transmission parameters
c rwrec = read/write record for asynchronous messages
c mqueue = message request queue
      common /mpisendrec/ monitor, curreq, header, rwrec, trash, mqueue
c function declarations
      integer OTSnd, OTLook
      external OTSnd, OTLook
c OT constants
      integer kOTFlowErr, kENOMEMErr, kOTOutStateErr, kOTLookErr
      parameter(kOTFlowErr=-3161,kENOMEMErr=-3211,kOTOutStateErr=-3155)
      parameter(kOTLookErr=-3158)
c local data
      integer i, np, nps, response
      i = request
      np = dest
c Obtain the current time stamp
      call OTGetTimeStamp(rwrec(11,i))
c Send data to a remote peer
   10 response = OTSnd(val4(rwrec(1,i)),val4(rwrec(3,i)),val4(rwrec(4,i)
     1),val4(rwrec(5,i)))
c Process data which has been sent
      if (response.ge.0) then
c Clear non-fatal error code
         rwrec(14,i) = 0
c Set actual length sent
         rwrec(8,i) = rwrec(8,i) + response
c Check for incomplete header
         if (rwrec(8,i).lt.0) then
            header(2,i) = header(2,i) + response
            header(3,i) = header(3,i) - response
c Check for incomplete data
         elseif (rwrec(7,i).gt.rwrec(8,i)) then
c Header is sent, readjust parameters to send data
            if (rwrec(2,i).eq.1) then
               rwrec(3,i) = rwrec(6,i)
               rwrec(4,i) = rwrec(7,i)
               response = response - header(3,i)
            endif
c Readjust buffer pointer
            rwrec(3,i) = rwrec(3,i) + response
            rwrec(4,i) = rwrec(4,i) - response
            rwrec(2,i) = rwrec(2,i) + 1
c Data is sent
         else
c Obtain the current time stamp
            call OTGetTimeStamp(rwrec(11,i))
            rwrec(2,i) = 0
c Get next message if messages are queued
            if (rwrec(13,i).gt.0) then
               i = rwrec(13,i)
               if (i.eq.mqueue(2,np)) mqueue(2,np) = 0
               ioc(4,np) = i
               go to 10
            else
c Clear pointer to current send record
               ioc(4,np) = 0
            endif
         endif
c Check for errors
      else
c Process non-fatal errors
         if ((response.eq.kOTFlowErr).or.(response.eq.kENOMEMErr)) then
            rwrec(14,i) = response
c Process fatal errors
         else
c Find actual destination
         nps = np - 1
         if (nps.eq.(MAXS+1)) nps = idproc
            write (2,*) 'OTSnd Error, ierr, dest, tag=', response, nps, 
     1header(8,i)
            if (response.eq.kOTOutStateErr) then
               write (2,*) 'Endpoint not in an appropriate state'
c Determine cause of a kOTLookErr
            elseif (response.eq.kOTLookErr) then
               write (2,*) 'OTLookErr, cause=', OTLook(val4(rwrec(1,i)))
            endif
c Set iocompletion flag to error
            rwrec(2,i) = response
c Clear pointer to current send record
            ioc(4,np) = 0
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine rcvmsgf(request,source)
c receive a message fragment
c request = request handle
c source = rank of source + 1
      implicit none
      integer request, source
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c ioc = array of context pointers for notifier function
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer MAXM, curreq, header, rwrec, monitor, trash, mqueue
      parameter(MAXM=2*MAXS)
      dimension curreq(5,MAXM), header(10,MAXM), rwrec(14,MAXM)
      dimension trash(256), mqueue(2,MAXS+1)
c curreq = request record for transmission parameters
c rwrec = read/write record for asynchronous messages
c trash = trash bin for unwanted data
c mqueue = message request queue
      common /mpisendrec/ monitor, curreq, header, rwrec, trash, mqueue
c function declarations
      integer ioresult, OTRcv, OTLook
      external ioresult, OTRcv, OTLook
c OT constants
      integer kOTFlowErr, kOTNoDataErr, kENOMEMErr
      integer kOTLookErr, T_GODATA, kOTOutStateErr
      parameter(kOTFlowErr=-3161,kOTNoDataErr=-3162,kENOMEMErr=-3211)
      parameter(kOTLookErr=-3158,T_GODATA=256,kOTOutStateErr=-3155)
c local data
      integer i, np, nps, response, j
      i = request
      np = source
c Obtain the current time stamp
      call OTGetTimeStamp(rwrec(11,i))
c Read data sent from a remote peer
   10 response = OTRcv(val4(rwrec(1,i)),val4(rwrec(3,i)),val4(rwrec(4,i)
     1),rwrec(5,i))
c Unexpected flag returned
      if (rwrec(5,i).gt.1) then
         rwrec(2,i) = -3
c Clear pointer to current receive record
         ioc(3,np) = 0
c Process data which has arrived
      elseif (response.ge.0) then
c Clear more flag
         rwrec(5,i) = 0
c Clear non-fatal error code
         rwrec(14,i) = 0
c Set actual length received
         rwrec(8,i) = rwrec(8,i) + response
c Check if all the data has arrived
         if (rwrec(8,i).lt.header(10,i)) then
c Incomplete data
            if (rwrec(4,i).gt.response) then
c Readjust buffer pointer
               rwrec(3,i) = rwrec(3,i) + response
               rwrec(4,i) = rwrec(4,i) - response
               if (rwrec(8,i).ge.0) rwrec(2,i) = rwrec(2,i) + 1
c Header is received, readjust parameters to receive data
            elseif (rwrec(2,i).eq.1) then
               rwrec(3,i) = rwrec(6,i)
               rwrec(4,i) = min(header(10,i),rwrec(7,i))
               rwrec(2,i) = 2
               go to 10
c Data is received and buffer is full
            else
               rwrec(3,i) = loc(trash)
               rwrec(4,i) = 1024
               go to 10
            endif
c Data is received
         else
c Obtain the current time stamp
            call OTGetTimeStamp(rwrec(11,i))
            rwrec(2,i) = 0
c Get next message if messages are queued
            if (rwrec(13,i).gt.0) then
               i = rwrec(13,i)
               if (i.eq.mqueue(1,np)) mqueue(1,np) = 0
               ioc(3,np) = i
c Clear pointer to current read record
            else
               ioc(3,np) = 0
            endif
         endif
c Check for errors
      else
c Process non-fatal errors
         if ((response.eq.kOTFlowErr).or.(response.eq.kOTNoDataErr).or.(
     1response.eq.kENOMEMErr)) then
            rwrec(14,i) = response
c Process potentially fatal errors
         else
c Determine cause of a kOTLookErr
            if (response.eq.kOTLookErr) then
               j = OTLook(val4(rwrec(1,i)))
               if (j.eq.T_GODATA) then
c Set pointer to current send record
                  j = ioc(4,np)
c Send next block of data
                  if ((j.gt.0).and.(j.le.MAXM)) then
                     write (2,*) 'rcvmsgf: OTLookErr for node=', np-1
                     call sndmsgf(j,np)
                     response = ioresult(rwrec(1,j))
c Report send error
                     if (response.lt.0) then
                        write (2,*) 'rcvmsgf: send error=', response 
                        ioc(4,np) = 0
                     endif
                  endif
c Clear more flag
                  rwrec(5,i) = 0
                  go to 10
               endif
            endif
c Process fatal errors
            nps = np - 1
            write (2,*) 'OTRcv Error, ierr, source, tag=', response, nps
     1, header(8,i)
            if (response.eq.kOTOutStateErr) then
               write (2,*) 'Endpoint not in an appropriate state'
            elseif (response.eq.kOTLookErr) then
               write (2,*) 'OTLookErr, unknown cause=', j
            endif
c Set iocompletion flag to error
            rwrec(2,i) = response
c Clear pointer to current receive record
            ioc(3,np) = 0
         endif
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_WAIT(request,status,ierror)
c wait for an MPI send or receive to complete
c request = request handle
c status = status object
c ierror = error indicator
c input: request
c output: request, status, ierror
      implicit none
      integer status(*)
      integer request, ierror
c local data
      logical flag
   10 call MPI_TEST(request,flag,status,ierror)
      if (.not.flag) go to 10
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_REQUEST_FREE(request,ierror)
c free a communication request object
c request = request handle
c ierror = error indicator
c input: request
c output: request, ierror
      implicit none
      integer request, ierror
c declare internal mpi common block
      integer nproc, idproc
      integer MAXS
      parameter(MAXS=32)
c nproc = number of real or virtual processors obtained
      common /mpiparms/ nproc, idproc
c declare common block for non-blocking messages
      integer MAXM, curreq, header, rwrec, monitor, trash, mqueue
      parameter(MAXM=2*MAXS)
      dimension curreq(5,MAXM), header(10,MAXM), rwrec(14,MAXM)
      dimension trash(256), mqueue(2,MAXS+1)
c curreq = request record for transmission parameters
c rwrec = read/write record for asynchronous messages
      common /mpisendrec/ monitor, curreq, header, rwrec, trash, mqueue
c function declarations
      integer ioresult
      external ioresult
c MPI constants
      integer MPI_REQUEST_NULL
      parameter(MPI_REQUEST_NULL=-1)
c local data
      integer i
      ierror = 0
c check for error conditions
      i = request
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c null request
      elseif (request.lt.0) then
         return
c invalid request handle
      elseif ((i.lt.1).or.(i.gt.MAXM)) then
         ierror = 16
      elseif (curreq(1,i).eq.0) then
         ierror = 16
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_REQUEST_FREE: ',ierror)
         return
      endif
c Check if data read or write has completed
      if (ioresult(rwrec(1,i)).le.0) then
c Nullify transmission mode
         curreq(1,i) = 0
c Nullify request handle
         request = MPI_REQUEST_NULL
      else
         write (2,*) 'MPI_REQUEST_FREE: Message not Completed'
c Write out readwrite record
         call rwstat(i,2)
         ierror = 32
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_SENDRECV (sendbuf,sendcount,sendtype,dest,sendtag,
     1recvbuf,recvcount,recvtype,source,recvtag,comm,status,ierror)
c blocking send and receive operation
c sendbuf = initial address of send buffer
c sendcount = number of entries to send
c sendtype = type of entries in send buffer
c dest = rank of destination
c sendtag = send tag
c recvbuf = initial address of receive buffer
c recvcount = max number of entries to receive
c recvtype = type of entries in receive buffer
c source = rank of source
c recvtag = receive tag
c comm = communicator
c status = return status
c ierror = error indicator
c input: sendbuf, sendcount, sendtype, dest, sendtag, recvbuf, recvcount
c        recvtype, source, recvtag, comm
c output: recvbuf, status, ierror
      implicit none
      integer sendbuf(*), recvbuf(*), status(*)
      integer sendcount, sendtype, dest, sendtag, recvcount, recvtype
      integer source, recvtag, comm, ierror
c local data
      integer recvreq, sendreq
c post non-blocking receive and send
      call MPI_IRECV(recvbuf,recvcount,recvtype,source,recvtag,comm,recv
     1req,ierror)
      call MPI_ISEND(sendbuf,sendcount,sendtype,dest,sendtag,comm,sendre
     1q,ierror)
c wait for send and receive
      call MPI_WAIT(sendreq,status,ierror)
      call MPI_WAIT(recvreq,status,ierror)
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_SSEND(buf,count,datatype,dest,tag,comm,ierror)
c blocking synchronous mode send
c buf = initial address of send buffer
c count = number of entries to send
c datatype = datatype of each entry
c dest = rank of destination
c tag = message tag
c comm = communicator
c ierror = error indicator
c input: buf, count, datatype, dest, tag, comm
c output: ierror
c comment: this is just a temporary stub
      implicit none
      integer buf(*)
      integer count, datatype, dest, tag, comm, ierror
      call MPI_SEND(buf,count,datatype,dest,tag,comm,ierror)
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_ISSEND(buf,count,datatype,dest,tag,comm,request,ier
     1ror)
c start a non-blocking synchronous mode send
c buf = initial address of send buffer
c count = number of entries to send
c datatype = datatype of each entry
c dest = rank of destination
c tag = message tag
c comm = communicator
c request = request handle
c ierror = error indicator
c input: buf, count, datatype, dest, tag, comm
c output: request, ierror
c comment: this is just a temporary stub
      implicit none
      integer buf(*)
      integer count, datatype, dest, tag, comm, request, ierror
      call MPI_ISEND(buf,count,datatype,dest,tag,comm,request,ierror)
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_WAITALL (count,array_of_requests,array_of_statuses,
     1ierror)
c wait for a collection of specified MPI sends or receives to complete
c count = list length
c array_of_requests = array of request handles
c array_of_statuses = array of status objects
c ierror = error indicator
c input: count, array_of_requests
c output: array_of_requests, array_of_statuses, ierror
      implicit none
c MPI constants
      integer MPI_STATUS_SIZE, MPI_ERR_IN_STATUS
      parameter(MPI_STATUS_SIZE=5,MPI_ERR_IN_STATUS=67)
      integer array_of_requests(*), array_of_statuses(MPI_STATUS_SIZE,*)
      integer count, ierror
c local data
      integer i, ierr
c invalid count
      if (count.lt.0) then
         write (2,*) 'Invalid list length = ', count
         ierror = 17
         call writerrs('MPI_WAITALL: ',ierror)
         return
      endif
      ierror = 0
      do 10 i = 1, count
      call MPI_WAIT(array_of_requests(i),array_of_statuses(1,i),ierr)
      if (ierr.ne.0) ierror = MPI_ERR_IN_STATUS
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_WAITANY(count,array_of_requests,index,status,ierror
     1)
c wait for any specified MPI send or receive to complete
c count = list length
c array_of_requests = array of request handles
c index = index of request handle that completed
c status = status object
c ierror = error indicator
c input: count, array_of_requests
c output: array_of_requests, index, status, ierror
      implicit none
      integer array_of_requests(*), status(*)
      integer count, index, ierror
c local data
      integer i, k
      logical flag
c invalid count
      if (count.lt.0) then
         write (2,*) 'Invalid list length = ', count
         ierror = 17
         call writerrs('MPI_WAITANY: ',ierror)
         return
      endif
c find number of requests already completed
      k = 0
      do 10 i = 1, count
      if (array_of_requests(i).lt.0) k = k + 1
   10 continue
      if (k.eq.count) then
         index = -1
         ierror = 0
         return
      endif
      i = 1
   20 flag = .false.
      if (array_of_requests(i).ge.0) then
         call MPI_TEST(array_of_requests(i),flag,status,ierror)
      endif
      if (flag) then
         index = i
      else
         i = i + 1
         if (i.gt.count) i = 1
         go to 20
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_GET_COUNT(status,datatype,count,ierror)
c get the number of "top level" elements
c status = return status of receive operation
c datatype = datatype of each receive buffer entry
c count = number of received entries
c ierror = error indicator
c input: status, datatype
c output: count, ierror
      implicit none
      integer status(*)
      integer datatype, count, ierror
c declare internal mpi common block
      integer nproc, idproc
c nproc = number of real or virtual processors obtained
      common /mpiparms/ nproc, idproc
c MPI constants
      integer MPI_INTEGER, MPI_REAL, MPI_DOUBLE_PRECISION
      parameter(MPI_INTEGER=18,MPI_REAL=19,MPI_DOUBLE_PRECISION=20)
      integer MPI_COMPLEX, MPI_DOUBLE_COMPLEX, MPI_BYTE
      parameter(MPI_COMPLEX=22,MPI_DOUBLE_COMPLEX=23,MPI_BYTE=2)
      integer MPI_2REAL, MPI_2DOUBLE_PRECISION, MPI_2INTEGER
      parameter(MPI_2REAL=35,MPI_2DOUBLE_PRECISION=36,MPI_2INTEGER=37)
      integer MPI_UNDEFINED
      parameter(MPI_UNDEFINED=-1)
      ierror = 0
      count = 0
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c mismatched datatype
      elseif (datatype.ne.status(5)) then
         ierror = 18
c calculate count
      elseif ((datatype.eq.MPI_INTEGER).or.(datatype.eq.MPI_REAL)) then
         count = status(4)/4
         if (4*count.ne.status(4)) count = MPI_UNDEFINED
      elseif ((datatype.eq.MPI_DOUBLE_PRECISION).or.(datatype.eq.MPI_COM
     1PLEX)) then
         count = status(4)/8
         if (8*count.ne.status(4)) count = MPI_UNDEFINED
      elseif (datatype.eq.MPI_DOUBLE_COMPLEX) then
         count = status(4)/16
         if (16*count.ne.status(4)) count = MPI_UNDEFINED
      elseif (datatype.eq.MPI_BYTE) then
         count = status(4)
      elseif ((datatype.eq.MPI_2INTEGER).or.(datatype.eq.MPI_2REAL)) the
     1n
         count = status(4)/8
         if (8*count.ne.status(4)) count = MPI_UNDEFINED
      elseif (datatype.eq.MPI_2DOUBLE_PRECISION) then
         count = status(4)/16
         if (16*count.ne.status(4)) count = MPI_UNDEFINED
c invalid datatype
      else
         ierror = 7
      endif
c handle errors
      if (ierror.ne.0) call writerrs('MPI_GET_COUNT: ',ierror)
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_INITIALIZED(flag,ierror)
c indicate whether MPI_init has been called
c flag = true if MPI_Init has been called, false otherwise
c ierror = error indicator
c output: flag, ierror
      implicit none
      logical flag
      integer ierror
c declare internal mpi common block
      integer nproc, idproc
c nproc = number of real or virtual processors obtained
      common /mpiparms/ nproc, idproc
      if (nproc.gt.0) then
         flag = .true.
      else
         flag = .false.
      endif
      ierror = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_COMM_SIZE(comm,size,ierror)
c determine the size of the group associated with a communicator
c comm = communicator
c size = number of processors in the group of comm
c ierror = error indicator
c input: comm
c output: size, ierror
      implicit none
      integer comm, size, ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c local data
      integer np
      ierror = 0
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c get size
      else
         np = mapcomm(MAXS+1,comm+1)
         if ((np.gt.0).and.(np.le.nproc)) then
            size = np
         else
            ierror = 2
         endif
      endif
c handle errors
      if (ierror.ne.0) call writerrs('MPI_COMM_SIZE: ',ierror)
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_COMM_RANK(comm,rank,ierror)
c determine the rank of the calling process in the communicator
c comm = communicator
c rank = rank of the calling process in group of comm
c ierror = error indicator
c input: comm
c output: rank, ierror
      implicit none
      integer comm, rank, ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c idproc = processor id
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c local data
      integer np
      ierror = 0
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c get rank
      else
         np = mapcomm(MAXS+1,comm+1)
         if ((np.gt.0).and.(np.le.nproc)) then
            rank = mapcomm(MAXS+2,comm+1)
            if ((rank.ge.0).and.(rank.lt.np)) then
               if (mapcomm(rank+1,comm+1).ne.idproc) then
                  ierror = 29
               endif
            else
               ierror = 29
            endif
         else
            ierror = 2
         endif
      endif
c handle errors
      if (ierror.ne.0) call writerrs('MPI_COMM_RANK: ',ierror)
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_COMM_DUP(comm,newcomm,ierror)
c duplicate existing communicator with all its cached information
c comm = communicator
c newcomm = new communicator
c ierror = error indicator
c input: comm
c output: newcomm, ierror
      implicit none
      integer comm, newcomm, ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c MPI constants
      integer MPI_COMM_NULL, MPI_MAX, MPI_MIN, MPI_INTEGER
      parameter(MPI_COMM_NULL=-1,MPI_MAX=0,MPI_MIN=1,MPI_INTEGER=18)
c local data
      integer np, i, j, k
      ierror = 0
      newcomm = MPI_COMM_NULL
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c communicator errors
      else
         np = mapcomm(MAXS+1,comm+1)
c communicator not mapped
         if ((np.le.0).or.(np.gt.nproc)) then
            ierror = 2
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_COMM_DUP: ',ierror)
         return
      endif
      if (monitor.eq.2) write (2,*) 'MPI_COMM_DUP started'
c find space for communication record
      i = 2
   10 i = i + 1
      if (i.gt.MAXC) then
         write (2,*) 'too many communicators'
         ierror = 25
         call writerrs('MPI_COMM_DUP: ',ierror)
         return
      elseif (mapcomm(MAXS+1,i).gt.0) then
         go to 10
      endif
c check if all nodes agree on new communicator
      call MPI_ALLREDUCE(i,j,1,MPI_INTEGER,MPI_MIN,comm,ierror)
      call MPI_ALLREDUCE(i,k,1,MPI_INTEGER,MPI_MAX,comm,ierror)
      if (j.ne.k) then
c try to find another communicator
         i = k - 1
         go to 10
      endif
c duplicate mapping
      do 20 j = 1, MAXS+MAXD+3
      mapcomm(j,i) = mapcomm(j,comm+1)
   20 continue
c assign communicator
      newcomm = i - 1
      if (monitor.eq.2) write (2,*) 'MPI_COMM_DUP complete'
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_COMM_SPLIT(comm,color,key,newcomm,ierror)
c create a new communicator based on color and key
c comm = communicator
c color = control of subset assignment
c key = control of rank assignment
c newcomm = new communicator
c ierror = error indicator
c input: comm, color, key
c output: newcomm, ierror
      implicit none
      integer comm, color, key, newcomm, ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c idproc = processor id
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c MPI constants
      integer MPI_MAX, MPI_MIN
      parameter(MPI_MAX=0,MPI_MIN=1)
      integer MPI_UNDEFINED, MPI_COMM_NULL, MPI_INTEGER
      parameter(MPI_UNDEFINED=-1,MPI_COMM_NULL=-1,MPI_INTEGER=18)
c local data
      integer np, i, mp, j, k, l, kmin
      ierror = 0
      newcomm = MPI_COMM_NULL
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c invalid color
      elseif (color.lt.(-1)) then
         write (2,*) 'Invalid color = ', color
         ierror = 23
c communicator errors
      else
         np = mapcomm(MAXS+1,comm+1)
c communicator not mapped
         if ((np.le.0).or.(np.gt.nproc)) then
            ierror = 2
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_COMM_SPLIT: ',ierror)
         return
      endif
      if (monitor.eq.2) write (2,*) 'MPI_COMM_SPLIT started'
c find space for communication record
      i = 2
   10 i = i + 1
      if (i.gt.MAXC) then
         write (2,*) 'too many communicators, color = ', color
         ierror = 25
         call writerrs('MPI_COMM_SPLIT: ',ierror)
         return
      elseif (mapcomm(MAXS+1,i).gt.0) then
         go to 10
      endif
c check if all nodes agree on new communicator
      call MPI_ALLREDUCE(i,j,1,MPI_INTEGER,MPI_MIN,comm,ierror)
      call MPI_ALLREDUCE(i,k,1,MPI_INTEGER,MPI_MAX,comm,ierror)
      if (j.ne.k) then
c try to find another communicator
         i = k - 1
         go to 10
      endif
c gather all the colors
      call MPI_ALLGATHER(color,1,MPI_INTEGER,mapcomm(1,i),1,MPI_INTEGER,
     1comm,ierror)
c this node does not participate
      if (color.eq.MPI_UNDEFINED) then
         if (monitor.eq.2) write (2,*) 'MPI_COMM_SPLIT complete'
         return
      endif
c find other processors with same color
      mp = 0
      mapcomm(MAXS+2,i) = -1
      do 20 j = 1, np
      if (mapcomm(j,i).eq.color) then
         mp = mp + 1
         k = mapcomm(j,comm+1)
         if ((k.ge.0).or.(k.lt.np)) then
            mapcomm(mp,i) = k
            if (k.eq.idproc) mapcomm(MAXS+2,i) = mp - 1
         else
            ierror = 2
         endif
      endif
   20 continue
c no nodes with found with color
      if (mp.eq.0) then
         write (2,*) 'Self color not found'
         ierror = 2
c no nodes found with idproc
      else if (mapcomm(MAXS+2,i).lt.0) then
         write (2,*) 'Self rank not found'
         ierror = 2
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_COMM_SPLIT: ',ierror)
         return
      endif
c finish mapping
      do 30 j = mp+1, MAXS
      mapcomm(j,i) = MPI_UNDEFINED
   30 continue
      mapcomm(MAXS+1,i) = mp
c assign communicator
      newcomm = i - 1
c gather all the keys into MPI_COMM_SELF mapping
      call MPI_ALLGATHER(key,1,MPI_INTEGER,mapcomm(1,2),1,MPI_INTEGER,ne
     1wcomm,ierror)
      k = 0
c find lowest remaining key
   40 kmin = mapcomm(k+1,2)
      do 50 j = k+2, mp
      if (mapcomm(j,2).lt.kmin) kmin = mapcomm(j,2)
   50 continue
c order all nodes with lowest remaining key
      do 70 j = k+1, mp
      if (mapcomm(j,2).eq.kmin) then
         k = k + 1
c right shift node and key order, if necessary
         if (j.gt.k) then
            np = mapcomm(j,i)
            do 60 l = 1, j-k
            mapcomm(j-l+1,i) = mapcomm(j-l,i)
            mapcomm(j-l+1,2) = mapcomm(j-l,2)
   60       continue
            mapcomm(k,i) = np
         endif
      endif
   70 continue
      if (k.lt.mp) go to 40
c find self rank
      mapcomm(MAXS+2,i) = -1
      do 80 j = 1, mp
      k = mapcomm(j,i)
      if (k.eq.idproc) mapcomm(MAXS+2,i) = j - 1
   80 continue
      mapcomm(MAXS+3,i) = 0
c restore MPI_COMM_SELF map
      mapcomm(1,2) = idproc
      do 90 j = 2, mp
      mapcomm(j,2) = MPI_UNDEFINED
   90 continue
      if (monitor.eq.2) write (2,*) 'MPI_COMM_SPLIT complete'
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_COMM_FREE(comm,ierror)
c mark the communicator object for deallocation
c comm = communicator
c ierror = error indicator
c input: comm
c output: comm, ierror
      implicit none
      integer comm, ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c MPI constants
      integer MPI_COMM_NULL, MPI_UNDEFINED
      parameter(MPI_COMM_NULL=-1,MPI_UNDEFINED=-1)
c local data
      integer np, i
      ierror = 0
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c release communicator
      elseif (comm.gt.1) then
         np = mapcomm(MAXS+1,comm+1)
         if ((np.gt.0).and.(np.le.nproc)) then
c clear mapping
            do 10 i = 1, MAXS
            mapcomm(i,comm+1) = 0
   10       continue
            mapcomm(MAXS+1,comm+1) = 0
            mapcomm(MAXS+2,comm+1) = MPI_UNDEFINED
            mapcomm(MAXS+3,comm+1) = 0
            comm = MPI_COMM_NULL
         else
            ierror = 2
         endif
      endif
c handle errors
      if (ierror.ne.0) call writerrs('MPI_COMM_FREE: ',ierror)
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_CART_CREATE(comm_old,ndims,dims,periods,reorder,com
     1m_cart,ierror)
c make new communicator to which topology information has been attached
c comm_old = input communicator
c ndims = number of dimensions of Cartesian grid
c dims = array specifying the number of processes in each dimension
c periods = array specifying whether grid is periodic or not
c reorder = specifies whether ranks may be reordered or not (ignored)
c comm_cart = communicator with new Carteisan topology
c ierror = error indicator
c input: comm_old, ndims, dims, periods, reorder
c output: comm_cart, ierror
      implicit none
      integer dims(*), comm_old, ndims, comm_cart, ierror
      logical periods(*), reorder
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c MPI constants
      integer MPI_MAX, MPI_MIN
      parameter(MPI_MAX=0,MPI_MIN=1)
      integer MPI_UNDEFINED, MPI_COMM_NULL, MPI_INTEGER
      parameter(MPI_UNDEFINED=-1,MPI_COMM_NULL=-1,MPI_INTEGER=18)
c local data
      integer np, mp, i, j, k
      ierror = 0
      comm_cart = MPI_COMM_NULL
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm_old.lt.0).or.(comm_old.ge.MAXC)) then
         ierror = 2
c invalid topology
      elseif ((ndims.lt.0).or.(ndims.gt.MAXD)) then
          write (2,*) 'Invalid topology dimension = ', ndims
          ierror = 26
c communicator errors
      else
         np = mapcomm(MAXS+1,comm_old+1)
c communicator not mapped
         if ((np.le.0).or.(np.gt.nproc)) then
            ierror = 2
c check topology length
         else
            mp = 1
            do 10 i = 1, ndims
            mp = mp*dims(i)
   10       continue
            if (ndims.eq.0) mp = 0
            if ((mp.lt.0).or.(mp.gt.np)) then
               write (2,*) 'Invalid Cartesian topology size = ', mp
               ierror = 27
            endif
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_CART_CREATE: ',ierror)
         return
      endif
      if (mp.eq.0) return
      if (monitor.eq.2) write (2,*) 'MPI_CART_CREATE started'
c find space for communication record
      i = 2
   20 i = i + 1
      if (i.gt.MAXC) then
         write (2,*) 'too many communicators'
         ierror = 25
         call writerrs('MPI_CART_CREATE: ',ierror)
         return
      elseif (mapcomm(MAXS+1,i).gt.0) then
         go to 20
      endif
c check if all nodes agree on new communicator
      call MPI_ALLREDUCE(i,j,1,MPI_INTEGER,MPI_MIN,comm_old,ierror)
      call MPI_ALLREDUCE(i,k,1,MPI_INTEGER,MPI_MAX,comm_old,ierror)
      if (j.ne.k) then
c try to find another communicator
         i = k - 1
         go to 20
      endif
c quit if processor is beyond range of topology
      if (mapcomm(MAXS+2,comm_old+1).ge.mp) then
         if (monitor.eq.2) write (2,*) 'MPI_CART_CREATE complete'
         return
      endif
c create mapping
      do 30 j = 1, mp
      mapcomm(j,i) = mapcomm(j,comm_old+1)
   30 continue
      do 40 j = mp+1, MAXS
      mapcomm(j,i) = MPI_UNDEFINED
   40 continue
      mapcomm(MAXS+1,i) = mp
      mapcomm(MAXS+2,i) = mapcomm(MAXS+2,comm_old+1)
c store topology
      mapcomm(MAXS+3,i) = ndims
      do 50 j = 1, ndims
      if (periods(j)) then
         mapcomm(MAXS+3+j,i) = dims(j)
      else
         mapcomm(MAXS+3+j,i) = -dims(j)
      endif
   50 continue
c assign communicator
      comm_cart = i - 1
      if (monitor.eq.2) write (2,*) 'MPI_CART_CREATE complete'
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_CART_COORDS(comm,rank,maxdims,coords,ierror)
c determine process coords in Cartesian topology, given rank in group
c comm = communicator with Cartesian structure
c rank = rank of a process within group of comm
c maxdims = length of vector coords in the calling program
c coords = array containing Cartesian coordinates of specified process
c ierror = error indicator
c input: comm, rank, maxdims
c output: coords, ierror
      implicit none
      integer coords(*), comm, rank, maxdims, ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c local data
      integer np, ndims, i, j
      ierror = 0
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c communicator errors
      else
         np = mapcomm(MAXS+1,comm+1)
         ndims = mapcomm(MAXS+3,comm+1)
c communicator not mapped
         if ((np.le.0).or.(np.gt.nproc)) then
            ierror = 2
c invalid topology
         elseif ((ndims.lt.1).or.(ndims.gt.MAXD)) then
            write (2,*) 'Invalid topology dimension = ', ndims
            ierror = 26
c invalid vector length
         elseif (maxdims.lt.ndims) then
             write (2,*) 'Vector length too small = ', maxdims
             ierror = 28
c invalid rank
         else
            if ((rank.lt.0).or.(rank.ge.np)) then
               write (2,*) 'Invalid rank = ', rank
               ierror = 29
            endif
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_CART_COORDS: ',ierror)
         return
      endif
c calculate cartesian coordinates
      j = rank
      do 10 i = 1, ndims
      np = np/abs(mapcomm(MAXS+3+i,comm+1))
      coords(i) = j/np
      j = j - coords(i)*np
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_CART_GET(comm,maxdims,dims,periods,coords,ierror)
c retrieve cartesian topology information associated with communicator
c comm = communicator with Cartesian structure
c maxdims = length of vectors dims, periods and coords
c dims = number of processes for each Cartesian dimension
c periods = periodicity (true/false) for each Cartesian dimension
c coords = array containing Cartesian coordinates of specified process
c ierror = error indicator
c input: comm, maxdims
c output: dims, periods, coords, ierror
      implicit none
      integer dims(*), coords(*), comm, maxdims, ierror
      logical periods(*)
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c idproc = processor id
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c local data
      integer np, ndims, rank, i, j, k
      ierror = 0
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c communicator errors
      else
         np = mapcomm(MAXS+1,comm+1)
         ndims = mapcomm(MAXS+3,comm+1)
         rank = mapcomm(MAXS+2,comm+1)
c communicator not mapped
         if ((np.le.0).or.(np.gt.nproc)) then
            ierror = 2
c invalid topology
         elseif ((ndims.lt.1).or.(ndims.gt.MAXD)) then
            write (2,*) 'Invalid topology dimension = ', ndims
            ierror = 26
c invalid vector length
         elseif (maxdims.lt.ndims) then
             write (2,*) 'Vector length too small = ', maxdims
             ierror = 28
c get rank
         elseif ((rank.ge.0).and.(rank.lt.np)) then
            if (mapcomm(rank+1,comm+1).ne.idproc) then
               ierror = 29
            endif
         else
            ierror = 29
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_CART_GET: ',ierror)
         return
      endif
c calculate dimension, periodicity, and cartesian coordinates
      j = rank
      do 10 i = 1, ndims
      k = mapcomm(MAXS+3+i,comm+1)
      if (k.gt.0) then
         periods(i) = .true.
      else
         periods(i) = .false.
         k = -k
      endif
      dims(i) = k
      np = np/k
      coords(i) = j/np
      j = j - coords(i)*np
   10 continue
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_CART_SHIFT(comm,direction,disp,rank_source,rank_des
     1t,ierror)
c return shifted source and destination ranks given shift direction and
c amount
c comm = communicator with Cartesian structure
c direction = coordinate dimension shift
c disp = displacement (> 0: upwards shift, < 0: downwards shift)
c rank_source = rank of source process
c rank_dest = rank of destination process
c ierror = error indicator
c input: comm, direction, disp
c output: rank_source, rank_dest, ierror
      implicit none
      integer comm, direction, disp, rank_source, rank_dest, ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c MPI constants
      integer MPI_PROC_NULL
      parameter(MPI_PROC_NULL=-3)
c local data
      integer np, ndims, rank, i, j, k, l, shift
      ierror = 0
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c communicator errors
      else
         np = mapcomm(MAXS+1,comm+1)
         ndims = mapcomm(MAXS+3,comm+1)
         rank = mapcomm(MAXS+2,comm+1)
c communicator not mapped
         if ((np.le.0).or.(np.gt.nproc)) then
            ierror = 2
c invalid topology
         elseif ((ndims.lt.1).or.(ndims.gt.MAXD)) then
            write (2,*) 'Invalid topology dimension = ', ndims
            ierror = 26
c invalid direction
         elseif ((direction.lt.0).or.(direction.ge.ndims)) then
             write (2,*) 'Invalid direction = ', direction
             ierror = 31
c get rank
         elseif ((rank.ge.0).and.(rank.lt.np)) then
            if (mapcomm(rank+1,comm+1).ne.idproc) then
               ierror = 29
            endif
         else
            ierror = 29
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_CART_SHIFT: ',ierror)
         return
      endif
c check if shift amount is zero
      if (disp.eq.0) then
         rank_source = rank
         rank_dest = rank
         return
      endif
c find coordinate for selected direction
      j = rank
      shift = np
      do 10 i = 1, direction+1
      shift = shift/abs(mapcomm(MAXS+3+i,comm+1))
      k = j/shift
      j = j - k*shift
   10 continue
c calculate size of shift
      shift = 1
      do 20 i = direction+2, ndims
      shift = shift*abs(mapcomm(MAXS+3+i,comm+1))
   20 continue
c size of selected direction
      l = mapcomm(MAXS+4+direction,comm+1)
c calculate new coordinate
c periodic boundary conditions
      if (l.gt.0) then
c make disp also periodic
         i = mod(abs(disp),l)
         if (disp.lt.0) i = -i
c right neighbor
         j = k + i
         if (j.lt.0) then
            j = j + l
         elseif (j.ge.l) then
            j = j - l
         endif
         rank_dest = rank + (j - k)*shift
c left neighbor
         j = k - i
         if (j.lt.0) then
            j = j + l
         elseif (j.ge.l) then
            j = j - l
         endif
         rank_source = rank + (j - k)*shift
c non-periodic boundary conditions
      elseif (l.lt.0) then
c right neighbor
         j = k + disp
         if ((j.lt.0).or.(j.ge.(-l))) then
            rank_dest = MPI_PROC_NULL
         else
            rank_dest = rank + (j - k)*shift
         endif
c left neighbor
         j = k - disp
         if ((j.lt.0).or.(j.ge.(-l))) then
            rank_source = MPI_PROC_NULL
         else
            rank_source = rank + (j - k)*shift
         endif
      endif
c verify ranks
      if (rank_source.ne.MPI_PROC_NULL) then
         if ((rank_source.lt.0).or.(rank_source.ge.np)) then
            write (2,*) 'rank_source = ', rank_source
            ierror = 29
         endif
      endif
      if (rank_dest.ne.MPI_PROC_NULL) then
         if ((rank_dest.lt.0).or.(rank_dest.ge.np)) then
            write (2,*) 'rank_dest = ', rank_dest
            ierror = 29
         endif
      endif
c process errors
      if (ierror.ne.0) call writerrs('MPI_CART_SHIFT: ',ierror)
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_CART_RANK(comm,coords,rank,ierror)
c determine process rank in communicator, given Cartesian location
c comm = communicator with Cartesian structure
c coords = array specifying Cartesian coordinates of a process
c rank = rank of specified process
c ierror = error indicator
c input: comm, coords
c output: rank, ierror
      implicit none
      integer coords(*), comm, rank, ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c local data
      integer np, ndims, i, j, k, l
      ierror = 0
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c communicator errors
      else
         np = mapcomm(MAXS+1,comm+1)
         ndims = mapcomm(MAXS+3,comm+1)
c communicator not mapped
         if ((np.le.0).or.(np.gt.nproc)) then
            ierror = 2
c invalid topology
         elseif ((ndims.lt.1).or.(ndims.gt.MAXD)) then
            write (2,*) 'Invalid topology dimension = ', ndims
            ierror = 26
c invalid coords
         else
            do 10 i = 1, ndims
            j = mapcomm(MAXS+3+i,comm+1)
            if (j.lt.0) then
               if ((coords(i).lt.0).or.(coords(i).ge.(-j))) then
                  write (2,*) 'Invalid ith coord = ', i, coords(i)
                  ierror = 30
               endif
            endif
   10       continue
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_CART_RANK: ',ierror)
         return
      endif
c calculate rank, wrap periodic coordinates
      l = 0
      do 20 i = 1, ndims
      j = mapcomm(MAXS+3+i,comm+1)
      k = coords(i)
      if (j.gt.0) then
         if (k.lt.0) then
            k = k + j
         elseif (k.ge.j) then
            k = k - j
         endif
      else
         j = -j
      endif
      l = k + j*l
   20 continue
c verify rank
      if ((l.lt.0).or.(l.ge.np)) then
         write (2,*) 'Invalid coords, resulting rank = ', l
         ierror = 30
         call writerrs('MPI_CART_RANK: ',ierror)
      else
         rank = l
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_CART_SUB(comm,remain_dims,newcomm,ierror)
c partition communicator into subgroups that form lower-dimensional
c Cartesian subgrids
c comm = communicator with Cartesian structure
c remain_dims = each entry of remain_dims specifies whether dimension is
c kept in the subgrid or not
c newcomm = communicator containing the subgrid
c ierror = error indicator
c input: comm, remain_dims
c output: newcomm, ierror
      implicit none
      integer comm, newcomm, ierror
      logical remain_dims(*)
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c idproc = processor id
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c MPI constants
      integer MPI_MAX, MPI_MIN
      parameter(MPI_MAX=0,MPI_MIN=1)
      integer MPI_UNDEFINED, MPI_COMM_NULL, MPI_INTEGER
      parameter(MPI_UNDEFINED=-1,MPI_COMM_NULL=-1,MPI_INTEGER=18)
c local data
      integer np, ndims, rank, i, j, k, l, m, mp, color, key
      ierror = 0
      newcomm = MPI_COMM_NULL
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c communicator errors
      else
         np = mapcomm(MAXS+1,comm+1)
         ndims = mapcomm(MAXS+3,comm+1)
         rank = mapcomm(MAXS+2,comm+1)
c communicator not mapped
         if ((np.le.0).or.(np.gt.nproc)) then
            ierror = 2
c invalid topology
         elseif ((ndims.lt.1).or.(ndims.gt.MAXD)) then
            write (2,*) 'Invalid topology dimension = ', ndims
            ierror = 26
c get rank
         elseif ((rank.ge.0).and.(rank.lt.np)) then
            if (mapcomm(rank+1,comm+1).ne.idproc) then
               ierror = 29
            endif
         else
            ierror = 29
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_CART_SUB: ',ierror)
         return
      endif
c find dimension of new topology
      k = 0
      do 10 j = 1, ndims
      if (remain_dims(j)) then
         k = k + 1
      endif
   10 continue
c empty topology
      if (k.eq.0) return
      if (monitor.eq.2) write (2,*) 'MPI_CART_SUB started'
c find space for communication record
      i = 2
   20 i = i + 1
      if (i.gt.MAXC) then
         write (2,*) 'too many communicators'
         ierror = 25
         call writerrs('MPI_CART_SUB: ',ierror)
         return
      elseif (mapcomm(MAXS+1,i).gt.0) then
         go to 20
      endif
c check if all nodes agree on new communicator
      call MPI_ALLREDUCE(i,j,1,MPI_INTEGER,MPI_MIN,comm,ierror)
      call MPI_ALLREDUCE(i,k,1,MPI_INTEGER,MPI_MAX,comm,ierror)
      if (j.ne.k) then
c try to find another communicator
         i = k - 1
         go to 20
      endif
c find color for missing dimension
      j = rank
      color = 0
      key = 0
      mp = np
      do 30 l = 1, ndims
      m = abs(mapcomm(MAXS+3+l,comm+1))
      mp = mp/m
      k = j/mp
      j = j - k*mp
      if (remain_dims(l)) then
         key = k + key*m
      else
         color = k + color*m
      endif
   30 continue
c gather all the colors
      call MPI_ALLGATHER(color,1,MPI_INTEGER,mapcomm(1,i),1,MPI_INTEGER,
     1comm,ierror)
c find other processors with same color
      mp = 0
      mapcomm(MAXS+2,i) = -1
      do 40 j = 1, np
      if (mapcomm(j,i).eq.color) then
         mp = mp + 1
         k = mapcomm(j,comm+1)
         if ((k.ge.0).or.(k.lt.np)) then
            mapcomm(mp,i) = k
            if (k.eq.idproc) mapcomm(MAXS+2,i) = mp - 1
         else
            ierror = 2
         endif
      endif
   40 continue
c no nodes with found with color
      if (mp.eq.0) then
         write (2,*) 'Self color not found'
         ierror = 2
c no nodes found with idproc
      else if (mapcomm(MAXS+2,i).lt.0) then
         write (2,*) 'Self rank not found'
         ierror = 2
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_CART_SUB: ',ierror)
         return
      endif
c finish mapping
      do 50 j = mp+1, MAXS
      mapcomm(j,i) = MPI_UNDEFINED
   50 continue
      mapcomm(MAXS+1,i) = mp
c assign communicator
      newcomm = i - 1
c create new topology
      k = 0
      do 60 j = 1, ndims
      if (remain_dims(j)) then
         k = k + 1
         mapcomm(MAXS+3+k,i) = mapcomm(MAXS+3+j,comm+1)
      endif
   60 continue
      mapcomm(MAXS+3,i) = k
      if (monitor.eq.2) write (2,*) 'MPI_CART_SUB complete'
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_DIMS_CREATE(nnodes,ndims,dims,ierror)
c create a division of processes in a Cartesian grid
c nnodes = number of nodes in a grid
c ndims = number of Cartesian dimensions
c dims = array specifying number of nodes in each dimension
c ierror = error indicator
c input: nnodes, ndims, dims
c output: dims, ierror
      implicit none
      integer dims(*), nnodes, ndims, ierror
c local data
      integer MAXS
      parameter(MAXS=32)
      integer i, j, nd, mp, md
      ierror = 0
c check for errors
      nd = 0
      mp = 1
      do 10 i = 1, ndims
      if (dims(i).eq.0) then
         nd = nd + 1
      elseif (dims(i).gt.0) then
         mp = mp*dims(i)
      elseif (dims(i).lt.0) then
         ierror = 26
      endif
   10 continue
      if (nd.eq.0) return
      if (ierror.ne.0) then
         write (2,*) 'Invalid topology dimension'
      elseif ((nnodes.lt.1).or.(nnodes.gt.MAXS)) then
         write (2,*) 'Invalid number of nodes = ', nnodes
         ierror = 24
      else
         md = nnodes/mp
         if ((md*mp).ne.nnodes) then
            write (2,*) 'Topology dimension, node number inconsistent'
            ierror = 26
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         write (2,*) 'MPI_DIMS_CREATE: Error code = ', ierror
         return
      endif
c look for nd factors of md
      do 30 i = 1, ndims
      if (dims(i).eq.0) then
         mp = exp(alog(real(md)/real(nd))) + .0001
         if ((mp**nd).lt.md) mp = mp + 1
   20    j = md/mp
         if ((j*mp).eq.md) then
            dims(i) = mp
            md = j
            nd = nd - 1
         else
            mp = mp + 1
            if (mp.le.md) go to 20
            write (2,*) 'MPI_DIMS_CREATE: factor not found'
            ierror = 26
         endif
      endif
   30 continue
c sanity check
      if ((md.ne.1).or.(nd.ne.0)) then
         write (2,*) 'MPI_DIMS_CREATE: missing factors'
         ierror = 26
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_BCAST(buffer,count,datatype,root,comm,ierror)
c broadcast a message from root to all processes in comm
c buffer = starting address of buffer
c count = number of entries in buffer
c datatype = datatype of buffer
c root = rank of broadcast root
c comm = communicator
c ierror = error indicator
c input: buffer, count, datatype, root, comm
c output: buffer, ierror
      implicit none
      integer buffer(*)
      integer count, datatype, root, comm, ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c idproc = processor id
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c MPI constants
      integer MPI_STATUS_SIZE
      parameter(MPI_STATUS_SIZE=5)
c local data
      integer i, np, id, rank, status
      dimension status(MPI_STATUS_SIZE)
      ierror = 0
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c invalid root
      elseif ((root.lt.0).or.(root.ge.nproc)) then
         ierror = 19
c communicator errors
      else
         np = mapcomm(MAXS+1,comm+1)
         id = mapcomm(root+1,comm+1)
c communicator not mapped
         if ((np.le.0).or.(np.gt.nproc)) then
            ierror = 2
c invalid root
         elseif ((root.lt.0).or.(root.ge.np)) then
            ierror = 19
c invalid mapping
         elseif ((id.lt.0).or.(id.ge.nproc)) then
            write (2,*) 'Invalid mapping, root, node = ', root, id
            ierror = 2
c get rank
         else
            rank = mapcomm(MAXS+2,comm+1)
            if ((rank.ge.0).and.(rank.lt.np)) then
               if (mapcomm(rank+1,comm+1).ne.idproc) then
                  ierror = 29
               endif
            else
               ierror = 29
            endif
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_BCAST: ',ierror)
         return
      endif
      if (monitor.eq.2) write (2,*) 'MPI_BCAST started'
c start broadcast
      if (rank.eq.root) then
         do 10 i = 1, np
         id = i - 1
         if (id.ne.root) call MPI_SEND(buffer,count,datatype,id,0,comm,i
     1error)
   10    continue
      else
         call MPI_RECV(buffer,count,datatype,root,0,comm,status,ierror)
      endif
      if (monitor.eq.2) write (2,*) 'MPI_BCAST complete'
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_BARRIER(comm,ierror)
c blocks each process in comm until all processes have called it.
c comm = communicator
c ierror = error indicator
c input: comm
c output: ierror
      implicit none
      integer comm, ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c idproc = processor id
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c MPI constants
      integer MPI_STATUS_SIZE, MPI_INTEGER
      parameter(MPI_STATUS_SIZE=5,MPI_INTEGER=18)
c local data
      integer np, rank, ntasks, isync, irync, i, status
      dimension status(MPI_STATUS_SIZE)
      ierror = 0
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c communicator errors
      else
         np = mapcomm(MAXS+1,comm+1)
c communicator not mapped
         if ((np.le.0).or.(np.gt.nproc)) then
            ierror = 2
c get rank
         else
            rank = mapcomm(MAXS+2,comm+1)
            if ((rank.ge.0).and.(rank.lt.np)) then
               if (mapcomm(rank+1,comm+1).ne.idproc) then
                  ierror = 29
               endif
            else
               ierror = 29
            endif
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_BARRIER: ',ierror)
         return
      endif
      if (monitor.eq.2) write (2,*) 'MPI_BARRIER started'
c begin synchronization
      ntasks = np - 1
      isync = -1
      if (rank.eq.0) then
c processor 0 receives a message from everyone else
         do 10 i = 1, ntasks
         call MPI_RECV(irync,1,MPI_INTEGER,i,0,comm,status,ierror)
         if (irync.ne.isync) write (2,*) 'sync error from proc', i
   10    continue
c then sends an acknowledgment back
         isync = 1
         call MPI_BCAST(isync,1,MPI_INTEGER,0,comm,ierror)
      else
c remaining processors send a message to processor 0
         call MPI_SEND(isync,1,MPI_INTEGER,0,0,comm,ierror)
c then receive an acknowledgement back
         isync = 1
         call MPI_BCAST(irync,1,MPI_INTEGER,0,comm,ierror)
         if (irync.ne.isync) write (2,*) 'rsync error at proc', rank
      endif
      if (monitor.eq.2) write (2,*) 'MPI_BARRIER complete'
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_REDUCE(sendbuf,recvbuf,count,datatype,op,root,comm,
     1ierror)
c applies a reduction operation to the vector sendbuf over the set of
c processes specified by comm and places the result in recvbuf on root
c sendbuf = address of send buffer
c recvbuf = address of receive buffer
c count = number of elements in send buffer
c datatype = datatype of elements in send buffer
c op = reduce operation
c (only max, min, sum, maxloc, and minloc currently supported)
c root = rank of root process
c comm = communicator
c ierror = error indicator
c input: sendbuf, count, datatype, op, root, comm
c output: recvbuf, ierror
      implicit none
      integer sendbuf(*), recvbuf(*)
      integer count, datatype, op, root, comm, ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c idproc = processor id
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c function declarations
      integer NewPtr
      external NewPtr
c MPI constants
      integer MPI_STATUS_SIZE, MPI_SUM, MPI_MAXLOC, MPI_MINLOC
      parameter(MPI_STATUS_SIZE=5,MPI_SUM=2,MPI_MAXLOC=4,MPI_MINLOC=5)
      integer MPI_INTEGER, MPI_REAL, MPI_DOUBLE_PRECISION
      parameter(MPI_INTEGER=18,MPI_REAL=19,MPI_DOUBLE_PRECISION=20)
      integer MPI_COMPLEX, MPI_DOUBLE_COMPLEX
      parameter(MPI_COMPLEX=22,MPI_DOUBLE_COMPLEX=23)
      integer MPI_2REAL, MPI_2DOUBLE_PRECISION, MPI_2INTEGER
      parameter(MPI_2REAL=35,MPI_2DOUBLE_PRECISION=36,MPI_2INTEGER=37)
c local data
      integer status, tmpbuf
      integer i, j, np, rank, id, ltmp, loct, nl, lcnt, lsize
      dimension status(MPI_STATUS_SIZE)
      ierror = 0
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c invalid root
      elseif ((root.lt.0).or.(root.ge.nproc)) then
         ierror = 19
c invalid op
      elseif ((op.lt.0).or.(op.gt.5).or.(op.eq.3)) then
         ierror = 20
c communicator errors
      else
         np = mapcomm(MAXS+1,comm+1)
         id = mapcomm(root+1,comm+1)
c communicator not mapped
         if ((np.le.0).or.(np.gt.nproc)) then
            ierror = 2
c invalid root
         elseif ((root.lt.0).or.(root.ge.np)) then
            ierror = 19
c invalid mapping
         elseif ((id.lt.0).or.(id.ge.nproc)) then
            write (2,*) 'Invalid mapping, root, node = ', root, id
            ierror = 2
c get rank
         else
            rank = mapcomm(MAXS+2,comm+1)
            if ((rank.ge.0).and.(rank.lt.np)) then
               if (mapcomm(rank+1,comm+1).ne.idproc) then
                  ierror = 29
               endif
            else
               ierror = 29
            endif
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_REDUCE: ',ierror)
         return
      endif
c determine size of temporary buffer
      if ((datatype.eq.MPI_INTEGER).or.(datatype.eq.MPI_REAL)) then
         lsize = 4
      elseif (datatype.eq.MPI_DOUBLE_PRECISION) then
         lsize = 8
      elseif ((datatype.eq.MPI_COMPLEX).and.(op.eq.MPI_SUM)) then
         lsize = 8
      elseif ((datatype.eq.MPI_DOUBLE_COMPLEX).and.(op.eq.MPI_SUM)) then
         lsize = 16
      elseif ((datatype.eq.MPI_2INTEGER).and.((op.eq.MPI_MAXLOC).or.(op.
     1eq.MPI_MINLOC))) then
         lsize = 8
      elseif ((datatype.eq.MPI_2REAL).and.((op.eq.MPI_MAXLOC).or.(op.eq.
     1MPI_MINLOC))) then
         lsize = 8
      elseif ((datatype.eq.MPI_2DOUBLE_PRECISION).and.((op.eq.MPI_MAXLOC
     1).or.(op.eq.MPI_MINLOC))) then
         lsize = 16
c invalid datatype
      else
         ierror = 7
         call writerrs('MPI_REDUCE: ',ierror)
         return
      endif
      if (rank.eq.root) then
         loct = 0
c initialize by copying from send to receive buffer
         if (datatype.eq.MPI_INTEGER) then
            call iredux(recvbuf,sendbuf,loct,count,-1)
         elseif (datatype.eq.MPI_REAL) then
            call fredux(recvbuf,sendbuf,loct,count,-1)
         elseif (datatype.eq.MPI_DOUBLE_PRECISION) then
            call dredux(recvbuf,sendbuf,loct,count,-1)
         elseif (datatype.eq.MPI_COMPLEX) then
            call fredux(recvbuf,sendbuf,loct,2*count,-1)
         elseif (datatype.eq.MPI_DOUBLE_COMPLEX) then
            call dredux(recvbuf,sendbuf,loct,2*count,-1)
         elseif (datatype.eq.MPI_2INTEGER) then
            call iredux(recvbuf,sendbuf,loct,2*count,-1)
         elseif (datatype.eq.MPI_2REAL) then
            call fredux(recvbuf,sendbuf,loct,2*count,-1)
         elseif (datatype.eq.MPI_2DOUBLE_PRECISION) then
            call dredux(recvbuf,sendbuf,loct,2*count,-1)
         endif
      else
         loct = 1
      endif
      ltmp = lsize*count
c allocate a nonrelocatable block of memory
      tmpbuf = NewPtr(val4(ltmp))
c memory not available
      if (tmpbuf.eq.0) then
         ierror = 21
         call writerrs('MPI_REDUCE: ',ierror)
         return
      endif
      if (monitor.eq.2) write (2,*) 'MPI_REDUCE started'
      ltmp = ltmp/lsize
c send messages in groups of ltmp
      nl = (count - 1)/ltmp + 1
      lcnt = ltmp
      lsize = lsize*ltmp/4
      do 20 j = 1, nl
      if (j.eq.nl) lcnt = count - ltmp*(nl - 1)
      if (rank.eq.root) then
c root receives data from everyone else
         do 10 i = 1, np
         id = i - 1
         if (id.ne.root) then
            call MPI_RECV(val4(tmpbuf),lcnt,datatype,id,j,comm,status,ie
     1rror)
c reduce data
            if (datatype.eq.MPI_INTEGER) then
               call iredux(recvbuf,val4(tmpbuf),loct,lcnt,op)
            elseif (datatype.eq.MPI_REAL) then
               call fredux(recvbuf,val4(tmpbuf),loct,lcnt,op)
            elseif (datatype.eq.MPI_DOUBLE_PRECISION) then
               call dredux(recvbuf,val4(tmpbuf),loct,lcnt,op)
            elseif (datatype.eq.MPI_COMPLEX) then
               call fredux(recvbuf,val4(tmpbuf),loct,2*lcnt,op)
            elseif (datatype.eq.MPI_DOUBLE_COMPLEX) then
               call dredux(recvbuf,val4(tmpbuf),loct,2*lcnt,op)
            elseif (datatype.eq.MPI_2INTEGER) then
               call iredux(recvbuf,val4(tmpbuf),loct,lcnt,op)
            elseif (datatype.eq.MPI_2REAL) then
               call fredux(recvbuf,val4(tmpbuf),loct,lcnt,op)
            elseif (datatype.eq.MPI_2DOUBLE_PRECISION) then
               call dredux(recvbuf,val4(tmpbuf),loct,lcnt,op)
            endif
         endif
   10    continue
         loct = loct + ltmp
      else
c remaining processors send data to root
         call MPI_SEND(sendbuf(loct),lcnt,datatype,root,j,comm,ierror)
         loct = loct + lsize
      endif
   20 continue
c release nonrelocatable memory block
      call DisposePtr(val4(tmpbuf))
      if (monitor.eq.2) write (2,*) 'MPI_REDUCE complete'
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_SCAN(sendbuf,recvbuf,count,datatype,op,comm,ierror)
c performs a parallel prefix reduction on data distributed across a
c group
c sendbuf = address of send buffer
c recvbuf = address of receive buffer
c count = number of elements in send buffer
c datatype = datatype of elements in send buffer
c op = reduce operation
c (only max, min, sum, maxloc, and minloc currently supported)
c comm = communicator
c ierror = error indicator
c input: sendbuf, count, datatype, op, comm
c output: recvbuf, ierror
      implicit none
      integer sendbuf(*), recvbuf(*)
      integer count, datatype, op, comm, ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c idproc = processor id
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c function declarations
      integer NewPtr
      external NewPtr
c MPI constants
      integer MPI_STATUS_SIZE, MPI_SUM, MPI_MAXLOC, MPI_MINLOC
      parameter(MPI_STATUS_SIZE=5,MPI_SUM=2,MPI_MAXLOC=4,MPI_MINLOC=5)
      integer MPI_INTEGER, MPI_REAL, MPI_DOUBLE_PRECISION
      parameter(MPI_INTEGER=18,MPI_REAL=19,MPI_DOUBLE_PRECISION=20)
      integer MPI_COMPLEX, MPI_DOUBLE_COMPLEX
      parameter(MPI_COMPLEX=22,MPI_DOUBLE_COMPLEX=23)
      integer MPI_2REAL, MPI_2DOUBLE_PRECISION, MPI_2INTEGER
      parameter(MPI_2REAL=35,MPI_2DOUBLE_PRECISION=36,MPI_2INTEGER=37
c local data
      integer status, tmpbuf
      integer i, j, np, rank, root, id, ltmp, loct, nl, lcnt, lsize
      dimension status(MPI_STATUS_SIZE)
      data root /0/
      ierror = 0
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c invalid root
      elseif ((root.lt.0).or.(root.ge.nproc)) then
         ierror = 19
c invalid op
      elseif ((op.lt.0).or.(op.gt.5).or.(op.eq.3)) then
         ierror = 20
c communicator errors
      else
         np = mapcomm(MAXS+1,comm+1)
         id = mapcomm(root+1,comm+1)
c communicator not mapped
         if ((np.le.0).or.(np.gt.nproc)) then
            ierror = 2
c invalid root
         elseif ((root.lt.0).or.(root.ge.np)) then
            ierror = 19
c invalid mapping
         elseif ((id.lt.0).or.(id.ge.nproc)) then
            write (2,*) 'Invalid mapping, root, node = ', root, id
            ierror = 2
c get rank
         else
            rank = mapcomm(MAXS+2,comm+1)
            if ((rank.ge.0).and.(rank.lt.np)) then
               if (mapcomm(rank+1,comm+1).ne.idproc) then
                  ierror = 29
               endif
            else
               ierror = 29
            endif
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_SCAN: ',ierror)
         return
      endif
c determine size of temporary buffer
      if ((datatype.eq.MPI_INTEGER).or.(datatype.eq.MPI_REAL)) then
         lsize = 4
      elseif (datatype.eq.MPI_DOUBLE_PRECISION) then
         lsize = 8
      elseif ((datatype.eq.MPI_COMPLEX).and.(op.eq.MPI_SUM)) then
         lsize = 8
      elseif ((datatype.eq.MPI_DOUBLE_COMPLEX).and.(op.eq.MPI_SUM)) then
         lsize = 16
      elseif ((datatype.eq.MPI_2INTEGER).and.((op.eq.MPI_MAXLOC).or.(op.
     1eq.MPI_MINLOC))) then
         lsize = 8
      elseif ((datatype.eq.MPI_2REAL).and.((op.eq.MPI_MAXLOC).or.(op.eq.
     1MPI_MINLOC))) then
         lsize = 8
      elseif ((datatype.eq.MPI_2DOUBLE_PRECISION).and.((op.eq.MPI_MAXLOC
     1).or.(op.eq.MPI_MINLOC))) then
         lsize = 16
c invalid datatype
      else
         ierror = 7
         call writerrs('MPI_SCAN: ',ierror)
         return
      endif
      if (rank.eq.root) then
         loct = 0
c initialize by copying from send to receive buffer
         if (datatype.eq.MPI_INTEGER) then
            call iredux(recvbuf,sendbuf,loct,count,-1)
         elseif (datatype.eq.MPI_REAL) then
            call fredux(recvbuf,sendbuf,loct,count,-1)
         elseif (datatype.eq.MPI_DOUBLE_PRECISION) then
            call dredux(recvbuf,sendbuf,loct,count,-1)
         elseif (datatype.eq.MPI_COMPLEX) then
            call fredux(recvbuf,sendbuf,loct,2*count,-1)
         elseif (datatype.eq.MPI_DOUBLE_COMPLEX) then
            call dredux(recvbuf,sendbuf,loct,2*count,-1)
         elseif (datatype.eq.MPI_2INTEGER) then
            call iredux(recvbuf,sendbuf,loct,2*count,-1)
         elseif (datatype.eq.MPI_2REAL) then
            call fredux(recvbuf,sendbuf,loct,2*count,-1)
         elseif (datatype.eq.MPI_2DOUBLE_PRECISION) then
            call dredux(recvbuf,sendbuf,loct,2*count,-1)
         endif
      else
         loct = 1
      endif
      ltmp = lsize*count
c allocate a nonrelocatable block of memory
      tmpbuf = NewPtr(val4(ltmp))
c memory not available
      if (tmpbuf.eq.0) then
         ierror = 21
         call writerrs('MPI_SCAN: ',ierror)
         return
      endif
      if (monitor.eq.2) write (2,*) 'MPI_SCAN started'
      ltmp = ltmp/lsize
c send messages in groups of ltmp
      nl = (count - 1)/ltmp + 1
      lcnt = ltmp
      lsize = lsize*ltmp/4
      do 20 j = 1, nl
      if (j.eq.nl) lcnt = count - ltmp*(nl - 1)
      if (rank.eq.root) then
c root receives data from everyone else
         do 10 i = 1, np
         id = i - 1
         if (id.ne.root) then
            call MPI_RECV(val4(tmpbuf),lcnt,datatype,id,j,comm,status,ie
     1rror)
c reduce data
            if (datatype.eq.MPI_INTEGER) then
               call iredux(recvbuf,val4(tmpbuf),loct,lcnt,op)
            elseif (datatype.eq.MPI_REAL) then
               call fredux(recvbuf,val4(tmpbuf),loct,lcnt,op)
            elseif (datatype.eq.MPI_DOUBLE_PRECISION) then
               call dredux(recvbuf,val4(tmpbuf),loct,lcnt,op)
            elseif (datatype.eq.MPI_COMPLEX) then
               call fredux(recvbuf,val4(tmpbuf),loct,2*lcnt,op)
            elseif (datatype.eq.MPI_DOUBLE_COMPLEX) then
               call dredux(recvbuf,val4(tmpbuf),loct,2*lcnt,op)
            elseif (datatype.eq.MPI_2INTEGER) then
               call iredux(recvbuf,val4(tmpbuf),loct,lcnt,op)
            elseif (datatype.eq.MPI_2REAL) then
               call fredux(recvbuf,val4(tmpbuf),loct,lcnt,op)
            elseif (datatype.eq.MPI_2DOUBLE_PRECISION) then
               call dredux(recvbuf,val4(tmpbuf),loct,lcnt,op)
            endif
c send partial result data to processor id
            call MPI_SEND(recvbuf(loct+1),lcnt,datatype,id,j+nproc,comm,
     1ierror)
         endif
   10    continue
         loct = loct + ltmp
      else
c remaining processors send data to root
         call MPI_SEND(sendbuf(loct),lcnt,datatype,root,j,comm,ierror)
c receive partial result data from root
         call MPI_RECV(recvbuf(loct),lcnt,datatype,root,j+nproc,comm,sta
     1tus,ierror)
         loct = loct + lsize
      endif
   20 continue
      if (rank.eq.root) then
c initialize by copying from send to receive buffer
         if (datatype.eq.MPI_INTEGER) then
            call iredux(recvbuf,sendbuf,0,count,-1)
         elseif (datatype.eq.MPI_REAL) then
            call fredux(recvbuf,sendbuf,0,count,-1)
         elseif (datatype.eq.MPI_DOUBLE_PRECISION) then
            call dredux(recvbuf,sendbuf,0,count,-1)
         elseif (datatype.eq.MPI_COMPLEX) then
            call fredux(recvbuf,sendbuf,0,2*count,-1)
         elseif (datatype.eq.MPI_DOUBLE_COMPLEX) then
            call dredux(recvbuf,sendbuf,0,2*count,-1)
         elseif (datatype.eq.MPI_2INTEGER) then
            call iredux(recvbuf,sendbuf,0,2*count,-1)
         elseif (datatype.eq.MPI_2REAL) then
            call fredux(recvbuf,sendbuf,0,2*count,-1)
         elseif (datatype.eq.MPI_2DOUBLE_PRECISION) then
            call dredux(recvbuf,sendbuf,0,2*count,-1)
         endif
      endif
c release nonrelocatable memory block
      call DisposePtr(val4(tmpbuf))
      if (monitor.eq.2) write (2,*) 'MPI_SCAN complete'
      return
      end
c-----------------------------------------------------------------------
      subroutine iredux(recvbuf,sendbuf,offset,count,op)
c perform reduction operation for integer types
c recvbuf = address of receive buffer
c sendbuf = address of send buffer
c offset = starting location minus one in receive buffer
c count = number of elements in send buffer
c op = reduce operation
c (only max, min, sum, maxloc, and minloc currently supported)
c input: recvbuf, sendbuf, offset, count, op
c output: recvbuf
      implicit none
      integer offset, count, op
      integer recvbuf(*), sendbuf(*)
c MPI constants
      integer MPI_MAX, MPI_MIN, MPI_SUM, MPI_MAXLOC, MPI_MINLOC
      parameter(MPI_MAX=0,MPI_MIN=1,MPI_SUM=2,MPI_MAXLOC=4,MPI_MINLOC=5)
c local data
      integer i, j, k
c perform reduction
      if (op.eq.MPI_MAX) then
         do 10 i = 1, count
         recvbuf(i+offset) = max(recvbuf(i+offset),sendbuf(i))
   10    continue
      elseif (op.eq.MPI_MIN) then
         do 20 i = 1, count
         recvbuf(i+offset) = min(recvbuf(i+offset),sendbuf(i))
   20    continue
      elseif (op.eq.MPI_SUM) then
         do 30 i = 1, count
         recvbuf(i+offset) = recvbuf(i+offset) + sendbuf(i)
   30    continue
c perform reduction and location
      elseif (op.eq.MPI_MAXLOC) then
         do 40 i = 1, count
         j = recvbuf(2*i+offset-1)
         k = sendbuf(2*i-1)
         if (j.lt.k) then
            recvbuf(2*i+offset-1) = sendbuf(2*i-1)
            recvbuf(2*i+offset) = sendbuf(2*i)
         elseif (j.eq.k) then
            recvbuf(2*i+offset) = min(recvbuf(2*i+offset),sendbuf(2*i))
         endif
   40    continue
      elseif (op.eq.MPI_MINLOC) then
         do 50 i = 1, count
         j = recvbuf(2*i+offset-1)
         k = sendbuf(2*i-1)
         if (j.gt.k) then
            recvbuf(2*i+offset-1) = sendbuf(2*i-1)
            recvbuf(2*i+offset) = sendbuf(2*i)
         elseif (j.eq.k) then
            recvbuf(2*i+offset) = max(recvbuf(2*i+offset),sendbuf(2*i))
         endif
   50    continue
c copy initial data
      elseif (op.eq.(-1)) then
         do 60 i = 1, count
         recvbuf(i+offset) = sendbuf(i)
   60    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine fredux(recvbuf,sendbuf,offset,count,op)
c perform reduction operation for real types
c recvbuf = address of receive buffer
c sendbuf = address of send buffer
c offset = starting location minus one in receive buffer
c count = number of elements in send buffer
c op = reduce operation
c (only max, min, sum, maxloc, and minloc currently supported)
c input: recvbuf, sendbuf, offset, count, op
c output: recvbuf
      implicit none
      integer offset, count, op
      real recvbuf(*), sendbuf(*)
c MPI constants
      integer MPI_MAX, MPI_MIN, MPI_SUM, MPI_MAXLOC, MPI_MINLOC
      parameter(MPI_MAX=0,MPI_MIN=1,MPI_SUM=2,MPI_MAXLOC=4,MPI_MINLOC=5)
c local data
      integer i
      real j, k
c perform reduction
      if (op.eq.MPI_MAX) then
         do 10 i = 1, count
         recvbuf(i+offset) = max(recvbuf(i+offset),sendbuf(i))
   10    continue
      elseif (op.eq.MPI_MIN) then
         do 20 i = 1, count
         recvbuf(i+offset) = min(recvbuf(i+offset),sendbuf(i))
   20    continue
      elseif (op.eq.MPI_SUM) then
         do 30 i = 1, count
         recvbuf(i+offset) = recvbuf(i+offset) + sendbuf(i)
   30    continue
c perform reduction and location
      elseif (op.eq.MPI_MAXLOC) then
         do 40 i = 1, count
         j = recvbuf(2*i+offset-1)
         k = sendbuf(2*i-1)
         if (j.lt.k) then
            recvbuf(2*i+offset-1) = sendbuf(2*i-1)
            recvbuf(2*i+offset) = sendbuf(2*i)
         elseif (j.eq.k) then
            recvbuf(2*i+offset) = min(recvbuf(2*i+offset),sendbuf(2*i))
         endif
   40    continue
      elseif (op.eq.MPI_MINLOC) then
         do 50 i = 1, count
         j = recvbuf(2*i+offset-1)
         k = sendbuf(2*i-1)
         if (j.gt.k) then
            recvbuf(2*i+offset-1) = sendbuf(2*i-1)
            recvbuf(2*i+offset) = sendbuf(2*i)
         elseif (j.eq.k) then
            recvbuf(2*i+offset) = max(recvbuf(2*i+offset),sendbuf(2*i))
         endif
   50    continue
c copy initial data
      elseif (op.eq.(-1)) then
         do 60 i = 1, count
         recvbuf(i+offset) = sendbuf(i)
   60    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine dredux(recvbuf,sendbuf,offset,count,op)
c perform reduction operation for double precision types
c recvbuf = address of receive buffer
c sendbuf = address of send buffer
c offset = starting location minus one in receive buffer
c count = number of elements in send buffer
c op = reduce operation
c (only max, min, sum, maxloc, and minloc currently supported)
c input: recvbuf, sendbuf, offset, count, op
c output: recvbuf
      implicit none
      integer offset, count, op
      double precision recvbuf(*), sendbuf(*)
c MPI constants
      integer MPI_MAX, MPI_MIN, MPI_SUM, MPI_MAXLOC, MPI_MINLOC
      parameter(MPI_MAX=0,MPI_MIN=1,MPI_SUM=2,MPI_MAXLOC=4,MPI_MINLOC=5)
c local data
      integer i
      double precision j, k
c perform reduction
      if (op.eq.MPI_MAX) then
         do 10 i = 1, count
         recvbuf(i+offset) = max(recvbuf(i+offset),sendbuf(i))
   10    continue
      elseif (op.eq.MPI_MIN) then
         do 20 i = 1, count
         recvbuf(i+offset) = min(recvbuf(i+offset),sendbuf(i))
   20    continue
      elseif (op.eq.MPI_SUM) then
         do 30 i = 1, count
         recvbuf(i+offset) = recvbuf(i+offset) + sendbuf(i)
   30    continue
c perform reduction and location
      elseif (op.eq.MPI_MAXLOC) then
         do 40 i = 1, count
         j = recvbuf(2*i+offset-1)
         k = sendbuf(2*i-1)
         if (j.lt.k) then
            recvbuf(2*i+offset-1) = sendbuf(2*i-1)
            recvbuf(2*i+offset) = sendbuf(2*i)
         elseif (j.eq.k) then
            recvbuf(2*i+offset) = min(recvbuf(2*i+offset),sendbuf(2*i))
         endif
   40    continue
      elseif (op.eq.MPI_MINLOC) then
         do 50 i = 1, count
         j = recvbuf(2*i+offset-1)
         k = sendbuf(2*i-1)
         if (j.gt.k) then
            recvbuf(2*i+offset-1) = sendbuf(2*i-1)
            recvbuf(2*i+offset) = sendbuf(2*i)
         elseif (j.eq.k) then
            recvbuf(2*i+offset) = max(recvbuf(2*i+offset),sendbuf(2*i))
         endif
   50    continue
c copy initial data
      elseif (op.eq.(-1)) then
         do 60 i = 1, count
         recvbuf(i+offset) = sendbuf(i)
   60    continue
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_ALLREDUCE(sendbuf,recvbuf,count,datatype,op,comm,ie
     1rror)
c applies a reduction operation to the vector sendbuf over the set of
c processes specified by comm and places result in recvbuf on all nodes
c sendbuf = address of send buffer
c recvbuf = address of receive buffer
c count = number of elements in send buffer
c datatype = datatype of elements in send buffer
c op = reduce operation
c (only max, min, sum, maxloc, and minloc currently supported)
c comm = communicator
c ierror = error indicator
c input: sendbuf, count, datatype, op, root, comm
c output: recvbuf, ierror
      implicit none
      integer sendbuf(*), recvbuf(*)
      integer count, datatype, op, comm, ierror
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c local data
      integer root, ierr
      if (monitor.eq.2) write (2,*) 'MPI_ALLREDUCE started'
      root = 0
      call MPI_REDUCE(sendbuf,recvbuf,count,datatype,op,root,comm,ierror
     1)
      call MPI_BCAST(recvbuf,count,datatype,root,comm,ierr)
      if (monitor.eq.2) write (2,*) 'MPI_ALLREDUCE complete'
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_GATHER(sendbuf,sendcount,sendtype,recvbuf,recvcount
     1,recvtype,root,comm,ierror)
c collect individual messages from each process in comm at root
c sendbuf = starting address of send buffer
c sendcount = number of elements in send buffer
c sendtype = datatype of send buffer elements
c recvbuf = address of receive buffer
c recvcount = number of elements for any single receive
c recvtype = datatype of recv buffer elements
c root = rank of receiving process
c comm = communicator
c ierror = error indicator
c input: sendbuf, sendcount, sendtype, recvcount, recvtype, root, comm
c output: recvbuf, ierror
      implicit none
      integer sendbuf(*), recvbuf(*)
      integer sendcount, sendtype, recvcount, recvtype, root, comm
      integer ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c idproc = processor id
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c MPI constants
      integer MPI_STATUS_SIZE
      parameter(MPI_STATUS_SIZE=5)
      integer MPI_INTEGER, MPI_REAL, MPI_DOUBLE_PRECISION
      parameter(MPI_INTEGER=18,MPI_REAL=19,MPI_DOUBLE_PRECISION=20)
      integer MPI_COMPLEX, MPI_DOUBLE_COMPLEX
      parameter(MPI_COMPLEX=22,MPI_DOUBLE_COMPLEX=23)
c local data
      integer np, loct, lsize, id, rank, i, j, status
      dimension status(MPI_STATUS_SIZE)
      ierror = 0
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c invalid root
      elseif ((root.lt.0).or.(root.ge.nproc)) then
         ierror = 19
c invalid count
      elseif (sendcount.lt.0) then
         ierror = 3
c communicator errors
      else
         np = mapcomm(MAXS+1,comm+1)
         id = mapcomm(root+1,comm+1)
c communicator not mapped
         if ((np.le.0).or.(np.gt.nproc)) then
            ierror = 2
c invalid root
         elseif ((root.lt.0).or.(root.ge.np)) then
            ierror = 19
c invalid mapping
         elseif ((id.lt.0).or.(id.ge.nproc)) then
            write (2,*) 'Invalid mapping, root, node = ', root, id
            ierror = 2
c get rank
         else
            rank = mapcomm(MAXS+2,comm+1)
            if ((rank.ge.0).and.(rank.lt.np)) then
               if (mapcomm(rank+1,comm+1).ne.idproc) then
                  ierror = 29
               endif
            else
               ierror = 29
            endif
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_GATHER: ',ierror)
         return
      endif
c root receives data
      if (rank.eq.root) then
c invalid count
         if (recvcount.lt.0) ierror = 3
c determine size of data to be sent
         if ((sendtype.eq.MPI_INTEGER).or.(sendtype.eq.MPI_REAL)) then
            loct = sendcount
         elseif ((sendtype.eq.MPI_DOUBLE_PRECISION).or.(sendtype.eq.MPI_
     1COMPLEX)) then
            loct = 2*sendcount
         elseif (sendtype.eq.MPI_DOUBLE_COMPLEX) then
            loct = 4*sendcount
c invalid datatype
         else
            loct = 0
            ierror = 7
         endif
c determine size of data to be received
         if ((recvtype.eq.MPI_INTEGER).or.(recvtype.eq.MPI_REAL)) then
            lsize = recvcount
         elseif ((recvtype.eq.MPI_DOUBLE_PRECISION).or.(recvtype.eq.MPI_
     1COMPLEX)) then
            lsize = 2*recvcount
         elseif (recvtype.eq.MPI_DOUBLE_COMPLEX) then
            lsize = 4*recvcount
c invalid datatype
         else
            lsize = 0
            ierror = 7
         endif
c unequal message length error
         if (loct.ne.lsize) then
            write (2,*) 'Unequal message length, send/receive bytes = ',
     1loct, lsize
            ierror = 22
         endif
c handle count, datatype and length errors
         if (ierror.ne.0) then
            call writerrs('MPI_GATHER: ',ierror)
            return
         endif
         if (monitor.eq.2) write (2,*) 'MPI_GATHER started'
         do 20 i = 1, np
         id = i - 1
         loct = lsize*id
c root copies its own data directly
         if (id.eq.root) then
            do 10 j = 1, lsize
            recvbuf(j+loct) = sendbuf(j)
   10       continue
c otherwise, root receives data from other processors
         else
            call MPI_RECV(recvbuf(loct+1),recvcount,recvtype,id,1,comm,s
     1tatus,ierror)
         endif
   20    continue
c processors other than root send data to root
      else
         if (monitor.eq.2) write (2,*) 'MPI_GATHER started'
         call MPI_SEND(sendbuf,sendcount,sendtype,root,1,comm,ierror)
      endif
      if (monitor.eq.2) write (2,*) 'MPI_GATHER complete'
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_ALLGATHER(sendbuf,sendcount,sendtype,recvbuf,recvco
     1unt,recvtype,comm,ierror)
c gather individual messages from each process in comm and distribute
c the resulting message to each process.
c sendbuf = starting address of send buffer
c sendcount = number of elements in send buffer
c sendtype = datatype of send buffer elements
c recvbuf = address of receive buffer
c recvcount = number of elements for any process
c recvtype = datatype of receive buffer elements
c comm = communicator
c ierror = error indicator
c input: sendbuf, sendcount, sendtype, recvcount, recvtype, comm
c output: recvbuf, ierror
      implicit none
      integer sendbuf(*), recvbuf(*)
      integer sendcount, sendtype, recvcount, recvtype, comm, ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c local data
      integer np, root, ierr
      if (monitor.eq.2) write (2,*) 'MPI_ALLGATHER started'
      root = 0
      call MPI_GATHER(sendbuf,sendcount,sendtype,recvbuf,recvcount,recvt
     1ype,root,comm,ierror)
      np = mapcomm(MAXS+1,comm+1)
      call MPI_BCAST(recvbuf,np*recvcount,recvtype,root,comm,ierr)
      if (monitor.eq.2) write (2,*) 'MPI_ALLGATHER complete'
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_SCATTER(sendbuf,sendcount,sendtype,recvbuf,recvcoun
     1t,recvtype,root,comm,ierror)
c distribute individual messages from root to each process in comm
c sendbuf = starting address of send buffer
c sendcount = number of elements in send buffer
c sendtype = datatype of send buffer elements
c recvbuf = address of receive buffer
c recvcount = number of elements for any single receive
c recvtype = datatype of recv buffer elements
c root = rank of sending process
c comm = communicator
c ierror = error indicator
c input: sendbuf, sendcount, sendtype, recvcount, recvtype, root, comm
c output: recvbuf, ierror
      implicit none
      integer sendbuf(*), recvbuf(*)
      integer sendcount, sendtype, recvcount, recvtype, root, comm
      integer ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c idproc = processor id
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c MPI constants
      integer MPI_STATUS_SIZE
      parameter(MPI_STATUS_SIZE=5)
      integer MPI_INTEGER, MPI_REAL, MPI_DOUBLE_PRECISION
      parameter(MPI_INTEGER=18,MPI_REAL=19,MPI_DOUBLE_PRECISION=20)
      integer MPI_COMPLEX, MPI_DOUBLE_COMPLEX
      parameter(MPI_COMPLEX=22,MPI_DOUBLE_COMPLEX=23)
c local data
      integer np, lsize, loct, id, rank, i, j, status
      dimension status(MPI_STATUS_SIZE)
      ierror = 0
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c invalid root
      elseif ((root.lt.0).or.(root.ge.nproc)) then
         ierror = 19
c invalid counts
      elseif (recvcount.lt.0) then
         ierror = 3
c communicator errors
      else
         np = mapcomm(MAXS+1,comm+1)
         id = mapcomm(root+1,comm+1)
c communicator not mapped
         if ((np.le.0).or.(np.gt.nproc)) then
            ierror = 2
c invalid root
         elseif ((root.lt.0).or.(root.ge.np)) then
            ierror = 19
c invalid mapping
         elseif ((id.lt.0).or.(id.ge.nproc)) then
            write (2,*) 'Invalid mapping, root, node = ', root, id
            ierror = 2
c get rank
         else
            rank = mapcomm(MAXS+2,comm+1)
            if ((rank.ge.0).and.(rank.lt.np)) then
               if (mapcomm(rank+1,comm+1).ne.idproc) then
                  ierror = 29
               endif
            else
               ierror = 29
            endif
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_SCATTER: ',ierror)
         return
      endif
c root sends data
      if (rank.eq.root) then
c invalid counts
         if (sendcount.lt.0) ierror = 3
c determine size of data to be sent
         if ((sendtype.eq.MPI_INTEGER).or.(sendtype.eq.MPI_REAL)) then
            lsize = sendcount
         elseif ((sendtype.eq.MPI_DOUBLE_PRECISION).or.(sendtype.eq.MPI_
     1COMPLEX)) then
            lsize = 2*sendcount
         elseif (sendtype.eq.MPI_DOUBLE_COMPLEX) then
            lsize = 4*sendcount
c invalid datatype
         else
            lsize = 0
            ierror = 7
         endif
c determine size of data to be received
         if ((recvtype.eq.MPI_INTEGER).or.(recvtype.eq.MPI_REAL)) then
            loct = recvcount
         elseif ((recvtype.eq.MPI_DOUBLE_PRECISION).or.(recvtype.eq.MPI_
     1COMPLEX)) then
            loct = 2*recvcount
         elseif (recvtype.eq.MPI_DOUBLE_COMPLEX) then
            loct = 4*recvcount
c invalid datatype
         else
            loct = 0
            ierror = 7
         endif
c unequal message length error
         if (loct.ne.lsize) then
            write (2,*) 'Unequal message length, send/receive bytes = ',
     1lsize, loct
            ierror = 22
         endif
c handle count, datatype and length errors
         if (ierror.ne.0) then
            call writerrs('MPI_SCATTER: ',ierror)
            return
         endif
         if (monitor.eq.2) write (2,*) 'MPI_SCATTER started'
         do 20 i = 1, np
         id = i - 1
         loct = lsize*id
c root copies its own data directly
         if (id.eq.root) then
            do 10 j = 1, lsize
            recvbuf(j) = sendbuf(j+loct)
   10       continue
c otherwise, root sends data to other processors
         else
            call MPI_SEND(sendbuf(loct+1),sendcount,sendtype,id,1,comm,i
     1error)
         endif
   20    continue
c processors other than root receive data from root
      else
         if (monitor.eq.2) write (2,*) 'MPI_SCATTER started'
         call MPI_RECV(recvbuf,recvcount,recvtype,root,1,comm,status,ier
     1ror)
      endif
      if (monitor.eq.2) write (2,*) 'MPI_SCATTER complete'
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_ALLTOALL(sendbuf,sendcount,sendtype,recvbuf,recvcou
     1nt,recvtype,comm,ierror)
c send a distinct message from each process to every process
c sendbuf = starting address of send buffer
c sendcount = number of elements in send buffer
c sendtype = datatype of send buffer elements
c recvbuf = address of receive buffer
c recvcount = number of elements for any single receive
c recvtype = datatype of recv buffer elements
c comm = communicator
c ierror = error indicator
c input: sendbuf, sendcount, sendtype, recvcount, recvtype, comm
c output: recvbuf, ierror
      implicit none
      integer sendbuf(*), recvbuf(*)
      integer sendcount, sendtype, recvcount, recvtype, comm
      integer ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c idproc = processor id
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c MPI constants
      integer MPI_STATUS_SIZE
      parameter(MPI_STATUS_SIZE=5)
      integer MPI_INTEGER, MPI_REAL, MPI_DOUBLE_PRECISION
      parameter(MPI_INTEGER=18,MPI_REAL=19,MPI_DOUBLE_PRECISION=20)
      integer MPI_COMPLEX, MPI_DOUBLE_COMPLEX
      parameter(MPI_COMPLEX=22,MPI_DOUBLE_COMPLEX=23)
c local data
      integer np, loct, lsize, id, rank, i, j, request, status
      dimension status(MPI_STATUS_SIZE)
      ierror = 0
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c invalid counts
      elseif ((sendcount.lt.0).or.(recvcount.lt.0)) then
         ierror = 3
c communicator errors
      else
         np = mapcomm(MAXS+1,comm+1)
c communicator not mapped
         if ((np.le.0).or.(np.gt.nproc)) then
            ierror = 2
c get rank
         else
            rank = mapcomm(MAXS+2,comm+1)
            if ((rank.ge.0).and.(rank.lt.np)) then
               if (mapcomm(rank+1,comm+1).ne.idproc) then
                  ierror = 29
               endif
            else
               ierror = 29
            endif
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_ALLTOALL: ',ierror)
         return
      endif
c determine size of data to be sent
      if ((sendtype.eq.MPI_INTEGER).or.(sendtype.eq.MPI_REAL)) then
         loct = sendcount
      elseif ((sendtype.eq.MPI_DOUBLE_PRECISION).or.(sendtype.eq.MPI_
     1COMPLEX)) then
         loct = 2*sendcount
      elseif (sendtype.eq.MPI_DOUBLE_COMPLEX) then
         loct = 4*sendcount
c invalid datatype
      else
         loct = 0
         ierror = 7
      endif
c determine size of data to be received
      if ((recvtype.eq.MPI_INTEGER).or.(recvtype.eq.MPI_REAL)) then
         lsize = recvcount
      elseif ((recvtype.eq.MPI_DOUBLE_PRECISION).or.(recvtype.eq.MPI_
     1COMPLEX)) then
         lsize = 2*recvcount
      elseif (recvtype.eq.MPI_DOUBLE_COMPLEX) then
         lsize = 4*recvcount
c invalid datatype
      else
         lsize = 0
         ierror = 7
      endif
c unequal message length error
      if (loct.ne.lsize) then
         write (2,*) 'Unequal message length, send/receive bytes = ',
     1loct, lsize
         ierror = 22
      endif
c handle count, datatype and length errors
      if (ierror.ne.0) then
         call writerrs('MPI_ALLTOALL: ',ierror)
         return
      endif
      if (monitor.eq.2) write (2,*) 'MPI_ALLTOALL started'
      do 20 i = 1, np
      id = i - rank - 1
      if (id.lt.0) id = id + np
      loct = lsize*id
c each node copies its own data directly
      if (rank.eq.id) then
         do 10 j = 1, lsize
         recvbuf(j+loct) = sendbuf(j+loct)
   10    continue
c otherwise, each node receives data from other nodes
      else
         call MPI_IRECV(recvbuf(loct+1),recvcount,recvtype,id,i,comm,req
     1uest,ierror)
         call MPI_SEND(sendbuf(loct+1),sendcount,sendtype,id,i,comm,ierr
     1or)
         call MPI_WAIT(request,status,ierror)
      endif
   20 continue
      if (monitor.eq.2) write (2,*) 'MPI_ALLTOALL complete'
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_GATHERV(sendbuf,sendcount,sendtype,recvbuf,recvcoun
     1ts,displs,recvtype,root,comm,ierror)
c collect individual messages from each process in comm at root
c messages can have different sizes and displacements
c sendbuf = starting address of send buffer
c sendcount = number of elements in send buffer
c sendtype = datatype of send buffer elements
c recvbuf = address of receive buffer
c recvcounts = integer array
c displs = integer array of displacements
c recvtype = datatype of recv buffer elements
c root = rank of receiving process
c comm = communicator
c ierror = error indicator
c input: sendbuf, sendcount, sendtype, recvcounts, displs, recvtype
c input: root, comm
c output: recvbuf, ierror
      implicit none
      integer sendbuf(*), recvbuf(*), recvcounts(*), displs(*)
      integer sendcount, sendtype, recvtype, root, comm
      integer ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c idproc = processor id
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c MPI constants
      integer MPI_STATUS_SIZE
      parameter(MPI_STATUS_SIZE=5)
      integer MPI_INTEGER, MPI_REAL, MPI_DOUBLE_PRECISION
      parameter(MPI_INTEGER=18,MPI_REAL=19,MPI_DOUBLE_PRECISION=20)
      integer MPI_COMPLEX, MPI_DOUBLE_COMPLEX
      parameter(MPI_COMPLEX=22,MPI_DOUBLE_COMPLEX=23)
c local data
      integer np, loct, lsize, id, rank, i, j, status
      dimension status(MPI_STATUS_SIZE)
      ierror = 0
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c invalid root
      elseif ((root.lt.0).or.(root.ge.nproc)) then
         ierror = 19
c invalid count
      elseif (sendcount.lt.0) then
         ierror = 3
c communicator errors
      else
         np = mapcomm(MAXS+1,comm+1)
         id = mapcomm(root+1,comm+1)
c communicator not mapped
         if ((np.le.0).or.(np.gt.nproc)) then
            ierror = 2
c invalid root
         elseif ((root.lt.0).or.(root.ge.np)) then
            ierror = 19
c invalid mapping
         elseif ((id.lt.0).or.(id.ge.nproc)) then
            write (2,*) 'Invalid mapping, root, node = ', root, id
            ierror = 2
c get rank
         else
            rank = mapcomm(MAXS+2,comm+1)
            if ((rank.ge.0).and.(rank.lt.np)) then
               if (mapcomm(rank+1,comm+1).ne.idproc) then
                  ierror = 29
               endif
            else
               ierror = 29
            endif
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_GATHERV: ',ierror)
         return
      endif
c root receives data
      if (rank.eq.root) then
c invalid counts
         do 10 i = 1, np
         if (recvcounts(i).lt.0) ierror = 3
   10    continue
c determine size of data to be sent
         if ((sendtype.eq.MPI_INTEGER).or.(sendtype.eq.MPI_REAL)) then
            loct = 1
         elseif ((sendtype.eq.MPI_DOUBLE_PRECISION).or.(sendtype.eq.MPI_
     1COMPLEX)) then
            loct = 2
         elseif (sendtype.eq.MPI_DOUBLE_COMPLEX) then
            loct = 4
c invalid datatype
         else
            ierror = 7
         endif
c determine size of data to be received
         if ((recvtype.eq.MPI_INTEGER).or.(recvtype.eq.MPI_REAL)) then
            lsize = 1
         elseif ((recvtype.eq.MPI_DOUBLE_PRECISION).or.(recvtype.eq.MPI_
     1COMPLEX)) then
            lsize = 2
         elseif (recvtype.eq.MPI_DOUBLE_COMPLEX) then
            lsize = 4
c invalid datatype
         else
            ierror = 7
         endif
c unequal message length error
         id = lsize*recvcounts(rank+1)
         if ((ierror.eq.0).and.(loct*sendcount.ne.id)) then
            write (2,*) 'Unequal self message, send/receive bytes = ',
     1loct*sendcount, id
            ierror = 22
         endif
c handle count and datatype errors
         if (ierror.ne.0) then
            call writerrs('MPI_GATHERV: ',ierror)
            return
         endif
         if (monitor.eq.2) write (2,*) 'MPI_GATHERV started'
         do 30 i = 1, np
         id = i - 1
         loct = lsize*displs(i)
c root copies its own data directly
         if (id.eq.root) then
            do 20 j = 1, lsize*recvcounts(i)
            recvbuf(j+loct) = sendbuf(j)
   20       continue
c otherwise, root receives data from other processors
         else
            call MPI_RECV(recvbuf(loct+1),recvcounts(i),recvtype,id,1,co
     1mm,status,ierror)
         endif
   30    continue
c processors other than root send data to root
      else
         if (monitor.eq.2) write (2,*) 'MPI_GATHERV started'
         call MPI_SEND(sendbuf,sendcount,sendtype,root,1,comm,ierror)
      endif
      if (monitor.eq.2) write (2,*) 'MPI_GATHERV complete'
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_ALLGATHERV(sendbuf,sendcount,sendtype,recvbuf,recvc
     1ounts,displs,recvtype,comm,ierror)
c gather individual messages from each process in comm and distribute
c the resulting message to each process.
c messages can have different sizes and displacements
c sendbuf = starting address of send buffer
c sendcount = number of elements in send buffer
c sendtype = datatype of send buffer elements
c recvbuf = address of receive buffer
c recvcounts = integer array
c displs = integer array of displacements
c recvtype = datatype of receive buffer elements
c comm = communicator
c ierror = error indicator
c input: sendbuf, sendcount, sendtype, recvcounts, displs, recvtype
c input: comm
c output: recvbuf, ierror
      implicit none
      integer sendbuf(*), recvbuf(*), recvcounts(*), displs(*)
      integer sendcount, sendtype, recvtype, comm, ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c local data
      integer np, root, i, ierr
      ierror = 0
c communicator errors
      if ((comm.ge.0).and.(comm.lt.MAXC)) then
         np = mapcomm(MAXS+1,comm+1)
         if ((np.le.0).or.(np.gt.nproc)) then
            ierror = 2
         endif
      else
         ierror = 2
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_ALLGATHERV: ',ierror)
         return
      endif
      if (monitor.eq.2) write (2,*) 'MPI_ALLGATHERV started'
      do 10 i = 1, np
      root = i - 1
      call MPI_GATHERV(sendbuf,sendcount,sendtype,recvbuf,recvcounts,dis
     1pls,recvtype,root,comm,ierr)
      if (ierr.ne.0) ierror = ierr
   10 continue
      if (monitor.eq.2) write (2,*) 'MPI_ALLGATHERV complete'
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_SCATTERV(sendbuf,sendcounts,displs,sendtype,recvbuf
     1,recvcount,recvtype,root,comm,ierror)
c distribute individual messages from root to each process in comm
c messages can have different sizes and displacements
c sendbuf = starting address of send buffer
c sendcounts = integer array
c displs = integer array of displacements
c sendtype = datatype of send buffer elements
c recvbuf = address of receive buffer
c recvcount = number of elements for any single receive
c recvtype = datatype of recv buffer elements
c root = rank of sending process
c comm = communicator
c ierror = error indicator
c input: sendbuf, sendcounts, displs, sendtype, recvcount, recvtype
c input: root, comm
c output: recvbuf, ierror
      implicit none
      integer sendbuf(*), recvbuf(*), sendcounts(*), displs(*)
      integer sendtype, recvcount, recvtype, root, comm
      integer ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c idproc = processor id
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c MPI constants
      integer MPI_STATUS_SIZE
      parameter(MPI_STATUS_SIZE=5)
      integer MPI_INTEGER, MPI_REAL, MPI_DOUBLE_PRECISION
      parameter(MPI_INTEGER=18,MPI_REAL=19,MPI_DOUBLE_PRECISION=20)
      integer MPI_COMPLEX, MPI_DOUBLE_COMPLEX
      parameter(MPI_COMPLEX=22,MPI_DOUBLE_COMPLEX=23)
c local data
      integer np, lsize, loct, id, rank, i, j, status
      dimension status(MPI_STATUS_SIZE)
      ierror = 0
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c invalid root
      elseif ((root.lt.0).or.(root.ge.nproc)) then
         ierror = 19
c invalid counts
      elseif (recvcount.lt.0) then
         ierror = 3
c communicator errors
      else
         np = mapcomm(MAXS+1,comm+1)
         id = mapcomm(root+1,comm+1)
c communicator not mapped
         if ((np.le.0).or.(np.gt.nproc)) then
            ierror = 2
c invalid root
         elseif ((root.lt.0).or.(root.ge.np)) then
            ierror = 19
c invalid mapping
         elseif ((id.lt.0).or.(id.ge.nproc)) then
            write (2,*) 'Invalid mapping, root, node = ', root, id
            ierror = 2
c get rank
         else
            rank = mapcomm(MAXS+2,comm+1)
            if ((rank.ge.0).and.(rank.lt.np)) then
               if (mapcomm(rank+1,comm+1).ne.idproc) then
                  ierror = 29
               endif
            else
               ierror = 29
            endif
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_SCATTERV: ',ierror)
         return
      endif
c root sends data
      if (rank.eq.root) then
c invalid counts
         do 10 i = 1, np
         if (sendcounts(i).lt.0) ierror = 3
   10    continue
c determine size of data to be sent
         if ((sendtype.eq.MPI_INTEGER).or.(sendtype.eq.MPI_REAL)) then
            lsize = 1
      elseif ((sendtype.eq.MPI_DOUBLE_PRECISION).or.(sendtype.eq.MPI_
     1COMPLEX)) then
            lsize = 2
         elseif (sendtype.eq.MPI_DOUBLE_COMPLEX) then
            lsize = 4
c invalid datatype
         else
            ierror = 7
         endif
c determine size of data to be received
         if ((recvtype.eq.MPI_INTEGER).or.(recvtype.eq.MPI_REAL)) then
            loct = 1
         elseif ((recvtype.eq.MPI_DOUBLE_PRECISION).or.(recvtype.eq.MPI_
     1COMPLEX)) then
            loct = 2
         elseif (recvtype.eq.MPI_DOUBLE_COMPLEX) then
            loct = 4
c invalid datatype
         else
            ierror = 7
         endif
c unequal message length error
         id = lsize*sendcounts(rank+1)
         if ((ierror.eq.0).and.(loct*recvcount.ne.id)) then
            write (2,*) 'Unequal self message, send/receive bytes = ',
     1id, loct*recvcount
            ierror = 22
         endif
c handle count and datatype errors
         if (ierror.ne.0) then
            call writerrs('MPI_SCATTERV: ',ierror)
            return
         endif
         if (monitor.eq.2) write (2,*) 'MPI_SCATTERV started'
         do 30 i = 1, np
         id = i - 1
         loct = lsize*displs(i)
c root copies its own data directly
         if (id.eq.root) then
            do 20 j = 1, lsize*sendcounts(i)
            recvbuf(j) = sendbuf(j+loct)
   20       continue
c otherwise, root sends data to other processors
         else
            call MPI_SEND(sendbuf(loct+1),sendcounts(i),sendtype,id,1,co
     1mm,ierror)
         endif
   30    continue
c processors other than root receive data from root
      else
         if (monitor.eq.2) write (2,*) 'MPI_SCATTERV started'
         call MPI_RECV(recvbuf,recvcount,recvtype,root,1,comm,status,ier
     1ror)
      endif
      if (monitor.eq.2) write (2,*) 'MPI_SCATTERV complete'
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_ALLTOALLV(sendbuf,sendcounts,sdispls,sendtype,recvb
     1uf,recvcounts,rdispls,recvtype,comm,ierror)
c send a distinct message from each process to every process
c messages can have different sizes and displacements
c sendbuf = starting address of send buffer
c sendcounts = integer array
c sdispls = integer array of send displacements
c sendtype = datatype of send buffer elements
c recvbuf = address of receive buffer
c recvcounts = integer array
c rdispls = integer array of receive displacements
c recvtype = datatype of recv buffer elements
c comm = communicator
c ierror = error indicator
c input: sendbuf, sendcount, sendtype, recvcount, recvtype, comm
c output: recvbuf, ierror
      implicit none
      integer sendbuf(*), recvbuf(*)
      integer sendcounts(*), sdispls(*), recvcounts(*), rdispls(*)
      integer sendtype, recvtype, comm
      integer ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c idproc = processor id
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c MPI constants
      integer MPI_STATUS_SIZE
      parameter(MPI_STATUS_SIZE=5)
      integer MPI_INTEGER, MPI_REAL, MPI_DOUBLE_PRECISION
      parameter(MPI_INTEGER=18,MPI_REAL=19,MPI_DOUBLE_PRECISION=20)
      integer MPI_COMPLEX, MPI_DOUBLE_COMPLEX
      parameter(MPI_COMPLEX=22,MPI_DOUBLE_COMPLEX=23)
c local data
      integer status
      integer np, locs, msize, loct, lsize, id, ld, rank, i, j, request
      dimension status(MPI_STATUS_SIZE)
      ierror = 0
c check for error conditions
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c invalid comm
      elseif ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c communicator errors
      else
         np = mapcomm(MAXS+1,comm+1)
c communicator not mapped
         if ((np.le.0).or.(np.gt.nproc)) then
            ierror = 2
c get rank
         else
            rank = mapcomm(MAXS+2,comm+1)
            if ((rank.ge.0).and.(rank.lt.np)) then
               if (mapcomm(rank+1,comm+1).ne.idproc) then
                  ierror = 29
               endif
            else
               ierror = 29
            endif
         endif
      endif
c invalid counts
      do 10 i = 1, np
      if ((sendcounts(i).lt.0).or.(recvcounts(i).lt.0)) ierror = 3
   10 continue
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_ALLTOALLV: ',ierror)
         return
      endif
c determine size of data to be sent
      if ((sendtype.eq.MPI_INTEGER).or.(sendtype.eq.MPI_REAL)) then
         msize = 1
      elseif ((sendtype.eq.MPI_DOUBLE_PRECISION).or.(sendtype.eq.MPI_
     1COMPLEX)) then
         msize = 2
      elseif (sendtype.eq.MPI_DOUBLE_COMPLEX) then
         msize = 4
c invalid datatype
      else
         ierror = 7
      endif
c determine size of data to be received
      if ((recvtype.eq.MPI_INTEGER).or.(recvtype.eq.MPI_REAL)) then
         lsize = 1
      elseif ((recvtype.eq.MPI_DOUBLE_PRECISION).or.(recvtype.eq.MPI_
     1COMPLEX)) then
         lsize = 2
      elseif (recvtype.eq.MPI_DOUBLE_COMPLEX) then
         lsize = 4
c invalid datatype
      else
         ierror = 7
      endif
c unequal message length error
      id = msize*sendcounts(rank+1)
      ld = lsize*recvcounts(rank+1)
      if ((ierror.eq.0).and.(id.ne.ld)) then
         write (2,*) 'Unequal self message length, send/receive bytes=',
     1id, ld
         ierror = 22
      endif
c handle count and datatype errors
      if (ierror.ne.0) then
         call writerrs('MPI_ALLTOALLV: ',ierror)
         return
      endif
      if (monitor.eq.2) write (2,*) 'MPI_ALLTOALLV started'
      do 30 i = 1, np
      id = i - rank - 1
      if (id.lt.0) id = id + np
      ld = id + 1
      locs = msize*sdispls(ld)
      loct = lsize*rdispls(ld)
c each node copies its own data directly
      if (rank.eq.id) then
         do 20 j = 1, lsize*recvcounts(ld)
         recvbuf(j+loct) = sendbuf(j+locs)
   20    continue
c otherwise, each node receives data from other nodes
      else
         call MPI_IRECV(recvbuf(loct+1),recvcounts(ld),recvtype,id,i,com
     1m,request,ierror)
         call MPI_SEND(sendbuf(locs+1),sendcounts(ld),sendtype,id,i,comm
     1,ierror)
         call MPI_WAIT(request,status,ierror)
      endif
   30 continue
      if (monitor.eq.2) write (2,*) 'MPI_ALLTOALLV complete'
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_REDUCE_SCATTER(sendbuf,recvbuf,recvcounts,datatype,
     1op,comm,ierror)
c applies a reduction operation to the vector sendbuf over the set of
c processes specified by comm and scatters the result according to the
c values in recvcounts
c sendbuf = starting address of send buffer
c recvbuf = starting address of receive buffer
c recvcounts = integer array
c datatype = datatype of elements in input buffer
c op = reduce operation
c (only max, min, sum, maxloc, and minloc currently supported)
c comm = communicator
c ierror = error indicator
c input: sendbuf, recvcounts, datatype, op, comm
c output: recvbuf, ierror
      implicit none
      integer sendbuf(*), recvbuf(*), recvcounts(*)
      integer datatype, op, comm, ierror
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c idproc = processor id
c mapcomm = communicator map
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c function declarations
      integer NewPtr
      external NewPtr
c MPI constants
      integer MPI_SUM, MPI_MAXLOC, MPI_MINLOC
      parameter(MPI_SUM=2,MPI_MAXLOC=4,MPI_MINLOC=5)
      integer MPI_INTEGER, MPI_REAL, MPI_DOUBLE_PRECISION
      parameter(MPI_INTEGER=18,MPI_REAL=19,MPI_DOUBLE_PRECISION=20)
      integer MPI_COMPLEX, MPI_DOUBLE_COMPLEX
      parameter(MPI_COMPLEX=22,MPI_DOUBLE_COMPLEX=23)
      integer MPI_2REAL, MPI_2DOUBLE_PRECISION, MPI_2INTEGER
      parameter(MPI_2REAL=35,MPI_2DOUBLE_PRECISION=36,MPI_2INTEGER=37)
c local data
      integer np, rank, root, count, lsize, ltmp, i, displs, tmpbuf
      data root /0/
      ierror = 0
c invalid comm
      if ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
c communicator errors
      else
         np = mapcomm(MAXS+1,comm+1)
c communicator not mapped
         if ((np.le.0).or.(np.gt.nproc)) then
            ierror = 2
c get rank
         else
            rank = mapcomm(MAXS+2,comm+1)
            if ((rank.ge.0).and.(rank.lt.np)) then
               if (mapcomm(rank+1,comm+1).ne.idproc) then
                  ierror = 29
               endif
            else
               ierror = 29
            endif
         endif
      endif
c handle errors
      if (ierror.ne.0) then
         call writerrs('MPI_REDUCE_SCATTER: ',ierror)
         return
      endif
      ltmp = 4*np
c allocate a nonrelocatable block of memory
      displs = NewPtr(val4(ltmp))
c memory not available
      if (displs.eq.0) then
         ierror = 21
         call writerrs('MPI_REDUCE_SCATTER: ',ierror)
         return
      endif
      count = 0
c determines overall count and displacements
      do 10 i = 1, np
      long(displs+4*(i-1)) = count
      count = count + recvcounts(i)
   10 continue
c determine size of temporary buffer
      if ((datatype.eq.MPI_INTEGER).or.(datatype.eq.MPI_REAL)) then
         lsize = 4
      elseif (datatype.eq.MPI_DOUBLE_PRECISION) then
         lsize = 8
      elseif ((datatype.eq.MPI_COMPLEX).and.(op.eq.MPI_SUM)) then
         lsize = 8
      elseif ((datatype.eq.MPI_DOUBLE_COMPLEX).and.(op.eq.MPI_SUM)) then
         lsize = 16
      elseif ((datatype.eq.MPI_2INTEGER).and.((op.eq.MPI_MAXLOC).or.(op.
     1eq.MPI_MINLOC))) then
         lsize = 8
      elseif ((datatype.eq.MPI_2REAL).and.((op.eq.MPI_MAXLOC).or.(op.eq.
     1MPI_MINLOC))) then
         lsize = 8
      elseif ((datatype.eq.MPI_2DOUBLE_PRECISION).and.((op.eq.MPI_MAXLOC
     1).or.(op.eq.MPI_MINLOC))) then
         lsize = 16
      else
         ierror = 7
         call writerrs('MPI_REDUCE_SCATTER: ',ierror)
         return
      endif
      if (monitor.eq.2) write (2,*) 'MPI_REDUCE_SCATTER started'
      ltmp = lsize*count
c allocate a nonrelocatable block of memory
      tmpbuf = NewPtr(val4(ltmp))
c memory not available
      if (tmpbuf.eq.0) then
         ierror = 21
         call writerrs('MPI_REDUCE_SCATTER: ',ierror)
         return
      endif
      call MPI_REDUCE(sendbuf,val4(tmpbuf),count,datatype,op,root,comm,i
     1error)
      call MPI_SCATTERV(val4(tmpbuf),recvcounts,val4(displs),datatype,re
     1cvbuf,recvcounts(rank+1),datatype,root,comm,ierror)
c release nonrelocatable memory block
      call DisposePtr(val4(tmpbuf))
c release nonrelocatable memory block
      call DisposePtr(val4(displs))
      if (monitor.eq.2) write (2,*) 'MPI_REDUCE_SCATTER complete'
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_ABORT(comm,errorcode,ierror)
c force all tasks on an MPI environment to terminate
c comm = communicator
c errorcode = error code to return to invoking environment
c ierror = error indicator
c input: comm, errorcode
c output: ierror
      implicit none
      integer comm, errorcode, ierror
c declare internal mpi common block
      integer nproc, idproc
c nproc = number of real or virtual processors obtained
      common /mpiparms/ nproc, idproc
c MPI not initialized
      if (nproc.le.0) then
         ierror = 1
c this is just a temporary patch, have not yet notified everyone else
      else
         call MPI_FINALIZE(ierror)
      endif
      return
      end
c-----------------------------------------------------------------------
      function MPI_WTIME()
c return an elapsed time on the calling processor in seconds
      implicit none
      double precision MPI_WTIME, tick
      parameter(tick=1.0d0/1000.0d0)
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c stime = first time stamp if MPI_Init successful
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c function declarations
      integer OTElapsedMilliseconds
      external OTElapsedMilliseconds
c calculate time elapsed in milliseconds
      MPI_WTIME = dble(OTElapsedMilliseconds(stime))*tick
      return
      end
c-----------------------------------------------------------------------
      function MPI_WTICK()
c return the resolution of MPI_WTIME in seconds
      implicit none
      double precision MPI_WTICK, tick
      parameter(tick=1.0d0/1000.0d0)
      MPI_WTICK = tick
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_TYPE_EXTENT(datatype,extent,ierror)
c returns the size of a datatype
c datatype = datatype
c extent = datatype extent
c ierror = error indicator
c input: dataype
c output: extent, ierror
      implicit none
      integer datatype, extent, ierror
c MPI constants
      integer MPI_INTEGER, MPI_REAL, MPI_DOUBLE_PRECISION
      parameter(MPI_INTEGER=18,MPI_REAL=19,MPI_DOUBLE_PRECISION=20)
      integer MPI_COMPLEX, MPI_DOUBLE_COMPLEX, MPI_BYTE
      parameter(MPI_COMPLEX=22,MPI_DOUBLE_COMPLEX=23,MPI_BYTE=2)
      integer MPI_2REAL, MPI_2DOUBLE_PRECISION, MPI_2INTEGER
      parameter(MPI_2REAL=35,MPI_2DOUBLE_PRECISION=36,MPI_2INTEGER=37)
      ierror = 0
c find size of datatype
      if ((datatype.eq.MPI_INTEGER).or.(datatype.eq.MPI_REAL)) then
         extent = 4
      elseif ((datatype.eq.MPI_DOUBLE_PRECISION).or.(datatype.eq.MPI_COM
     1PLEX)) then
         extent = 8
      elseif (datatype.eq.MPI_DOUBLE_COMPLEX) then
         extent = 16
      elseif (datatype.eq.MPI_BYTE) then
         extent = 1
      elseif ((datatype.eq.MPI_2INTEGER).or.(datatype.eq.MPI_2REAL)) the
     1n
         extent = 8
      elseif (datatype.eq.MPI_2DOUBLE_PRECISION) then
         extent = 16
c invalid datatype
      else
         ierror = 7
         write (2,*) 'MPI_TYPE_EXTENT: Invalid datatype'
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_ERRHANDLER_SET(comm,errhandler,ierror)
c associates error handler with communicator comm
c only MPI_COMM_WORLD and built-in error handlers are supported
c comm = communicator
c errhandler = new MPI error handler for communicator
c ierror = error indicator
c input: comm, errhandler
c output: ierror
      implicit none
      integer comm, errhandler, ierror
c declare error handler common block
      integer MAXC
      parameter(MAXC=10)
      integer errh
      dimension errh(MAXC)
c errh = error handler
      common /mpierrh/ errh
c MPI constants
      integer MPI_ERRORS_RETURN,MPI_ERRORS_ARE_FATAL
      parameter(MPI_ERRORS_RETURN=0,MPI_ERRORS_ARE_FATAL=1)
      ierror = 0
c invalid comm
      if ((comm.lt.0).or.(comm.ge.MAXC)) then
         ierror = 2
         write (2,*) 'MPI_ERRHANDLER_SET: Invalid communicator'
         return
      endif
      if ((errhandler.eq.MPI_ERRORS_RETURN).or.(errhandler.eq.MPI_ERRORS
     1_ARE_FATAL)) then
         errh(comm+1) = errhandler
c set MPI_COMM_WORLD error handler as well
         errh(1) = errhandler
      else
         ierror = 33
         write (2,*) 'MPI_ERRHANDLER_SET: Invalid errhandler'
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MPI_GET_PROCESSOR_NAME(name,resultlen,ierror)
c returns the name of the processor on which it was called
c name = a unique specifier for the current physical node
c resultlen = length of the result returned in name
c ierror = error indicator
c input: none
c output: name, resultlen, ierror
      implicit none
      character*(*) name
      integer resultlen, ierror
c declare internal mpi common block
      integer nproc, idproc
c idproc = processor id
      common /mpiparms/ nproc, idproc
c function declarations
      integer OTInetGetInterfaceInfo
      integer lencstr, lenstr, instr
      external OTInetGetInterfaceInfo
      external lencstr, lenstr, instr
c MPI constants
      integer kDefaultInetInterface
      parameter(kDefaultInetInterface=-1)
c local data
      integer*2 intw(2), info(276)
      integer longw, iv, nv, lv, oss
      character*16 myself, ename
c longw and intw are used to convert types between integer*2 and integer
      equivalence(longw,intw)
      ierror = 0
c convert processor id to printable integer
      write (ename,'(i5)') idproc
      iv = instr(ename)
      nv = lenstr(ename(iv:))
c Obtain information about the Internet environment
      oss = OTInetGetInterfaceInfo(info,val4(kDefaultInetInterface))
      if (oss.ne.0) then
         ierror = oss
         write (2,*) 'MPI_GET_PROCESSOR_NAME: GetInterfaceInfo Error'
         return
      endif
      intw(1) = info(1)
      intw(2) = info(2)
      call OTInetHostToString(val4(longw),myself)
      lv = lencstr(myself)
      name = myself(1:lv)//':'//ename(iv:iv+nv-1)
      resultlen = min(lv+nv+1,len(name))
      return
      end
c-----------------------------------------------------------------------
      subroutine writerrs(source,ierror)
c this subroutine writes out error descriptions from error codes
c source = source subroutine of error message
c ierror = error indicator
c input: source, ierror
      implicit none
      character*(*) source
      integer ierror
c declare error handler common block
      integer MAXC
      parameter(MAXC=10)
      integer errh
      dimension errh(MAXC)
c errh = error handler
      common /mpierrh/ errh
c MPI constants
      integer MPI_COMM_WORLD
      parameter(MPI_COMM_WORLD=0)
c local data
      integer ierr
c check error code and print corresponding message
      if (ierror.eq.1) then
         write (2,*) source, 'MPI not initialized'
      elseif (ierror.eq.2) then
         write (2,*) source, 'Invalid Communicator'
      elseif (ierror.eq.3) then
         write (2,*) source, 'Invalid count'
      elseif (ierror.eq.4) then
         write (2,*) source, 'Invalid destination'
      elseif (ierror.eq.5) then
         write (2,*) source, 'Invalid source'
      elseif (ierror.eq.6) then
         write (2,*) source, 'Invalid tag'
      elseif (ierror.eq.7) then
         write (2,*) source, 'Invalid datatype'
      elseif (ierror.eq.12) then
         write (2,*) source, 'Incomplete read'
      elseif (ierror.eq.16) then
         write (2,*) source, 'Invalid request handle'
      elseif (ierror.eq.18) then
         write (2,*) source, 'Mismatched datatype'
      elseif (ierror.eq.19) then
         write (2,*) source, 'Invalid root'
      elseif (ierror.eq.20) then
         write (2,*) source, 'Invalid operation'
      elseif (ierror.eq.21) then
         write (2,*) source, 'Unable to allocate memory'
      elseif (ierror.eq.29) then
         write (2,*) source, 'Invalid rank or communicator'
      elseif (ierror.eq.(-1)) then
         write (2,*) source, 'Orderly disconnect request received'
      elseif (ierror.eq.(-2)) then
         write (2,*) source, 'Disorderly disconnect request received'
      elseif (ierror.eq.(-3)) then
         write (2,*) source, 'Unexpected flag returned on OTRcv'
      elseif (ierror.eq.(-9)) then
         write (2,*) source, 'Escape requested'
c unlisted error code
      else
         write (2,*) source, 'Error code = ', ierror
      endif
c abort if error is fatal
c only MPI_COMM_WORLD is currently supported
      if (errh(1).eq.1) then
         end file 2
         backspace 2
         call MPI_ABORT(MPI_COMM_WORLD,ierror,ierr)
         write (2,*) 'MPI_ABORT complete'
         close(unit=2)
         stop
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine rwstat(request,iunit)
c this subroutine writes out the status of a read-write record
c request = request handle
c iunit = Fortran output unit number
c input: request, iunit
      implicit none
      integer request, iunit
c declare common block for non-blocking messages
      integer MAXS, MAXM, curreq, header, rwrec, monitor, trash, mqueue
      parameter(MAXS=32,MAXM=2*MAXS)
      dimension curreq(5,MAXM), header(10,MAXM), rwrec(14,MAXM)
      dimension trash(256), mqueue(2,MAXS+1)
c rwrec = read/write record for asynchronous messages
      common /mpisendrec/ monitor, curreq, header, rwrec, trash, mqueue
c invalid request handle
      if ((request.lt.1).or.(request.gt.MAXM)) return
      if (curreq(1,request).eq.(-1)) then
         write (iunit,*) 'transmission mode = send'
      elseif (curreq(1,request).eq.1) then
         write (iunit,*) 'transmission mode = receive'
      endif
      write (iunit,*) 'destination/source = ', curreq(2,request)
      write (iunit,*) 'communicator = ', curreq(3,request)
      write (iunit,*) 'tag = ', curreq(4,request)
      write (iunit,*) 'datatype = ', curreq(5,request)
      write (iunit,*) 'request handle = ', request
      write (iunit,*) 'Endpoint reference pointer = ', rwrec(1,request)
      write (iunit,*) 'iocompletion flag = ', rwrec(2,request)
      write (iunit,*) 'current buffer pointer = ', rwrec(3,request)
      write (iunit,*) 'current buffer length = ', rwrec(4,request)
      write (iunit,*) 'T_MORE flag = ', rwrec(5,request)
      write (iunit,*) 'data pointer = ', rwrec(6,request)
      write (iunit,*) 'data length = ', rwrec(7,request)
      write (iunit,*) 'actual length sent/received = ', rwrec(8,request)
      write (iunit,*) 'first time stamp = ', rwrec(9,request), rwrec(10,
     1request)
      write (iunit,*) 'second time stamp = ', rwrec(11,request), rwrec(1
     12,request)
      write (iunit,*) 'next message pointer = ', rwrec(13,request)
      write (iunit,*) 'non-fatal error code = ', rwrec(14,request)
      write (iunit,*)
      return
      end
c-----------------------------------------------------------------------
      subroutine wqueue(iunit)
c this subroutine writes queue of pending messages
c iunit = Fortran output unit number
c input: iunit
      implicit none
      integer iunit
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c idproc = processor id
c ioc = array of context pointers for notifier function
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer MAXM, curreq, header, rwrec, monitor, trash, mqueue
      parameter(MAXM=2*MAXS)
      dimension curreq(5,MAXM), header(10,MAXM), rwrec(14,MAXM)
      dimension trash(256), mqueue(2,MAXS+1)
c mqueue = message request queue
      common /mpisendrec/ monitor, curreq, header, rwrec, trash, mqueue
c local data
      integer i, id, is
      is = 0
      do 10 i = 1, nproc
      id = i - 1
      if (ioc(3,i).ne.0) then
         write (iunit,*) 'Incomplete receive: node, handle =', id,
     1ioc(3,i)
         is = is + 1
      endif
      if (mqueue(1,i).ne.0) then
         write (iunit,*) 'Non-empty receive queue end: node, handle =',
     1id, mqueue(1,i)
         is = is + 1
      endif
      if (ioc(4,i).ne.0) then
         write (iunit,*) 'Incomplete send: node, handle =', id,
     1ioc(4,i)
         is = is + 1
      endif
      if (mqueue(2,i).ne.0) then
         write (iunit,*) 'Non-empty send queue end: node, handle =', id,
     1mqueue(2,i)
         is = is + 1
      endif
   10 continue
      i = MAXS + 1
      if (ioc(4,i).ne.0) then
         write (iunit,*) 'Incomplete selfsend: node, handle =', idproc,
     1ioc(4,i)
         is = is + 1
      endif
      if (mqueue(2,i).ne.0) then
         write (iunit,*) 'Non-empty selfsend queue end: node, handle =',
     1idproc, mqueue(2,i)
         is = is + 1
      endif
      if (is.gt.0) write (iunit,*)
      return
      end
c-----------------------------------------------------------------------
      subroutine messwin(nvp)
c this subroutine creates a window for showing MPI message status
c nvp = number of real or virtual processors
c input argument: nvp
      implicit none
      integer nvp
c function declarations
      integer*2 TextWidth
      integer GetMainDevice, NewCWindow, GetWindowPort
      external GetMainDevice, NewCWindow, GetWindowPort, TextWidth
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c common block for message window
c cpptr = pointer to window structure
c crect = current drag region
c nsp = amount of space between boxes
c nbx = size of box
c nds = number of message sizes monitored
c mbs = maximum maximum bandwidth in MB/sec expected
      integer cpptr
      integer*2 crect(4), nsp, nbx, nds, mbs
      common /winmess/ cpptr, crect, nsp, nbx, nds, mbs
c local data
      integer*2 wrect(4), trect(4), iw, is, ix, iy, it
      integer handle, wptr, n, i
      character*36 name
c get handle to main graphics device that carries a menu bar
      handle = GetMainDevice()
      wptr = long(handle)
c get size of screen
      crect(3) = word(wptr+38)
      crect(4) = word(wptr+40)
      wrect(3) = crect(3) - 40
      wrect(4) = crect(4) - 12
c find which grafPort is currently active
      call GetPort(wptr)
c calculate size of window
      n = min(nvp,nproc)
      wrect(1) = wrect(3) - 2*(nbx + nsp) - nbx
      wrect(2) = wrect(4) - ((nbx + nsp)*max(n,8) + nsp)
c add more space for distribution function
      wrect(1) = wrect(1) - (6*nbx + nsp)
c add more space for speedometer
      wrect(1) = wrect(1) - (4*nbx + nsp)
c add more space for user-defined label
      wrect(1) = wrect(1) - nbx
c name = label for window
      name = char(6)//'MacMPI'
c create a window
      if (cpptr.eq.0) then
         cpptr = NewCWindow(val4(0),wrect,name,val1(.true.),val2(4),val4
     1(-1),val1(.false.),val4(0))
      endif
c activate a GrafPort
      call SetPort(val4(GetWindowPort(val4(cpptr))))
c keep a rectangular area from being updated
      call GetPortBounds(GetWindowPort(val4(cpptr)),trect)
      call ValidWindowRect(val4(cpptr),trect)
c calculate clipping region
      wrect(3) = wrect(3) - wrect(1)
      wrect(4) = wrect(4) - wrect(2)
      wrect(1) = 0
      wrect(2) = 0
c select color to use in foreground drawing to black
      call BackColor(val4(33))
c fill rectangle with background pattern
      call EraseRect(wrect)
c calculate drag region
      crect(1) = 4
      crect(2) = 4
      crect(3) = crect(3) - 4
      crect(4) = crect(4) - 4
c select color to use in foreground drawing to white
      call ForeColor(val4(30))
c set dimensions of pen for current Grafport
      call PenSize(val2(1),val2(1))
c set point size for subsequent text drawing to 12 points
      call TextSize(val2(12))
c set initial x and y coordinates for label
      iw = nbx + nsp
      ix = nsp
      iy = iw + nbx
c write label for processor id
      do 10 i = 1, n
c convert node number to string
      write (name,'(i2)') i-1
c get width of unformatted text
      is = TextWidth(name,val2(0),val2(2))
c set pen location without drawing
      call MoveTo(val2(ix+(nbx-is)/2),val2(iy))
c draw text from any arbitrary buffer
      call DrawText(name,val2(0),val2(2))
c update x coordinate
      ix = ix + iw
   10 continue
c write label for message state display
      name = 'Send=Green,Receive=Red,Both=Yellow'
c set point size for subsequent text drawing to 10 points
      call TextSize(val2(10))
      ix = nsp
      iy = nsp + 3*nbx
c set pen location without drawing
      call MoveTo(val2(ix),val2(iy))
c select color to use in foreground drawing to green
      call ForeColor(val4(341))
c draw text from any arbitrary buffer
      call DrawText(name,val2(0),val2(11))
c get width of unformatted text
      is = TextWidth(name,val2(0),val2(11))
      ix = ix + is
c set pen location without drawing
      call MoveTo(val2(ix),val2(iy))
c select color to use in foreground drawing to red
      call ForeColor(val4(205))
c draw text from any arbitrary buffer
      call DrawText(name,val2(11),val2(12))
c get width of unformatted text
      is = TextWidth(name,val2(11),val2(12))
      ix = ix + is
c set pen location without drawing
      call MoveTo(val2(ix),val2(iy))
c select color to use in foreground drawing to yellow
      call ForeColor(val4(69))
c draw text from any arbitrary buffer
      call DrawText(name,val2(23),val2(11))
c select color to use in foreground drawing to white
      call ForeColor(val4(30))
c set point size for subsequent text drawing to 9 points
      call TextSize(val2(9))
c set initial x and y coordinates for label
      it = nbx/4
      iw = nbx/2 + nsp
      ix = nsp
      iy = 8*nbx + 2*nsp
c write label for message size, odd numbers only
      do 20 i = 1, nds, 2
c convert node number to string
      write (name,'(i2)') i-1
c get width of unformatted text
      is = TextWidth(name,val2(0),val2(2))
c set pen location without drawing
      call MoveTo(val2(ix+(it-is)/2),val2(iy))
c draw text from any arbitrary buffer
      call DrawText(name,val2(0),val2(2))
c underline location of half maximum
      if (i.eq.13) then
         call MoveTo(val2(ix+(it-is)/2),val2(iy+3))
         call DrawText('__',val2(0),val2(2))
      endif
c update x coordinate
      ix = ix + iw
   20 continue
c write second label
      name = 'Log2 Number vs. Log2 Message Size'
c set point size for subsequent text drawing to 10 points
      call TextSize(val2(10))
c set pen location without drawing
      call MoveTo(val2(nsp),val2(9*nbx+2*nsp))
c draw text from any arbitrary buffer
      call DrawText(name,val2(0),val2(33))
c set rectangle for speedometer
      wrect(1) = 9*nbx + 3*nsp
      wrect(2) = 2*nsp
      wrect(3) = wrect(1) + 4*nbx + 6
      iy = 11*nbx + 3*nsp + 2
c set point size for subsequent text drawing to 9 points
      call TextSize(val2(9))
      do 30 i = 1, 2
      if (i.eq.2) wrect(2) = wrect(2) + 4*(nbx + nsp) + (nsp - 2)
      wrect(4) = wrect(2) + 4*nbx + 1
c set dimensions of pen for current Grafport
      call PenSize(val2(2),val2(2))
c fill a wedge with current pen pattern and mode
      call PaintArc(wrect,val2(-90),val2(180))
      if (i.eq.1) then
c select color to use in foreground drawing to cyan
         call ForeColor(val4(273))
c write third label
         name = 'Communication % '
      elseif (i.eq.2) then
c select color to use in foreground drawing to red
         call ForeColor(val4(205))
c write third label
         name = 'Receiving (MB/s)'
      endif
c draw an arc
      call FrameArc(wrect,val2(-90),val2(180))
c set pen location without drawing
      call MoveTo(val2(wrect(2)),val2(iy+1))
c draw a line to specified coordinates
      call LineTo(val2(wrect(4)-2),val2(iy+1))
c set dimensions of pen for current Grafport
      call PenSize(val2(1),val2(1))
c select color to use in foreground drawing to white
      call ForeColor(val4(30))
c set pen location without drawing
      call MoveTo(val2(wrect(2)-nsp+2),val2(iy+nbx))
c draw text from any arbitrary buffer
      call DrawText(name,val2(0),val2(16))
c set pen location without drawing
      call MoveTo(val2(wrect(2)-nsp),val2(iy))
c draw text from any arbitrary buffer
      call DrawText('0',val2(0),val2(1))
c set pen location without drawing
      call MoveTo(val2(wrect(4)+2),val2(iy))
c convert node number to string
      if (i.eq.1) then
         write (name,'(i3)') 100
c draw text from any arbitrary buffer
         call DrawText(name,val2(0),val2(3))
      else
         write (name,'(i2)') mbs
c draw text from any arbitrary buffer
         call DrawText(name,val2(0),val2(2))
      endif
   30 continue
c write third label
      name = ' Current and Average Message-Passing'
c set point size for subsequent text drawing to 10 points
      call TextSize(val2(10))
c set pen location without drawing
      call MoveTo(val2(nsp),val2(13*nbx+3*nsp))
c draw text from any arbitrary buffer
      call DrawText(name,val2(0),val2(36))
c activate the GrafPort that was originally active
      if (wptr.ne.0) call SetPort(val4(wptr))
c display status
      call logmess(0,0,0,0,0)
      return
      end
c-----------------------------------------------------------------------
      subroutine logmess(idp,lstat,lsize,mticks,tag)
c this subroutine logs MPI message state change and displays status
c idp = remote processor id
c lstat = (-1,1,-2,2) = (clear send,add send,clear receive,add receive)
c lstat = 0 means display current status for all processors
c lsize = size of message (in bytes)
c lsize = -1 means print out message size distribution function
c mticks = wait time in microseconds
c tag = message tag
c input argument: idp, lstat, lsize, mticks, tag
      implicit none
      integer idp, lstat, lsize, mticks, tag
c function declarations
      integer GetWindowPort
      external GetWindowPort
c declare internal mpi common block
      integer nproc, idproc, cfig0
      integer MAXS, MAXC, MAXD, epref, ioc, nevents, stime, mapcomm
      parameter(MAXS=32,MAXC=10,MAXD=6)
      dimension epref(MAXS+1), ioc(4,MAXS+1), nevents(MAXS+1), stime(2)
      dimension mapcomm(MAXS+MAXD+3,MAXC)
c nproc = number of real or virtual processors obtained
c idproc = processor id
      common /mpiparms/ nproc, idproc, cfig0, epref, ioc, nevents, stime
     1, mapcomm
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c common block for message window
c cpptr = pointer to window structure
      integer cpptr
      integer*2 crect(4), nsp, nbx, nds, mbs
      common /winmess/ cpptr, crect, nsp, nbx, nds, mbs
c function declaration
      double precision MPI_WTIME
      external MPI_WTIME
c local data
      integer NSIZE, NDSIZE, NDSIZE2
      parameter(NSIZE=2*MAXS,NDSIZE=24,NDSIZE2=2*NDSIZE)
      integer mstat, msize, istat, istyle, lmax, nmax, i, i1, i2
      integer nums, numr
      real ar, cr, at, ct
      double precision ssize, times, rsize, timer, dt, ptime
c mstat = number of outstanding sends and receives for each node
c msize = message size distribution function
c nmax = maximum number of messages in display
c lmax = log2 of nmax
c ssize/times = total bytes/time sent
c rsize/timer = total bytes/time received
c ptime = time since last short term average
      dimension mstat(2,MAXS), msize(2*NDSIZE)
      dimension ssize(2), times(2), rsize(2), timer(2)
      save mstat, msize, lmax, nmax, nums, numr
      save ssize, times, rsize, timer, ptime
      data mstat /NSIZE*0/
      data msize /NDSIZE2*0/
      data lmax, nmax, nums, numr /8,256,0,0/
      data ssize, rsize, times, timer /2*0.0d0,2*0.0d0,2*0.0d0,2*0.0d0/
      data ptime /0.0d0/
c check for errors
      if ((idp.lt.0).or.(idp.ge.nproc)) return
c print out message size distribution function
      if (lsize.lt.0) then
         write (2,*) ' Message Size Distribution Function'
         do 10 i = 1, NDSIZE
         i1 = i + NDSIZE
         write (2,*) ' Size(bytes) = ', 2**(i-1), ' Sends = ', msize(i1)
     1, ' Receives = ', msize(i)-msize(i1)
   10    continue
         dt = 1.0d0/MPI_WTIME()
         i1 = 100.0*real(times(1)*dt) + .5
         ar = real(ssize(1)/times(1))
         write (2,*) ' Sending Time = ',i1,'%',', Speed = ',ar,'MB/s'
         i1 = 100.0*real(timer(1)*dt) + .5
         ar = real(rsize(1)/timer(1))
         write (2,*) ' Receiving Time = ',i1,'%',', Speed = ',ar,'MB/s'
         return
c process message size and time data
      elseif ((lstat.lt.0).and.(lsize.gt.0)) then
c accumulate sending data
         if (lstat.eq.(-1)) then
            dt = 1.0d-6*dble(lsize)
            ssize(1) = ssize(1) + dt
            ssize(2) = ssize(2) + dt
            dt = 1.0d-6*dble(mticks)
            times(1) = times(1) + dt
            times(2) = times(2) + dt
c calculate short term average
            nums = nums + 1
            if (nums.eq.4) then
               nums = 0
               ssize(2) = 0.0d
               times(2) = 0.0d0
            endif
c accumulate receiving data
         elseif (lstat.eq.(-2)) then
            dt = 1.0d-6*dble(lsize)
            rsize(1) = rsize(1) + dt
            rsize(2) = rsize(2) + dt
            dt = 1.0d-6*dble(mticks)
            timer(1) = timer(1) + dt
            timer(2) = timer(2) + dt
c calculate short term averge
            numr = numr + 1
            if (numr.eq.4) then
               dt = MPI_WTIME()
               ct = real(int(100.0*real((times(1)+timer(1))/dt) + .5))
               dt = dt - ptime
               at = real(int(100.0*real((times(2)+timer(2))/dt) + .5))
               ptime = MPI_WTIME()
               cr = real(rsize(1)/timer(1))
               ar = real(rsize(2)/timer(2))
               call shospeed(at,ct,ar,cr)
               numr = 0
               rsize(2) = 0.0d0
               timer(2) = 0.0d0
            endif
         endif
c increment message size distribution function
         i1 = min(int(alog(real(lsize))/alog(2.) + 1.5),NDSIZE)
         msize(i1) = msize(i1) + 1
         i2 = i1 + NDSIZE
c increment message size distribution function for sends
         if (lstat.eq.(-1)) msize(i2) = msize(i2) + 1
         if (msize(i1).gt.nmax) then
c erase display of distribution function
            call showdism(1,nmax,0,lmax,0)
            lmax = lmax + 1
            nmax = nmax + nmax
c redisplay distribution function
            do 20 i = 1, NDSIZE
            call showdism(i,msize(i),msize(i+NDSIZE),lmax,1)
   20       continue 
c display distribution function
         else
            call showdism(i1,msize(i1),msize(i2),lmax,1)
         endif
      endif
      i1 = idp + 1
c calculate all current statuses
      if (lstat.eq.0) then
         do 30 i = 1, nproc
         i1 = i - 1
         istat = 0
         if (mstat(1,i).ge.1) istat = istat + 1
         if (mstat(2,i).ge.1) istat = istat + 2
c differentiate single from multiple sends/receives
c        if ((mstat(1,i).gt.1).or.(mstat(2,i).gt.1)) then
c           istat = istat + 3
c        endif
c display status, outline local node
         if (i1.eq.idproc) then
            istyle = 1
         else
            istyle = 0
         endif
         call showmess(i1,istat,istyle)
   30    continue
c display scale
         call showdism(1,msize(1),msize(NDSIZE+1),lmax,0)
c display current distribution function
         do 40 i = 1, NDSIZE
         call showdism(i,msize(i),msize(i+NDSIZE),lmax,1)
   40    continue
         return
c add state change to log
      elseif (lstat.eq.1) then
         mstat(1,i1) = mstat(1,i1) + 1
         if (monitor.eq.2) then
            write (2,*) 'send posted: destination=', idp, ' size=' 
     1, lsize, 'tag=', tag
         endif
      elseif (lstat.eq.(-1)) then
         mstat(1,i1) = mstat(1,i1) - 1
         if (monitor.eq.2) then
            write (2,*) 'sent: destination=', idp, ' size=', lsize
     1, ' time=', mticks, 'tag=', tag
         endif
      elseif (lstat.eq.2) then
         mstat(2,i1) = mstat(2,i1) + 1
         if (monitor.eq.2) then
            write (2,*) 'receive posted: source=', idp, ' size=', lsize
     1, 'tag=', tag
         endif
      elseif (lstat.eq.(-2)) then
         mstat(2,i1) = mstat(2,i1) - 1
         if (monitor.eq.2) then
            write (2,*) 'received: source=', idp, ' size=', lsize
     1, ' time=', mticks, 'tag=', tag
         endif
      endif
c calculate current status
      istat = 0
      if (mstat(1,i1).ge.1) istat = istat + 1
      if (mstat(2,i1).ge.1) istat = istat + 2
c differentiate single from multiple sends/receives
c     if ((mstat(1,i1).gt.1).or.(mstat(2,i1).gt.1)) then
c        istat = istat + 3
c     endif
c display status, outline local node
      if (idp.eq.idproc) then
         istyle = 1
      else
         istyle = 0
      endif
      call showmess(idp,istat,istyle)
      if (cpptr.ne.0) then
c flush window buffer to screen
        call QDFlushPortBuffer(val4(GetWindowPort(val4(cpptr))),val4(0))
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine showmess(idp,istat,istyle)
c this subroutine shows MPI message status
c idp = remote processor id
c istat = message status = (0,1,2,3) = (none,sending,receiving,both)
c istyle = (0,1) = (no,yes) outline rectangle
c input argument: idp, istat, istyle
      implicit none
      integer idp, istat, istyle
c function declarations
      integer*2 TextWidth
      integer GetWindowPort
      external GetWindowPort, TextWidth
c common block for message window
c cpptr = pointer to window structure
c nsp = amount of space between boxes
c nbx = size of box
      integer cpptr
      integer*2 crect(4), nsp, nbx, nds, mbs
      common /winmess/ cpptr, crect, nsp, nbx, nds, mbs
c local data
      integer*2 wrect(4)
      integer wptr, icolor(8)
c icolor = white, green, red, yellow, blue, magenta, cyan, black
      data icolor /30,341,205,69,409,137,273,33/
c check for errors
      if ((istat.lt.0).or.(istat.gt.7)) return
      if (cpptr.eq.0) return
c find which grafPort is currently active
      call GetPort(wptr)
c activate a GrafPort
      call SetPort(val4(GetWindowPort(val4(cpptr))))
c set rectangle
      wrect(1) = nsp
      wrect(2) = (nbx + nsp)*idp + nsp
      wrect(3) = wrect(1) + nbx
      wrect(4) = wrect(2) + nbx
c select color to use in foreground drawing
      call ForeColor(val4(icolor(istat+1)))
c draw the outline of a rectangle
      call FrameRect(wrect)
c fill rectangle with current pen pattern and mode
      call PaintRect(wrect)
c outline rectangle
      if (istyle.ne.0) then
c shrink or expand a rectangle
         call InsetRect(wrect,val2(-2),val2(-2))
c draw the outline of a rectangle
         call FrameRect(wrect)
      endif
c activate the GrafPort that was originally active
      if (wptr.ne.0) call SetPort(val4(wptr))
      return
      end
c-----------------------------------------------------------------------
      subroutine showdism(ibin,nbin,mbin,lmax,istyle)
c this subroutine shows distribution of MPI messages
c ibin = bin number for distribution function
c nbin = number of total messages in ibin
c mbin = number of send messages in ibin
c lmax = log2 of maximum number of messages in display
c istyle = (0,1) = (erase and rescale,draw) display 
c input argument: ibin, nbin, nmax
      implicit none
      integer ibin, nbin, mbin, lmax, istyle
c function declarations
      integer*2 TextWidth
      integer GetWindowPort
      external TextWidth, GetWindowPort
c common block for message window
c cpptr = pointer to window structure
c nsp = amount of space between boxes
c nbx = size of box
c nds = number of message sizes monitored
      integer cpptr
      integer*2 crect(4), nsp, nbx, nds, mbs
      common /winmess/ cpptr, crect, nsp, nbx, nds, mbs
c local data
      integer*2 wrect(4), iw, nscale, nsub
      integer wptr
      real scale
      character*8 name
c check for errors
      if (cpptr.eq.0) return
c find which grafPort is currently active
      call GetPort(wptr)
c activate a GrafPort
      call SetPort(val4(GetWindowPort(val4(cpptr))))
c set rectangle for entire distribution
      iw = nbx/4
      wrect(3) = 7*nbx + 2*nsp
c erase and rescale display
      if (istyle.eq.0) then
c select color to use in foreground drawing to black
         call ForeColor(val4(33))
         wrect(1) = wrect(3) - 4*nbx
         wrect(4) = (iw + nsp/2)*nds + nsp/2
         wrect(2) = 0
c fill rectangle with background pattern
         call EraseRect(wrect)
c select color to use in foreground drawing to white
         call ForeColor(val4(30))
c set point size for subsequent text drawing to 9 points
         call TextSize(val2(9))
c convert node number to string
         write (name,'(i2)') lmax
c set pen location without drawing
         call MoveTo(val2(0),val2(wrect(1)+9))
c draw text from any arbitrary buffer
         call DrawText(name,val2(0),val2(2))
c draw display
      else
c calculate scale of total image
         scale = real(4*nbx)/(real(lmax)*alog(2.))
c calculate size of image
         nscale = min(int(alog(real(nbin+1))*scale),4*nbx)
         if (nbin.eq.0) go to 10
c set width of individual bin
         wrect(4) = (iw + nsp/2)*ibin + nsp/2
         wrect(2) = wrect(4) - iw
c calculate size of image for sends
         nsub = (real(mbin)/real(nbin))*real(nscale) + .5
c set height of individual bin for sends
         wrect(1) = wrect(3) - nsub
c select color to use in foreground drawing to green
         call ForeColor(val4(341))
c draw the outline of a rectangle
         call FrameRect(wrect)
c fill rectangle with current pen pattern and mode
         call PaintRect(wrect)
c set height of individual bin for receives
         wrect(1) = wrect(3) - nscale
         wrect(3) = wrect(3) - nsub
c select color to use in foreground drawing to red
         call ForeColor(val4(205))
c draw the outline of a rectangle
         call FrameRect(wrect)
c fill rectangle with current pen pattern and mode
         call PaintRect(wrect)
      endif
c activate the GrafPort that was originally active
   10 if (wptr.ne.0) call SetPort(val4(wptr))
      return
      end
c-----------------------------------------------------------------------
      subroutine shospeed(atime,ctime,arate,crate)
c this subroutine shows communication rates for MPI messages
c atime = short term average communication time (% of total)
c ctime = long term average communication time (% of total)
c arate = short term average reception rate (MB/sec)
c crate = long term average reception rate (MB/sec)
c input argument: atime, ctime, arate, crate
      implicit none
      real atime, ctime, arate, crate
c function declarations
      integer GetWindowPort
      external GetWindowPort
c common block for message window
c cpptr = pointer to window structure
c nsp = amount of space between boxes
c nbx = size of box
c mbs = maximum maximum bandwidth in MB/sec expected
      integer cpptr
      integer*2 crect(4), nsp, nbx, nds, mbs
      common /winmess/ cpptr, crect, nsp, nbx, nds, mbs
c local data
      real scale, pi
      parameter(pi=.5*6.28318530717959)
      integer*2 ix0, iy0, ix(4), iy(4)
      integer wptr, i
      real angle, dx
      save ix, iy
      data ix, iy /4*0,4*0/
c check for errors
      if (cpptr.eq.0) return
c find which grafPort is currently active
      call GetPort(wptr)
c activate a GrafPort
      call SetPort(val4(GetWindowPort(val4(cpptr))))
c set dimensions of pen for current Grafport
      call PenSize(val2(2),val2(2))
c find center of speedometer
      iy0 = 11*nbx + 3*nsp + 1
c parameters for communication plot
      ix0 = 2*(nbx + nsp) + 1
      scale = .01
c plot speedometer lines
      do 10 i = 1, 4
c parameters for reception plot
      if (i.eq.3) then
         ix0 = ix0 + 4*(nbx + nsp) + (nsp - 2)
         scale = 1./real(mbs)
      endif
c erase old speedometer line
c select color to use in foreground drawing to white
      call ForeColor(val4(30))
c set pen location without drawing
      call MoveTo(val2(ix0),val2(iy0))
c draw a line a specified distance
      call Line(val2(ix(i)),val2(-iy(i)))
c draw new speedometer line
      if (i.eq.1) then
c select color to use in foreground drawing to cyan
         call ForeColor(val4(273))
         angle = max(0.,(1. - atime*scale))*pi
      elseif (i.eq.2) then
c select color to use in foreground drawing to black
         call ForeColor(val4(33))
         angle = max(0.,(1. - ctime*scale))*pi
      elseif (i.eq.3) then
c select color to use in foreground drawing to red
         call ForeColor(val4(205))
         angle = max(0.,(1. - arate*scale))*pi
      elseif (i.eq.4) then
c select color to use in foreground drawing to black
         call ForeColor(val4(33))
         angle = max(0.,(1. - crate*scale))*pi
      endif
      dx = real(2*nbx-3)
      ix(i) = dx*cos(angle) + .5
      iy(i) = dx*sin(angle) + .5
c set pen location without drawing
      call MoveTo(val2(ix0),val2(iy0))
c draw a line a specified distance
      call Line(val2(ix(i)),val2(-iy(i)))
   10 continue
c set dimensions of pen for current Grafport
      call PenSize(val2(1),val2(1))
c activate the GrafPort that was originally active
      if (wptr.ne.0) call SetPort(val4(wptr))
      return
      end
c-----------------------------------------------------------------------
      subroutine LOGNAME(name)
c this subroutine records and displays user-defined label
c name = label to display
      implicit none
      character*(*) name
      integer MAXS
      parameter(MAXS=32)
c function declarations
      integer GetWindowPort
      external GetWindowPort
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c common block for message window
c cpptr = pointer to window structure
c nsp = amount of space between boxes
c nbx = size of box
      integer cpptr
      integer*2 crect(4), nsp, nbx, nds, mbs
      common /winmess/ cpptr, crect, nsp, nbx, nds, mbs
c local data
      integer*2 wrect(4), nl
      integer wptr
      if (monitor.eq.0) return
      if (monitor.eq.2) write (2,*) name
c check for errors
      if (cpptr.eq.0) return
c find which grafPort is currently active
      call GetPort(wptr)
c activate a GrafPort
      call SetPort(val4(GetWindowPort(val4(cpptr))))
c set rectangle
      wrect(1) = 13*nbx + 3*nsp + 2
      wrect(2) = nsp
      wrect(3) = wrect(1) + nbx
      wrect(4) = wrect(2) + 8*(nbx + nsp)
c select color to use in foreground drawing to black
      call BackColor(val4(33))
c fill rectangle with background pattern
      call EraseRect(wrect)
c select color to use in foreground drawing to white
      call ForeColor(val4(30))
c set pen location without drawing
      call MoveTo(val2(wrect(2)),val2(wrect(3)-2))
      nl = min(36,len(name))
c draw text from any arbitrary buffer
      call DrawText(name,val2(0),val2(nl))
c activate the GrafPort that was originally active
      if (wptr.ne.0) call SetPort(val4(wptr))
      return
      end
c-----------------------------------------------------------------------
      subroutine SET_MON(monval)
c this subroutine sets new monitor value and corresponding window
c monval = new monitor value
      implicit none
      integer monval
c declare internal mpi common block
      integer nproc
c nproc = number of real or virtual processors obtained
      common /mpiparms/ nproc
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
c create or destroy window if MPI has been initialized
      if (nproc.gt.0) then
c open window
         if ((monitor.eq.0).and.(monval.gt.0)) call messwin(nproc)
c close window
         if ((monitor.gt.0).and.(monval.lt.1)) call delmess()
      endif
c reset monitor value
      if (monval.gt.1) then
         monitor = 2
      elseif (monval.lt.1) then
         monitor = 0
      else
         monitor = 1
      endif
      return
      end
c-----------------------------------------------------------------------
      function GET_MON()
c this function gets current monitor value
      implicit none
      integer GET_MON
c declare common block for non-blocking messages
      integer monitor
c monitor = (0,1,2) = (suppress,display,display & log) monitor messages
      common /mpisendrec/ monitor
      GET_MON = monitor
      return
      end
c-----------------------------------------------------------------------
      subroutine delmess()
c this subroutine closes a window for showing MPI message status
c common block for message window
c cpptr = pointer to window structure
      implicit none
      integer cpptr
      integer*2 crect(4), nsp, nbx, nds, mbs
      common /winmess/ cpptr, crect, nsp, nbx, nds, mbs
c remove window from screen; keep WindowRecord
      if (cpptr.ne.0) call DisposeWindow(val4(cpptr))
      cpptr = 0
      return
      end
