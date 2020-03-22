c-----------------------------------------------------------------------
c Multitasking library for Apple Macintosh
c based on the Multiprocessing Library in the MacOS, v. 2.1
c Fortran unit 2 is used throughout for error messages
c written by viktor k. decyk, ucla
c copyright 2000, regents of the university of california.
c all rights reserved.
c no warranty for proper operation of this software is given or implied.
c software or information may be copied, distributed, and used at own
c risk; it may not be distributed without this notice included verbatim
c with each file. 
c update: may 5, 2003
      block data
      implicit none
c common block for tasks
      integer MAXARGS, MAXTASKS
      parameter(MAXARGS=24,MAXTASKS=16)
      integer notifyq, ntask, itask, params, mqueues, stacksize
      dimension ntask(MAXTASKS), itask(MAXTASKS)
      dimension params(MAXARGS+4,MAXTASKS), mqueues(2,MAXTASKS)
      common /tasks/ notifyq, ntask, itask, params, mqueues, stacksize
      save /tasks/
      data notifyq, stacksize /0,0/
      end
c-----------------------------------------------------------------------
      subroutine MP_INIT(nproc)
c initialize multitasking environment
c nproc = number of processors, 0 if Multitasking is not available
c output: nproc
      implicit none
      integer nproc
c function declarations
      integer*2 GetSharedLibrary, FindSymbol
      integer exfunc, MPProcessors, MPCreateQueue
      external GetSharedLibrary, FindSymbol
      external exfunc, MPProcessors, MPCreateQueue
c common block for tasks
c MAXARGS = maximum number of arguments in task procedure
c MAXTASKS = maximum number of tasks supported
c notifyq = notifyqueue record for tasks
c ntask = task record
c itask = MPTaskID record
c mqueues = message queues for tasks (1=request,2=result)
c stacksize = size of stack for task in bytes (0 = default)
      integer MAXARGS, MAXTASKS
      parameter(MAXARGS=24,MAXTASKS=16)
      integer notifyq, ntask, itask, params, mqueues, stacksize
      dimension ntask(MAXTASKS), itask(MAXTASKS)
      dimension params(MAXARGS+4,MAXTASKS), mqueues(2,MAXTASKS)
      common /tasks/ notifyq, ntask, itask, params, mqueues, stacksize
c local data
      integer kCompiledCFragArch, kReferenceCFrag
      parameter(kCompiledCFragArch=1886875747,kReferenceCFrag=1)
      integer*2 ierr
      integer i, response, connid, mpinit, oss
c Check if MP is already initialized
      if (notifyq.ne.0) then
c Return the number of processors on the host computer
         nproc = MPProcessors()
         return
      endif
      nproc = 0
      mpinit = 0
      notifyq = 0
c Get information about the operating environment
      call Gestalt(VAL4('sysv'),response)
c Check if MacOS installed is Classic
      if (response.lt.4096) then
c Determine if MPLibrary is available
         ierr = GetSharedLibrary(char(9)//'MPLibrary',
     1val4(kCompiledCFragArch),val4(kReferenceCFrag),
     2connid,val4(0),val4(0))
         if (ierr.ne.0) then
            write (2,*) 'GetSharedLibrary Error, ierr = ', ierr
            return
         endif
c Determine if MP API Services are available
         ierr = FindSymbol(val4(connid),char(21)//'_MPIsFullyInitialized
     1',mpinit,val4(0))
         if (ierr.ne.0) then
            write (2,*) 'FindSymbol Error, ierr = ', ierr
            return
         endif
c Execute Boolean procedure
         nproc = exfunc(val4(mpinit))
c Assume that MPIsFullyInitialized is always true with Mac OS X
      else
c        nproc = MPIsFullyInitialized()
         nproc = 1
      endif
c Exit if MP API Services are not available
      if (nproc.eq.0) return
c Return the number of processors on the host computer
      nproc = MPProcessors()
c Create a message queue
      oss = MPCreateQueue(notifyq)
      if (oss.ne.0) then
         write (2,*) 'MPCreateQueue error, oss = ', oss
         if (oss.eq.(-4)) then
            write (2,*) 'Unimplemented Core Routine'
         endif
         nproc = 0
         notifyq = 0
         return
      endif
c Clear task records
      do 10 i = 1, MAXTASKS
      ntask(i) = 0
      itask(i) = 0
      mqueues(1,i) = 0
      mqueues(2,i) = 0
   10 continue
c set current stacksize (in bytes) for tasks
      stacksize = 16384
      return
      end
c-----------------------------------------------------------------------
      subroutine MP_CR8TASK(entrypoint,ttype,proc,parameter,nargs,taskid
     1)
c start task
c overhead is about 130 microseconds on Macintosh G4/450
c entrypoint = pointer to task procedure
c ttype = task type (0=non-reusable,1=reusable(signaling),2=non-running)
c proc = pointer to actual procedure
c parameter = pointer to task procedure argument
c nargs = number of arguments in actual procedure
c taskid = index to task record (0 if error)
c input: all except taskid
c output: taskid
      implicit none
      integer entrypoint, ttype, proc, parameter, taskid, nargs
      dimension parameter(*)
c function declarations
      integer MPCreateTask, MPCreateSemaphore
      external MPCreateTask, MPCreateSemaphore
c common block for tasks
c MAXARGS = maximum number of arguments in task procedure
c MAXTASKS = maximum number of tasks supported
c notifyq = notifyqueue record for tasks
c ntask = task record
c itask = MPTaskID record
c params = arguments in task procedures
c mqueues = message queues for tasks (1=request,2=result)
c stacksize = size of stack for task in bytes (0 = default)
      integer MAXARGS, MAXTASKS
      parameter(MAXARGS=24,MAXTASKS=16)
      integer notifyq, ntask, itask, params, mqueues, stacksize
      dimension ntask(MAXTASKS), itask(MAXTASKS)
      dimension params(MAXARGS+4,MAXTASKS), mqueues(2,MAXTASKS)
      common /tasks/ notifyq, ntask, itask, params, mqueues, stacksize
c local data
      integer i, j, oss, tp, options
      external tproc, entrypoint
      save tp
      data tp, options /0,0/
      taskid = 0
c Check for errors
      if ((notifyq.eq.0).and.(ttype.lt.2)) then
         write (2,*) 'MP_CR8TASK: MP not initialized'
         return
      else if (nargs.gt.MAXARGS) then
         write (2,*) 'Too many arguments in task, nargs = ', nargs,
     1', MAXARGS = ', MAXARGS
         return
      endif
c Find space for record
      i = 0
   10 i = i + 1
      if (i.gt.MAXTASKS) then
         write (2,*) 'Exceeded maximum number of tasks = ', MAXTASKS
         return
      else if (ntask(i).ne.0) then
         go to 10
      endif
c Create semaphores
      if (ttype.eq.1) then
         do 30 j = 1, 2
         oss = MPCreateSemaphore(val4(1),val4(0),mqueues(j,i))
         if (oss.ne.0) then
            write (2,*) 'MPCreateSemaphore error, oss = ', oss
            mqueues(j,i) = 0
            return
         endif
   30    continue
      endif
      ntask(i) = i
c Copy number of arguments
      params(1,i) = nargs
c Copy arguments
      do 40 j = 1, nargs
      params(j+1,i) = parameter(j)
   40 continue
c Copy procedure name and message queue pointers
      params(nargs+2,i) = loc(proc)
      params(nargs+3,i) = mqueues(1,i)
      params(nargs+4,i) = mqueues(2,i)
c Create a preemptive task
      if (ttype.lt.2) then
         oss = MPCreateTask(entrypoint,params(1,i),val4(stacksize),
     1val4(notifyq),ntask(i),tp,val4(options),itask(i))
         if (oss.ne.0) then
            write (2,*) 'MPCreateTask error, oss = ', oss
            ntask(i) = 0
            itask(i) = 0
            return
         endif
      endif
      taskid = i
      return
      end
c-----------------------------------------------------------------------
      subroutine MP_TASKWAIT(taskid)
c wait for task to complete
c taskid = index to task record (0 if task completed successfully)
c input and output: taskid
      implicit none
      integer taskid
c function declarations
      integer MPWaitOnQueue, MPDeleteSemaphore
      external MPWaitOnQueue, MPDeleteSemaphore
c common block for tasks
c MAXARGS = maximum number of arguments in task procedure
c MAXTASKS = maximum number of tasks supported
c notifyq = notifyqueue record for tasks
c ntask = task record
c itask = MPTaskID record
c mqueues = message queues for tasks (1=request,2=result)
      integer MAXARGS, MAXTASKS
      parameter(MAXARGS=24,MAXTASKS=16)
      integer notifyq, ntask, itask, params, mqueues, stacksize
      dimension ntask(MAXTASKS), itask(MAXTASKS)
      dimension params(MAXARGS+4,MAXTASKS), mqueues(2,MAXTASKS)
      common /tasks/ notifyq, ntask, itask, params, mqueues, stacksize
c local data
      integer kDurationForever
      parameter(kDurationForever=2147483647)
      integer i, j, oss, message, tp, result, id
      i = taskid
c Check for errors
      if (notifyq.eq.0) then
         write (2,*) 'MP_TASKWAIT: MP not initialized'
         return
      else if ((i.lt.1).or.(i.gt.MAXTASKS)) then
         write (2,*) 'MP_TASKWAIT: Invalid taskid = ', taskid
         return
      else if ((mqueues(2,i).ne.0).and.(itask(i).ne.0)) then
         write (2,*) 'MP_TASKWAIT: Invalid type for task = ', taskid
         return
      endif
c Look for message from taskid
   10 if (ntask(i).gt.0) then
c Obtain a message from a specified message queue
         oss = MPWaitOnQueue(val4(notifyq),message,tp,result,
     1val4(kDurationForever))
         if (oss.ne.0) then
            write (2,*) 'MPWaitOnQueue error for taskid, oss = ', i, oss
            return
         endif
         id = long(message)
         if ((id.lt.1).or.(id.gt.MAXTASKS)) then
            write (2,*) 'Invalid taskid from notifyq = ', id
         else if (ntask(id).eq.0) then
            write (2,*) 'Task from notifyq already ended, taskid = ', id
         else
c Mark task as completed, may not be the one we are waiting for
            ntask(id) = -id
         endif
         go to 10
c Task marked as completed earlier
      else if (ntask(i).lt.0) then
         ntask(i) = 0
      endif
c Remove semaphores
      do 20 j = 1, 2
      if (mqueues(j,i).ne.0) then
         oss = MPDeleteSemaphore(val4(mqueues(j,i)))
         if (oss.ne.0) then
            write (2,*) 'MPDeleteSemaphore error, oss = ', oss
         endif
         mqueues(j,i) = 0
      endif
   20 continue
      itask(i) = 0
      taskid = 0
      return
      end
c-----------------------------------------------------------------------
      integer function MP_SNDSIG(taskid)
c send signal to task
c overhead for sending a signal to a task and receiving a reply 
c is about 60 microseconds on Macintosh G4/450
c taskid = index to task record
c input: all
      implicit none
      integer taskid
c function declarations
      integer MPSignalSemaphore
      external MPSignalSemaphore
c common block for tasks
c MAXARGS = maximum number of arguments in task procedure
c MAXTASKS = maximum number of tasks supported
c notifyq = notifyqueue record for tasks
c ntask = task record
c mqueues = message queues for tasks (1=request,2=result)
      integer MAXARGS, MAXTASKS
      parameter(MAXARGS=24,MAXTASKS=16)
      integer notifyq, ntask, itask, params, mqueues, stacksize
      dimension ntask(MAXTASKS), itask(MAXTASKS)
      dimension params(MAXARGS+4,MAXTASKS), mqueues(2,MAXTASKS)
      common /tasks/ notifyq, ntask, itask, params, mqueues, stacksize
c local data
      integer i, oss
      i = taskid
      MP_SNDSIG = 0
c Check for errors
      if (notifyq.eq.0) then
         MP_SNDSIG = 1
         write (2,*) 'MP_SNDSIG: MP not initialized'
         return
      else if ((i.lt.1).or.(i.gt.MAXTASKS)) then
         MP_SNDSIG = 2
         write (2,*) 'MP_SNDSIG: Invalid taskid = ', taskid
         return
      else if (ntask(i).eq.0) then
         MP_SNDSIG = 3
         write (2,*) 'MP_SNDSIG: Task already ended, taskid = ', taskid
         return
      else if (mqueues(1,i).eq.0) then
         MP_SNDSIG = 4
         write (2,*) 'MP_SNDSIG: Invalid type for task = ', taskid
         return
      endif
c Signal a semaphore
      oss = MPSignalSemaphore(val4(mqueues(1,i)))
      if (oss.ne.0) then
         MP_SNDSIG = oss
         write (2,*) 'MPSignalSemaphore error for taskid, oss = ', i,oss
      endif
      return
      end
c-----------------------------------------------------------------------
      integer function MP_WAITSIG(taskid)
c receive signal from task
c overhead for sending a signal to a task and receiving a reply 
c is about 60 microseconds on Macintosh G4/450
c taskid = index to task record
c input: taskid
      implicit none
      integer taskid
c function declarations
      integer MPWaitOnSemaphore
      external MPWaitOnSemaphore
c common block for tasks
c MAXARGS = maximum number of arguments in task procedure
c MAXTASKS = maximum number of tasks supported
c notifyq = notifyqueue record for tasks
c ntask = task record
c mqueues = message queues for tasks (1=request,2=result)
      integer MAXARGS, MAXTASKS
      parameter(MAXARGS=24,MAXTASKS=16)
      integer notifyq, ntask, itask, params, mqueues, stacksize
      dimension ntask(MAXTASKS), itask(MAXTASKS)
      dimension params(MAXARGS+4,MAXTASKS), mqueues(2,MAXTASKS)
      common /tasks/ notifyq, ntask, itask, params, mqueues, stacksize
c local data
      integer kDurationForever
      parameter(kDurationForever=2147483647)
      integer i, oss
      i = taskid
      MP_WAITSIG = 0
c Check for errors
      if (notifyq.eq.0) then
         MP_WAITSIG = 1
         write (2,*) 'MP_WAITSIG: MP not initialized'
         return
      else if ((i.lt.1).or.(i.gt.MAXTASKS)) then
         MP_WAITSIG = 2
         write (2,*) 'MP_WAITSIG: Invalid taskid = ', taskid
         return
      else if (ntask(i).eq.0) then
         MP_WAITSIG = 3
         write (2,*) 'MP_WAITSIG: Task already ended, taskid = ', taskid
         return
      else if (mqueues(2,i).eq.0) then
         MP_WAITSIG = 4
         write (2,*) 'MP_WAITSIG: Invalid type for task = ', taskid
         return
      endif
c Wait on semaphore
      oss = MPWaitOnSemaphore(val4(mqueues(2,i)),val4(kDurationForever))
      if (oss.ne.0) then
          MP_WAITSIG = 5
         write (2,*) 'MPWaitOnSemaphore error for taskid, oss = ', i,oss
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MP_KILLTASK(taskid)
c terminate multitasking environment
c taskid = index to task record (0 if task killed successfully)
      implicit none
      integer taskid
c function declarations
      integer MPTerminateTask
      external MPTerminateTask
c common block for tasks
c MAXARGS = maximum number of arguments in task procedure
c MAXTASKS = maximum number of tasks supported
c notifyq = notifyqueue record for tasks
c itask = MPTaskID record
      integer MAXARGS, MAXTASKS
      parameter(MAXARGS=24,MAXTASKS=16)
      integer notifyq, ntask, itask, params, mqueues, stacksize
      dimension ntask(MAXTASKS), itask(MAXTASKS)
      dimension params(MAXARGS+4,MAXTASKS), mqueues(2,MAXTASKS)
      common /tasks/ notifyq, ntask, itask, params, mqueues, stacksize
c local data
      integer i, oss, result
      data result /1/
      i = taskid
c Check for errors
      if (notifyq.eq.0) then
         write (2,*) 'MP_KILLTASK: MP not initialized'
         return
      else if ((i.lt.1).or.(i.gt.MAXTASKS)) then
         write (2,*) 'MP_KILLTASK: Invalid taskid = ', taskid
         return
      else if (ntask(i).eq.0) then
         write (2,*) 'Task already killed, taskid = ', taskid
         return
      endif
c Terminate an existing task
      oss = MPTerminateTask(val4(itask(i)),val4(result))
      if (oss.ne.0) then
         write (2,*) 'MPTerminateTask error for taskid, oss = ', i, oss
      else
         itask(i) = 0
         call MP_TASKWAIT(taskid)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MP_END()
c terminate multitasking environment
      implicit none
c function declarations
      integer MPDeleteQueue
      external MPDeleteQueue
c common block for tasks
c MAXARGS = maximum number of arguments in task procedure
c MAXTASKS = maximum number of tasks supported
c notifyq = notifyqueue record for tasks
      integer MAXARGS, MAXTASKS
      parameter(MAXARGS=24,MAXTASKS=16)
      integer notifyq, ntask, itask, params, mqueues, stacksize
      dimension ntask(MAXTASKS), itask(MAXTASKS)
      dimension params(MAXARGS+4,MAXTASKS), mqueues(2,MAXTASKS)
      common /tasks/ notifyq, ntask, itask, params, mqueues, stacksize
c local data
      integer i, oss, taskid
c Check for errors
      if (notifyq.eq.0) then
         write (2,*) 'MP_END: MP not initialized'
         return
      endif
c Check if any tasks are outstanding
      do 10 i = 1, MAXTASKS
      if (ntask(i).gt.0) then
         taskid = i
         write (2,*) 'MP_END: Task still outstanding, taskid = ', i
         call MP_KILLTASK(taskid)
      endif
   10 continue
c Delete a message queue
      oss = MPDeleteQueue(val4(notifyq))
      if (oss.ne.0) then
         write (2,*) 'MPDeleteQueue error, oss = ', oss
      endif
      notifyq = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine MP_TASKSTART(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5
     1,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16
     2,arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24)
c create a task by packing arguments into a single structure
c taskid = index to notify queue for task (0 if error)
c proc = pointer to actual procedure
c nargs = number of arguments in actual procedure
c arg1-arg24 = arguments for actual procedure proc
c restrictions: actual arguments should not be temporary variables, such
c as return values of functions.  Also, local variables in procedure
c which called MP_TASKSTART should not be used if MP_TASKWAIT is not
c called in the same procedure.
c character variables should also be avoided.
      implicit none
      integer taskid, nargs
      external proc
      integer arg1(*), arg2(*), arg3(*), arg4(*), arg5(*), arg6(*)
      integer arg7(*), arg8(*), arg9(*), arg10(*), arg11(*), arg12(*)
      integer arg13(*), arg14(*), arg15(*), arg16(*), arg17(*)
      integer arg18(*), arg19(*), arg20(*), arg21(*), arg22(*)
      integer arg23(*), arg24(*)
c local data
      integer MAXARGS
c MAXARGS = maximum number of arguments in task procedure
      parameter(MAXARGS=24)
      integer param
      dimension param(MAXARGS)
      external tproc
c Create argument list for task
      if (nargs.eq.0) go to 10
      param(1) = loc(arg1)
      if (nargs.eq.1) go to 10
      param(2) = loc(arg2)
      if (nargs.eq.2) go to 10
      param(3) = loc(arg3)
      if (nargs.eq.3) go to 10
      param(4) = loc(arg4)
      if (nargs.eq.4) go to 10
      param(5) = loc(arg5)
      if (nargs.eq.5) go to 10
      param(6) = loc(arg6)
      if (nargs.eq.6) go to 10
      param(7) = loc(arg7)
      if (nargs.eq.7) go to 10
      param(8) = loc(arg8)
      if (nargs.eq.8) go to 10
      param(9) = loc(arg9)
      if (nargs.eq.9) go to 10
      param(10) = loc(arg10)
      if (nargs.eq.10) go to 10
      param(11) = loc(arg11)
      if (nargs.eq.11) go to 10
      param(12) = loc(arg12)
      if (nargs.eq.12) go to 10
      param(13) = loc(arg13)
      if (nargs.eq.13) go to 10
      param(14) = loc(arg14)
      if (nargs.eq.14) go to 10
      param(15) = loc(arg15)
      if (nargs.eq.15) go to 10
      param(16) = loc(arg16)
      if (nargs.eq.16) go to 10
      param(17) = loc(arg17)
      if (nargs.eq.17) go to 10
      param(18) = loc(arg18)
      if (nargs.eq.18) go to 10
      param(19) = loc(arg19)
      if (nargs.eq.19) go to 10
      param(20) = loc(arg20)
      if (nargs.eq.20) go to 10
      param(21) = loc(arg21)
      if (nargs.eq.21) go to 10
      param(22) = loc(arg22)
      if (nargs.eq.22) go to 10
      param(23) = loc(arg23)
      if (nargs.eq.23) go to 10
      param(24) = loc(arg24)
c Start task
   10 call MP_CR8TASK(tproc,0,proc,param,nargs,taskid)
      if (taskid.eq.0) then
         write (2,*) 'MP_TASKSTART failed'
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MP_TASKINIT(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,
     1arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
     2arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24)
c create a reusable task by packing arguments into a single structure
c taskid = index to notify queue for task (0 if error)
c proc = pointer to actual procedure
c nargs = number of arguments in actual procedure
c arg1-arg24 = arguments for actual procedure proc
c restrictions: actual arguments should not be temporary variables, such
c as return values of functions.  Also, local variables in procedure
c which called MP_TASKINIT should not be used if MP_SNDSIG and
c MP_WAITSIG are not called in the same procedure.
c character variables should also be avoided.
      implicit none
      integer taskid, nargs
      external proc
      integer arg1(*), arg2(*), arg3(*), arg4(*), arg5(*), arg6(*)
      integer arg7(*), arg8(*), arg9(*), arg10(*), arg11(*), arg12(*)
      integer arg13(*), arg14(*), arg15(*), arg16(*), arg17(*)
      integer arg18(*), arg19(*), arg20(*), arg21(*), arg22(*)
      integer arg23(*), arg24(*)
c local data
      integer MAXARGS
c MAXARGS = maximum number of arguments in task procedure
      parameter(MAXARGS=24)
      integer param
      dimension param(MAXARGS)
      external tsproc
c Create argument list for task
      if (nargs.eq.0) go to 10
      param(1) = loc(arg1)
      if (nargs.eq.1) go to 10
      param(2) = loc(arg2)
      if (nargs.eq.2) go to 10
      param(3) = loc(arg3)
      if (nargs.eq.3) go to 10
      param(4) = loc(arg4)
      if (nargs.eq.4) go to 10
      param(5) = loc(arg5)
      if (nargs.eq.5) go to 10
      param(6) = loc(arg6)
      if (nargs.eq.6) go to 10
      param(7) = loc(arg7)
      if (nargs.eq.7) go to 10
      param(8) = loc(arg8)
      if (nargs.eq.8) go to 10
      param(9) = loc(arg9)
      if (nargs.eq.9) go to 10
      param(10) = loc(arg10)
      if (nargs.eq.10) go to 10
      param(11) = loc(arg11)
      if (nargs.eq.11) go to 10
      param(12) = loc(arg12)
      if (nargs.eq.12) go to 10
      param(13) = loc(arg13)
      if (nargs.eq.13) go to 10
      param(14) = loc(arg14)
      if (nargs.eq.14) go to 10
      param(15) = loc(arg15)
      if (nargs.eq.15) go to 10
      param(16) = loc(arg16)
      if (nargs.eq.16) go to 10
      param(17) = loc(arg17)
      if (nargs.eq.17) go to 10
      param(18) = loc(arg18)
      if (nargs.eq.18) go to 10
      param(19) = loc(arg19)
      if (nargs.eq.19) go to 10
      param(20) = loc(arg20)
      if (nargs.eq.20) go to 10
      param(21) = loc(arg21)
      if (nargs.eq.21) go to 10
      param(22) = loc(arg22)
      if (nargs.eq.22) go to 10
      param(23) = loc(arg23)
      if (nargs.eq.23) go to 10
      param(24) = loc(arg24)
c Start task
   10 call MP_CR8TASK(tsproc,1,proc,param,nargs,taskid)
      if (taskid.eq.0) then
         write (2,*) 'MP_TASKINIT failed'
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MP_SETSTACK(stackval)
c set current stacksize (in bytes) for tasks
c input: stackval
      implicit none
      integer stackval
c common block for tasks
c MAXARGS = maximum number of arguments in task procedure
c MAXTASKS = maximum number of tasks supported
c stacksize = size of stack for task in bytes (0 = default)
      integer MAXARGS, MAXTASKS
      parameter(MAXARGS=24,MAXTASKS=16)
      integer notifyq, ntask, itask, params, mqueues, stacksize
      dimension ntask(MAXTASKS), itask(MAXTASKS)
      dimension params(MAXARGS+4,MAXTASKS), mqueues(2,MAXTASKS)
      common /tasks/ notifyq, ntask, itask, params, mqueues, stacksize
      if (stackval.ge.0) stacksize = stackval
      return
      end
c-----------------------------------------------------------------------
      subroutine MP_TASKBUILD(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5
     1,arg6,arg7,arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16
     2,arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24)
c create a non-running task by packing arguments into a single structure
c generally used for debugging in conjunction with MP_RUNTASK
c taskid = index to notify queue for task (0 if error)
c proc = pointer to actual procedure
c nargs = number of arguments in actual procedure
c arg1-arg24 = arguments for actual procedure proc
c restrictions: actual arguments should not be temporary variables, such
c as return values of functions.  Also, local variables in procedure
c which called MP_TASKBUILD should not be used if MP_RUNTASK is not
c called in the same procedure.
c character variables should also be avoided.
      implicit none
      integer taskid, nargs
      external proc
      integer arg1(*), arg2(*), arg3(*), arg4(*), arg5(*), arg6(*)
      integer arg7(*), arg8(*), arg9(*), arg10(*), arg11(*), arg12(*)
      integer arg13(*), arg14(*), arg15(*), arg16(*), arg17(*)
      integer arg18(*), arg19(*), arg20(*), arg21(*), arg22(*)
      integer arg23(*), arg24(*)
c local data
      integer MAXARGS
c MAXARGS = maximum number of arguments in task procedure
      parameter(MAXARGS=24)
      integer param
      dimension param(MAXARGS)
      external tproc
c Create argument list for task
      if (nargs.eq.0) go to 10
      param(1) = loc(arg1)
      if (nargs.eq.1) go to 10
      param(2) = loc(arg2)
      if (nargs.eq.2) go to 10
      param(3) = loc(arg3)
      if (nargs.eq.3) go to 10
      param(4) = loc(arg4)
      if (nargs.eq.4) go to 10
      param(5) = loc(arg5)
      if (nargs.eq.5) go to 10
      param(6) = loc(arg6)
      if (nargs.eq.6) go to 10
      param(7) = loc(arg7)
      if (nargs.eq.7) go to 10
      param(8) = loc(arg8)
      if (nargs.eq.8) go to 10
      param(9) = loc(arg9)
      if (nargs.eq.9) go to 10
      param(10) = loc(arg10)
      if (nargs.eq.10) go to 10
      param(11) = loc(arg11)
      if (nargs.eq.11) go to 10
      param(12) = loc(arg12)
      if (nargs.eq.12) go to 10
      param(13) = loc(arg13)
      if (nargs.eq.13) go to 10
      param(14) = loc(arg14)
      if (nargs.eq.14) go to 10
      param(15) = loc(arg15)
      if (nargs.eq.15) go to 10
      param(16) = loc(arg16)
      if (nargs.eq.16) go to 10
      param(17) = loc(arg17)
      if (nargs.eq.17) go to 10
      param(18) = loc(arg18)
      if (nargs.eq.18) go to 10
      param(19) = loc(arg19)
      if (nargs.eq.19) go to 10
      param(20) = loc(arg20)
      if (nargs.eq.20) go to 10
      param(21) = loc(arg21)
      if (nargs.eq.21) go to 10
      param(22) = loc(arg22)
      if (nargs.eq.22) go to 10
      param(23) = loc(arg23)
      if (nargs.eq.23) go to 10
      param(24) = loc(arg24)
c Start task
   10 call MP_CR8TASK(tproc,2,proc,param,nargs,taskid)
      if (taskid.eq.0) then
         write (2,*) 'MP_TASKBUILD failed'
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine MP_RUNTASK(taskid)
c executes a non-running task manually
c generally used for debugging in conjunction with MP_TASKBUILD
c taskid = index to task record (0 if task completed successfully)
c input: taskid
      implicit none
      integer taskid
c common block for tasks
c MAXARGS = maximum number of arguments in task procedure
c MAXTASKS = maximum number of tasks supported
c params = arguments in task procedures
      integer MAXARGS, MAXTASKS
      parameter(MAXARGS=24,MAXTASKS=16)
      integer notifyq, ntask, itask, params, mqueues, stacksize
      dimension ntask(MAXTASKS), itask(MAXTASKS)
      dimension params(MAXARGS+4,MAXTASKS), mqueues(2,MAXTASKS)
      common /tasks/ notifyq, ntask, itask, params, mqueues, stacksize
c local data
      integer i, ierr, tproc
      external tproc
      i = taskid
c Check for errors
      if ((i.lt.1).or.(i.gt.MAXTASKS)) then
         write (2,*) 'MP_RUNTASK: Invalid taskid = ', taskid
         return
      endif
      ierr = tproc(params(1,i))
      ntask(i) = 0
      taskid = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine MP_INITIALIZED(flag)
c indicate whether MP_INIT has been called
c flag = 1 if MP_INIT has been called, 0 otherwise
c output: flag
      implicit none
      integer flag
c common block for tasks
c MAXARGS = maximum number of arguments in task procedure
c MAXTASKS = maximum number of tasks supported
c notifyq = notifyqueue record for tasks
      integer MAXARGS, MAXTASKS
      parameter(MAXARGS=24,MAXTASKS=16)
      integer notifyq, ntask, itask, params, mqueues, stacksize
      dimension ntask(MAXTASKS), itask(MAXTASKS)
      dimension params(MAXARGS+4,MAXTASKS), mqueues(2,MAXTASKS)
      common /tasks/ notifyq, ntask, itask, params, mqueues, stacksize
      if (notifyq.ne.0) then
         flag = 1
      else
         flag = 0
      endif
      return
      end
c-----------------------------------------------------------------------
c Internal functions for Multitasking library
c-----------------------------------------------------------------------
      integer function exfunc(mpinit)
c this utility executes a dummy Boolean procedure
      implicit none
      integer*1 mpinit
      external mpinit
      exfunc = mpinit()
      return
      end
c-----------------------------------------------------------------------
      integer function tproc(parameter)
c generic task procedure
      implicit none
      integer parameter(*)
c local data
      integer nargs
      nargs = parameter(1)
      call mproc(val4(parameter(nargs+2)),nargs,parameter(2))
      tproc = 0
      return
      end
c-----------------------------------------------------------------------
      integer function tsproc(parameter)
c generic task procedure which waits on semaphores
      implicit none
      integer parameter(*)
c function declarations
      integer MPWaitOnSemaphore, MPSignalSemaphore
      external MPWaitOnSemaphore, MPSignalSemaphore
c local data
      integer kDurationForever
      parameter(kDurationForever=2147483647)
      integer nargs, oss
      nargs = parameter(1)
c Wait on semaphore
   10 oss = MPWaitOnSemaphore(val4(parameter(nargs+3)),
     1val4(kDurationForever))
c Signal arrived
      if (oss.eq.0) then
         call mproc(val4(parameter(nargs+2)),nargs,parameter(2))
c Signal a semaphore
         oss = MPSignalSemaphore(val4(parameter(nargs+4)))
c Wait for next signal
         if (oss.eq.0) go to 10
      endif
      tsproc = oss
      return
      end
c-----------------------------------------------------------------------
      subroutine mproc(proc,nargs,parameter)
c task subroutine with multiple possible arguments
c proc = pointer to actual procedure
c nargs = number of arguments in task procedure
c parameter = pointer to task procedure argument
      implicit none
      external proc
      integer nargs, parameter(*)
      select case (nargs)
      case (0)
         call proc()
      case (1)
         call proc(val4(parameter(1)))
      case (2)
         call proc(val4(parameter(1)),val4(parameter(2)))
      case (3)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)))
      case (4)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)))
      case (5)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)),val4(parameter(5)))
      case (6)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)),val4(parameter(5)),
     2val4(parameter(6)))
      case (7)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)),val4(parameter(5)),
     2val4(parameter(6)),val4(parameter(7)))
      case (8)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)),val4(parameter(5)),
     2val4(parameter(6)),val4(parameter(7)),val4(parameter(8)))
      case (9)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)),val4(parameter(5)),
     2val4(parameter(6)),val4(parameter(7)),val4(parameter(8)),
     3val4(parameter(9)))
      case (10)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)),val4(parameter(5)),
     2val4(parameter(6)),val4(parameter(7)),val4(parameter(8)),
     3val4(parameter(9)),val4(parameter(10)))
      case (11)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)),val4(parameter(5)),
     2val4(parameter(6)),val4(parameter(7)),val4(parameter(8)),
     3val4(parameter(9)),val4(parameter(10)),val4(parameter(11)))
      case (12)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)),val4(parameter(5)),
     2val4(parameter(6)),val4(parameter(7)),val4(parameter(8)),
     3val4(parameter(9)),val4(parameter(10)),val4(parameter(11)),
     4val4(parameter(12)))
      case (13)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)),val4(parameter(5)),
     2val4(parameter(6)),val4(parameter(7)),val4(parameter(8)),
     3val4(parameter(9)),val4(parameter(10)),val4(parameter(11)),
     4val4(parameter(12)),val4(parameter(13)))
      case (14)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)),val4(parameter(5)),
     2val4(parameter(6)),val4(parameter(7)),val4(parameter(8)),
     3val4(parameter(9)),val4(parameter(10)),val4(parameter(11)),
     4val4(parameter(12)),val4(parameter(13)),val4(parameter(14)))
      case (15)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)),val4(parameter(5)),
     2val4(parameter(6)),val4(parameter(7)),val4(parameter(8)),
     3val4(parameter(9)),val4(parameter(10)),val4(parameter(11)),
     4val4(parameter(12)),val4(parameter(13)),val4(parameter(14)),
     5val4(parameter(15)))
      case (16)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)),val4(parameter(5)),
     2val4(parameter(6)),val4(parameter(7)),val4(parameter(8)),
     3val4(parameter(9)),val4(parameter(10)),val4(parameter(11)),
     4val4(parameter(12)),val4(parameter(13)),val4(parameter(14)),
     5val4(parameter(15)),val4(parameter(16)))
      case (17)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)),val4(parameter(5)),
     2val4(parameter(6)),val4(parameter(7)),val4(parameter(8)),
     3val4(parameter(9)),val4(parameter(10)),val4(parameter(11)),
     4val4(parameter(12)),val4(parameter(13)),val4(parameter(14)),
     5val4(parameter(15)),val4(parameter(16)),val4(parameter(17)))
      case (18)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)),val4(parameter(5)),
     2val4(parameter(6)),val4(parameter(7)),val4(parameter(8)),
     3val4(parameter(9)),val4(parameter(10)),val4(parameter(11)),
     4val4(parameter(12)),val4(parameter(13)),val4(parameter(14)),
     5val4(parameter(15)),val4(parameter(16)),val4(parameter(17)),
     6val4(parameter(18)))
      case (19)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)),val4(parameter(5)),
     2val4(parameter(6)),val4(parameter(7)),val4(parameter(8)),
     3val4(parameter(9)),val4(parameter(10)),val4(parameter(11)),
     4val4(parameter(12)),val4(parameter(13)),val4(parameter(14)),
     5val4(parameter(15)),val4(parameter(16)),val4(parameter(17)),
     6val4(parameter(18)),val4(parameter(19)))
      case (20)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)),val4(parameter(5)),
     2val4(parameter(6)),val4(parameter(7)),val4(parameter(8)),
     3val4(parameter(9)),val4(parameter(10)),val4(parameter(11)),
     4val4(parameter(12)),val4(parameter(13)),val4(parameter(14)),
     5val4(parameter(15)),val4(parameter(16)),val4(parameter(17)),
     6val4(parameter(18)),val4(parameter(19)),val4(parameter(20)))
      case (21)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)),val4(parameter(5)),
     2val4(parameter(6)),val4(parameter(7)),val4(parameter(8)),
     3val4(parameter(9)),val4(parameter(10)),val4(parameter(11)),
     4val4(parameter(12)),val4(parameter(13)),val4(parameter(14)),
     5val4(parameter(15)),val4(parameter(16)),val4(parameter(17)),
     6val4(parameter(18)),val4(parameter(19)),val4(parameter(20)),
     7val4(parameter(21)))
      case (22)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)),val4(parameter(5)),
     2val4(parameter(6)),val4(parameter(7)),val4(parameter(8)),
     3val4(parameter(9)),val4(parameter(10)),val4(parameter(11)),
     4val4(parameter(12)),val4(parameter(13)),val4(parameter(14)),
     5val4(parameter(15)),val4(parameter(16)),val4(parameter(17)),
     6val4(parameter(18)),val4(parameter(19)),val4(parameter(20)),
     7val4(parameter(21)),val4(parameter(22)))
      case (23)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)),val4(parameter(5)),
     2val4(parameter(6)),val4(parameter(7)),val4(parameter(8)),
     3val4(parameter(9)),val4(parameter(10)),val4(parameter(11)),
     4val4(parameter(12)),val4(parameter(13)),val4(parameter(14)),
     5val4(parameter(15)),val4(parameter(16)),val4(parameter(17)),
     6val4(parameter(18)),val4(parameter(19)),val4(parameter(20)),
     7val4(parameter(21)),val4(parameter(22)),val4(parameter(23)))
      case (24)
         call proc(val4(parameter(1)),val4(parameter(2)),
     1val4(parameter(3)),val4(parameter(4)),val4(parameter(5)),
     2val4(parameter(6)),val4(parameter(7)),val4(parameter(8)),
     3val4(parameter(9)),val4(parameter(10)),val4(parameter(11)),
     4val4(parameter(12)),val4(parameter(13)),val4(parameter(14)),
     5val4(parameter(15)),val4(parameter(16)),val4(parameter(17)),
     6val4(parameter(18)),val4(parameter(19)),val4(parameter(20)),
     7val4(parameter(21)),val4(parameter(22)),val4(parameter(23)),
     8val4(parameter(24)))
      end select
      return
      end
c-----------------------------------------------------------------------
      subroutine prparms(taskid)
c debugging subroutine for printing task arguments
c used for debugging in conjunction with MP_TASKBUILD
c taskid = index to task record (0 if task completed successfully)
c input: taskid
      implicit none
      integer taskid
c common block for tasks
c MAXARGS = maximum number of arguments in task procedure
c MAXTASKS = maximum number of tasks supported
c params = arguments in task procedures
      integer MAXARGS, MAXTASKS
      parameter(MAXARGS=24,MAXTASKS=16)
      integer notifyq, ntask, itask, params, mqueues, stacksize
      dimension ntask(MAXTASKS), itask(MAXTASKS)
      dimension params(MAXARGS+4,MAXTASKS), mqueues(2,MAXTASKS)
      common /tasks/ notifyq, ntask, itask, params, mqueues, stacksize
c local data
      integer i, n, nargs, larg, iarg(2)
      real arg
      double precision darg
      equivalence(darg,iarg)
      i = taskid
c Check for errors
      if ((i.lt.1).or.(i.gt.MAXTASKS)) then
         return
      endif
c extract number of arguments
      nargs =  params(1,i)
      write (2,*) 'taskid, nargs = ', taskid, nargs
c write location and values of arguments (as integer, real and double)
      do 10 n = 1, nargs
      larg = params(n+1,i)
      iarg(1) = long(larg)
      iarg(2) = long(larg+4)
      arg = long(larg)
      write (2,*) 'n, loc, int, real, double = ',n,larg,iarg(1),arg,darg
   10 continue
      return
      end
