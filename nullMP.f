c-----------------------------------------------------------------------
c Null Multitasking library for Apple Macintosh
c includes a few null routines from MacMPI
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
c-----------------------------------------------------------------------
      subroutine MP_INIT(nproc)
c initialize multitasking environment
c nproc = number of processors, 0 if Multitasking is not available
c output: nproc
      implicit none
      integer nproc
      nproc = 0
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
      taskid = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine MP_TASKWAIT(taskid)
c wait for task to complete
c taskid = index to task record (0 if task completed successfully)
c input and output: taskid
      implicit none
      integer taskid
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
      MP_SNDSIG = 0
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
      MP_WAITSIG = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine MP_KILLTASK(taskid)
c terminate multitasking environment
c taskid = index to task record (0 if task killed successfully)
      implicit none
      integer taskid
      taskid = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine MP_END()
c terminate multitasking environment
      implicit none
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
      taskid = 0
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
      taskid = 0
      return
      end
c-----------------------------------------------------------------------
      subroutine MP_SETSTACK(stackval)
c set current stacksize (in bytes) for tasks
c input: stackval
      implicit none
      integer stackval
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
      taskid = 0
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
      flag = 0
      return
      end
