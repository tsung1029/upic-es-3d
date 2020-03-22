/* This is a Fortran callable interface library to the C version of a
   Multitasking library for Apple Macintosh
   based on the Multiprocessing Library in the MacOS, v. 2.1
   The file MPerrs is used throughout for error messages
   written by viktor k. decyk, ucla
   copyright 2000, regents of the university of california.
   all rights reserved.
   no warranty for proper operation of this software is given or implied.
   software or information may be copied, distributed, and used at own
   risk; it may not be distributed without this notice included verbatim
   with each file. 
   update: october 25, 2004                                             */

#include "MacMP.h"

void mp_init(int *nproc) {
   MP_Init(nproc);
   return;
}

void mp_cr8task(long (*entryPoint)(void *parameter), int *ttype,
                void (*proc)(), void *parameter, int *nargs,
                int *taskid) {
   MP_Cr8task(entryPoint,ttype,proc,parameter,nargs,taskid);
   return;
}

void mp_taskwait(int *taskid) {
   MP_Taskwait(taskid);
   return;
}

int mp_sndsig(int *taskid) {
   return MP_Sndsig(taskid);
}

int mp_waitsig(int *taskid) {
   return MP_Waitsig(taskid);
}

void mp_killtask(int *taskid) {
   MP_Killtask(taskid);
   return;
}

void mp_end() {
   MP_End();
   return;
}

void mp_taskstart(int *taskid, void (*proc)(), int *nargs, int *arg1,
                  int *arg2, int *arg3, int *arg4, int *arg5, int *arg6,
                  int *arg7, int *arg8, int *arg9, int *arg10, int *arg11,
                  int *arg12, int *arg13, int *arg14, int *arg15,
                  int *arg16, int *arg17, int *arg18, int *arg19,
                  int *arg20, int *arg21, int *arg22, int *arg23,
                  int *arg24) {
   switch (*nargs) {
   case 0:
      MP_Taskstart(taskid,proc,nargs);
      return;
   case 1:
      MP_Taskstart(taskid,proc,nargs,arg1);   
      return;
   case 2:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2);
      return;
   case 3:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3);
      return;
   case 4:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4);
      return;
   case 5:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5);
      return;
   case 6:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6);
      return;
   case 7:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7);
      return;
   case 8:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8);
      return;
   case 9:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9);
      return;
   case 10:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10);
      return;
   case 11:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11);
      return;
   case 12:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12);
      return;
   case 13:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13);
      return;
   case 14:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14);
      return;
   case 15:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15);
      return;
   case 16:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16);
      return;
   case 17:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17);
      return;
   case 18:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18);
      return;
   case 19:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19);
      return;
   case 20:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20);
      return;
   case 21:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21);
      return;
   case 22:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22);
      return;
   case 23:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22,arg23);
      return;
   case 24:
      MP_Taskstart(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24);
      return;
   }
   return;
}

void mp_taskinit(int *taskid, void (*proc)(), int *nargs, int *arg1,
                  int *arg2, int *arg3, int *arg4, int *arg5, int *arg6,
                  int *arg7, int *arg8, int *arg9, int *arg10, int *arg11,
                  int *arg12, int *arg13, int *arg14, int *arg15,
                  int *arg16, int *arg17, int *arg18, int *arg19,
                  int *arg20, int *arg21, int *arg22, int *arg23,
                  int *arg24) {
   switch (*nargs) {
   case 0:
      MP_Taskinit(taskid,proc,nargs);
      return;
   case 1:
      MP_Taskinit(taskid,proc,nargs,arg1);   
      return;
   case 2:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2);
      return;
   case 3:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3);
      return;
   case 4:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4);
      return;
   case 5:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5);
      return;
   case 6:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6);
      return;
   case 7:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7);
      return;
   case 8:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8);
      return;
   case 9:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9);
      return;
   case 10:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10);
      return;
   case 11:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11);
      return;
   case 12:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12);
      return;
   case 13:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13);
      return;
   case 14:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14);
      return;
   case 15:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15);
      return;
   case 16:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16);
      return;
   case 17:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17);
      return;
   case 18:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18);
      return;
   case 19:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19);
      return;
   case 20:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20);
      return;
   case 21:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21);
      return;
   case 22:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22);
      return;
   case 23:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22,arg23);
      return;
   case 24:
      MP_Taskinit(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24);
      return;
   }
   return;
}

void mp_setstack(int *stackval) {
   int lstackval;
   lstackval = *stackval;
   MP_Setstack(lstackval);
   return;
}

void mp_taskbuild(int *taskid, void (*proc)(), int *nargs, int *arg1,
                  int *arg2, int *arg3, int *arg4, int *arg5, int *arg6,
                  int *arg7, int *arg8, int *arg9, int *arg10, int *arg11,
                  int *arg12, int *arg13, int *arg14, int *arg15,
                  int *arg16, int *arg17, int *arg18, int *arg19,
                  int *arg20, int *arg21, int *arg22, int *arg23,
                  int *arg24) {
   switch (*nargs) {
   case 0:
      MP_Taskbuild(taskid,proc,nargs);
      return;
   case 1:
      MP_Taskbuild(taskid,proc,nargs,arg1);   
      return;
   case 2:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2);
      return;
   case 3:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3);
      return;
   case 4:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4);
      return;
   case 5:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5);
      return;
   case 6:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6);
      return;
   case 7:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7);
      return;
   case 8:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8);
      return;
   case 9:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9);
      return;
   case 10:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10);
      return;
   case 11:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11);
      return;
   case 12:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12);
      return;
   case 13:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13);
      return;
   case 14:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14);
      return;
   case 15:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15);
      return;
   case 16:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16);
      return;
   case 17:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17);
      return;
   case 18:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18);
      return;
   case 19:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19);
      return;
   case 20:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20);
      return;
   case 21:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21);
      return;
   case 22:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22);
      return;
   case 23:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22,arg23);
      return;
   case 24:
      MP_Taskbuild(taskid,proc,nargs,arg1,arg2,arg3,arg4,arg5,arg6,arg7,
                   arg8,arg9,arg10,arg11,arg12,arg13,arg14,arg15,arg16,
                   arg17,arg18,arg19,arg20,arg21,arg22,arg23,arg24);
      return;
   }
   return;
}

void mp_runtask(int *taskid) {
   MP_Runtask(taskid);
   return;
}

void mp_initialized(int *flag) {
   MP_Initialized(flag);
   return;
}
