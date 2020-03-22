/* Multitasking library for Apple Macintosh
   based on the Multiprocessing Library in the MacOS, v. 2.1
   The file MPerrs is used throughout for error messages
   written by viktor k. decyk, ucla
   copyright 2000, regents of the university of california.
   all rights reserved.
   no warranty for proper operation of this software is given or implied.
   software or information may be copied, distributed, and used at own
   risk; it may not be distributed without this notice included verbatim
   with each file. 
   update: june 14, 2003                                              */

#include <stdarg.h>
#include <stdio.h>
#include "MacMP.h"

#include <Multiprocessing.h>

/* MAXARGS = maximum number of arguments in task procedure */
#define MAXARGS                 24
/* MAXTASKS = maximum number of tasks supported */
#define MAXTASKS                16

/* internal common block for tasks
   notifyq = notifyqueue record for tasks
   ntask = task record
   itask = MPTaskID record
   params = arguments in task procedures
   mqueues = message queues for tasks (1=request,2=result)
   stackSize = size of stack for task in bytes (0 = default)  */

   static MPQueueID notifyq = 0;
   static int ntask[MAXTASKS];
   static MPTaskID itask[MAXTASKS];
   static void *params[MAXTASKS][MAXARGS+4];
   static MPSemaphoreID mqueues[MAXTASKS][2];
   static ByteCount stackSize = 0;

static FILE *unit2;

/* prototypes for internal procedures */

void MP_Cr8task(TaskProc entryPoint, int *ttype, void (*proc)(),
               void *parameter, int *nargs, int *taskid);

OSStatus tproc(void *parameter);

OSStatus tsproc(void *parameter);

void  mproc(void (*proc)(), int nargs, void *parameter);

void MP_Init(int *nproc) {
/* initialize multitasking environment
   nproc = number of processors, 0 if Multitasking is not available
   output: nproc
local data                             */
   int i;
   long oss;
/* internal common block for tasks
   notifyq = notifyqueue record for tasks
   ntask = task record
   itask = MPTaskID record
   mqueues = message queues for tasks (1=request,2=result)
   stackSize = size of stack for task in bytes (0 = default)  */
   *nproc = 0;
   notifyq = 0;
/* Determine if MP API Services are available */
   *nproc = MPLibraryIsLoaded();
/* Exit if MP API Services are not available */
   if (!(*nproc))
      return;
/* Open error file */
   unit2 = fopen("MPerrs","w");
/* Return the number of processors on the host computer */
   *nproc = MPProcessors();
/* Create a message queue */
   oss = MPCreateQueue(&notifyq);
   if (oss) {
      fprintf(unit2,"MPCreateQueue error, oss = %d\n",oss);
      *nproc = 0;
      notifyq = 0;
      return;
   }
/* Clear task records */
   for (i = 0; i < MAXTASKS; i++) {
      ntask[i] = 0;
      itask[i] = 0;
      mqueues[i][0] = 0;
      mqueues[i][1] = 0;
   }
/* set current stacksize (in bytes) for tasks */
   stackSize = 16384;
   return;
}

void MP_Cr8task(TaskProc entryPoint, int *ttype, void (*proc)(),
               void *parameter, int *nargs, int *taskid) {
/* start task
   overhead is about 130 microseconds on Macintosh G4/450
   entrypoint = pointer to task procedure
   ttype = task type (0=non-reusable,1=reusable(signaling),2=non-running)
   proc = pointer to actual procedure
   parameter = pointer to task procedure argument
   nargs = number of arguments in task procedure
   taskid = index to task record (0 if error)
   input: all except taskid
   output: taskid
local data                                              */
   void **param;
   int i, j;
   long oss;
   static int tp = 0;
   MPTaskOptions options = 0;
/* internal common block for tasks
   notifyq = notifyqueue record for tasks
   ntask = task record
   itask = MPTaskID record
   params = arguments in task procedures
   mqueues = message queues for tasks (1=request,2=result)
   stackSize = size of stack for task in bytes (0 = default)  */
   param = (void **)parameter;
   *taskid = 0;
/* Check for errors */
   if ((!notifyq) && (*ttype < 2)) {
      fprintf(unit2,"MP_Cr8task: MP not initialized\n");
      return;
   }
   else if (*nargs > MAXARGS) {
      fprintf(unit2,"Too many arguments in task, nargs = %d, MAXARGS = %d\n",
              *nargs,MAXARGS);
      return;
   }
/* Find space for record */
   i = -1;
L10: i += 1;
   if (i >= MAXTASKS) {
      fprintf(unit2,"Exceeded maximum number of tasks = %d\n",MAXTASKS);
      return;
   }
   else if (ntask[i])
      goto L10;
/* Create semaphores */
   if (*ttype==1) {
      for (j = 0; j < 2; j++) {
         oss = MPCreateSemaphore(1,0,&mqueues[i][j]);
         if (oss) {
            fprintf(unit2,"MPCreateSemaphore error, oss = %d\n",oss);
            mqueues[i][j] = 0;
            return;
         }
      }
   }
   ntask[i] = i + 1;
/* Copy number of arguments */
   params[i][0] = (void *)(*nargs);
/* Copy arguments */
   for (j = 0; j < *nargs; j++) {
      params[i][j+1] = param[j];
   }
/* Copy procedure name and message queue pointers */
   params[i][*nargs+1] = (void *)proc;
   params[i][*nargs+2] = (void *)mqueues[i][0];
   params[i][*nargs+3] = (void *)mqueues[i][1];
/* Create a preemptive task */
   if (*ttype < 2) {
      oss = MPCreateTask(entryPoint,&params[i][0],stackSize,notifyq,
            &ntask[i],&tp,options,&itask[i]);
      if (oss) {
         fprintf(unit2,"MPCreateTask error, oss = %d\n",oss);
         ntask[i] = 0;
         itask[i] = 0;
         return;
      }
   }
   *taskid = i + 1;
   return;
}

void MP_Taskwait(int *taskid) {
/* wait for task to complete
   taskid = index to task record (0 if task completed successfully)
   input and output: taskid
local data                                             */
   int i, j, id;
   long oss;
   int *message = 0, *tp = 0, *result = 0;
/* internal common block for tasks
   notifyq = notifyqueue record for tasks
   ntask = task record
   itask = MPTaskID record
   mqueues = message queues for tasks (1=request,2=result) */
   i = *taskid - 1;
/* Check for errors */
   if (!notifyq) {
      fprintf(unit2,"MP_Taskwait: MP not initialized\n");
      return;
   }
   else if ((i < 0) || (i >=MAXTASKS)) {
      fprintf(unit2,"MP_Taskwait: Invalid taskid = %d\n",*taskid);
      return;
   }
   else if (mqueues[i][1] && itask[i]) {
      fprintf(unit2,"MP_Taskwait: Invalid type for task = %d\n",*taskid);
      return;
   }
/* Look for message from taskid */
L10: if (ntask[i] > 0) {
/* Obtain a message from a specified message queue */
      oss = MPWaitOnQueue(notifyq,(void *)&message,(void *)&tp,
                         (void *)&result,kDurationForever);
      if (oss) {
         fprintf(unit2,"MPWaitOnQueue error for taskid, oss = %d,%d\n",
                 i+1,oss);
         return;
      }
      id = *message - 1;
      if ((id < 0) || (id >=MAXTASKS))
         fprintf(unit2,"Invalid taskid from notifyq = %d\n",id+1);
      else if (ntask[id]==0)
         fprintf(unit2,"task from notifyq already ended, taskid = %d\n",
                 id+1);
      else
/* Mark task as completed, may not be the one we are waiting for */
        ntask[id] = -(id + 1);
      goto L10;
   }
/* Task marked as completed earlier */
   else if (ntask[i] < 0)
      ntask[i] = 0;
/* Remove semaphores */
   for (j = 0; j < 2; j++) {
      if (mqueues[i][j]) {
         oss = MPDeleteSemaphore(mqueues[i][j]);
         if (oss)
            fprintf(unit2,"MPDeleteSemaphore error, oss = %d\n",oss);
         mqueues[i][j] = 0;
      }
   }
   itask[i] = 0;
   *taskid = 0;
   return;
}

int MP_Sndsig(int *taskid) {
/* send message to task
   overhead for sending a message to a task and receiving a reply 
   is about 60 microseconds on Macintosh G4/450
   taskid = index to task record
   input: all
local data                                          */
   int i;
   long oss;
/* internal common block for tasks
   notifyq = notifyqueue record for tasks
   ntask = task record
   mqueues = message queues for tasks (1=request,2=result) */
   i = *taskid - 1;
/* Check for errors */
   if (!notifyq) {
      fprintf(unit2,"MP_Sndsig: MP not initialized\n");
      return 1;
   }
   else if ((i < 0) || (i >=MAXTASKS)) {
      fprintf(unit2,"MP_Sndsig: Invalid taskid = %d\n",*taskid);
      return 2;
   }
   else if (!(ntask[i])) {
      fprintf(unit2,"MP_Sndsig: Task already ended, taskid = %d\n",
              *taskid);
      return 3;
   }
   else if (!mqueues[i][0]) {
      fprintf(unit2,"MP_Sndsig: Invalid type for task = %d\n",*taskid);
      return 4;
   }
/* Signal a semaphore */
   oss = MPSignalSemaphore(mqueues[i][0]);
   if (oss) {
      fprintf(unit2,"MPSignalSemaphore error for taskid, oss = %d,%d\n",
              i+1,oss);
      return 5;
   }
   return 0;
}

int MP_Waitsig(int *taskid) {
/* receive signal from task
   overhead for sending a message to a task and receiving a reply 
   is about 60 microseconds on Macintosh G4/450
   taskid = index to task record
   input: taskid
   input: all
local data                                          */
   int i;
   long oss;
/* internal common block for tasks
   notifyq = notifyqueue record for tasks
   ntask = task record
   mqueues = message queues for tasks (1=request,2=result) */
   i = *taskid - 1;
/* Check for errors */
   if (!notifyq) {
      fprintf(unit2,"MP_Waitsig: MP not initialized\n");
      return 1;
   }
   else if ((i < 0) || (i >=MAXTASKS)) {
      fprintf(unit2,"MP_Waitsig: Invalid taskid = %d\n",*taskid);
      return 2;
   }
   else if (!(ntask[i])) {
      fprintf(unit2,"MP_Waitsig: Task already ended, taskid = %d\n",
              *taskid);
      return 3;
   }
   else if (!mqueues[i][1]) {
      fprintf(unit2,"MP_Waitsig: Invalid type for task = %d\n",*taskid);
      return 4;
   }
/* Wait on semaphore */
   oss = MPWaitOnSemaphore(mqueues[i][1],kDurationForever);
   if (oss) {
      fprintf(unit2,"MPWaitOnSemaphore for taskid, error, oss = %d\n",
              i+1,oss);
      return 5;
   }
   return 0;
}

void MP_Killtask(int *taskid) {
/* terminate a task
   taskid = index to task record (0 if task killed successfully)
local data                             */
   int i, result = 1;
   long oss;
/* internal common block for tasks
   notifyq = notifyqueue record for tasks
   itask = MPTaskID record                */
   i = *taskid - 1;
/* Check for errors */
   if (!notifyq) {
      fprintf(unit2,"MP_Killtask: MP not initialized\n");
      return;
   }
   else if ((i < 0) || (i >=MAXTASKS)) {
      fprintf(unit2,"MP_Killtask: Invalid taskid = %d\n",*taskid);
      return;
   }
   else if (!(ntask[i])) {
      fprintf(unit2,"Task already killed, taskid = %d\n",*taskid);
      return;
   }
/* Terminate an existing task */
   oss = MPTerminateTask(itask[i],result);
   if (oss)
      fprintf(unit2,"MPTerminateTask error for taskid, oss = %d,%d\n",
              i+1,oss);
   else {
      itask[i] = 0;
      MP_Taskwait(taskid);
   }
   return;
}

void MP_End() {
/* terminate multitasking environment
local data                             */
   int i, taskid;
   long oss;
/* internal common block for tasks
   notifyq = notifyqueue record for tasks */
/* Check for errors */
   if (!notifyq) {
      fprintf(unit2,"MP_End: MP not initialized\n");
      return;
   }
/* Check if any tasks are outstanding */
   for (i = 0; i < MAXTASKS; i++) {
      if (ntask[i] > 0) {
         taskid = i + 1;
         fprintf(unit2,"MP_End: Task still outstanding, taskid=%d\n",i+1);
         MP_Killtask(&taskid);
      }
   }
/* Delete a message queue */
   oss = MPDeleteQueue(notifyq);
   if (oss)
      fprintf(unit2,"MPDeleteQueue error, oss = %d\n",oss);
/* Delete file if empty */
   if (!fseek(unit2,0,SEEK_END)) {
      oss = ftell(unit2);
      fclose(unit2);
      if (!oss)
         remove("MPerrs");
   }
   notifyq = 0;
   return;
}

void MP_Taskstart(int *taskid, void (*proc)(), int *nargs, ...) {
/* create a task by packing arguments into a single structure
   taskid = index to notify queue for task (0 if error)
   proc = pointer to actual procedure
   nargs = number of arguments in actual procedure
   restrictions: actual arguments should be pointers and not temporary
   variables, such as return values of functions.  Also, local variables
   in procedure which called MP_Taskstart should not be used if
   MP_Taskwait is not called in the same procedure.
local data                                                           */
   void *param[MAXARGS];
   int i, ttype = 0;
   va_list argptr;
   va_start(argptr,nargs);
/* Create argument list for task */
   for (i = 0; i < *nargs; i++) {
      param[i] = va_arg(argptr,void*);
   }
   va_end(argptr);
/* Start task */
   MP_Cr8task(&tproc,&ttype,proc,&param[0],nargs,taskid);
   if (!(*taskid))
      fprintf(unit2,"MP_Taskstart failed\n");
   return;
}

void MP_Taskinit(int *taskid, void (*proc)(), int *nargs, ...) {
/* create a reusable task by packing arguments into a single structure
   taskid = index to notify queue for task (0 if error)
   proc = pointer to actual procedure
   nargs = number of arguments in actual procedure
   restrictions: actual arguments should be pointers and not temporary
   variables, such as return values of functions.  Also, local variables
   in procedure which called MP_Taskinit should not be used if
   MP_Sndsig and MP_Waitsig are not called in the same procedure.
local data                                                           */
   void *param[MAXARGS];
   int i, ttype = 1;
   va_list argptr;
   va_start(argptr,nargs);
/* Create argument list for task */
   for (i = 0; i < *nargs; i++) {
      param[i] = va_arg(argptr,void*);
   }
   va_end(argptr);
/* Start task */
   MP_Cr8task(&tsproc,&ttype,proc,&param[0],nargs,taskid);
   if (!(*taskid))
      fprintf(unit2,"MP_Taskinit failed\n");
   return;
}

void MP_Setstack(int stackval) {
/* set current stacksize (in bytes) for tasks 
   input: stackval                            */
/* internal common block for tasks
   stackSize = size of stack for task in bytes (0 = default) */
   if (stackval >= 0) 
      stackSize = stackval;
   return;
}

void MP_Taskbuild(int *taskid, void (*proc)(), int *nargs, ...) {
/* create a non-running task by packing arguments into a single structure
   generally used for debugging in conjunction with MP_Runtask
   taskid = index to notify queue for task (0 if error)
   proc = pointer to actual procedure
   nargs = number of arguments in actual procedure
   restrictions: actual arguments should be pointers and not temporary
   variables, such as return values of functions.  Also, local variables
   in procedure which called MP_Taskbuild should not be used if
   MP_Runtask is not called in the same procedure.
local data                                                           */
   void *param[MAXARGS];
   int i, ttype = 2;
   va_list argptr;
   va_start(argptr,nargs);
/* Create argument list for task */
   for (i = 0; i < *nargs; i++) {
      param[i] = va_arg(argptr,void*);
   }
   va_end(argptr);
/* Start task */
   MP_Cr8task(&tproc,&ttype,proc,&param[0],nargs,taskid);
   if (!(*taskid))
      fprintf(unit2,"MP_Taskbuild failed\n");
   return;
}

void MP_Runtask(int *taskid) {
/* executes a non-running task manually
   generally used for debugging in conjunction with MP_Taskbuild
   taskid = index to task record (0 if task completed successfully)
   input: taskid
local data  */
   int i;
   OSStatus ierr;
/* internal common block for tasks
   params = arguments in task procedures */
   i = *taskid - 1;
/* Check for errors */
   if ((i < 0) || (i >=MAXTASKS)) {
      fprintf(unit2,"MP_Runtask: Invalid taskid = %d\n",*taskid);
      return;
   }
   ierr = tproc(&params[i][0]);
   ntask[i] = 0;
   *taskid = 0;
   return;
}

void MP_Initialized(int *flag) {
/* indicate whether MP_Init has been called
   flag = true if MP_Init has been called, false otherwise
   output: flag                                             */
/* internal common block for tasks
   notifyq = notifyqueue record for tasks                   */
   if (notifyq)
      *flag = 1;
   else
      *flag = 0;
   return;
}

/*--------------------------------------------------------------------*/
/* Internal functions for Multitasking library                        */
/*--------------------------------------------------------------------*/

void mproc(void (*proc)(), int nargs, void *parameter);

OSStatus tproc(void *parameter) {
/* generic task procedure
local data                */
   void **param;
   int nargs;
   param = (void **)parameter;
   nargs = (int )param[0];
   mproc(param[nargs+1],nargs,&param[1]);
   return 0;
}

OSStatus tsproc(void *parameter) {
/* generic task procedure which waits on semaphores
local data                */
   void **param;
   int nargs; 
   long oss;
   param = (void **)parameter;
   nargs = (int )param[0];
/* Wait on semaphore */
L10: oss = MPWaitOnSemaphore(param[nargs+2],kDurationForever);
/* Signal arrived */
   if (!oss) {
      mproc(param[nargs+1],nargs,&param[1]);
/* Signal a semaphore */
      oss = MPSignalSemaphore(param[nargs+3]);
/* Wait for next signal */
      if (!oss)
         goto L10;
   }
   return oss;
}

void mproc(void (*proc)(), int nargs, void *parameter) {
/* task subroutine with multiple possible arguments
   proc = pointer to actual procedure
   nargs = number of arguments in task procedure
   parameter = pointer to task procedure argument */
   void **param;
   param = (void **)parameter;
   switch (nargs) {
      case 0:
         proc();
         break;
      case 1:
         proc((void *)param[0]);
         break;
      case 2:
         proc((void *)param[0],(void *)param[1]);
         break;
      case 3:
         proc((void *)param[0],(void *)param[1],(void *)param[2]);
         break;
      case 4:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3]);
         break;
      case 5:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4]);
         break;
      case 6:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5]);
         break;
      case 7:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6]);
         break;
      case 8:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7]);
         break;
      case 9:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8]);
         break;
      case 10:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9]);
         break;
      case 11:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10]);
         break;
      case 12:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11]);
         break;
      case 13:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12]);
         break;
      case 14:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13]);
         break;
      case 15:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14]);
         break;
      case 16:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14],
              (void *)param[15]);
         break;
      case 17:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14],
              (void *)param[15],(void *)param[16]);
         break;
      case 18:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14],
              (void *)param[15],(void *)param[16],(void *)param[17]);
         break;
      case 19:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14],
              (void *)param[15],(void *)param[16],(void *)param[17],
              (void *)param[18]);
         break;
      case 20:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14],
              (void *)param[15],(void *)param[16],(void *)param[17],
              (void *)param[18],(void *)param[19]);
         break;
      case 21:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14],
              (void *)param[15],(void *)param[16],(void *)param[17],
              (void *)param[18],(void *)param[19],(void *)param[20]);
         break;
      case 22:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14],
              (void *)param[15],(void *)param[16],(void *)param[17],
              (void *)param[18],(void *)param[19],(void *)param[20],
              (void *)param[21]);
         break;
      case 23:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14],
              (void *)param[15],(void *)param[16],(void *)param[17],
              (void *)param[18],(void *)param[19],(void *)param[20],
              (void *)param[21],(void *)param[22]);
         break;
      case 24:
         proc((void *)param[0],(void *)param[1],(void *)param[2],
              (void *)param[3],(void *)param[4],(void *)param[5],
              (void *)param[6],(void *)param[7],(void *)param[8],
              (void *)param[9],(void *)param[10],(void *)param[11],
              (void *)param[12],(void *)param[13],(void *)param[14],
              (void *)param[15],(void *)param[16],(void *)param[17],
              (void *)param[18],(void *)param[19],(void *)param[20],
              (void *)param[21],(void *)param[22],(void *)param[23]);
   }
   return;
}

void prparms(int taskid) {
/* debugging subroutine for printing task arguments
   used for debugging in conjunction with MP_Taskbuild
   taskid = index to task record (0 if task completed successfully)
   input: taskid
local data    */
   int i, n, nargs, iarg;
   void *larg;
   float arg;
   double darg;
/* internal common block for tasks
   params = arguments in task procedures */
   i = taskid - 1;
/* Check for errors */
   if ((i < 0) || (i >=MAXTASKS)) {
      return;
   }
/* extract number of arguments */
   nargs = (int )params[i][0];
   fprintf(unit2,"taskid, nargs = %d,%d\n",taskid,nargs);
/* write location and values of arguments (as integer, real and double) */
   for (n = 0; n < nargs; n++) {
      larg = params[i][n+1];
      iarg = *(int *) larg;
      arg = *(float *) larg;
      darg = *(double *) larg;
      fprintf(unit2,"n, loc, int, float, double = %d,%p,%d,%f,%f\n",
	         n,larg,iarg,arg,darg);
   }
   return;
}
