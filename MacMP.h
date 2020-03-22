/* header file for MacMP.c   */

void MP_Init(long *nproc);
	       
void MP_Taskwait(long *taskid);

int MP_Sndsig(long *taskid);

int MP_Waitsig(long *taskid);

void MP_Killtask(long *taskid);

void MP_End();

void MP_Taskstart(long *taskid, void (*proc)(), long *nargs, ...);

void MP_Taskinit(long *taskid, void (*proc)(), long *nargs, ...);

void MP_Setstack(long stackval);

void MP_Taskbuild(long *taskid, void (*proc)(), long *nargs, ...);

void MP_Runtask(long *taskid);

void MP_Initialized(long *flag);

void prparms(long taskid);
