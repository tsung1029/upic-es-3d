/*

 $URL: svn+ssh://exppmaster/svn_repositories/osiris/trunk/source/os-sys-macosx-c.c $
 $Id: os-sys-macosx-c.c 55 2005-11-30 17:31:32Z zamb $

 Interfacing routines for POSIX system calls
 (should work under any POSIX compliant system)

*/

#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include <sys/resource.h>

#include <time.h>

/* Macros for fortran functions */
#include "fortran.h"

/* Force stdout and stderr to be line buffered */

void FNAME(linebuf_stdio_f)
(void)
{
  setvbuf(stdout, (char *)NULL, _IOLBF, 0);
  setvbuf(stderr, (char *)NULL, _IOLBF, 0);
}

/* Remove a file */

void FNAME(remove_f)
(const char *path, int *ierr)
{
 if ((*ierr = remove(path))) *ierr = errno;
}


/* Create a directory */

void FNAME(mkdir_f)
(const char *path, int *ierr)
{
 char uppath[256], *p;
    
 if (mkdir(path,S_IRWXU | (S_IRGRP | S_IXGRP ) | (S_IROTH | S_IXOTH) ))
   switch (errno) {
     case ENOENT : /* A component of the path does not exist */
        
        /* get upper path */
        
        strcpy(uppath, path);
        p = uppath + strlen(uppath);
        while(*p!='/') p--;
        *p=0;
        
        /* recursively build the path */
        
        FNAME(mkdir_f) (uppath, ierr);
        if (!*ierr) FNAME(mkdir_f) (path, ierr);
     
     case EEXIST : /* if directory already exists ignore the error */
        *ierr = 0;
        break;
     default: *ierr = errno;
   }
 else *ierr = 0;
 
}

/* change working a directory */

void FNAME(chdir_f)
(const char *path, int *ierr)
{
  if ((*ierr = chdir(path))) *ierr=errno;
}

/* get current working a directory */

void FNAME(getcwd_f)
(char *path, int *len, int *ierr)
{
  char *buf = (char*) malloc( (*len+1) );
  
  if (getcwd(buf, *len)) {
     strcpy(path, buf);
     *ierr = 0;
   } else {
    strcpy(path, "\0");
    *ierr=errno;
   }
  
  free(buf);
}

/* Create a symbolic link */

void FNAME(symlink_f)
(const char *name1, const char *name2, int *ierr)
{
  if ((*ierr = symlink(name1,name2))) *ierr=errno;
}

/* Get the error string */

void FNAME(strerror_f)
(const int *errnum, char *err)
{
 strcpy(err, strerror(*errnum));
}

/* Get the environment variable value */

void FNAME(getenv_f)
(const char *name, char *value, int *ierr)
{
 char *lvalue;
  
 lvalue = (char *) getenv(name);
  
 if (lvalue == NULL) *ierr = -1;
 else 
 {
   *ierr = 0;
   strcpy(value,lvalue);
 }
 
}

/* Get the hostname */

void FNAME(gethostname_f)
(char *hostname, int *ierr)
{
 char lhostname[256];
  
 if (gethostname(lhostname, 255)) {
   strcpy(hostname,"\0");
   *ierr = errno; 
 } else {
   strcpy(hostname,lhostname);
   *ierr = 0;
 }
}

/* Set the umask */
void FNAME(umask_f)
(int *numask)
{
  umask( *numask );
}

/* check if floating point is infinity */

void FNAME(isinf_f)
( double *x, int *res )
{
  *res = isinf( *x);
}

/* check if floating point is nan */

void FNAME(isnan_double_f)
( double *x, int *res )
{
  *res = isnan( *x);
}

void FNAME(isnan_single_f)
( float *x, int *res )
{
  *res = isnan( *x);
}


/*
void FNAME(getrusage_f)
( double *utime, double *stime, long *maxrss, long *ixrss, long *idrss, long *isrss )
{
   struct rusage selfUsage;
   
   getrusage( RUSAGE_SELF, &selfUsage );

   *utime = selfUsage.ru_utime.tv_sec + 1.0e-6*selfUsage.ru_utime.tv_usec; 
   *stime = selfUsage.ru_stime.tv_sec + 1.0e-6*selfUsage.ru_stime.tv_usec;
   
   *maxrss = selfUsage.ru_maxrss;
   *ixrss  = selfUsage.ru_ixrss;
   *idrss  = selfUsage.ru_idrss;
   *isrss  = selfUsage.ru_isrss;

}
*/

#include <unistd.h>

/* return maximum memory usage of process in kBytes */
void FNAME(used_memory_f)
( long *used_kb )
{
  struct rusage usage;
  if ( getrusage(RUSAGE_SELF, &usage) ) *used_kb = 0;
  else *used_kb = usage.ru_maxrss;

#ifdef __APPLE__
  /* OS X returns this as bytes, not kb */
  *used_kb /= 1024;
#endif

}

/******************************************************************************
		To get the Unix time.  Ian needs this for his automatic restarts to
		operate efficiently.
*******************************************************************************/
int64_t FNAME(unix_time)
(void)
{
  int64_t utime = 0;
	utime = (int64_t) time(NULL);
	return utime;
}


/****************************************************************************** 
   
   System dependent high resolution timers.
   
   Note: There is a significant overhead on converting the time from ticks
         to seconds, so for high precision measurements a tick based approach
         is preferrable.
         
         Resolution is defined in terms of the minimum time difference between
         calls to the system clock, without accounting the time for converting
         the ticks to seconds.
   
*******************************************************************************/

#ifdef __AIX_TIMER__
#define HAVE_TIMER

/* 
  Use high resolution timers from the Standard C Library (libc.a) on IBM AIX 
  systems running on PowerPC architecture. The resolution is system dependent
  and is on the order of tens of nanoseconds (~ 70 ns on a JS21 PPC5+ blade).    
*/

#include <sys/time.h>

int tb_flag;

void FNAME(timer_init)
(void)
{
  /* store real time flag for future use */
  
  timebasestruct_t tb;
  read_real_time(&tb, TIMEBASE_SZ);
  tb_flag = tb.flag;
}

uint64_t FNAME(timer_ticks)
(void)
{
  timebasestruct_t tb;
  read_real_time(&tb, TIMEBASE_SZ);
 
  return ((uint64_t)tb.tb_high)<<32 | (uint64_t)tb.tb_low;

}

double FNAME(timer_interval_seconds)
(uint64_t *start, uint64_t *end)
{
  int sec, nsec;
  timebasestruct_t tb_start, tb_end;

  /* Convert start/end values to real time */
  
  tb_start.tb_high = (unsigned int)((*start) >> 32);
  tb_start.tb_low  = (unsigned int)((*start) & 0xFFFFFFFF);
  tb_start.flag = tb_flag;
  time_base_to_time(&tb_start, TIMEBASE_SZ);

  tb_end.tb_high = (unsigned int)((*end) >> 32);
  tb_end.tb_low  = (unsigned int)((*end) & 0xFFFFFFFF);
  tb_end.flag = tb_flag;
  time_base_to_time(&tb_end, TIMEBASE_SZ);
    
  sec  = tb_end.tb_high - tb_start.tb_high;
  nsec = tb_end.tb_low  - tb_start.tb_low;
  if ( nsec < 0 ) {
    sec --;
    nsec += 1000000000;
  } 
  return (double) sec + 1.0e-9 * ((double) nsec);
}

double FNAME(timer_cpu_seconds)
( void )
{
    timebasestruct_t tb;
    double wtime;
    
    read_real_time(&tb, TIMEBASE_SZ);
    time_base_to_time(&tb, TIMEBASE_SZ);

    wtime  =  (double)tb.tb_low * 1.0e-9;
    wtime +=  (double)tb.tb_high;
    return wtime;
}

double FNAME(timer_resolution)
( void )
{
  timebasestruct_t tb1,tb2;
  int res;
   
  read_real_time(&tb1, TIMEBASE_SZ);
  do {
       read_real_time(&tb2, TIMEBASE_SZ);
  } while (tb1.tb_low == tb2.tb_low);
  
  time_base_to_time(&tb1, TIMEBASE_SZ);
  time_base_to_time(&tb2, TIMEBASE_SZ);
  
  res = tb2.tb_low - tb1.tb_low;
  if ( res < 0 ) res += 1000000000;
    
  return (double)(res) * 1e-9;

}

#endif /* __AIX_TIMER__ */


#ifdef __MACH_TIMER__
#define HAVE_TIMER

/* 
  Use high resolution timers from the Mach kernel. The resolution is system
  dependent and is on the order of tens of nanoseconds.
  
  - Resolution in a G4 800Mhz: 60.22 ns
   
  There is a (unlikely) possibility of overflow that is not checked. 
  (this could be fixed by recording mach_absolute_time() when timer_init 
  is called and then checking for overflow on timer_cpu_seconds)
  
*/

#include <mach/mach_time.h>

static mach_timebase_info_data_t    sTimebaseInfo;
static double tick_conv;

void FNAME(timer_init)
(void)
{
  (void) mach_timebase_info(&sTimebaseInfo);
  tick_conv = (double) sTimebaseInfo.numer / (double) sTimebaseInfo.denom * 1e-9;
}

uint64_t FNAME(timer_ticks)
(void)
{
  return mach_absolute_time();
}

double FNAME(timer_interval_seconds)
(uint64_t *start, uint64_t *end)
{
  /* we should check for overflow */
  
  return (*end - *start) * tick_conv;
}

double FNAME(timer_cpu_seconds)
( void )
{
  return (double)mach_absolute_time() * tick_conv;
  
}

double FNAME(timer_resolution)
( void )
{
  uint64_t        start, end;
  start = mach_absolute_time();
  end = mach_absolute_time();
  
  return (end-start)*tick_conv;
}

#endif /* __MACH_TIMER__ */

#ifdef __HRTIME_TIMER__
#define HAVE_TIMER

/* 
  Use high resolution timers based on the gethrtime function,
  available on most linux distros.  
*/

#include <sys/time.h>

void FNAME(timer_init)
(void)
{
 /* No initialization required */
}

hrtime_t FNAME(timer_ticks)
(void)
{
  return gethrtime();
}

double FNAME(timer_interval_seconds)
(hrtime_t *start, hrtime_t *end)
{
  /* we should check for overflow */
  
  return (*end - *start) * 1.0e-9;
}

double FNAME(timer_cpu_seconds)
( void )
{
  return (double)gethrtime() * 1.0e-9;
  
}

double FNAME(timer_resolution)
( void )
{
  hrtime_t        start, end;
  start = gethrtime();
  end = gethrtime();
  
  return (end-start)*1.0e-9;
}

#endif /* __HRTIME_TIMER__ */

/***************************************************************************************************
  BlueGene/P specific code (includes timers)
***************************************************************************************************/

#ifdef __bgp__         

#include "mpi.h"
#include "mpix.h"

/* returns BlueGene hardware topology has an MPI topology ordered has [T,Z,Y,X] */

void FNAME(bgp_mpix_cart_comm_create_f)
(MPI_Fint *comm, int *ierr)
{
  MPI_Comm cart_comm_bg;

  *ierr = MPIX_Cart_comm_create (& cart_comm_bg);
  if (*ierr != MPI_SUCCESS ) cart_comm_bg = MPI_COMM_NULL;

  *comm = MPI_Comm_c2f( cart_comm_bg );
}


/* returns BlueGene memory per core in MBytes */

#include <sys/resource.h>
#include <common/bgp_personality.h>
#include <common/bgp_personality_inlines.h>
#include <spi/kernel_interface.h>
static _BGP_Personality_t mybgp;

void FNAME( bgp_core_memory_f )
( long *bg_coreMB )
{
  unsigned procMB, coreMB;
  Kernel_GetPersonality(&mybgp, sizeof(_BGP_Personality_t));
  procMB = BGP_Personality_DDRSizeMB(&mybgp);
  coreMB = procMB/Kernel_ProcessCount();
  *bg_coreMB = coreMB;
}

#include <unistd.h>

/* returns BlueGene process memory in kbytes */

void FNAME(bgp_used_memory_f)
( long *used_kb )
{
  struct rusage usage;
  if ( getrusage(RUSAGE_SELF, &usage) ) *used_kb = 0;
  else *used_kb = usage.ru_maxrss;

}


/* BG specific timers */
#define HAVE_TIMER

static double tick_conv;

void FNAME(timer_init)
(void)
{
   unsigned clockMHz;
   Kernel_GetPersonality(&mybgp, sizeof(_BGP_Personality_t));
   clockMHz = BGP_Personality_clockMHz(&mybgp);
   tick_conv = 1.0e-6/(double)clockMHz;
}

uint64_t FNAME(timer_ticks)
( void )
{
   return _bgp_GetTimeBase();
}

double FNAME(timer_interval_seconds)
(uint64_t  *start, uint64_t *end)
{
  /* we should check for overflow */
  
  return (*end - *start) * tick_conv;
}

double FNAME(timer_resolution)
( void )
{
  uint64_t        start, end;
  start = _bgp_GetTimeBase();
  end   = _bgp_GetTimeBase();
  
  return (end-start)*tick_conv;
}

double FNAME(timer_cpu_seconds)
( void )
{
    return (double)_bgp_GetTimeBase() * tick_conv;
}


#endif


/***************************************************************************************************
  BlueGene/Q specific code
***************************************************************************************************/

#ifdef __bgq__

#include <unistd.h>

/* returns BlueGene/Q process memory in kbytes */

#include <spi/include/kernel/memory.h>

void FNAME(bgq_used_memory_f)
( long *used_kb )
{
  uint64_t allocated_memory = 0;
  
  /* Size in bytes of the heap size*/
  Kernel_GetMemorySize(KERNEL_MEMSIZE_HEAP, &allocated_memory);  
  *used_kb = allocated_memory / 1024;
}


/* BG/Q specific timers */
#define HAVE_TIMER

#include <hwi/include/bqc/testint_inlines.h>
static double tick_conv;

void FNAME(timer_init)
(void)
{
   /* The core frequency should always be 1600 MHz */
   
   /* This causes a signal 11 */
   /*tick_conv = 1.0e-6 / (double) TI_CoreFrequency();*/

   tick_conv = 1.0e-6 / 1600.0;
}


uint64_t FNAME(timer_ticks)
( void )
{
   return GetTimeBase();
}

double FNAME(timer_interval_seconds)
(uint64_t  *start, uint64_t *end)
{
  /* we should check for overflow */
  
  return (*end - *start) * tick_conv;
}

double FNAME(timer_resolution)
( void )
{
  uint64_t        start, end;
  start = GetTimeBase();
  end   = GetTimeBase();
  
  return (end-start)*tick_conv;
}

double FNAME(timer_cpu_seconds)
( void )
{
    return (double) GetTimeBase() * tick_conv;
}

#endif


/***************************************************************************************************
  MPI Timers
***************************************************************************************************/
#ifdef __MPI_TIMER__
#define HAVE_TIMER

#include <stdint.h>

void FNAME(timer_init)
(void)
{
  /* No initialization required for these timers */
}

uint64_t FNAME(timer_ticks)
( void )
{
  /* Return number of microsseconds */
  return 1.0e6*MPI_Wtime();
}

double FNAME(timer_interval_seconds)
(uint64_t  *start, uint64_t *end)
{
  return (*end - *start) * 1.0e-6;
}

double FNAME(timer_resolution)
( void )
{
  return MPI_Wtick();
}

double FNAME(timer_cpu_seconds)
( void )
{
    return MPI_Wtime();
}

#endif /* __MPI_TIMER__ */


/***************************************************************************************************
  Default timers (use POSIX timers)
***************************************************************************************************/

#ifndef HAVE_TIMER

/*
  Use standard POSIX timers based on gettimeofday. (This is actually how
  LAM and OpenMPI implement MPI_Wtime). Minimum resolution is 1 microssecond.
  
  The actual resolution is measured by finding the minimum difference >0 
  between succsessive gettimeofday calls. 
*/

#include <sys/time.h>
#include <stdint.h>

void FNAME(timer_init)
(void)
{
  /*  No initialization is required for this type of timers */
}

uint64_t FNAME(timer_ticks)
(void)
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  
  return ((uint64_t)tv.tv_sec)*1000000 + (uint64_t)tv.tv_usec;
}

double FNAME(timer_interval_seconds)
(uint64_t *start, uint64_t *end)
{
  return (*end - *start) * 1.0e-6;
}

double FNAME(timer_cpu_seconds)
( void )
{
    struct timeval tv;
    double wtime;
    gettimeofday(&tv, NULL);
    wtime = tv.tv_sec;
    wtime += (double)tv.tv_usec * 1.0e-6;
    return wtime;
}

double FNAME(timer_resolution)
( void )
{
  struct timeval tv1, tv2;
   
  gettimeofday(&tv1, NULL);
  do {
       gettimeofday(&tv2, NULL);
  } while (tv1.tv_usec == tv2.tv_usec);
  
  return (double)(tv2.tv_usec - tv1.tv_usec) * 1.0e-6;
}


#endif 
