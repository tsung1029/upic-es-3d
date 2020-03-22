#include <string.h>
#include <stdio.h>

#ifdef __bgp__
#include <spi/kernel_interface.h>
#elif __bgq__
#include <spi/include/kernel/process.h>
#endif

int MPProcessors() {
/* Return the number of processors/hardware threads on the
   host computer on Linux platforms
local data                                              */
   int nproc = 0;

#ifdef __bg__

   /* simply the number of processes (Virtual Nodes) running 
      on this Physical Node times the number of Processors 
      (cores) running in this Process (Virtual Node). */

   nproc = Kernel_ProcessCount()*Kernel_ProcessorCount();

#else

   char *cnerr = 0;
   char cline[82];
   FILE *unit;
   unit = fopen("/proc/cpuinfo","r");
/* Quit if file does not exist */
   if (!unit)
      return nproc;
/* Read next line */
   while (cnerr = fgets(cline,81,unit)) {
      cline[9] = '\0';
      if (!strcmp(cline,"processor"))
         nproc += 1;
   }
   fclose(unit);

#endif

   return nproc;
}
