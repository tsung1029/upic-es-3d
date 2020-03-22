/*

 Interfacing routines for POSIX system calls
 under macosx (should work under any POSIX
 compliant system)

*/


#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <string.h>


/* Remove a file */

void remove_f_(const char *path, int *ierr)
{
 if (*ierr = remove(path)) *ierr = errno;
}


/* Create a directory */

void mkdir_f_(const char *path, int *ierr)
{
 char uppath[256], *p;
    
 if (mkdir(path,0755))
   switch (errno) {
     case ENOENT : // A component of the path does not exist
        
        // get upper path
        
        strcpy(uppath, path);
        p = uppath + strlen(uppath);
        while(*p!='/') p--;
        *p=0;
        
        // recursively build the path
        
        mkdir_f_(uppath, ierr);
        if (!*ierr) mkdir_f_(path, ierr);
     
     case EEXIST : // if directory already exists ignore the error
        *ierr = 0;
        break;
     default: *ierr = errno;
   }
 else *ierr = 0;
 
}

/* Create a symbolic link */

void symlink_f_(const char *name1, const char *name2, int *ierr)
{
  if (*ierr = symlink(name1,name2)) *ierr=errno;
}

/* Get the error string */

void strerror_f_(const int *errnum, char *err)
{
 strcpy(err, strerror(*errnum));
}

/* Get the environment variable value */

void getenv_f_(const char *name, char *value, int *ierr)
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

/* check if floating point is infinity */

void isinf_f_( double *x, int *res )
{
  *res = ! finite( *x);
}

/* check if floating point is nan */

void isnan_f_( double *x, int *res )
{
  *res = isnan( *x);
}
