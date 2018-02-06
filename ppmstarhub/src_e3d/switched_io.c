/* switched_io.c : Routines for switched Input and Output
 * David Porter
 * Laboratory for Computational Science & Engineering
 * University of Minnesota
 * August 2004
 *
 * Purpose: provide a unified application program interface for FORTRAN
 * and C routines to write blocks of binary data to local or remote disk or
 * memory via the most efficient method available.
 *
 * Routines provided:
 *    integer sw_write(host, directory, filename, offset, nbytes, buffer)
 *    integer sw_read (host, directory, filename, offset, nbytes, buffer)
 *
 * Arguments:
 *    host, directory, and filename are character strings
 *    offset is a 8-byte integer (use integer(kind=8) in Intel Fortran 8.0)
 *    nbytes is a 4-byte integer
 *    buffer is a pointer to memory to start writting from
 *
 * Return Values: integer
 *    If succesful     : number of bytes read or written
 *    If NOT succesful : error code < 0
 *
 * Both routines are blocking.
 */

#ifdef _WIN32
#include <windows.h>
#include <process.h>
#include <io.h>
#include <winbase.h>
typedef __int64 Offset_t;
#else
#include <unistd.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
typedef long long Offset_t;
// #include <pthread.h>
#endif

#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include "mnt.h"
// #include <errmsg.h>

int sw_write( char *host, char *dir, char *name, Offset_t offset, int nbytes, void *buffer)
{
#ifdef _WIN32
    HANDLE hFile,fd2;
#else
    static int fd = -1;
#endif
    Offset_t shift32, shift8;
    long dwHigh, dwLow;
    long nwrote;
    char cDirName[256]; char * buf = NULL;
    int itriedsomething;
    int * iostatus;
    char path[256];
    int err,index;
    char cHostName[256];
    unsigned long nSize[256];
    int *result, igot;
    char interface[256];

    itriedsomething = 0;

    shift8  = (Offset_t)256;
    shift32 = shift8 * shift8 * shift8 * shift8;
    dwHigh  = (long)(offset / shift32);
    dwLow   = (long)(offset % shift32);

    if(strncmp(host, "localdisk",9) == 0) {
       itriedsomething = 1;

      // Generage full path to local file and open it
      sprintf(cDirName, "%s/%s", dir, name);
#ifdef _WIN32
      hFile = CreateFile( cDirName, GENERIC_WRITE, 
			0, NULL, OPEN_ALWAYS, 
                        FILE_ATTRIBUTE_NORMAL,
			// ( FILE_FLAG_OVERLAPPED | FILE_FLAG_SEQUENTIAL_SCAN),
                        NULL);
      if (hFile == INVALID_HANDLE_VALUE) {
        printf("Error -1010:  SW_WRITE: Local disk write - CreateFile generated an invalid handle\n");
        return -1010;
      }

      // Write to the file
      SetFilePointer(hFile,dwLow, &dwHigh, FILE_BEGIN);
      WriteFile(hFile, buffer, nbytes, &nwrote, NULL);

      // Close the file
      CloseHandle(hFile);
#else
      if(nbytes == 0) { // Just toggle open or closed
          if(fd < 0) { fd = open(cDirName, O_RDWR|O_CREAT, 0644); }
          else       { close(fd); fd = -1; }
          return 0;
      }

      if(fd < 0) {
          fd = open(cDirName, O_RDWR|O_CREAT, 0644);
          if ( fd == -1 ) {
            printf("Error -1010:  SW_WRITE: Local disk write - file open failed\n");
            return -1010;
          }
          lseek( fd, offset, SEEK_SET );
          nwrote= write( fd, buffer, nbytes );
          close(fd);
          fd = -1;
      } else {
          lseek( fd, offset, SEEK_SET );
          nwrote= write( fd, buffer, nbytes );
      }
#endif

      if (nwrote != nbytes) {
        printf("Error -1020:  SW_WRITE: Local disk write - short write\n");
        return -1020;
      }

    } else {
      /***else use RPC ****/
     itriedsomething = 1;

      int af, port, bytes;
/* port = (host[0]-'0')*10000 + (host[1]-'0')*1000 + (host[2]-'0')*100 +
          (host[3]-'0')*10 + (host[4]-'0'); */
      igot = sscanf(host, "%d %d %s", &af, &port, interface);

      if(igot != 3  ||  !(af == 2 || af == 26)) {
        printf("sw_write: problem with specfied host: %s\n", host);
        printf("   host should either be \"localdisk\" or \"af port interface\"\n");
        printf("    af        = address family (2 or 26)\n");
        printf("    port      = port number\n");
        printf("    interface = network address of the interface an memservt is attached to\n");
        return -1;
      }
      
      if(strcmp(dir,"memory")==0) {
	err = -1;
		mn_sockets_init();
		/* err = mn_put(offset,buffer,nbytes,2,19950,host+5); */
		err = mn_put(offset,buffer,nbytes,af,port,interface);
		mn_sockets_end();
		if (err != MN_SUCCESS) { /*err = -1;*/ printf("put err = %d\n", err); usleep(100000);}
		else                   { return(err = nbytes); }
      } else {
	sprintf(cDirName, "%s/%s", dir, name);
                mn_sockets_init();
                /* err = mn_put_mem2disk(cDirName,offset,buffer,nbytes,2,19950,host+5); */
                err = mn_put_mem2disk(cDirName,offset,buffer,nbytes,af,port,interface);
                mn_sockets_end();
                if (err != MN_SUCCESS) { /*err = -1;*/  
			printf("put2disk error %d\n", err);
			usleep(100000); 
		}
                else                   { return nbytes; }


/*
          int mn_put_mem2disk(          char    *rem_disk,
                                        size_t  rem_offset,
                                        void    *loc_mem,
                                        size_t  len,
                                        u_short AF,
                                        u_short rem_port,
                                        char *rem_host    )

*/
      }
    }

    if(itriedsomething != 1) {
      printf("Error -1900:  SW_WRITE: Local disk write - IO method not supported\n");
      return -1900;
    }

    return (int)nwrote;
}

int sw_read( char *host, char *dir, char *name, Offset_t offset, int nbytes, void *buffer)
{
#ifdef _WIN32
    HANDLE hFile,fd2;
#else
    static int fd = -1;
#endif
    Offset_t shift32, shift8;
    long dwHigh, dwLow;
    long ngot;
    char cDirName[256];
    int itriedsomething; char * buf = NULL;
    int * iostatus; int err, index;
    char cHostName[256], path[256];
    unsigned long nSize[256];

    itriedsomething = 0;

    shift8  = (Offset_t)256;
    shift32 = shift8 * shift8 * shift8 * shift8;
    dwHigh  = (long)(offset / shift32);
    dwLow   = (long)(offset % shift32);

    if(strcmp(host, "localdisk") == 0) {
       itriedsomething = 1;

      // Generage full path to local file and open it
      sprintf(cDirName, "%s/%s", dir, name);

#ifdef _WIN32
      hFile = CreateFile( cDirName, GENERIC_READ, 
			FILE_SHARE_READ, NULL, OPEN_EXISTING,
                        FILE_ATTRIBUTE_NORMAL,
			// ( FILE_FLAG_OVERLAPPED | FILE_FLAG_SEQUENTIAL_SCAN),
                        NULL);
      if (hFile == INVALID_HANDLE_VALUE) {
        printf("Error -2010:  SW_READ: Local disk read - CreateFile generated an invalid handle\n");
        return -2010;
      }

      // read from the file
      SetFilePointer(hFile,dwLow, &dwHigh, FILE_BEGIN);
      ReadFile(hFile, buffer, nbytes, &ngot, NULL);

      // Close the file
      CloseHandle(hFile);
#else
      if(nbytes == 0) { // Just toggle open or closed
          if(fd < 0) { fd = open(cDirName, O_RDONLY); }
          else       { close(fd);  fd = -1; }
          return 0;
      }

      if(fd < 0) {
          fd = open(cDirName, O_RDONLY);
          if ( fd == -1 ) {
            printf("Error -2010:  SW_READ: Local disk read - file open failed\n");
            printf("SW_READ failed to open: -->%s<--\n", cDirName);
            return -2010;
          }
          lseek( fd, offset, SEEK_SET );
          ngot= read( fd, buffer, nbytes );
          close(fd);
          fd = -1;
      } else {
          lseek( fd, offset, SEEK_SET );
          ngot= read( fd, buffer, nbytes );
      }
#endif

      if (ngot != nbytes) {
        printf("Error -2020:  SW_READ: Local disk read - short read\n");
        printf("File path:=%s\n", cDirName);
        printf("Offset = %lld\n", offset);
        printf("nbytes got=%d   nbytes requested=%d\n", ngot, nbytes);
        return -2020;
      } else {
        return (int)ngot;
      }
   }

   /***else use RPC ****/
   if (strcmp(dir,"memory")==0) {
		int port;
		err = -1;
		 	//printf("began reading\n");
			mn_sockets_init();
			port = (host[0]-'0')*10000 + (host[1]-'0')*1000 +
			       (host[2]-'0')*100 + (host[3]-'0')*10 + (host[4]-'0');
			//printf("SW_READ:  hostname is %s bytes = %d port = %d\n",host+6,nbytes,port);
			err = mn_get(offset,buffer,nbytes,2,19950,host+6);//$$$$$
			mn_sockets_end();
			if (err != MN_SUCCESS) {
        			printf("SW_READ:  SW_READ: Local disk read - short read %d %d\n",err,nbytes); 
				err = -1;
				usleep(100000);
      			} else {
				return (err = nbytes);
			}
	} else {	
	int port; int bytes;
        port = (host[0]-'0')*10000 + (host[1]-'0')*1000 + (host[2]-'0')*100 +
             (host[3]-'0')*10 + (host[4]-'0');

	sprintf(cDirName, "%s/%s", dir, name);
                mn_sockets_init();
                err = mn_get_disk2mem(cDirName,offset,buffer,nbytes,2,port,host+6);//$$$$$
                mn_sockets_end();
                if (err != MN_SUCCESS) { /*err = -1;*/  
			printf("getdisk error %d\n", err); 
			usleep(100000);
			return (int)err;
		} else { 
			return nbytes; 
		}
 	}

   if(itriedsomething != 1) {
      printf("Error -2900:  SW_READ: Local disk read - IO method not supported\n");
      return -2900;
   }

    return (int)ngot;
}

/*******************************************************************
 *       WRAPPERS FOR CALLS FROM FORTRAN
 */

#ifdef _WIN32
#define sw_write_ SW_WRITE
#define sw_read_ SW_READ
#endif

#define CARG(t,s,n) i = 0; while(i<n  &&  s[i] != ' ') { t[i] = s[i]; i++; } t[i] = 0

int sw_write_( char *host, char *dir, char *name, Offset_t *offset, int *nbytes, void *buffer,
              int n1,      int n2,    int n3)
{
  int i;
  char str1[1024],str2[1024],str3[1024];
  //CARG(str1,host,n1);
  for(i=0; i<n1; i++) { str1[i] = host[i]; }
  i=n1; while(i > 0 && str1[i] == ' ') { str1[i] = 0;  i--; }
  CARG(str2,dir,n2);
  CARG(str3,name,n3);
  return sw_write(str1, str2, str3, *offset, *nbytes, buffer);
}

int sw_read_( char *host, char *dir, char *name, Offset_t *offset, int *nbytes, void *buffer,
              int n1,      int n2,    int n3)
{
  int i;
  char str1[1024],str2[1024],str3[1024];
  CARG(str1,host,n1);
  CARG(str2,dir,n2);
  CARG(str3,name,n3);
  return sw_read(str1, str2, str3, *offset, *nbytes, buffer);
}


/*****************************************************************
 * return endian-ness:  returns 1 for littleend (intel, sun, etc)
 *                               0 for bigend (mips,ppc, ...)
 */
#ifdef _WIN32
int ENDIAN()
#else
int endian_()
#endif
{
  unsigned int i=1;
  if ( *((unsigned char*)&i) ) return 1;
  else                         return 0;
}

