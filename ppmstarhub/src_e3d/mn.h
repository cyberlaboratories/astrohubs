/* ********************************************************************************* */
/* ********************* Minimal Networking API ************************************ */
/* ********************************************************************************* */

/* mn.h  November 12, 2004  Ted Wetherbee ted@fdltcc.edu

	rev. Nov. 22, 2004: completed disk<->mem xfer functions, cleaned up redundancies
	rev. Jan. 24, 2005: completed C++ and blocking/default API for windows
	rev. Mar. 25, 2005: mem file system, put & get only mem2mem, renamed mn from ib
	rev. Mar. 28, 2005: removed __int64, int & void* will be OK in emt64
	                    shared memory tests good
						     port to Linux tests good
	rev. Mar. 30, 2005: separated shmem from sockets server
						     added rpc() features
	rev. Mar. 31, 2005: added atomic_rpc() 
	rev. Apr.  6, 2005: added diskIO and integrated server(disk, mem, rpc)
	rev. Apr. 25, 2005: ported to Linux, aligned data types and structs,
						     successfully tested sockets code on Mac OS X
	rev. Apr. 27, 2005: touch file page memory ahead threading
	rev. May   3, 2005: removed touch page memory.  Overlapped IO is better for windows.
	                    sequential IO on Linux is just as good as touching page memory
	rev. May   3, 2005: added put_disk2mem, get_mem2disk, and overloaded put & get,
						and local put and get permutations of mem & disk xfer
	rev. May  17, 2005: added security, 3 levels
	rev. July 30, 2005: removed high security, used size_t instead of __int64 throughout
	                    when delaing with memory and file sizes and offsets to simplify
						arithmetic, use of fstat() verses fstat64() in Linux since
						fstat64() is not found but should go to 64 bit stat using
						HUGEFILE options if not yet enabled here
	rev. Aug. 19, 2005: 1) Combined all source files into mn.c, all headers into mn.h
	                    This fits the programming needs of FORTRAN programmers.
						2) Removed multi-threaded server disk and memory operations, left
						rpc multi-threaded
						3) A write to a disk, local or remote, with zero length deletes
						the target disk if it exists.
	rev. Aug. 20, 2005: mn_swap() exchanges values between local memory and remote memory,
						a "small" bit of memory which can fit within a control packet,
						certainly 128 bytes which is less than SERV_DATA_LEN.  Currently,
						SERV_DATA_LEN is 200 bytes.  mn_swap() is atomic on a particular port
						so it can be used to manage concurrency, say, by sharing a token
						between processes.
						 
*************************************************************************************** */

/* ******** programmer options *********************************
   1.  set system LINUX or WINDOWS in "Building Options"
   2.  to print or not print BUG() and SBUG() lines in "Building Options"
   3.  set Security level in mn_build.h; must be same for client and server sides.
   4.  set server packet buffer length  SERV_BUF_LEN
   5a. define memory chunk size MEM_CHUNK_SIZE. 16 or 32MB should be fine
   5b. define DISK_CHUNK_SIZE.  16 or 32MB should be fine
   6.  define SERV_BUF_LEN as needed to tune for pattern of sending memory or as
       needed for big rpc arguments, remex executable paths, or remote file names.
	   256 should be fine for most applications.


NOTES
   1. Porting between 64bit and 32bit
Windows 64 bit code uses 64 bit pointers but int is 32 bits.  Thus, casting pointers
for some arithmetic requires INT_PTR.  Linux might be different in this matter.  
In any case, it would probably be safe to do pointer arithmetic entirely in 64bit
casting with __int64 in windows and __int64_t in Linux.  mn_build.h replaces
__int64 with __int64_t when LINUX is defined, so the __int64 casting of pointers
appears to be a good scheme.  This was changed to using size_t later.

*/

// **********************************************************
// **** Building Options ************************************
// **********************************************************

// $$$$$$$$ comment out one of these ************************
//#define WINDOWS
#define LINUX

// $$$$$$$$ comment out one of these out ********************
#define BUG(report, arg) printf(#report, arg)
//#define BUG(report, arg) //for no user-level debugging

// $$$$$$$$ comment one of these out ************************
#define SBUG(report) perror(#report)
//#define SBUG(report) //for no system perror messages 



// **********************************************************
// **** Include Files ***************************************
// **********************************************************
// ***** common header for everything lower-level ***********
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include<time.h>

#ifdef WINDOWS
#include<winsock2.h> // add wsock32.lib to project settings
#include <io.h>
//#include<process.h>
#include<stddef.h>
#else // Linux
#include <sys/mman.h>
#include <sys/ipc.h>
#include <sys/shm.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <unistd.h> // close
#include<pthread.h> 
#endif

// **********************************************************
// **** Definitions  ****************************************
// **********************************************************
/*
#ifndef WINDOWS // defined in stddef.h for windows

#ifndef u_char
  #define u_char unsigned char // 1 byte
#endif

#ifndef u_short
  #define u_short unsigned short // 2 byte
#endif

#ifndef u_long
  #define u_long unsigned long // 4 byte
#endif

#ifndef u_int
  #define u_int unsigned int
#endif

#endif
*/

// protocol numbers
#ifndef AF_IBT
  #define AF_IBT 26 
#endif

#ifndef AF_INET
	#define AF_INET 2 
#endif


// **********************************************************
// ***** Security *******************************************
// **********************************************************
#define MAGIC_NUMBER 1201762731
/* The magic number is arbitrary, but it must be the same between 
   client and server.  If this is commented out, there is no security
   other than the server feature which drops unrecognized packets
   without reply.  This should always be defined whether security is
   enabled or not.                                           */

//#define SECURITY
/* Turns on security check within server loop.  Comment out
   to turn off magic number checking.                        */


// **********************************************************
// **** Definitions  ****************************************
// **********************************************************


// ******** shmem ID ****************************************
#ifdef LINUX
#define MN_IPC_KEY 1856253913 // for remote acces to shmem
#endif
/* This value is arbitrary.  Windows uses a character string
   formed from concatenation of "mn" and the server port number. 
   But, it has to be the same for all communicating shmem codes */

// These two chunk sizes better be the same now
#define MEM_CHUNK_SIZE 0x1000000 // 2^24 ~ 16MBytes, mem2mem
#define DISK_CHUNK_SIZE 0x1000000 // 2^24 ~ 16MBytes, disk side of anything
/* Twice this size should be available to malloc and Valloc (windows),
   so it should be bumped down if necessary, or bumped up if desired.
   These must be powers of 2 else the algorithms in MN break.

   Keep both these chunk sizes exactly the same.      */


// ******** porting between Linux and Windows **********************
#ifndef WINDOWS // LINUX
#define SOCKET int// convenience
#define __int64 __int64_t//linux 64 bit signed integer
#define __int32 long// long is an ANSI 32 bit signed integer
#define closesocket(arg) close(arg)// convenience, save typing & #ifdef'ing
#define _endthread() pthread_exit(NULL)// very different threads, here convenience
#define _O_RDONLY O_RDONLY// lots of differences in file IO
#define _O_CREAT O_CREAT
#define _O_WRONLY O_WRONLY
#define _O_BINARY 0
#define _S_IWRITE 00777
#define _lseeki64 lseek64
#define _read read
#define _write write
#define _close close
#define _open open
//#define _stati64 stat64
//#define _fstati64 fstat64
//#define _stat stat
#endif

// *******************************************************************
// ***** Communication Buffers and Structures ************************
// *******************************************************************

// ***** flags for server_init
#define MN_SHMEM 1
#define MN_RPC 2
#define MN_DISK 4
#define MN_REMEX 8
#define MN_HANG 16

#define SERV_BUF_LEN 256
// *******************************************************************
// Common packets for communication between client and server
/* This is the common size for communication packets between client and
   server functions.  The beginning of each packet contains necessary
	data for functioning of the client and server functions.  The remander
	is available for data: memory for small put and get, file names, arguments
	for rpc, executable pathname for remex.

   Currently, at most 40 bytes are used for packet code, so SERV_BUF_LEN - 40
	is available for carrying data within a packet.

   If one is sending and receiving small amounts of memory, say 460 bytes using
	put() and get(), it would be a good idea to raise SER_BUF_LEN to 512 so
	that put() and get() could do these small data transfers more quickly without
	the overhead of setting up a memory transfer for data which cannot fit within
	a packet.

	This should be a power of 2 so that the algorithms work.       */ 

// ****** structure alignment ****************************************
/* These packet structures must be identical across the network on a 
   variety of systems, so there is this messy looking alignment with 
   unions and a ".x" after every variable aligned with char s[8]. */
// **** memory get & put packet **************************************
#define SERV_CODE_LEN 56 // max length of code fields, 8 bite
#define SERV_DATA_LEN  SERV_BUF_LEN - SERV_CODE_LEN  // data len
typedef struct{
	union{__int64 x; char s[8];} code; // action 
	union{__int64 x; char s[8];} magic_number;
	union{size_t x; char s[8];} rem_offset;// from base address, if shmem
	union{size_t x; char s[8];} len;// length of data
	union{SOCKET x; char s[8];} sock; // socket at server, internal
	union{size_t x; char s[8];} chunk;// 0 if mem2mem, otherwise a client disk
	union{void  *x; char s[8];} addr;// base address if shmem, internal 
	char data[SERV_DATA_LEN];// data if carried, carries whatever
} packet;

// **** rpc packet ***************************************************
#define RPC_SERV_CODE_LEN 40 // max needed for passed variables
#define RPC_SERV_DATA_LEN   SERV_BUF_LEN - RPC_SERV_CODE_LEN  // data len
typedef struct{
	union{__int64 x; char s[8];} code;//code
	union{__int64 x; char s[8];} magic_number;
	union{size_t x; char s[8];} len; // length of data within packet
	union{void *x; char s[8];} addr;// the function address to be passed, internally
	union{SOCKET x; char s[8];} sock;// socket if rpc handler is spawned, internal
	char data[RPC_SERV_DATA_LEN];// data if carried
} rpc_packet;

// **** disk get & put packet ****************************************
#define DISK_SERV_CODE_LEN 48  // needed for passed variables 
#define DISK_SERV_DATA_LEN  SERV_BUF_LEN - DISK_SERV_CODE_LEN // data
typedef struct{
	union{__int64 x; char s[8];} code; // action to perform by code
	union{__int64 x; char s[8];} magic_number;
	union{size_t x; char s[8];} offset; // offset on remote disk
	union{size_t x; char s[8];} len; // length of data to xchange
	union{SOCKET x; char s[8];} sock;// to pass from server to sub thread, internal
	union{size_t x; char s[8];} chunk;// mediated size of transfer buffers 
	char rem_file[DISK_SERV_DATA_LEN];// rem disk file name
} disk_packet;


// ***** Codes between client and server, carried in packets ******
// put and get for memory on server side
#define CODE_PUT 11				
#define CODE_PUT_DONE 12	
#define CODE_PUT_BIG 101		
#define CODE_PUT_BIG_OK 102		
#define CODE_PUT_BIG_DONE 103	
#define CODE_GET 21					
#define CODE_GET_BIG 201		

// rpc
#define CODE_RPC 31	
#define CODE_RPC_DONE 32		
#define CODE_ATOMIC_RPC 301
#define CODE_ATOMIC_RPC_DONE 302

// swap
#define CODE_SWAP 41
#define CODE_SWAP_DONE 401

// disk put and get for disk on server side
#define CODE_DPUT 51
#define CODE_DPUT_DEL 510 // delete target if len=0
#define CODE_DPUT_OK 52
#define CODE_DPUT_DONE 53
#define CODE_DGET 61

// remote execution
#define CODE_REMEX 71
#define CODE_REMEX_DONE 72

// server control
#define CODE_STOP 91			
#define CODE_DIE 92			

// disk server abort client
#define CODE_ABORT 666

// **************************************************************
// ***** MN Function Return Codes *******************************
// ************************************************************** 
// for all except when returning an address in which failure = NULL

// general errors
#define MN_SUCCESS 0 // success for all except getting shmem address
#define MN_FAIL -1	// used by make_sock for failure, also remex
#define MN_LEN_UNKNOWN -2 // only an error if not a disk2disk 
#define MN_LEN -3 // too big

// socket errors 
#define MN_MSOCK_FAIL -10
#define MN_SOCK -11

// packet communications errors
#define MN_SEND -31
#define MN_RECV -32
#define MN_CODE -33
#define MN_DONE -33

// memory and shmem errors
#define MN_MALLOC -40
#define MN_UNMAP -41
#define MN_SHMDT -42
#define MN_SHMGET -43
#define MN_SHMCTL -44
#define MN_2BIG -45

// disk errors
#define MN_FSTAT -50
#define MN_OPEN -51
#define MN_READ -52
#define MN_LSEEK -53
#define MN_WRITE -54

// server errors
#define MN_SHMEM_RPC -65
#define MN_BEGTHREAD -66
#define MN_TIMEOUT -67
#define MN_BIND -68
#define MN_LISTEN -69


// **************************************************************
// **************************************************************
// ***** Internal Functions *************************************
// **************************************************************
// **************************************************************

// ***** makesock ***********************************************
SOCKET make_sock(u_short AF, char *rem_host, u_short rem_port);
// returns -1 on failure, the socket (!= -1) on success


// ***** server calls *******************************************
/* These are called by the server thread only to handle communications
   in a separate thread so that the server listening port remains
   available for new incoming packets.     */
#ifdef WINDOWS
void put_disk_serv(void *in_pack);// called by _beginthread
#else // LINUX
void *put_disk_serv(void *in_pack);// called by pthread
#endif

#ifdef WINDOWS
void get_disk_serv(void *in_pack);// called by _beginthread
#else // LINUX
void *get_disk_serv(void *in_pack);// called by pthread
#endif

#ifdef WINDOWS
void put_big_serve(void *in_pack);// called by _beginthread
#else
void *put_big_serve(void *in_pack);// called by pthread 
#endif

#ifdef WINDOWS
void get_big_serve(void *in_pack);// called by _beginthread 
#else
void *get_big_serve(void *in_pack);// called by pthread 
#endif

#ifdef WINDOWS
void rpcserv(void *in_rpc_pack);// called by _beginthread
#else
void *rpcserv(void *in_rpc_pack);// called by pthread 
#endif

// This is called directly, not as a thread, and it creates a new
// process, not a thread.
int remex(char *command);



/* ********************************************************************************* */
/* ********************************************************************************* */
/* ***** MN API ******************************************************************** */
/* make_shmem() and get_shmem() return NULL on failure, an address on success
   All other functions return 0 on failure, 1 on success */

/* ********************************************************************************* */
/* ***** winsock initialization **************************************************** */ 
/* needed for all Windows applications, does nothing for UNIX
*/
int mn_sockets_init();
int mn_sockets_end();

/* ********************************************************************************* */
/* ***** server functions ********************************************************** */

// ***** shared memory ************************************************
/* make_shmem returns an address on success, NULL (0) on failure.
*/
void *mn_make_shmem(u_short port, size_t mem_file_size);
int mn_remove_shmem(u_short port, void * shared_memory_start_addr);

// ***** memory socket server *****************************************
/* None of the client functions, except shmem, do anything without a target server 
   running at a particular port, hostname, and address family.
*/
int mn_start_serv(u_short AF, u_short port, char *host, 
		void* mem_block_addr, // memory for shmem, function for rpc, NULL otherwise 
		u_char flag // flag to set allowable calls
		);
// ***** flag values for start_serv
// MN_SHMEM 1
// MN_RPC   2
// MN_DISK  4
// MN_REMEX 8
// MN_HANG 16
// Cannot use both MN_SHMEM and NM_RPC together since the server binds either to a memory 
// base address or a function address, not both.  All other combinations are allowed.
// E.g. MN_SMEM | MN_DISK | MN_REMEX 

int mn_stop_sock_serv(u_short AF, u_short rem_port, char *host);
int mn_kill_sock_serv_process(u_short AF, u_short rem_port, char *host);

/* ********************************************************************************* */
/* ***** client functions ********************************************************** */

// ***** shared memory ***********************************************
/* If shmem was created by a process on an SMP, another process can use it
   with get_shmem.  Several shmems can be created each with its own port
   number.  The port number is an identifier only, but it is presumed to be
   a valid port used by any server which would attach to a shmem.

   detach_shmem only removes the client access and a system reference count, 
   not the shmem.
*/
void *mn_get_shmem(u_short port);
int mn_detach_shmem(void *shared_memory_start_addr);

// ***** RPC services ************************************************
/* rpc spawns the remote function then returns immediately.  atomic_rpc runs the
   remote function and returns the argument, possibly modified, and only one
   atomic rpc is processed at a time.  Rpc functions must be like this..
   Windows:  void (*)(void *)  , e.g. void func(void *arg)
   Linux:    void* (*)(void *) , e.g. void *func(void *arg)
*/
int mn_rpc(void *arg, 
			size_t len, 
			u_short AF, u_short port, char *host);
int mn_atomic_rpc(void *arg, 
			size_t len, 
			u_short AF, u_short port, char *host);

// ***** remote execution *********************************************
int mn_remex(char *command, 
			u_short AF, u_short rem_port, char *rem_host);

// ***** swap (atomic) ************************************************
int mn_swap(size_t rem_offset,
			void *loc_mem,
			size_t len,
			u_short AF, u_short port, char *host);

// ***** put and get **************************************************
int mn_put(size_t rem_offset, 
		   void *loc_mem, 
		   size_t len, 
		   u_short AF, u_short port, char *host);

int mn_get(size_t rem_offset, 
		   void *loc_mem, 
		   size_t len, 
		   u_short AF, u_short port, char *host);

int mn_put_disk2disk(char *rem_disk, size_t rem_offset, 
					 char *loc_disk, size_t loc_offset,
					size_t len, 
					u_short AF, u_short rem_port, char *rem_host);

int mn_put_mem2disk(char *rem_disk, size_t rem_offset, 
					void *loc_mem,
					size_t len,
					u_short AF, u_short rem_port, char *rem_host);

int mn_put_disk2mem(size_t rem_offset,// from offset of server base memory 
					 char *loc_disk, size_t loc_offset, 
					 size_t len, 
					 u_short AF, u_short rem_port, char *rem_host);

int mn_get_disk2disk(char *rem_disk, size_t rem_offset, 
					 char *loc_disk, size_t loc_offset, 
					size_t len, 
					u_short AF, u_short rem_port, char *rem_host);

int mn_get_disk2mem(char *rem_disk, size_t rem_offset, 
					void *loc_mem, 
					size_t len,
					u_short AF, u_short rem_port, char *rem_host);

int mn_get_mem2disk(size_t rem_offset, // offset from server base memory 
					char *loc_disk, size_t loc_offset, 
					size_t len, 
					u_short AF, u_short rem_port, char *rem_host);

