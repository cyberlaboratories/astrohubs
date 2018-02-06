// mn_th.c

#include "mnt.h"

// *************************************************************************
// *********** mn_error() **************************************************
// *************************************************************************
char *mn_error(int error_num){// returns a character string from an error code number 
	switch(error_num){
	case 0: return "MN_SUCCESS"; // success for all except getting shmem address
	case -1: return "MN_FAIL";	// used by make_sock for failure, also remex
	case -2: return "MN_LEN_UNKNOWN"; // only an error if not a disk2disk 
	case -3: return "MN_LEN"; // too big
	case -10: return "MN_MSOCK_FAIL";
	case -11: return "MN_SOCK";
	case -31: return "MN_SEND";
	case -32: return "MN_RECV";
	case -33: return "MN_CODE";
	case -34: return "MN_DONE";
	case -40: return "MN_MALLOC";
	case -41: return "MN_UNMAP";
	case -42: return "MN_SHMDT";
	case -43: return "MN_SHMGET";
	case -44: return "MN_SHMCTL";
	case -45: return "MN_2BIG";
	case -50: return "MN_FSTAT";
	case -51: return "MN_OPEN";
	case -52: return "MN_READ";
	case -53: return "MN_LSEEK";
	case -54: return "MN_WRITE";
	case -65: return "MN_SHMEM_RPC";
	case -66: return "MN_BEGTHREAD";
	case -67: return "MN_TIMEOUT";
	case -68: return "MN_BIND";
	case -69: return "MN_LISTEN";
	default: return "unknown";
	}
}

// *************************************************************************
// *************************************************************************
// ****** Sockets **********************************************************
// *************************************************************************
// *************************************************************************

// *************************************************************************
// *********** mn_sockets_init() *******************************************
// *************************************************************************
int mn_sockets_init(){
#ifdef WINDOWS // only needed for Windows
	int ret;
	WORD ver_req = MAKEWORD(1,1);
	WSADATA winsock_data;

	ret = WSAStartup(ver_req, &winsock_data);// init winsock & check version
	if (ret != 0 || winsock_data.wVersion != ver_req) return MN_FAIL;
#endif
	return MN_SUCCESS;// 1
}

// *************************************************************************
// ********** mn_sockets_end() *********************************************
// *************************************************************************
int mn_sockets_end(){
#ifdef WINDOWS // only needed for Windows
		int ret;
		ret = WSACleanup();
		if(0 == ret) return MN_SUCCESS;
		else return MN_FAIL;
#endif
	return MN_SUCCESS;// 1
}

// *************************************************************************
// ********** make_sock() **************************************************
// *************************************************************************
SOCKET make_sock(u_short AF, char *rem_host, u_short rem_port){
	SOCKET	c_sock; // client socket
#ifdef WINDOWS
	LPHOSTENT host_id;
	SOCKADDR_IN s_sock_addr; // server and client socket addresses
#else
	struct hostent *host_id;
	struct sockaddr_in s_sock_addr;
#endif
	int n;

#ifdef WINDOWS
	c_sock = socket(AF, SOCK_STREAM, IPPROTO_TCP);
	if (c_sock == INVALID_SOCKET) return INVALID_SOCKET;
#else
	c_sock = socket(AF, SOCK_STREAM, 0);
	if (c_sock == -1) return -1;
#endif

	s_sock_addr.sin_family = AF_INET;
	host_id = gethostbyname(rem_host);

	if (host_id == NULL){
#ifdef WINDOWS
		return INVALID_SOCKET;// failure
#else
		return -1;
#endif
	}

#ifdef WINDOWS
	s_sock_addr.sin_addr = *((LPIN_ADDR)*host_id->h_addr_list); // Server addr
#else
	memcpy(&s_sock_addr.sin_addr.s_addr, host_id->h_addr_list[0], host_id->h_length );
#endif

	s_sock_addr.sin_port = htons((u_short)(rem_port));	// Port number 

#ifdef WINDOWS
	n = connect(c_sock, (LPSOCKADDR)&s_sock_addr, sizeof(struct sockaddr));
	if(SOCKET_ERROR == n) return INVALID_SOCKET;//failure, 0==connect() if success
#else
	n = connect(c_sock, (struct sockaddr *)&s_sock_addr, sizeof(struct sockaddr));
	if(n<0) return -1;//failure
#endif
	
	return c_sock;// success, a socket handle
}

// *************************************************************************
// *************************************************************************
// ****** Shared Memory ****************************************************
// *************************************************************************
// *************************************************************************

/* make_shmem() is a server call, and remove_shmem() prevents others
   from attaching to the shmem.  get_shmem() and detach_shmem() are client
   calls to access the existing shmem and detach.  Behavior is not the
   same as malloc() and free().  make_shmem() creates the shmem and attaches
   the base address to the server calling it.  remove_shmem() does not actually
   remove the shmem; instead, it prevents additional attachment and removes
   the servers attachment to the shmem.  The shmem is not freed untill all
   client get_shmem() calls are matched by detach_shmem() calls.          */ 
//**************************************************************************

// *************************************************************************
// ********** mn_make_shmem() **********************************************
// *************************************************************************
void  *mn_make_shmem(u_short port, __int64 mem_file_size){// returns address, NULL for failure
	void *addr;
#ifdef WINDOWS 
	char key[10] = "mn", str_port[7];
	HANDLE fh;

	itoa(port,str_port, 10); strcat(key, str_port);	// form key from port
	fh = CreateFileMapping(
		(void *)0xFFFFFFFF, NULL, PAGE_READWRITE, 0, (u_int)mem_file_size, key);
	if(NULL == fh) return NULL;

	addr = MapViewOfFile(
		fh, FILE_MAP_ALL_ACCESS, 0, 0, 0);
	return addr;
#else // Linux
	int shmid;

	shmid = shmget(MN_IPC_KEY+port, (u_int)mem_file_size, IPC_CREAT | SHM_R | SHM_W);
	if(-1 == shmid) return NULL;
	
	addr = shmat(shmid, 0, 0);// NULL if failure
	return addr;
#endif
}

// *************************************************************************
// ******** mn_remove_shmem() **********************************************
// *************************************************************************
int mn_remove_shmem(u_short port, void * mbr){
	int ret;
#ifdef WINDOWS
	ret = UnmapViewOfFile(mbr);// detach
	if(ret == 0) return MN_UNMAP;
#else
	int shmid;
	ret = shmdt(mbr);// detach
	if(ret == -1) return MN_SHMDT;

	shmid = shmget(MN_IPC_KEY+port, 0, SHM_R | SHM_W);// get key
	if(shmid == -1) return MN_SHMGET;

	ret = shmctl(shmid, IPC_RMID, 0);// remove from system space
	if(-1 == ret) return MN_SHMCTL;
#endif
	return MN_SUCCESS;
}

// *************************************************************************
// ***** mn_get_shmem() ****************************************************
// *************************************************************************
void *mn_get_shmem(u_short port){// address on success, NULL on failure
	void *addr;
#ifdef WINDOWS
	HANDLE fh;
	char key[10]="mn", str_port[7];

	itoa(port, str_port, 10); strcat(key, str_port);// form key
	
	fh = OpenFileMapping(FILE_MAP_ALL_ACCESS, FALSE, key);// NULL on failure
	if(NULL == fh) return NULL;
	
	addr = MapViewOfFile(fh, FILE_MAP_ALL_ACCESS, 0, 0, 0);// NULL on failure
#else
	int shmid;
	
	shmid = shmget(MN_IPC_KEY+port, 0, SHM_R | SHM_W);// get key
	if(-1==shmid) return NULL;
	
	addr = shmat(shmid, 0, 0);
	if(-1==(__int64)addr) return NULL;
#endif
	return addr;
}

// *************************************************************************
// *********** mn_detach_shmem *********************************************
// *************************************************************************
int mn_detach_shmem(void * mbr){
	int ret;
#ifdef WINDOWS
	ret = UnmapViewOfFile(mbr);// detach
	if(ret == 0) return MN_UNMAP;
	else return MN_SUCCESS;// 0
#else // Linux
	ret = shmdt(mbr);// detach
	if(-1==ret) return MN_SHMDT;
	else return MN_SUCCESS;// 0
#endif
}

// *************************************************************************
// *************************************************************************
// ****** mn_remex() <<client>> ********************************************
// *************************************************************************
// *************************************************************************
int mn_remex(char *command, 
			 u_short AF, u_short rem_port, char *rem_host){
	rpc_packet p;
	SOCKET	c_sock; // client socket 
	__int64 sent, got;
	int n;
	
	if(strlen(command)+1 > RPC_SERV_DATA_LEN) return MN_LEN;// error

	// prepare packet
	p.magic_number.x = MAGIC_NUMBER;
	p.code.x = CODE_REMEX;
	memcpy(p.data, command, strlen(command) + 1);
   
	c_sock = make_sock(AF, rem_host, rem_port);
#ifdef WINDOWS
	if(c_sock == INVALID_SOCKET) return MN_SOCK;
#else
	if(c_sock<0) return MN_SOCK;
#endif

	sent = 0; // remex request to server
	do{n = send(c_sock, (char *)(sent+(__int64)&p), SERV_BUF_LEN-sent, 0); 
		sent += n;
		if(n<=0){closesocket(c_sock);return MN_SEND;}
	}while(sent<SERV_BUF_LEN);

	got=0; // get packet with result from remex
	do{n = recv(c_sock, (char *)(got+(__int64)&p), SERV_BUF_LEN-got, 0); 
		got += n;
		if(n<=0){closesocket(c_sock); return MN_DONE;}
   	}while(got<SERV_BUF_LEN);

	closesocket(c_sock);
	return (int)p.code.x;// success or error
}

// *************************************************************************
// *************************************************************************
// ****** RPC **************************************************************
// *************************************************************************
// *************************************************************************
/* The remote function called must be attached to a server with the flag
   MN_RPC in order for it to work.  rpc() spawns a function in a separate 
   thread with the passed argument, and the function must free() the 
   argument else there will be memory leak.  
   Something like this should be done:

   void *function(void *arg){ // Linux server
   void function(void *arg){ // Windows server
		struct my_arg_type mat;

		mat = *((struct my_arg_type *)void);
		free(arg);
		...
   atomic_rpc() runs the function without spawning a thread.  Thus, atomic_rpc()
   is blocking.  Also, the argument cannot not be freed in atomic_rpc(), and
   the argument is returned and modified by the called function to modify the
   argument passed by the client. 
*/

// *************************************************************************
// ******** mn_rpc() *******************************************************
// *************************************************************************
int mn_rpc(void *arg, __int64 len, 
		   u_short AF, u_short port, char *host){
	rpc_packet p;
	SOCKET	c_sock; // client socket 
	__int64 sent, got, bytes_sent;
	int n;
	__int64 chunk;
	void *arg_ptr;

	// prepare packet
	p.magic_number.x = MAGIC_NUMBER;
	p.len.x = len;
	p.code.x = CODE_RPC;

	if(len > RPC_SERV_DATA_LEN) memcpy(p.data, arg, RPC_SERV_DATA_LEN);
	else memcpy(p.data, arg, len);
   
	c_sock = make_sock(AF, host, port);
#ifdef WINDOWS
	if(c_sock == INVALID_SOCKET) return MN_SOCK;
#else
	if(c_sock<0) return MN_SOCK;
#endif 

	sent = 0; // rpc request to server
	do{n = send(c_sock, (char *)(sent+(__int64)&p), SERV_BUF_LEN-sent, 0); 
		sent += n;
		if(n<=0){ closesocket(c_sock); return MN_SEND;}
	}while(sent<SERV_BUF_LEN);

	// here if data needs additional sends
	if(len > RPC_SERV_DATA_LEN){// send additional data to server spawned rpc handler
		arg_ptr = (void *)((__int64)arg + RPC_SERV_DATA_LEN);
		bytes_sent = RPC_SERV_DATA_LEN;
		while(bytes_sent < len){
			if(len - bytes_sent < CHUNK_SIZE) chunk = (__int64)(len - bytes_sent);
			else chunk = CHUNK_SIZE;

   		sent = 0; // send data to server
   		do{n = send(c_sock, (char *)(sent + (__int64)bytes_sent + (__int64)arg_ptr), 
				chunk - sent, 0);
				sent += n;
				if(n<=0){ closesocket(c_sock);return MN_SEND;}
   		}while(sent < chunk);

			bytes_sent += chunk;
		}
	}

	got=0; // get packet with OK that rpc was spawned by server
	do{	n = recv(c_sock, (char *)(got+(__int64)&p), SERV_BUF_LEN-got, 0); 
		got += n;
		if(n<=0){ closesocket(c_sock); return MN_RECV;}
   }while(got<SERV_BUF_LEN);

	closesocket(c_sock);
   if(p.code.x == CODE_RPC_DONE){ return MN_SUCCESS;}
	return MN_CODE;// error
}

// *************************************************************************
// ********** mn_atomic_rpc ************************************************
// *************************************************************************
int mn_atomic_rpc(void *arg, __int64 len,
		u_short AF, u_short rem_port, char *rem_host){
	rpc_packet p;
	SOCKET	c_sock; // client socket
	__int64 sent, got;
	int n;
	
	if(len>RPC_SERV_DATA_LEN) return MN_LEN;// arg too big

	// prepare packet
	p.magic_number.x = MAGIC_NUMBER;
	p.len.x = len; 
	p.code.x = CODE_ATOMIC_RPC; 
	memcpy(p.data, arg, len);
   
	c_sock = make_sock(AF, rem_host, rem_port);
#ifdef WINDOWS
	if(c_sock == INVALID_SOCKET) return MN_SOCK;
#else
	if(c_sock<0) return MN_SOCK;
#endif 

	sent = 0; // send rpc request to server
	do{n = send(c_sock, (char *)(sent+(__int64)&p), SERV_BUF_LEN-sent, 0); 
		sent += n;
		if(n<=0){ closesocket(c_sock); return MN_SEND;}
	}while(sent<SERV_BUF_LEN);

	got=0; // get packet with data from server
	do{n = recv(c_sock, (char *)(got+(__int64)&p), SERV_BUF_LEN-got, 0); 
		got += n;
		if(n<=0){ closesocket(c_sock); return MN_RECV;}
   }while(got<SERV_BUF_LEN);

	closesocket(c_sock);
   if(p.code.x == CODE_ATOMIC_RPC_DONE){
		memcpy(arg, p.data, len);// copy back argument
		return MN_SUCCESS;
	}
	return MN_CODE; // error
}

// *************************************************************************
// ********* mn_swap() *****************************************************
// *************************************************************************
int mn_swap(__int64 rem_offset,
			void *loc_mem,
			__int64 len,
			u_short AF, u_short rem_port, char *rem_host){
	packet p;// the packet within mn_put
	SOCKET	c_sock; // client socket
	int n;// return values from sockets send and recv
	__int64 sent, got;// byte counters for socket send and recv

	if(len <= 0) return MN_LEN_UNKNOWN;// must have a length

	// load packet
	p.magic_number.x = MAGIC_NUMBER;
	p.rem_offset.x = rem_offset; 
	p.len.x = len; 
	p.chunk.x = CHUNK_SIZE;// set in mn_build.h

	// choose type of put
	if(len>SERV_DATA_LEN) return MN_2BIG;
	else{ // data is within packet
		p.code.x = CODE_SWAP; 
		memcpy(p.data, loc_mem, (__int64)p.len.x);
	}

	// make socket
	c_sock = make_sock(AF, rem_host, rem_port);
#ifdef WINDOWS
	if(c_sock == INVALID_SOCKET) return MN_SOCK;
#else
	if(c_sock<0) return MN_SOCK;
#endif 

	// send packet to server 
	sent = 0;
	do{n = send(c_sock, (char *)((__int64)sent+(__int64)&p), SERV_BUF_LEN-sent, 0);
		sent += n;
		if(n<=0){ closesocket(c_sock); return MN_SEND;}
	}while(sent<SERV_BUF_LEN);

	got = 0; 
	do{n = recv(c_sock, (char *)((__int64)got+(__int64)&p), SERV_BUF_LEN-got, 0);
		got += n;
		if(n<=0){ closesocket(c_sock); return MN_RECV;}
	}while(got < SERV_BUF_LEN);

	if(p.code.x == CODE_SWAP_DONE){	// copy back memory from server now
		memcpy(loc_mem, p.data, (__int64)p.len.x);
		closesocket(c_sock);
		return MN_SUCCESS;
	}

	closesocket(c_sock); return MN_CODE;// failure
}

// *************************************************************************
// ********* mn_ping() *****************************************************
// *************************************************************************
int mn_ping(u_short AF, u_short rem_port, char *rem_host){
	packet p;// the packet within mn_put
	SOCKET	c_sock; // client socket
	int n;// return values from sockets send and recv
	__int64 sent, got;// byte counters for socket send and recv

	// load packet
	p.magic_number.x = MAGIC_NUMBER;
	p.code.x = CODE_PING;

	// make socket
	c_sock = make_sock(AF, rem_host, rem_port);
#ifdef WINDOWS
	if(c_sock == INVALID_SOCKET) return MN_SOCK;
#else
	if(c_sock<0) return MN_SOCK;
#endif 

	// send packet to server 
	sent = 0;
	do{	n = send(c_sock, (char *)((__int64)sent+(__int64)&p), SERV_BUF_LEN-sent, 0);
		sent += n;
		if(n<=0){ closesocket(c_sock); return MN_SEND;}
	}while(sent<SERV_BUF_LEN);

	got = 0; 
	do{n = recv(c_sock, (char *)((__int64)got+(__int64)&p), SERV_BUF_LEN-got, 0);
		got += n;
		if(n<=0){ closesocket(c_sock); return MN_RECV;}
	}while(got < SERV_BUF_LEN);

	closesocket(c_sock);
	if(p.code.x == CODE_PING_OK) return MN_SUCCESS;
	return MN_FAIL;// failure
}


// *************************************************************************
// ****** mn_put() *********************************************************
// *************************************************************************
int mn_put(__int64 rem_offset, void *loc_mem, __int64 len,
			u_short AF, u_short rem_port, char *rem_host){
	packet p;// the packet within mn_put
	SOCKET	c_sock; // client socket
	int n;// return values from sockets send and recv
	__int64 sent, got, chunk;// byte counters for socket send and recv, max size data
	__int64 bytes_sent;// counter for sending bytes

	if(len <= 0) return MN_LEN_UNKNOWN;// must have a length

	// load packet
	p.magic_number.x = MAGIC_NUMBER;
	p.rem_offset.x = rem_offset; 
	p.len.x = len; 
	p.chunk.x = CHUNK_SIZE;// set in mn_build.h

	// choose type of put
	if(len>SERV_DATA_LEN) p.code.x = CODE_PUT_BIG;// send data after packet
	else{ // data is within packet
		p.code.x = CODE_PUT; 
		memcpy(p.data, loc_mem, (__int64)p.len.x);
	}

	// make socket
	c_sock = make_sock(AF, rem_host, rem_port);
#ifdef WINDOWS
	if(c_sock == INVALID_SOCKET) return MN_SOCK;
#else
	if(c_sock<0) return MN_SOCK;
#endif 

	// send packet to server 
	sent = 0;
	do{	n = send(c_sock, (char *)((__int64)sent+(__int64)&p), SERV_BUF_LEN-sent, 0);
		sent += n;
		if(n<=0){ closesocket(c_sock); return MN_SEND;}
	}while(sent<SERV_BUF_LEN);

	if(p.code.x == CODE_PUT){ // data was in packet, get response
		got = 0; 
   	do{n = recv(c_sock, (char *)((__int64)got+(__int64)&p), SERV_BUF_LEN-got, 0);
			got += n;
			if(n<=0){ closesocket(c_sock); return MN_RECV;}
   	}while(got < SERV_BUF_LEN);

		closesocket(c_sock);
   	if(p.code.x == CODE_PUT_DONE) return MN_SUCCESS;
		return MN_CODE;// failure
	}
	else if(p.code.x == CODE_PUT_BIG){// send data in chunks to server

		// send data chunks to server
		bytes_sent = 0;
		while(bytes_sent < len){
			if(len - bytes_sent < CHUNK_SIZE) chunk = (__int64)(len - bytes_sent);
			else chunk = CHUNK_SIZE;

   		sent = 0; // send data to server
   		do{n = send(c_sock, (char *)((__int64)sent + (__int64)bytes_sent + (__int64)loc_mem), 
				chunk - sent, 0);
				sent += n;
				if(n<=0){ closesocket(c_sock);return MN_SEND;}
   		}while(sent < chunk);

			bytes_sent += chunk;
		}

		// get reply from server that data was received
   	got = 0;// get reply from server that data was received
   	do{n = recv(c_sock, (char *)((__int64)got+(__int64)&p), SERV_BUF_LEN-got, 0);
			got += n;
			if(n<=0){ closesocket(c_sock);return MN_RECV;}
   	}while(got < SERV_BUF_LEN);

		closesocket(c_sock);
   	if(p.code.x == CODE_PUT_BIG_DONE) return MN_SUCCESS;
		return MN_DONE;// failure
	}
	else{ closesocket(c_sock); return MN_CODE;}// failure
}

// *************************************************************************
// ****** mn_get() *********************************************************
// *************************************************************************
int mn_get(__int64 rem_offset, 
		   void *loc_mem, 
		   __int64 len,
			u_short AF, u_short rem_port, char *rem_host){
	packet p;
	SOCKET	c_sock; // client socket
	__int64 sent, got;
	int n;
	__int64 chunk;
	__int64 bytes_recvd;
	
	if(len <= 0) return MN_LEN_UNKNOWN;

	// prepare packet
	p.magic_number.x = MAGIC_NUMBER;
	p.rem_offset.x = rem_offset; 
	p.len.x = len; 
	p.chunk.x = CHUNK_SIZE;// this should be set in mn_build.h

	if(len>SERV_DATA_LEN) p.code.x = CODE_GET_BIG; 
	else p.code.x = CODE_GET;
   
	c_sock = make_sock(AF, rem_host, rem_port);
#ifdef WINDOWS
	if(c_sock == INVALID_SOCKET) return MN_SOCK;
#else
	if(c_sock<0) return MN_SOCK;
#endif

	sent = 0; // send to server to xfer data
	do{ n = send(c_sock, (char *)(sent+(__int64)&p), SERV_BUF_LEN-sent, 0); 
		sent += n;
		if(n<=0){ closesocket(c_sock);return MN_SEND;}
	}while(sent<SERV_BUF_LEN);

	if(p.code.x == CODE_GET_BIG){// data coming in

		// recv the data now
		bytes_recvd = 0;
		while(bytes_recvd < len){
			if(len - bytes_recvd < CHUNK_SIZE) chunk = (__int64)(len - bytes_recvd);
			else chunk = CHUNK_SIZE;

			got=0;
			do{ n=recv(c_sock,(char *)(got + (__int64)bytes_recvd + (__int64)loc_mem), 
				chunk - got, 0);
				got += n;
				if(n<=0){ closesocket(c_sock); return MN_RECV;}
   		}while(got < chunk);

			bytes_recvd += chunk;
		}
				
		closesocket(c_sock);return MN_SUCCESS;
	}
	else{// GET : a packet back with data inside
   		
		got=0;
		do{ n = recv(c_sock, (char *)(got+(__int64)&p), SERV_BUF_LEN-got, 0); 
			got += n;
			if(n<=0){ closesocket(c_sock); return MN_RECV;}
   	}while(got<SERV_BUF_LEN);

		memcpy(loc_mem, p.data, (__int64)p.len.x);// copy data over
     	closesocket(c_sock); return MN_SUCCESS;
	}
}


// *************************************************************************
// *************************************************************************
// *************************************************************************
// ************ Client Disk Side Calls *************************************
// *************************************************************************
// *************************************************************************
// *************************************************************************


// ******* buffer management: client disk-side calls ***********************
/* This buffer is shared by disk-side client routines.  If it is busy, then
   a disk-side call (put_d2d, put_d2m, get_d2d, get_m2d)
	allocates a bufer using malloc. 

   The mutex lock is for linux. Windows uses InterlockedIncrement() and
	InterlockedDecrement(), simple & handy calls for multiple threads.
*/
	static char buf_stat[CHUNK_SIZE];//  buffer for reading and sending data
	static int buf_busy = 0;// 0 if not busy, >= 1 if busy; 
#ifdef LINUX
	pthread_mutex_t mutex1 = PTHREAD_MUTEX_INITIALIZER;
#endif

// *************************************************************************
// **********  mn_put_disk2disk() ******************************************
// *************************************************************************
int mn_put_disk2disk(char *rem_disk, __int64 rem_offset, 
					 char *loc_disk, __int64 loc_offset, 
					 __int64 len, 
					 u_short AF, u_short rem_port, char *rem_host){
#ifdef WINDOWS
	struct _stati64 statbuf;
#else
	struct stat statbuf;
#endif
	__int64 left_to_send, chunk_size, sent, got;
	int	fin, ret, n;
	disk_packet p;
	SOCKET c_sock;
	char *buf;//  buffer for reading and sending data
	int error_code, dynamic_buf = 0;// ,0 if using static buffer, 1 if malloc
	
// setup buffers
#ifdef WINDOWS
	if(InterlockedIncrement(&buf_busy) > 1){
		InterlockedDecrement(&buf_busy);
		dynamic_buf = 1;
		buf = (char *)malloc(CHUNK_SIZE);
	}
	else buf = buf_stat;// use static buffer
#else // linux
	pthread_mutex_lock(&mutex1);// ******* lock ********
	++buf_busy;
	if(buf_busy > 1){
		--buf_busy;
		dynamic_buf = 1;
		buf = (char *)malloc(CHUNK_SIZE);
	}
	else buf = buf_stat;
	pthread_mutex_unlock(&mutex1);// ****** unlock ******
#endif
	
	// open the file and find size if necessary
	fin = _open(loc_disk, _O_RDONLY | _O_BINARY );
	if(fin < 0){ error_code = MN_OPEN;goto exit_put_d2d2;}

	p.code.x = CODE_DPUT;
	if(len == 0){
#ifdef WINDOWS
		ret = _fstati64(fin, &statbuf);
#else
		ret = fstat64(fin, &statbuf);
#endif
		if(ret < 0){ error_code = MN_FSTAT;goto exit_put_d2d1;}
		len = (__int64)statbuf.st_size;
		p.code.x = CODE_DPUT_DEL;// delete target disk
	}

	// load packet
	p.magic_number.x = MAGIC_NUMBER; 
	p.len.x = len;
	strcpy(p.rem_file, rem_disk); 
	p.offset.x = rem_offset; 
	p.chunk.x = CHUNK_SIZE;

	c_sock = make_sock(AF, rem_host, rem_port);
	if(c_sock < 0){ error_code = MN_SOCK;goto exit_put_d2d1;}

	sent = 0;// request to put data
	do{n = send(c_sock, (char *)(sent+(__int64)&p), SERV_BUF_LEN-sent, 0);
		sent += n;
		if(n<=0){ error_code = MN_SEND;goto exit_put_d2d0;}
	}while(sent<SERV_BUF_LEN);

	ret = (int)_lseeki64(fin, loc_offset, SEEK_SET);
	if(ret<0){ error_code = MN_LSEEK;goto exit_put_d2d0;}

	got = 0;// get response (all set to recv at rem host)
	do{n = recv(c_sock, (char *)(got+(__int64)&p), SERV_BUF_LEN-got, 0);
		got += n;
		if(n<=0){ error_code = MN_RECV;goto exit_put_d2d0;}
   }while(got < SERV_BUF_LEN);

	left_to_send = len; // init
	while(left_to_send > 0){ // ***** data send loop *******************

		if(left_to_send < CHUNK_SIZE ) chunk_size = (__int64)left_to_send;
		else chunk_size = CHUNK_SIZE;

		n = _read(fin, buf, chunk_size);
		if(n != (int)chunk_size){ error_code = MN_READ;goto exit_put_d2d0;}

		// send data
		sent = 0;// request to put data
		do{n = send(c_sock, (char *)(sent+(__int64)buf), chunk_size-sent, 0);
			sent += n;
			if(n<=0){ error_code = MN_SEND; goto exit_put_d2d0;}
		}while(sent<chunk_size);

		left_to_send -= chunk_size;// this is sending, or is gone
	}// ********************************** data send loop *****************
	_close(fin);
	if(dynamic_buf) free(buf); 
	else {
#ifdef WINDOWS
		InterlockedDecrement(&buf_busy);
#else // Linux
		pthread_mutex_lock(&mutex1);// ******* lock ********
		--buf_busy;
		pthread_mutex_unlock(&mutex1);// ****** unlock ******
#endif
	}

	got = 0;// get response done
	do{n = recv(c_sock, (char *)(got+(__int64)&p), SERV_BUF_LEN-got, 0);
		got += n;
		if(n<=0){ closesocket(c_sock); return MN_DONE;}
   	}while(got < SERV_BUF_LEN);

	if(p.code.x != CODE_DPUT_DONE){ closesocket(c_sock); return MN_CODE;}
	closesocket(c_sock);  
	return MN_SUCCESS;

// *** cleanup routines ***
exit_put_d2d0:
	closesocket(c_sock);
exit_put_d2d1:
	_close(fin);
exit_put_d2d2:
	if(dynamic_buf) free(buf); 
	else {
#ifdef WINDOWS
		InterlockedDecrement(&buf_busy);
#else // Linux
		pthread_mutex_lock(&mutex1);// ******* lock ********
		--buf_busy;
		pthread_mutex_unlock(&mutex1);// ****** unlock ******
#endif
	}
	return error_code;
}

// *************************************************************************
// *********** mn_put_disk2mem() *******************************************
// *************************************************************************
int mn_put_disk2mem(__int64 rem_offset,// from base of server memory 
					 char *loc_disk, __int64 loc_offset, 
					 __int64 len, 
					 u_short AF, u_short rem_port, char *rem_host){
#ifdef WINDOWS
	struct _stati64 statbuf;
#else
	struct stat statbuf;
#endif
	__int64 left_to_send, chunk_size, sent, got;
	int n, fin, ret, error_code;
	packet p;
	SOCKET c_sock;
	char *buf;//  buffer for reading and sending data
	int dynamic_buf = 0;// 0 if using static buffer, 1 if malloc
	
	// *********************
	// setup buffers
	#ifdef WINDOWS
	if(InterlockedIncrement(&buf_busy) > 1){
		InterlockedDecrement(&buf_busy);
		dynamic_buf = 1;
		buf = (char *)malloc(CHUNK_SIZE);
	}
	else buf = buf_stat;// use static buffer
#else // linux
	pthread_mutex_lock(&mutex1);// ******* lock ********
	++buf_busy;
	if(buf_busy > 1){
		--buf_busy;
		dynamic_buf = 1;
		buf = (char *)malloc(CHUNK_SIZE);
	}
	else buf = buf_stat;
	pthread_mutex_unlock(&mutex1);// ****** unlock ******
#endif

	if(len <= 0) return MN_LEN_UNKNOWN;

	// open the file and find size if necessary
	fin = _open(loc_disk, _O_RDONLY | _O_BINARY );
	if(fin < 0){ error_code = MN_OPEN;goto exit_put_d2m2;} 
	if(len == 0){
#ifdef WINDOWS
		ret = _fstati64(fin, &statbuf);
#else
		ret = fstat64(fin, &statbuf);
#endif
		if(ret < 0){ error_code = MN_FSTAT;goto exit_put_d2m1;}
		len = (__int64)statbuf.st_size;
	}

	// load packet
	p.magic_number.x = MAGIC_NUMBER;
	p.code.x = CODE_PUT_BIG; 
	p.len.x = len; 
	p.rem_offset.x = rem_offset; 
	p.chunk.x = CHUNK_SIZE;

	c_sock = make_sock(AF, rem_host, rem_port);
	if(c_sock < 0){ error_code = MN_SOCK;goto exit_put_d2m1;}

	sent = 0;// request to put data
	do{n = send(c_sock, (char *)(sent+(__int64)&p), SERV_BUF_LEN-sent, 0);
		sent += n;
		if(n<=0){ error_code = MN_SEND;goto exit_put_d2m0;}
	}while(sent<SERV_BUF_LEN);

	ret = (int)_lseeki64(fin, loc_offset, SEEK_SET);
	if(ret<0){ error_code = MN_LSEEK;goto exit_put_d2m0;}

	left_to_send = len; // init
	// *****************************************************************
	while(left_to_send > 0){ // ***** data send loop *******************

		// adjust packet size to remote mapped memory needs
		if(left_to_send < CHUNK_SIZE ) chunk_size = (__int64)left_to_send;
		else chunk_size = CHUNK_SIZE;

		n = _read(fin, buf, chunk_size);
		if(n != (int)chunk_size){ error_code = MN_READ;goto exit_put_d2m0;}

		// send data
		sent = 0;// request to put data
		do{n = send(c_sock, (char *)(sent+(__int64)buf), chunk_size-sent, 0);
			sent += n;
			if(n<=0){ error_code = MN_SEND;goto exit_put_d2m0;}
		}while(sent<chunk_size);

		left_to_send -= chunk_size;// this is sending, or is gone
	}// ********************************** data send loop *****************
	_close(fin);
	if(dynamic_buf) free(buf); 
	else {
#ifdef WINDOWS
		InterlockedDecrement(&buf_busy);
#else // Linux
		pthread_mutex_lock(&mutex1);// ******* lock ********
		--buf_busy;
		pthread_mutex_unlock(&mutex1);// ****** unlock ******
#endif
	}
	got = 0;// get response done
	do{n = recv(c_sock, (char *)(got+(__int64)&p), SERV_BUF_LEN-got, 0);
		got += n;
		if(n<=0){ closesocket(c_sock); return MN_DONE;}
   	}while(got < SERV_BUF_LEN);

	if(p.code.x != CODE_PUT_BIG_DONE){ closesocket(c_sock);return MN_CODE;}
	closesocket(c_sock);  
	return MN_SUCCESS;

	// *** cleanup routines ***
exit_put_d2m0:
	closesocket(c_sock);
exit_put_d2m1:
	_close(fin);
exit_put_d2m2:
	if(dynamic_buf) free(buf); 
	else {
#ifdef WINDOWS
		InterlockedDecrement(&buf_busy);
#else // Linux
		pthread_mutex_lock(&mutex1);// ******* lock ********
		--buf_busy;
		pthread_mutex_unlock(&mutex1);// ****** unlock ******
#endif
	}
	return error_code;
}

// *************************************************************************
// ****** mn_put_mem2disk() ************************************************
// *************************************************************************
int mn_put_mem2disk(char *rem_disk, __int64 rem_offset, 
					void *loc_mem, 
					__int64 len, 
					u_short AF, u_short rem_port, char *rem_host){
	
	disk_packet p;
	__int64 left_to_send; 
	int n;
	__int64 sent, got, chunk_size, bytes_sent;
	SOCKET c_sock;

	if(len <= 0) return MN_LEN_UNKNOWN;

	// load packet for remote server
	p.magic_number.x = MAGIC_NUMBER;
	p.code.x = CODE_DPUT; p.len.x = len;
	strcpy(p.rem_file, rem_disk);
	p.offset.x = rem_offset;
	p.chunk.x = CHUNK_SIZE;

	c_sock = make_sock(AF, rem_host, rem_port);
#ifdef WINDOWS
	if(c_sock == INVALID_SOCKET) return MN_SOCK;
#else
	if(c_sock<0) return MN_SOCK;
#endif  

	sent = 0;// request to put data 
	do{n = send(c_sock, (char *)(sent+(__int64)&p), SERV_BUF_LEN-sent, 0);
		sent += n;
		if(n<=0){closesocket(c_sock); return MN_SEND;}
	}while(sent<SERV_BUF_LEN);

	got = 0;// ack
	do{n = recv(c_sock, (char *)(got+(__int64)&p), SERV_BUF_LEN-got, 0);
		got += n;
		if(n<=0){closesocket(c_sock); return MN_RECV; }
	}while(got < SERV_BUF_LEN);

	bytes_sent = 0;
	left_to_send = len;// the entire length of data to send
	// *****************************************************************
	while(left_to_send > 0){// ********* send data loop ********************************

		// adjust packet size to remote mapped memory needs
		if(left_to_send < CHUNK_SIZE ) chunk_size = (__int64)left_to_send;
		else chunk_size = CHUNK_SIZE;

		sent = 0; // send data 
		do{n = send(c_sock,(char *)(sent + (__int64)bytes_sent + (__int64)loc_mem), 
				chunk_size - sent, 0);
			sent += n;
			if(n<=0){closesocket(c_sock); return MN_SEND;}
		}while(sent<chunk_size);
		bytes_sent += chunk_size;
		
		left_to_send -= chunk_size;// decrement by data sent
	}// *********** end send loop ****************************************

	got = 0;// get response done
	do{n = recv(c_sock, (char *)(got+(__int64)&p), SERV_BUF_LEN-got, 0);
		got += n;
   		if(n<=0){closesocket(c_sock); return MN_DONE;}
   	}while(got < SERV_BUF_LEN);
	if(p.code.x != CODE_DPUT_DONE){closesocket(c_sock); return MN_CODE;}

	closesocket(c_sock); return MN_SUCCESS;
}

// *************************************************************************
// *********** mn_get_disk2disk() ******************************************
// *************************************************************************
int mn_get_disk2disk(char *rem_disk, __int64 rem_offset, 
					char *loc_disk, __int64 loc_offset, 
					__int64 len, 
					u_short AF, u_short rem_port, char *rem_host){

	int fout, ret, n, error_code;
	__int64 chunk_size;
	__int64 sent, got, left_to_recv;
	disk_packet p;			// info packet to send
	SOCKET	c_sock;			// client socket
	
	char *buf;//  buffer for reading and sending data
	int dynamic_buf = 0;// 0 if using static buffer, 1 if malloc
	
	// *********************
	// setup buffers
	#ifdef WINDOWS
	if(InterlockedIncrement(&buf_busy) > 1){
		InterlockedDecrement(&buf_busy);
		dynamic_buf = 1;
		buf = (char *)malloc(CHUNK_SIZE);
	}
	else buf = buf_stat;// use static buffer
#else // linux
	pthread_mutex_lock(&mutex1);// ******* lock ********
	++buf_busy;
	if(buf_busy > 1){
		--buf_busy;
		dynamic_buf = 1;
		buf = (char *)malloc(CHUNK_SIZE);
	}
	else buf = buf_stat;
	pthread_mutex_unlock(&mutex1);// ****** unlock ******
#endif
	
	// remove disk if len = 0
	if(len == 0) remove(loc_disk);

	// prepare packet
	p.magic_number.x = MAGIC_NUMBER;
	p.offset.x = rem_offset; 
	p.len.x = len;
	strcpy(p.rem_file, rem_disk);
	p.code.x = CODE_DGET;
	p.chunk.x = CHUNK_SIZE;

	// make socket for connection and data xfer
	c_sock = make_sock(AF, rem_host, rem_port);
	if(c_sock<0) return MN_SOCK;

	// send to server to start a thread
	// if sent len = 0, we get back the length of the file

	sent = 0; 
	do{n = send(c_sock, (char *)(sent+(__int64)&p), SERV_BUF_LEN-sent, 0); 
		sent += n;
		if(n<=0){  closesocket(c_sock);error_code = MN_SEND;goto exit_get_d2d2;}
	}while(sent<SERV_BUF_LEN);

	got = 0;// get response, len value in common 
	do{n = recv(c_sock, (char *)(got+(__int64)&p), SERV_BUF_LEN - got, 0);
		got += n;
		if(n<=0){closesocket(c_sock); error_code = MN_RECV;goto exit_get_d2d2;}
	}while(got < SERV_BUF_LEN );

	// check for server abort
	if(p.code.x == CODE_ABORT){ closesocket(c_sock); return MN_OPEN;}

	len = p.len.x;// now set len in common if len was 0 initially

	// open file and position pointer
	fout = _open(loc_disk, _O_WRONLY | _O_CREAT | _O_BINARY, _S_IWRITE);	
	if(fout < 0){closesocket(c_sock); error_code = MN_OPEN;goto exit_get_d2d2;}
	ret = (int)_lseeki64(fout, loc_offset, SEEK_SET);
	if(ret<0){error_code = MN_LSEEK;goto exit_get_d2d0;}

	left_to_recv = len;// amount of data to send, and decrement

	while(left_to_recv > 0){// **** data send loop ***************************

		// adjust data chunk, trim if under default chunk size
		if(left_to_recv < CHUNK_SIZE) chunk_size = (int)left_to_recv;
		else chunk_size = CHUNK_SIZE;

		got=0; // receive the remote file data
		do{n = recv(c_sock,(char *)(got + (__int64)buf),chunk_size - got, 0); 
			got += n;
			if(n<=0){error_code = MN_RECV;goto exit_get_d2d0;}
		}while(got < chunk_size);

		// write the data to disk
		ret = _write(fout, buf, chunk_size);
		if(ret < (int)chunk_size){error_code = MN_WRITE;goto exit_get_d2d0;}

		left_to_recv -= chunk_size;// decrement
	}// *************************** end data send loop ************************
	
	error_code = MN_SUCCESS;

// *** error handling or closure on success ***
exit_get_d2d0:
	closesocket(c_sock);
//exit_get_d2d1:
	_close(fout);
exit_get_d2d2:
	if(dynamic_buf) free(buf); 
	else {
#ifdef WINDOWS
		InterlockedDecrement(&buf_busy);
#else // Linux
		pthread_mutex_lock(&mutex1);// ******* lock ********
		--buf_busy;
		pthread_mutex_unlock(&mutex1);// ****** unlock ******
#endif
	}
	return error_code;
}

// *************************************************************************
// ***** mn_get_mem2disk() *************************************************
// *************************************************************************
int mn_get_mem2disk(__int64 rem_offset, // offset from memory base 
					char *loc_disk, __int64 loc_offset, 
					__int64 len, 
					u_short AF, u_short rem_port, char *rem_host){

	int fout, ret, n, error_code;
	__int64 chunk_size;
	__int64 sent, got, left_to_recv;	
	packet p;			// info packet to send
	SOCKET	c_sock;			// client socket
	char *buf;//  buffer for reading and sending data
	int dynamic_buf = 0;// 0 if using static buffer, 1 if malloc
	
	// setup buffers
	#ifdef WINDOWS
	if(InterlockedIncrement(&buf_busy) > 1){
		InterlockedDecrement(&buf_busy);
		dynamic_buf = 1;
		buf = (char *)malloc(CHUNK_SIZE);
	}
	else buf = buf_stat;// use static buffer
#else // linux
	pthread_mutex_lock(&mutex1);// ******* lock ********
	++buf_busy;
	if(buf_busy > 1){
		--buf_busy;
		dynamic_buf = 1;
		buf = (char *)malloc(CHUNK_SIZE);
	}
	else buf = buf_stat;
	pthread_mutex_unlock(&mutex1);// ****** unlock ******
#endif

	if(len <= 0) return MN_LEN_UNKNOWN;// must know length!

	// prepare packet
	p.magic_number.x = MAGIC_NUMBER;
	p.rem_offset.x = rem_offset; 
	p.len.x = len;
	p.code.x = CODE_GET_BIG;
	p.chunk.x = CHUNK_SIZE;

	// make socket for connection and xfers
	c_sock = make_sock(AF, rem_host, rem_port);
	if(c_sock<0) return MN_SOCK;

	// send to server to start a thread
	sent = 0; 
	do{n = send(c_sock, (char *)(sent+(__int64)&p), SERV_BUF_LEN-sent, 0); 
		sent += n;
		if(n<=0){closesocket(c_sock);error_code = MN_SEND;goto exit_get_m2d2;}
	}while(sent<SERV_BUF_LEN);

	// open the file and position to the offset
	fout = _open(loc_disk, _O_WRONLY | _O_CREAT | _O_BINARY, _S_IWRITE);	
	if(fout < 0){closesocket(c_sock); error_code = MN_OPEN;goto exit_get_m2d2;}
	ret = (int)_lseeki64(fout, loc_offset, SEEK_SET);
	if(ret<0){error_code = MN_LSEEK;goto exit_get_m2d0;}

	left_to_recv = len;// amount of data to send, and decrement
	while(left_to_recv > 0){// **** data send loop ***************************

		// adjust data chunk, trim if under default chunk size
		if(left_to_recv < CHUNK_SIZE) chunk_size = (__int64)left_to_recv;
		else chunk_size = CHUNK_SIZE;

		got=0; // receive the remote file data
		do{n = recv(c_sock,(char *)(got + (__int64)buf),chunk_size - got, 0); 
			got += n;
			if(n<=0){error_code = MN_RECV;goto exit_get_m2d0;}
		}while(got < chunk_size);

		// write the data to disk
		ret = _write(fout, buf, chunk_size);
		if(ret < (int)chunk_size){error_code = MN_WRITE;goto exit_get_m2d0;}

		left_to_recv -= chunk_size;// decrement
	}// *************************** end data send loop ************************
	
	error_code = MN_SUCCESS;

	// *** error handling or closure on success ***
exit_get_m2d0:
	closesocket(c_sock);
//exit_get_m2d1:
	_close(fout);
exit_get_m2d2:
	if(dynamic_buf) free(buf); 
	else {
#ifdef WINDOWS
		InterlockedDecrement(&buf_busy);
#else // Linux
		pthread_mutex_lock(&mutex1);// ******* lock ********
		--buf_busy;
		pthread_mutex_unlock(&mutex1);// ****** unlock ******
#endif
	}
	return error_code;
}

// *************************************************************************
// ****** mn_get_disk2mem() ************************************************
// *************************************************************************
int mn_get_disk2mem(char *rem_disk, __int64 rem_offset, 
					void *loc_mem, 
					__int64 len,
					u_short AF, u_short rem_port, char *rem_host){

	// local variables
	disk_packet p;
	SOCKET	c_sock; // client socket
	__int64 sent, got, left_to_send;
	int n;
	__int64 chunk_size, received;

	// prepare packet
	p.magic_number.x = MAGIC_NUMBER;
	p.offset.x = rem_offset; p.len.x = len;
	strcpy(p.rem_file, rem_disk);
	p.code.x = CODE_DGET;
	p.chunk.x = CHUNK_SIZE;

	if(len <= 0) return MN_LEN_UNKNOWN;// we have to know this

	c_sock = make_sock(AF, rem_host, rem_port);
#ifdef WINDOWS
	if(c_sock == INVALID_SOCKET) return MN_SOCK;
#else
	if(c_sock<0) return MN_SOCK;
#endif

	sent = 0; // send to server to get file data
	do{n = send(c_sock, (char *)(sent+(__int64)&p), SERV_BUF_LEN-sent, 0); 
		sent += n;
		if(n<=0){ closesocket(c_sock);return MN_SEND;}
	}while(sent<SERV_BUF_LEN);

	got = 0;// get response, chunk values in common 
	do{n = recv(c_sock, (char *)(got+(__int64)&p), SERV_BUF_LEN - got, 0);
		got += n;
		if(n<=0){closesocket(c_sock);return MN_RECV;}
	}while(got < SERV_BUF_LEN );

	// check for server abort
	if(p.code.x == CODE_ABORT){ closesocket(c_sock); return MN_OPEN;}

	len = p.len.x;

	left_to_send = len;
	received = 0;
	while(left_to_send){
		// adjust size based on remainder to recv
		if(left_to_send < CHUNK_SIZE) chunk_size = (__int64)left_to_send;
		else chunk_size = CHUNK_SIZE;		

		got=0; // receive the data as fast as it can be sent
		do{n=recv(c_sock,(char *)(got + (__int64)received + (__int64)loc_mem), chunk_size - got, 0); 
			got += n;
			if(n<=0){ closesocket(c_sock); return MN_RECV; }
   		}while(got < chunk_size);
		received += chunk_size;

		left_to_send = left_to_send - chunk_size;
	} // **********************************************************
	closesocket(c_sock);
	return MN_SUCCESS;
}


/* ****** for use in server
  The structure is used by all threads, but only on startup for
  access to s_sock and addr.  The thread_id and protocol fields are only
  for information and obtained on startup.
*/

typedef struct{
	SOCKET *s_sock;
	void *addr;
	int thread_id;
	int protocol;
} serv_pack; 

// *************************************************************************
// *************************************************************************
// *************************************************************************
// ******       mn_serv_thread  <<called by start_serv>>       *************
// *************************************************************************
// *************************************************************************
// *************************************************************************
#ifdef WINDOWS
void serv_thread(void *inpack){
#else
void *serv_thread(void *inpack){
#endif

#ifdef WINDOWS
	void (*proc)(void *);// address of function for rpc
#else
	pthread_t thread;
	void *(*proc)(void *);// address of function for rpc
#endif
	void *addr;	
	int thread_id;// just a number
	int protocol;// 2 or 26, for debugging display
	int run = 1;// flag to keep looping server

	SOCKET s_sock, r_sock;// server and remote sockets
	int n, ret, fin, fout;
	__int64 sent, got, bytes_got, bytes_sent, chunk_size, left_to_send;
	char buf[SERV_BUF_LEN], swap_buf[SERV_DATA_LEN];
	packet *pp;// put and get and stop and kill
	rpc_packet *rp; // rpc
	disk_packet *dp; // disk
	void *rtdata;// for rpc data
	char *DB;// for all disk server functions
#ifdef WINDOWS
	struct _stati64 statbuf;
#else
	struct stat statbuf;
#endif
	serv_pack *sp;// info is passed to this thread with this

	sp = (serv_pack *)inpack;// access data
	s_sock = *(sp->s_sock);// copy socket over 
	addr = sp->addr;
	thread_id = sp->thread_id;
	protocol = sp->protocol;
	free(sp);// now get rid of it.  It is not shared.

	pp = (packet *)buf;// pointer to packet buffer
	DB = (char *)malloc(CHUNK_SIZE);// for disk side transfer buffer
	
	// Wait for an incoming request 
	// **********************************
	// ***** main loop ******************
	// **********************************
	while(run){
		serv_recv_loop: // goto destination ***** GOTO GOTO serv_recv_loop 
		printf("server # %d (AF = %d) waiting on accept:\n", thread_id, protocol);//PPPPPPPPP
		r_sock = accept(s_sock, NULL, NULL);// r_sock created
#ifdef WINDOWS
		if(r_sock == INVALID_SOCKET) goto serv_recv_loop;
#else
		if(r_sock<0) goto serv_recv_loop;
#endif

//		fcntl(r_sock, F_SETFD, FD_CLOEXEC);// to carry socket over

		n = recv(r_sock, buf, SERV_BUF_LEN, 0); 
		
		if(n != SERV_BUF_LEN){  closesocket(r_sock); goto serv_recv_loop; }
#ifdef SECURITY
		if(pp->magic_number.x != MAGIC_NUMBER){ closesocket(r_sock); goto serv_recv_loop; }
#endif


		switch (pp->code.x){// ********** Handle received packet *******************
		case CODE_STOP:// stop sockets server ************************************
			if(thread_id == 0){// only stop if not a spawned thread
				closesocket(r_sock);closesocket(s_sock);free(DB);run = 0;}
			break;
		case CODE_DIE: // terminate the process **********************************
			// without these commented out, the server "rolls" in a funky manner
//			closesocket(r_sock);
//			closesocket(s_sock);
//			mn_sockets_end();// try to close something gracefully
			free(addr);// free up the memory allocated
			free(DB);// free the buffer allocated
			exit(0);// die all
//			break;
		case CODE_PING:// check on server
			pp->code.x = CODE_PING_OK;
			sent = 0;
			do{n = send(r_sock, (char *)(sent + (__int64)buf), 
						SERV_BUF_LEN-sent, 0); 
				sent += n;
				if(n<=0){closesocket(r_sock); goto serv_recv_loop;}
			}while(sent<SERV_BUF_LEN);
			closesocket(r_sock);
			break;	
		case CODE_PUT: // *** mem 2 mem ******************************************
			if(NULL == addr){closesocket(r_sock); break;} 

			pp->addr.x = (char *)( (__int64)addr + pp->rem_offset.x );

			// put the data to memeory
			memcpy(pp->addr.x, pp->data, pp->len.x);
			
			pp->code.x = CODE_PUT_DONE;
			sent = 0;
			do{n = send(r_sock, (char *)(sent + (__int64)buf), 
						SERV_BUF_LEN-sent, 0); 
				sent += n;
				if(n<=0){closesocket(r_sock); goto serv_recv_loop;}
			}while(sent<SERV_BUF_LEN);
			closesocket(r_sock);

        	break;
		case CODE_SWAP:
			if(NULL == addr){closesocket(r_sock); break;} 
			pp->addr.x = (char *)( (__int64)addr + pp->rem_offset.x );

			// put the data to memory
			memcpy(swap_buf, pp->addr.x, pp->len.x);//move mem to temp buf
			memcpy(pp->addr.x, pp->data, pp->len.x);// move packet data to mem
			memcpy(pp->data, swap_buf, pp->len.x);//move mem(in swap_buf) to packet data
			
			pp->code.x = CODE_SWAP_DONE;
			sent = 0;
			do{n = send(r_sock, (char *)(sent + (__int64)buf), 
						SERV_BUF_LEN-sent, 0); 
				sent += n;
				if(n<=0){closesocket(r_sock); goto serv_recv_loop;}
			}while(sent<SERV_BUF_LEN);
			closesocket(r_sock);
			break;
		case CODE_GET: // *** mem 2 mem *****************************************
			if(NULL == addr){closesocket(r_sock); break;} 
			pp->addr.x = (char *)( (__int64)addr + pp->rem_offset.x ) ;

			// get data from memory
			memcpy(pp->data, pp->addr.x, pp->len.x);
			
			sent = 0;
			do{n = send(r_sock, (char *)(sent+(__int64)buf), 
						SERV_BUF_LEN-sent, 0); sent += n;
				if(n<=0){ closesocket(r_sock); goto serv_recv_loop;}
			}while(sent<SERV_BUF_LEN);
			closesocket(r_sock);
			break;
		case CODE_PUT_BIG: // *** mem 2 mem *************************************
			if(NULL == addr){closesocket(r_sock); break;} 
			
			//turn offset into address for put() and get() ops
			pp->addr.x = (char *)( (__int64)addr + pp->rem_offset.x);

			bytes_got = 0;
			while(bytes_got < pp->len.x){ // get remote data in chunks as defined remotely
				if(pp->len.x - bytes_got < CHUNK_SIZE) chunk_size = (__int64)pp->len.x - bytes_got;
				else chunk_size = CHUNK_SIZE;

				got = 0;// recv the data from client 
				do{n = recv(r_sock, (char *)(bytes_got + (__int64)got + (__int64)pp->addr.x ), 
					chunk_size - got, 0);
					if(n<=0){  closesocket(r_sock); goto serv_recv_loop;}
					got += n;
				}while(got < chunk_size);
				bytes_got += chunk_size;
			}

			pp->code.x = CODE_PUT_BIG_DONE;// tell client all was put in mem
			sent = 0;  
			do{	n = send(r_sock,(char *)(sent + (__int64)pp), SERV_BUF_LEN-sent, 0);
				if(n<=0){  closesocket(r_sock); goto serv_recv_loop;}
				sent += n;
			}while(sent<SERV_BUF_LEN);
			closesocket(r_sock); 
        	break;
		case CODE_GET_BIG:// *** mem 2 mem **************************************
			if(NULL == addr){closesocket(r_sock); break;}
			pp->addr.x = (char *)( (__int64)addr + pp->rem_offset.x );

			bytes_sent = 0;
			while(bytes_sent < pp->len.x){
				if(pp->len.x - bytes_sent < CHUNK_SIZE) chunk_size = (__int64)pp->len.x - bytes_sent;
				else chunk_size = CHUNK_SIZE;
		
				sent = 0; // send data to client
				do{n = send(r_sock, (char *)(bytes_sent + (__int64)sent + (__int64)pp->addr.x ), 
					chunk_size - sent, 0);
					sent += n;
					if(n<=0){  closesocket(r_sock);goto serv_recv_loop;}
				}while(sent < chunk_size);
				bytes_sent += chunk_size;
			}
			closesocket(r_sock); 
        	break;
		case CODE_RPC:// *** rpc ************************************************
			// rpc called function must handle memory free() itself
			if(NULL == addr){closesocket(r_sock); break;} 
			rp = (rpc_packet *)buf;
			rp->addr.x = addr;

			rtdata = (void *)malloc((__int64)rp->len.x);// rpc function must free() it
			if(rp->len.x <= RPC_SERV_DATA_LEN){
				memcpy(rtdata, rp->data, rp->len.x);// take the data in packet
			}
			else{
				memcpy(rtdata, rp->data, RPC_SERV_DATA_LEN);// take the data in packet

				// now rest of data comes in
				bytes_got = RPC_SERV_DATA_LEN;
				while(bytes_got < rp->len.x){ // get remote data in chunks as defined remotely
					if(rp->len.x - bytes_got < CHUNK_SIZE) chunk_size = (__int64)rp->len.x - bytes_got;
					else chunk_size = CHUNK_SIZE;

					got = 0;// recv the data from client 
					do{n = recv(r_sock, (char *)(bytes_got + (__int64)got + (__int64)rtdata), 
						chunk_size - got, 0);
						if(n<=0){ closesocket(r_sock); goto serv_recv_loop;}
						got += n;
					}while(got < chunk_size);
					bytes_got += chunk_size;
				}
			}

#ifdef WINDOWS
			proc = (void (*)(void *)) rp->addr.x;
			_beginthread(proc, 4096, rtdata);
#else
			proc = (void* (*)(void *)) rp->addr.x;
			pthread_create(&thread, NULL, proc, rtdata);
			pthread_detach(thread);
			sched_yield();
#endif

			rp->code.x = CODE_RPC_DONE;// send done to client
			sent = 0;
			do{	n = send(r_sock , (char *)(sent + (__int64)rp), SERV_BUF_LEN - sent, 0);
				sent += n;
				if(n<=0){ closesocket(r_sock); goto serv_recv_loop;}
			} while ( sent < SERV_BUF_LEN );
			closesocket(r_sock);
			break;
		case CODE_ATOMIC_RPC:// *** atomic rpc ********************************
			if(NULL == addr){closesocket(r_sock); break;} 
			rp = (rpc_packet *)buf;

#ifdef WINDOWS
			proc = (void (*)(void *)) addr;
#else
			proc = (void* (*)(void *)) addr;
#endif
			proc(rp->data); // run process

			rp->code.x = CODE_ATOMIC_RPC_DONE;// send done to client
			sent = 0;
			do{n = send(r_sock, (char *)(sent+(__int64)buf), SERV_BUF_LEN-sent, 0); sent += n;
				if(n<=0){closesocket(r_sock); goto serv_recv_loop;}
			}while(sent<SERV_BUF_LEN);
			closesocket(r_sock);

			break;
		case CODE_DPUT:// *** disk or mem 2 disk **********************************
		case CODE_DPUT_DEL://delete target disk
			dp = (disk_packet *)buf;// create mem for xfer

			// open the file
			if(dp->code.x == CODE_DPUT_DEL) remove(dp->rem_file);// remove if len=0 in client
			fout = _open(dp->rem_file, _O_WRONLY | _O_CREAT | _O_BINARY, _S_IWRITE);
			if(fout < 0) {  closesocket(r_sock);goto serv_recv_loop;}
	
			sent = 0;// send to client to fix chunk size 
			do{	n = send(r_sock,(char *)(sent+(__int64)dp), SERV_BUF_LEN - sent, 0);
				sent += n;
				if(n<=0){ closesocket(r_sock);_close(fout);goto serv_recv_loop;}
			}while(sent<SERV_BUF_LEN);
			ret = (int)_lseeki64(fout, dp->offset.x, SEEK_SET);
			if(ret < 0){  closesocket(r_sock);_close(fout);goto serv_recv_loop;}
			left_to_send = dp->len.x;// set for page_ahead thread
			// **********************************************************************
			while(left_to_send > 0){// ***** data recv loop *************************

				if(left_to_send < CHUNK_SIZE) chunk_size = (__int64)left_to_send;
				else chunk_size = CHUNK_SIZE;
				sent = 0;// recv the data from client 
				do{n = recv(r_sock, (char *)(sent + (__int64)DB), 
					(__int64)chunk_size - sent, 0);
					sent += n;
					if(n<=0){ closesocket(r_sock);_close(fout);goto serv_recv_loop;}
				}while(sent < chunk_size);
				n = _write(fout, DB, chunk_size);
				if(n < (int)chunk_size){  closesocket(r_sock);_close(fout);goto serv_recv_loop;}

				left_to_send = left_to_send - chunk_size;
			}// ************ end recv loop *****************************************
			_close(fout); 
			dp->code.x = CODE_DPUT_DONE;
			sent = 0; // tell client all was put in mem 
			do{n = send(r_sock,(char *)(sent + (__int64)dp), SERV_BUF_LEN-sent, 0);
				sent += n;
				if(n<=0){ closesocket(r_sock); goto serv_recv_loop;}
			}while(sent<SERV_BUF_LEN);
			closesocket(r_sock);	
        	break;
		case CODE_DGET:// get disk 2 disk or mem **********************************
			dp = (disk_packet *)buf;// create mem for xfer

			fin = _open(dp->rem_file, _O_RDONLY | _O_BINARY );
			if(fin < 0){ 
				// abort client since remote file does not exist
				dp->code.x = CODE_ABORT;
				sent = 0;
				do{n = send(r_sock,(char *)(sent+(__int64)dp), SERV_BUF_LEN - sent, 0);
					sent += n;
					if(n<=0){ closesocket(r_sock);_close(fin);goto serv_recv_loop;}
				}while(sent < SERV_BUF_LEN);
				closesocket(r_sock);
				goto serv_recv_loop;;
			} 

			if(dp->len.x == 0){
#ifdef WINDOWS
				ret = _fstati64(fin, &statbuf);
#else
				ret = fstat64(fin, &statbuf);
#endif
				if(ret < 0){  closesocket(r_sock);_close(fin);goto serv_recv_loop;} 
				dp->len.x = (__int64)statbuf.st_size;
			}

			left_to_send = dp->len.x;

			ret = (int)_lseeki64(fin, dp->offset.x, SEEK_SET);
			if(ret < 0){ closesocket(r_sock);_close(fin);goto serv_recv_loop;}

			sent = 0;
			do{ n = send(r_sock,(char *)(sent+(__int64)dp), SERV_BUF_LEN - sent, 0);
				sent += n;
				if(n<=0){ closesocket(r_sock);_close(fin);goto serv_recv_loop;}
			}while(sent < SERV_BUF_LEN);

			// *******************************************************************
			while(left_to_send > 0){// ***** data send loop **********************
				// adjust packet size to remote mapped memory needs
				if(left_to_send < CHUNK_SIZE ) chunk_size = (u_int)left_to_send;
				else chunk_size = CHUNK_SIZE;

				n = _read(fin, DB, chunk_size);
				if(n != (int)chunk_size){  closesocket(r_sock);_close(fin);goto serv_recv_loop;}

				sent = 0;
				do{n = send(r_sock,(char *)(sent+(__int64)DB), chunk_size - sent, 0);
					sent += n;
					if(n<=0){ closesocket(r_sock);_close(fin);goto serv_recv_loop;
					}
				}while(sent < chunk_size);

				left_to_send = left_to_send - chunk_size;
			}// ************ end data send loop ************************************
			closesocket(r_sock);
			_close(fin);

        	break;
		case CODE_REMEX:
			rp = (rpc_packet *)buf;
			ret = remex(rp->data); // run process
			
			rp->code.x = ret;// send return code back to client
			sent = 0;
			do{n = send(r_sock, sent+buf, SERV_BUF_LEN - sent, 0); 
				sent += n;
				if(n<=0){closesocket(r_sock); goto serv_recv_loop;}
			}while(sent<SERV_BUF_LEN);
			closesocket(r_sock);
			break;
		default: // bad packet ****************************************************
			closesocket(r_sock); 
			break; // goes back to accept() at top of receiver loop
		}
	}// ******* recv loop *********************************************************
}// >>>>>>>>>>>>>>>>>> threaded >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>


// *************************************************************************
// *************************************************************************
// ******                  mn_server   threaded            *****************
// *************************************************************************
// *************************************************************************
/* This is the multi-purpose server function which is spawned through
   the mn_start_serv() call.  It can do the following:

   1. Serve a block of memory through put() and get() calls
	   OR
      Serve a function for rpc() and atomic_rpc() remote execution
   2. Serve local disks for remote reading and writing
   3. Perform remote execution of local executables                       */
// *************************************************************************
int mn_server_th(u_short AF, u_short port, char *host, void *addr, u_char num_threads){
#ifdef WINDOWS
	SOCKADDR_IN s_sock_addr;
	LPHOSTENT host_id;
//	DWORD optval = 1;// for setsockopt()
#else
	struct hostent *host_id;
	struct sockaddr_in s_sock_addr;
	pthread_t thread;
#endif
	SOCKET s_sock;//, r_sock;// server and remote sockets
	int n, ret;
	serv_pack *sp;


	// Create sock to listen: socket( addr fam, sock type, protocol)
#ifdef WINDOWS
	s_sock = socket(AF, SOCK_STREAM, IPPROTO_TCP);
#else
	s_sock = socket(AF, SOCK_STREAM, 0);
#endif

#ifdef WINDOWS
	if(s_sock == INVALID_SOCKET) return MN_SOCK;
#else
	if(s_sock<0) return MN_SOCK;
#endif 

//	fcntl(s_sock, F_SETFD, FD_CLOEXEC);// keep open during remex or thread

	if(num_threads > 0){// spawn threads so set sock opt
#ifdef WINDOWS
		n = 1;
		ret = setsockopt(s_sock, SOL_SOCKET, SO_REUSEADDR, (char *)&n, sizeof(n));
		if(ret != 0) return 0;
#else
		n = 1;
		ret = setsockopt(s_sock, SOL_SOCKET, SO_REUSEADDR, &n, sizeof(n));
		if(ret < 0) return 0;
#endif
	}

	// set address family always to AF_INET; AF in socket() is enough
	s_sock_addr.sin_family = AF_INET;
	s_sock_addr.sin_port = htons(port);

	// find host
	if(0 == strcmp("0", host) || 0 == strcmp("NULL", host) ){
		s_sock_addr.sin_addr.s_addr = INADDR_ANY;
	}
	else{
		host_id = gethostbyname(host);
		if (host_id == NULL){closesocket(s_sock);s_sock_addr.sin_addr.s_addr = INADDR_ANY;}
		else{
#ifdef WINDOWS
			s_sock_addr.sin_addr = *((LPIN_ADDR)*host_id->h_addr_list);
#else
			memcpy(&s_sock_addr.sin_addr.s_addr, host_id->h_addr_list[0], host_id->h_length );	
#endif
		}
	}


	// bind the name to the socket bind(sock, addr, size addr)
#ifdef WINDOWS
	n = bind(s_sock, (LPSOCKADDR)&s_sock_addr, sizeof(struct sockaddr));
#else
	n = bind(s_sock, (struct sockaddr *)&s_sock_addr, sizeof(struct sockaddr));
#endif
	if (n<0){  closesocket(s_sock); return MN_BIND;}

	// Set the socket to listen 
	n = listen(s_sock, SOMAXCONN);
	if (n<0){ closesocket(s_sock); return MN_LISTEN;}

	

	if(num_threads > 0){	// spawn all threads, then exit
		for(n=1; n <= num_threads; n++){
			sp = (serv_pack *)malloc(sizeof(serv_pack));// freed by server thread
			sp->s_sock = &s_sock;
			sp->addr = addr;
			sp->thread_id = n;
			sp->protocol = AF;
#ifdef WINDOWS
			ret = _beginthread(serv_thread, 4096, (void *)sp);
			Sleep(1);// give up a slice
#else
			ret = pthread_create(&thread, NULL, serv_thread, (void *)sp);
			if(ret != 0) printf("Failed to create thread # %d\n", n);
			pthread_detach(thread);// lose them
			sched_yield();//let thread go a bit
#endif
		}
		printf("%d server threads started ...\n", num_threads);
	}
	else { // single hanging server
		sp = (serv_pack *)malloc(sizeof(serv_pack));// freed by server thread
		sp->s_sock = &s_sock;
		sp->addr = addr;
		sp->thread_id = 0;
		sp->protocol = AF;
		printf("Un-threaded server hanging here ...\n");
		serv_thread((void *)sp);// hang here, no threads
	}

	return MN_SUCCESS;
} // end mn_start_serv_th



// *************************************************************************
// ******** mn_start_server ************************************************
// ************************************************************************* 
int mn_start_serv_th(u_short AF, u_short port, char *host, void* addr, u_char num_threads){
	if(num_threads < 0) return MN_FAIL;
	return mn_server_th(AF, port, host, addr, num_threads);
}

// starts then hangs, no threads spawned
int mn_start_serv(u_short AF, u_short port, char *host, void* addr, u_char flag){
	return mn_start_serv_th(AF, port, host, addr, 0);
}

// *************************************************************************
// **** mn_kill_sock_serv_process ******************************************
// *************************************************************************
// This kills it dead.  Kill a server with threads.
int mn_kill_sock_serv_process(u_short AF, u_short rem_port, char *rem_host){
	packet p;
	SOCKET	c_sock; // client socket
	int n;
	__int64 sent;

	// prepare packet
	p.magic_number.x = MAGIC_NUMBER;
	p.code.x = CODE_DIE;

	c_sock = make_sock(AF, rem_host, rem_port);
#ifdef WINDOWS
	if(c_sock == INVALID_SOCKET) return MN_SOCK;
#else
	if(c_sock<0) return MN_SOCK;
#endif 

	sent = 0;// send the killer packet
	do{n = send(c_sock, (char *)(sent+(__int64)&p), SERV_BUF_LEN-sent, 0);
		sent += n;
		if(n<=0){ closesocket(c_sock); return MN_SEND;}
	}while(sent<SERV_BUF_LEN);

	closesocket(c_sock); return MN_SUCCESS;
}

// stops a single server hanging, no threads
// *************************************************************************
// **** mn_stop_sock_serv **************************************************
// *************************************************************************
int mn_stop_sock_serv(u_short AF,  u_short rem_port, char *rem_host){
		packet p;
	SOCKET	c_sock; // client socket
	int n;
	__int64 sent;

	// prepare packet
	p.magic_number.x = MAGIC_NUMBER;
	p.code.x = CODE_STOP;

	c_sock = make_sock(AF, rem_host, rem_port);
#ifdef WINDOWS
	if(c_sock == INVALID_SOCKET) return MN_SOCK;
#else
	if(c_sock<0) return MN_SOCK;
#endif 

	sent = 0;// send the stop packet
	do{n = send(c_sock, (char *)(sent+(__int64)&p), SERV_BUF_LEN-sent, 0);
		sent += n;
		if(n<=0){ closesocket(c_sock); return MN_SEND;}
	}while(sent<SERV_BUF_LEN);
	closesocket(c_sock); return MN_SUCCESS;	
}


// *************************************************************************
// ***** remex <<server>> **************************************************
// *************************************************************************
// This is called up by a server
int remex(char *command_in){
	char command_str[SERV_BUF_LEN], *command;
#ifdef WINDOWS
// The windows version is straightforward, CreateProcess() does the job.
	PROCESS_INFORMATION procInfo;
	STARTUPINFO	startupInfo;
	typedef DWORD processId;

	strcpy(command_str, command_in); command = command_str;// this is copy from memory reused by server
	GetStartupInfo(&startupInfo);
	if( CreateProcess( NULL, 
		command, // command line to execute including path & arguments
		NULL, 
		NULL, 
		TRUE,
		DETACHED_PROCESS | NORMAL_PRIORITY_CLASS, 
		NULL, 
		NULL,
		&startupInfo, 
		&procInfo) ){
		return MN_SUCCESS;
	}
	return MN_FAIL;
#else // Linux
#define WR_MAX_ARG 20
// For Linux (UNIX) there is fork and some version of exec; vfork might be OK here too.
	char exec_path_name[80], *pArg, *pPtr, *argv[WR_MAX_ARG + 1];// max
	int argc;
	pid_t pid;
	char *environment[] = { "USER=user", "PATH=/tmp", NULL};

	strcpy(command_str, command_in); command = command_str;// this is copy from memory reused by server
	// get the executable name with path from command, 1st str	
	sscanf(command, "%s", exec_path_name);
	// extract the executable name less the path
	// start of last slash, then inc past
	if( (pArg = strrchr(exec_path_name, '/') ) != NULL) pArg++;
	else pArg = exec_path_name;// no slashes case
	argv[0] = pArg;// this is the executable name

	argc = 1;
	command += strlen(exec_path_name) + 1;// go to start of arguments
	if(command != NULL && *command != '\0'){
		pArg = strtok_r(command, " ", &pPtr);
		while(pArg != NULL){
			argv[argc] = pArg;
			argc++;
			if(argc >= WR_MAX_ARG) break;
			pArg = strtok_r(NULL, " ", &pPtr);
		}
	}
	argv[argc] = NULL;//terminate argument list

	if( (pid=vfork()) == 0){//child
	// execute it, here THIS program is replaced by what is execv'ed
		execve(exec_path_name, argv, environment);
		exit(1);// here only if exec fails
	}
	else if(pid < 0){
		return MN_FAIL;
	}
	else{// this is the parent process, and we print the spawned child's pid
		return MN_SUCCESS;
	}
#endif
}
