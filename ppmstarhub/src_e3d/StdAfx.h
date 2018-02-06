
#if !defined(AFX_STDAFX_H__3FAEA897_5B8B_4437_B161_44351198DFA0__INCLUDED_)
#define AFX_STDAFX_H__3FAEA897_5B8B_4437_B161_44351198DFA0__INCLUDED_

#ifdef LINUX
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>

#define __int64 long long
#define DWORD unsigned int
#define LPVOID void *
#define LPCTSTR const char *
#define _O_WRONLY O_WRONLY
#define _O_CREAT O_CREAT
typedef size_t SIZE_T;
#endif

#if _MSC_VER > 1000
#pragma once
#endif


#include <stdio.h>
#endif
