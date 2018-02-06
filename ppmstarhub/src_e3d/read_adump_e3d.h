      integer MAXDVAR
      parameter(MAXDVAR=100)               ! Maximum number of dump variables

      character*256 cDumpFile              ! Path to dump file
      integer ncube                        ! cube size in dump format (0: this dim is ignored)
      integer ndvar                        ! Number of dump variables
      integer nbricksperfile               ! Max # of bricks in each file
      integer nx,ny,nz                     ! Number of zones in XYZ directions (per brick)
      integer ncx,ncy,ncz                  ! Number of cubes in XYZ directions
      integer nbx,nby,nbz                  ! Number of bricks in XYZ directions

      integer idvar(MAXDVAR)               ! Variable number in ddbrick array
      integer i2max                        ! range of 2-byte integers
      integer i2off                        ! offset of 2-byte integers
      integer nheadsize                    ! Size of header info in bytes
      real dvarmin(MAXDVAR)                ! Minimum of dump variable range
      real dvarmax(MAXDVAR)                ! Maximum of dump variable range
      character*256 cVarSymb(MAXDVAR)      ! Variable Symbols
      character*256 cVarName(MAXDVAR)      ! Variable Names/description
      character*256 cDumpMap(MAXDVAR)      ! Map from original to dump format
      character*256 cLoXbndry,cHiXbndry    ! X boundary types
      character*256 cLoYbndry,cHiYbndry    ! Y boundary types
      character*256 cLoZbndry,cHiZbndry    ! Z boundary types
      integer istep                        ! Current time step
      real time, dtime                     ! Current time and time step of dump
      real xbndyL,xbndyR                   ! X Left & Right boundaries of full domain
      real ybndyL,ybndyR                   ! Y Left & Right boundaries of full domain
      real zbndyL,zbndyR                   ! Z Left & Right boundaries of full domain

      common /dumpinfo/ ncube, ndvar, idvar, dvarmin, dvarmax
      common /dumpinfo/ cVarSymb, cVarName, cDumpMap, nbricksperfile
      common /dumpinfo/ nx, ny, nz, nbx, nby, nbz, ncx, ncy, ncz
      common /dumpinfo/ time, dtime, istep, i2max, i2off, nheadsize
      common /dumpinfo/ xbndyL,xbndyR,ybndyL,ybndyR,zbndyL,zbndyR
      common /dumpinfo/ cDumpFile
      common /dumpinfo/ cLoXbndry,cLoYbndry,cLoZbndry
      common /dumpinfo/ cHiXbndry,cHiYbndry,cHiZbndry

      include "amr_e3d.h"
