      integer MAXDVAR                      ! Maximum number of dump variables
      parameter(MAXDVAR=500)

      character*256 cDumpFile              ! Path to pointer file
      character*256 cDumpFiles(MAXDVAR)    ! Path to all files for current dump
      integer ncube                        ! cube size in dump format (0: this dim is ignored)
      integer ndvar                        ! Number of dump variables
      integer nbricksperfile               ! Max # of bricks in each file
      integer nx,ny,nz                     ! Number of zones in XYZ directions (per brick)
      integer multfac                      ! mesh mult fac for low res data
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
      character*16  cxcoord,cycoord,czcoord ! XYZ coordinate symbols
      integer istep                        ! Current time step

      integer nthdump                      ! Current dump number
      integer nstepperdump                 ! Number of steps per dump
      real    dtdump                       ! Time per dump
      real time, dtime                     ! Current time and time step of dump
      real xbndyL,xbndyR                   ! X Left & Right boundaries of full domain
      real ybndyL,ybndyR                   ! Y Left & Right boundaries of full domain
      real zbndyL,zbndyR                   ! Z Left & Right boundaries of full domain
      common /dumpinfo/ ncube, ndvar, idvar, dvarmin, dvarmax
      common /dumpinfo/ cVarSymb, cVarName, cDumpMap, nbricksperfile
      common /dumpinfo/ nx,ny,nz, nbx,nby,nbz, ncx,ncy,ncz, multfac
      common /dumpinfo/ time, dtime, istep, i2max, i2off, nheadsize
      common /dumpinfo/ xbndyL,xbndyR,ybndyL,ybndyR,zbndyL,zbndyR
      common /dumpinfo/ nthdump, nstepperdump, dtdump
      common /dumpinfo/ cDumpFile, cDumpFiles
      common /dumpinfo/ cLoXbndry,cLoYbndry,cLoZbndry
      common /dumpinfo/ cHiXbndry,cHiYbndry,cHiZbndry
      common /dumpinfo/ cxcoord,cycoord,czcoord

