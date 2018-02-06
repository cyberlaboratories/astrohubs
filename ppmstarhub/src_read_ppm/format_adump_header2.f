
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   Routine: adump_format_header
c   Format data into bufffer "abuf" suitable for use in AIO output
c   
c   * global mappings: one for each output variable, specified in input
c     + 4-byte format (original floating point numbers)
c     + 2-byte formats:
c       - clipped linear
c       - clipped log scale
c   * flexible hierarchical cube layout (for optimal postprocesing)
c     + options for format specified in input
c     + one format file, includes mesh coordinates, & time between dumps
c     + nzonespercube[xyz] (usually 1 or 4 or full brick size)
c     + variables
c     + ncubesperbrick[xyz] (number of "cubes" in each direction in brick)
c     + nbrickstotal[xyz]
c     + nbricksperfile
c     + serial order in file(s): (izone[xyz],ivar,icube[xyz],ibrick[xyz])
c       q[xyz] is short for (qx,qy,qz) in FORTRAN order
c
c       ndcubesx, ndcubesy, ndcubesz     Number of cubes across brick in [xyz]
c       nbricksx, nbricksy, nbricksz     Number of bricks across whole domain in [xyz]
c       d(ndcube,ndcube,ndcube,ndvar,ndcubesx,ndcubesy,ndcubesz)
c
c       nbricksperfile = 1 + (nxbricks*nybricks*nzbricks-1) / nioservers
c       ibrick = ibrickx + nbricksx(ibricky-1 + nbricksy*(ibrickz-1))
c       ifile  = 1 + (ibrick-1) / nbricksperfile
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C TODO: Modify set_dump_parameters for your list of variables
c   OR: just call set_dump_var for each variable

      subroutine set_dump_parameters(ntx,nty,ntz)
      implicit none

      integer ntx, nty, ntz

      ! dump format and content parameters
      !
      integer MAXDVAR                      ! Maximum number of dump variables
      parameter(MAXDVAR=100)
      integer ndcube                       ! cube size in dump format
      integer ndvar                        ! Number of dump variables
      integer nbricksperfile               ! Max # of bricks in each file
      integer idvar(MAXDVAR)               ! Variable number in ddd array
      real dvarmin(MAXDVAR)                ! Minimum of dump variable range
      real dvarmax(MAXDVAR)                ! Maximum of dump variable range
      character*256 cVarSymb(MAXDVAR)      ! Variable Symble
      character*256 cVarName(MAXDVAR)      ! Variable Name and description
      character*256 cDumpMap(MAXDVAR)      ! Map from original to dump format
      character*256 cLoXbndry, cHiXbndry   ! Low and High boundary conditions in X
      character*256 cLoYbndry, cHiYbndry   ! Low and High boundary conditions in Y
      character*256 cLoZbndry, cHiZbndry   ! Low and High boundary conditions in Z
      common /dumpformat/ ndcube, ndvar, idvar, dvarmin, dvarmax
      common /dumpformat/ cVarSymb, cVarName, cDumpMap, nbricksperfile
      common /dumpformat/ cLoXbndry,cLoYbndry,cLoZbndry
      common /dumpformat/ cHiXbndry,cHiYbndry,cHiZbndry

      ! Temporary variables
      !
      integer      i
      character*52 st(100)
      character*40 bdylo, bdyhi


      ! TODO: Set list of variables for dump
      ! In this example:
      !     8 raw fields in ddd: {rho,xyz-momenta, B[xyz], total energy}
      !     However, want the 8 fields {rho,prs,V[xyz], B[xyz]}
      ! idvar: select field to be written out
      !   <  90     raw field in ddd array
      !   >= 90     derived field calculated in a user provided routine
      !   99        pressure (see subroutine getprs)
      !   [101-103] velocity component (see subroutine getvel)
      ! cVarSymb: Symbol to refer to field in dump
      ! cVarName: Description of field
      ! cDumpMap: "Log_e", "Linear", "real4"
      !   Log_e and Linear map field to a 16-bit integer,
      !   vaues are clipped outside of [dvarmin,dvarmax]
      !   real4 writes field as a 4-byte (real*4) floating point number
      ! dvarmin & dvarmax: min and max limits of Linear or Log_e map.
      !                    not used with real4 format
      !
      st(1) = "  1  RHO  Mass_Density      real4        0.0     0.0"
      st(2) = "  2  Vx   X-Velocity        real4        0.0     0.0"
      st(3) = "  3  Vy   Y-Velocity        real4        0.0     0.0"
      st(4) = "  4  Vz   Z-Velocity        real4        0.0     0.0"
      st(5) = "  5  Bx   Z-MagneticFeild   real4        0.0     0.0"
      st(6) = "  6  By   Z-MagneticFeild   real4        0.0     0.0"
      st(7) = "  7  Bz   Z-MagneticFeild   real4        0.0     0.0"
      st(8) = "  8  Etot Total_Energy      real4        0.0     0.0"
      ndvar = 8
      do i=1,ndvar
c        write (6,*) 'I will read for i = ', i
        read (st(i),  *) idvar(i), cVarSymb(i), cVarName(i),
     1                   cDumpMap(i), dvarmin(i), dvarmax(i)
c999     format(i3, 2x, a5, a18, a8, 2f10.1)
c        write (6,*) 'reading done' 
      enddo


      ! TODO: Set number of bricks per file
      ! In this example, want one output file ==> all bricks in one file
      !
      nbricksperfile = ntx*nty*ntz


      ! TODO: Set global boundary condition flags
      ! Values: "continuation", "periodic", "mirror"
      ! Helps document run.
      ! Will someday be used in f3d for spatial derivatives & filters
      !
      bdylo = "continuation  continuation  continuation"
      bdyhi = "continuation  continuation  continuation"
      read (bdylo,*) cLoXbndry,cLoYbndry,cLoZbndry
      read (bdyhi,*) cHiXbndry,cHiYbndry,cHiZbndry


      return
      end

C=====================================================================72
      subroutine set_dump_var(iv,csym,cnam,cmap,amin,amax)
      implicit none

      integer iv
      character*(*) csym, cnam, cmap
      real amin, amax

      ! dump format and content parameters
      !
      integer MAXDVAR                      ! Maximum number of dump variables
      parameter(MAXDVAR=100)
      integer ndcube                       ! cube size in dump format
      integer ndvar                        ! Number of dump variables
      integer nbricksperfile               ! Max # of bricks in each file
      integer idvar(MAXDVAR)               ! Variable number in ddd array
      real dvarmin(MAXDVAR)                ! Minimum of dump variable range
      real dvarmax(MAXDVAR)                ! Maximum of dump variable range
      character*256 cVarSymb(MAXDVAR)      ! Variable Symble
      character*256 cVarName(MAXDVAR)      ! Variable Name and description
      character*256 cDumpMap(MAXDVAR)      ! Map from original to dump format
      character*256 cLoXbndry, cHiXbndry   ! Low and High boundary conditions in X
      character*256 cLoYbndry, cHiYbndry   ! Low and High boundary conditions in Y
      character*256 cLoZbndry, cHiZbndry   ! Low and High boundary conditions in Z
      common /dumpformat/ ndcube, ndvar, idvar, dvarmin, dvarmax
      common /dumpformat/ cVarSymb, cVarName, cDumpMap, nbricksperfile
      common /dumpformat/ cLoXbndry,cLoYbndry,cLoZbndry
      common /dumpformat/ cHiXbndry,cHiYbndry,cHiZbndry

       integer i

      write (6,*) "iv=", iv, "     name=", trim(cnam)

      do i=1,256
        cVarSymb(iv)(i:i) = ' '
        cVarName(iv)(i:i) = ' '
        cDumpMap(iv)(i:i) = ' '
      enddo

      read (csym,*) cVarSymb(iv)
      read (cnam,*) cVarName(iv)
      read (cmap,*) cDumpMap(iv)

      ndvar       = iv
      idvar(iv)   = iv
      dvarmin(iv) = amin
      dvarmax(iv) = amax

      return
      end
C=====================================================================72
      subroutine format_adump_header(abuf,ndump,time,dtime,istep,
     1         xyzranges,nbytes,nvar,nx,ny,nz,itx,ity,itz,ntx,nty,ntz)
      implicit none

      ! format_adump argument list
      !
      integer ndump,itx,ity,itz,ntx,nty,ntz
      integer istep,nvar,nx,ny,nz
      real dtime, time, xyzranges(6)
      integer*2 abuf(nbytes/2)

      ! dump format and content parameters
      !
      integer MAXDVAR                      ! Maximum number of dump variables
      parameter(MAXDVAR=100)
      integer ndcube                       ! cube size in dump format
      integer ndvar                        ! Number of dump variables
      integer nbricksperfile               ! Max # of bricks in each file
      integer idvar(MAXDVAR)               ! Variable number in ddd array
      real dvarmin(MAXDVAR)                ! Minimum of dump variable range
      real dvarmax(MAXDVAR)                ! Maximum of dump variable range
      character*256 cVarSymb(MAXDVAR)      ! Variable Symble
      character*256 cVarName(MAXDVAR)      ! Variable Name and description
      character*256 cDumpMap(MAXDVAR)      ! Map from original to dump format
      character*256 cLoXbndry, cHiXbndry   ! Low and High boundary conditions in X
      character*256 cLoYbndry, cHiYbndry   ! Low and High boundary conditions in Y
      character*256 cLoZbndry, cHiZbndry   ! Low and High boundary conditions in Z
      common /dumpformat/ ndcube, ndvar, idvar, dvarmin, dvarmax
      common /dumpformat/ cVarSymb, cVarName, cDumpMap, nbricksperfile
      common /dumpformat/ cLoXbndry,cLoYbndry,cLoZbndry
      common /dumpformat/ cHiXbndry,cHiYbndry,cHiZbndry

      ! Scratch space
      !
c      integer      MAXZONES
c      parameter   (MAXZONES=200*200*200)
c      real d1field(MAXZONES)
c      common /dscratch/ d1field

      ! Temporary and equivalenced variables
      !
      integer MAXHEAD, MAXBUF, MAXVAR
      parameter (MAXHEAD=64*1024)
c      parameter (MAXVAR=8*MAXZONES)
c      parameter (MAXBUF =MAXHEAD + MAXVAR)
      integer   nhead
      integer nexpand(6)
      real fexpand(6)
      real dx, dy, dz, xbdyL,xbdyR,ybdyL,ybdyR,zbdyL,zbdyR
      character*80 cline
      character*256 cstr, cfile
      integer joff, mbytes, ierr, sw_compress2, nfill
      integer i, j, islot, ndcsx, ndcsy, ndcsz, iloge, i2max, nbytes
      integer ncsymb, ncname, ncdmap, nheadk, ninbuf, ivalue,i2off,idv
      real fmin, fmax, fmapmin, fmapmax, fi2max, fmapscale
      real value
      integer idash, length
      character*256 cctype


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      ! Exit if insufficient scratch space
c      !
c      if(nx*ny*nz .gt. MAXZONES) then
c        write (6,*) 'Insufficient space allocated for scratch array'
c        write (6,*) '   MAXZONES = ', MAXZONES
c        write (6,*) '   nx*ny*nz = ', nx*ny*nz
c        write (6,*) '   nx,ny,nz = ', nx,ny,nz
c        stop
c      endif
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      write (6,*) 'nbytes = ', nbytes

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c              Pack header info into the buffer: abuf                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! Scaled integer*2 values will range from 0 to 65530 inclusive
      ! This is an odd number of bins so that if dvarmin and dvarmax
      ! are symmetical around 0, then one bin will also be centered
      ! on zero.
      ! 
      i2off = 32750
      i2max = 65501
      fi2max = float(i2max)

      ! Pad header region of buffer with carriage returns
      nheadk = nbytes
      call padheader(abuf, nheadk)


      ! TODO: Set XYZ ranges of mesh.
      ! 
      dx = (xyzranges(2)-xyzranges(1)) / float(ntx*nx)
      dy = (xyzranges(4)-xyzranges(3)) / float(nty*ny)
      dz = (xyzranges(6)-xyzranges(5)) / float(ntz*nz)
      xbdyL = xyzranges(1)
      ybdyL = xyzranges(2)
      zbdyL = xyzranges(3)
      xbdyR = xyzranges(4)
      ybdyR = xyzranges(5)
      zbdyR = xyzranges(6)

      ! TODO: If needed, set parameters for exponentially expanding mesh here.
      ! nexpand(1:6)=0 means a mesh is uniform everywhere.
      do i = 1, 6
        nexpand(i) = 0
        fexpand(i) = 1.0
      enddo


      ndcube = 0  ! ndcube>0 not supported yet
      if(ndcube .gt. 0) then
        ! Data will be cubed up, ndcube on a side
        ! validate dimensions
        if(mod(nx,ndcube) .ne. 0  .or.
     1     mod(ny,ndcube) .ne. 0  .or.
     1     mod(nz,ndcube) .ne. 0      ) then
          write (6,*) "ndcube does not divine evenly into brick"
          write (6,*) "nx,ny,nz=", nx, ny, nz
          write (6,*) "ndcube  =", ndcube
          stop
        endif
 
        ndcsx = nx / ndcube
        ndcsy = ny / ndcube
        ndcsz = nz / ndcube
      else
        ! If ndcube = 0, it is ignored as a dimension
        ndcsx = nx
        ndcsy = ny
        ndcsz = nz
      endif


      ! Get max string lengths (padded by 2 characters)
      ncsymb = 2
      ncname = 2
      ncdmap = 2
      do i=1,ndvar
        do j=1,250
          if(cVarSymb(i)(j:j) .ne. ' ') ncsymb = max(ncsymb, j+2)
          if(cVarName(i)(j:j) .ne. ' ') ncname = max(ncname, j+2)
          if(cDumpMap(i)(j:j) .ne. ' ') ncdmap = max(ncdmap, j+2)
        enddo
      enddo


      nhead = 0
      call addline_1f("time      ",   time,            nhead, abuf)
      call addline_1f("dtime     ",  dtime,            nhead, abuf)
      call addline_1i("step      ",  istep,            nhead, abuf)
      call addline_1i("ncube     ", ndcube,            nhead, abuf)
      call addline_1i("i2max     ",  i2max,            nhead, abuf)
      call addline_1i("i2off     ",  i2off,            nhead, abuf)
      call addline_3i("ncubesxyz ", ndcsx,ndcsy,ndcsz, nhead, abuf)
      call addline_3i("nbrickxyz ",   ntx,  nty,  ntz, nhead, abuf)
      call addline_1i("ndvar     ", ndvar,             nhead, abuf)
      call addline_1i("nbpf      ", nbricksperfile,    nhead, abuf)
      call addline_2f("xlimits   ", xbdyL, xbdyR,      nhead, abuf)
      call addline_2f("ylimits   ", ybdyL, ybdyR,      nhead, abuf)
      call addline_2f("zlimits   ", zbdyL, zbdyR,      nhead, abuf)

      call addline_1i1f("expandxl  ",nexpand(1),fexpand(1),nhead,abuf)
      call addline_1i1f("expandxr  ",nexpand(2),fexpand(2),nhead,abuf)
      call addline_1i1f("expandyl  ",nexpand(3),fexpand(3),nhead,abuf)
      call addline_1i1f("expandyr  ",nexpand(4),fexpand(4),nhead,abuf)
      call addline_1i1f("expandzl  ",nexpand(5),fexpand(5),nhead,abuf)
      call addline_1i1f("expandzr  ",nexpand(6),fexpand(6),nhead,abuf)

      cstr(1:10) = "lobndry   "
      call add_string(     cstr,  10, nhead, abuf)
      call add_string(cLoXbndry,  16, nhead, abuf)
      call add_string(cLoYbndry,  16, nhead, abuf)
      call add_string(cLoZbndry, -16, nhead, abuf)  ! nadd < 0: add a charage return

      cstr(1:10) = "hibndry   "
      call add_string(     cstr,  10, nhead, abuf)
      call add_string(cHiXbndry,  16, nhead, abuf)
      call add_string(cHiYbndry,  16, nhead, abuf)
      call add_string(cHiZbndry, -16, nhead, abuf)  ! nadd < 0: add a charage return

      do i=1,ndvar
        cstr(1:10) = "field3d   "
        call add_string(       cstr,     10, nhead, abuf)
        call add_string(cVarSymb(i), ncsymb, nhead, abuf)
        call add_string(cVarName(i), ncname, nhead, abuf)
        call add_string(cDumpMap(i), ncdmap, nhead, abuf)
        call addline_2f("          ",dvarmin(i),dvarmax(i),nhead,abuf)
      enddo

      call addline_1i("headsize  ", nheadk, nhead, abuf)
      ! headsize is always the last line in the header info
      ! Header info is followd by MANY charage returns

      if(nhead .gt. nbytes) then
        write (6,*) 'PROBLEM in adump_format_header:'
        write (6,*) '  size of header buffer (in bytes):', nbytes
        write (6,*) '  size of header data   (in bytes):', nheadk
        stop
      endif

      return
      end
c=======================================================================
c
c      subroutine fill_fbuf(fbuf, d1field, joff, nfill, MAXBUF)
c      parameter   (MAXZONES=128*128*128)
c      real d1field(MAXZONES)
c      real*4 fbuf(MAXBUF/2)
c
c      do j=1, nfill
c        fbuf(joff+j) = d1field(j)
c      enddo
c      return
c      end
c
c=======================================================================

      subroutine padheader(cbuf, nheadk)
c      parameter (MAXHEAD=64*1024)
      character cbuf(nheadk)
      do i = 1, nheadk
        cbuf(i:i) = char(10)
      enddo
      return
      end

c=======================================================================

      subroutine uniformmeshlr(nx,xl,nexp1,nexp2,bdyL,bdyR)
      real xl(-6:4096+1)

      n1 = 1 + nexp1              ! Left edge of 1st uniform cell
      n2 = 1 + nx - nexp2         ! Left edge of 1st non-uniform cell on right
      dx = (xl(n2) - xl(n1)) / float(n2-n1)
      bdyL = xl(n1) - dx*nexp1    ! Left  edge of non-expanding mesh
      bdyR = xl(n2) + dx*nexp2    ! Right edge of non-expanding mesh

      return
      end

c=======================================================================
      subroutine add_string(cstr, madd, nhead, cbuf)
      parameter (MAXHEAD=64*1024)
      character cbuf(MAXHEAD)
      character*256 cstr
      
      nadd = iabs(madd)
      do i=1,nadd
        nhead = nhead + 1
        cbuf(nhead:nhead) = cstr(i:i)
      enddo
      if(madd .le. 0) then
        nhead = nhead + 1
        cbuf(nhead:nhead) = char(10)
      endif

      return
      end

c=======================================================================
      subroutine addline_1i(ckey, ivar, nhead, cbuf)
      parameter (MAXHEAD=64*1024)
      character*10 ckey
      character cbuf(MAXHEAD)
      character*80 cline

      write (cline,999) ckey, ivar
999   format(a10, i12)
      nadd = 23
      cline(nadd:nadd) = char(10)
      
      do i=1,nadd
        nhead = nhead + 1
        cbuf(nhead:nhead) = cline(i:i)
      enddo

      return
      end

c=======================================================================
      subroutine addline_3i(ckey, iv1, iv2, iv3, nhead, cbuf)
      parameter (MAXHEAD=64*1024)
      character*10 ckey
      character cbuf(MAXHEAD)
      character*80 cline

      ! Construct new line
      write (cline,999) ckey, iv1, iv2, iv3
999   format(a10, 3i8)
      nadd = 35
      cline(nadd:nadd) = char(10)

      ! Append line to buffer
      do i=1,nadd
        nhead = nhead + 1
        cbuf(nhead:nhead) = cline(i:i)
      enddo

      return
      end

c=======================================================================
      subroutine addline_1f(ckey, val, nhead, cbuf)
      parameter (MAXHEAD=64*1024)
      character*10 ckey
      character cbuf(MAXHEAD)
      character*80 cline

      ! Construct new line
      if(ckey(1:1) .ne. ' ') then
        write (cline,999) ckey, val
999     format(a10, 1pe15.6)
c        nadd = 28
        nadd = 26
      else
        write (cline,998) val
998     format(1pe15.6)
        nadd = 16
      endif
      cline(nadd:nadd) = char(10)

      ! Append line to buffer
      do i=1,nadd
        nhead = nhead + 1
        cbuf(nhead:nhead) = cline(i:i)
      enddo

      return
      end


c=======================================================================
      subroutine addline_1i1f(ckey, ival, val, nhead, cbuf)
      parameter (MAXHEAD=64*1024)
      character*10 ckey
      character cbuf(MAXHEAD)
      character*80 cline

      ! Construct new line
      if(ckey(1:1) .ne. ' ') then
        write (cline,999) ckey, ival, val
999     format(a10, i10, 1pe15.6)
        nadd = 36
      else
        write (cline,998) ival, val
998     format(i10, 1pe15.6)
        nadd = 26
      endif
      cline(nadd:nadd) = char(10)

      ! Append line to buffer
      do i=1,nadd
        nhead = nhead + 1
        cbuf(nhead:nhead) = cline(i:i)
      enddo

      return
      end


c=======================================================================
      subroutine addline_2f(ckey, v1, v2, nhead, cbuf)
      parameter (MAXHEAD=64*1024)
      character*10 ckey
      character cbuf(MAXHEAD)
      character*80 cline

      ! Construct new line
      if(ckey(1:1) .ne. ' ') then
        write (cline,999) ckey, v1, v2
999     format(a10, 1p2e15.6)
        nadd = 41
      else
        write (cline,998) v1, v2
998     format(1p2e15.6)
        nadd = 31
      endif
      cline(nadd:nadd) = char(10)

      ! Append line to buffer
      do i=1,nadd
        nhead = nhead + 1
        cbuf(nhead:nhead) = cline(i:i)
      enddo

      return
      end

C=====================================================================72
C=====================================================================72

C TODO: Modify routines getone, getprs, & getvel for your data structures
C Examples here assume variables are interleaved, and no ghost zones
C See getflds_cube for examples with "cubed" data and using common.inc

C=====================================================================72
      subroutine getone(ivar, nvar, nx, ny, nz,  dd, d1field)
      real dd(nvar, nx, ny, nz)
      real d1field(nx,ny,nz)

      do iz=1,nz
      do iy=1,ny
      do ix=1,nx
         d1field(ix,iy,iz) = dd(ivar,ix,iy,iz)
      enddo
      enddo
      enddo

      return
      end

C=====================================================================72
      subroutine getprs(nvar, nx, ny, nz,  dd, d1field)
      real dd(nvar, nx, ny, nz)
      real d1field(nx,ny,nz)

      gam = 5.0 / 3.0
      gam1 = gam - 1.0

      do iz=1,nz
      do iy=1,ny
      do ix=1,nx
         rho  = dd(1,ix,iy,iz)
         rvx  = dd(2,ix,iy,iz)
         rvy  = dd(3,ix,iy,iz)
         rvz  = dd(4,ix,iy,iz)
         bx   = dd(5,ix,iy,iz)
         by   = dd(6,ix,iy,iz)
         bz   = dd(7,ix,iy,iz)
         etot = dd(8,ix,iy,iz)

         ekin = 0.5 * (rvx**2 + rvy**2 + rvz**2) / rho
         emag = 0.5 * ( bx**2 +  by**2 +  bz**2)
         prs  = gam1*(etot-ekin-emag)

         d1field(ix,iy,iz) = prs
      enddo
      enddo
      enddo

      return
      end

C=====================================================================72
      subroutine getvel(ivar, nvar, nx, ny, nz,  dd, d1field)
      real dd(nvar, nx, ny, nz)
      real d1field(nx,ny,nz)

      do iz=1,nz
      do iy=1,ny
      do ix=1,nx
         rho = dd(1,ix,iy,iz)
         if(ivar .eq. 101) d1field(ix,iy,iz) = dd(2,ix,iy,iz) / rho
         if(ivar .eq. 102) d1field(ix,iy,iz) = dd(3,ix,iy,iz) / rho
         if(ivar .eq. 103) d1field(ix,iy,iz) = dd(4,ix,iy,iz) / rho
      enddo
      enddo
      enddo

      return
      end

C=====================================================================72
