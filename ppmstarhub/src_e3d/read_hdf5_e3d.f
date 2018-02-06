c  Read Simple HDF5 files
c  David H. Porter
c  Versions:
c  Date    2013 June 25    Adapte to hdf5
c
c-------------------------------------------------------------
      subroutine read_header_info(cFile)
      implicit none
      include 'read_hdf5_e3d.h'
      character*256 cFile, cline, ckey, cPDTfile
      integer i, j, iField
      data nx/-1/

      if(cDumpFile .eq. cfile) return   ! Return if already read
      cDumpFile = cFile

      call set_read_header_defaults()   ! set/reset defaults

      do i = 1, ndvar                   ! for now: real*4 data type
        cDumpMap(i)(1:5) = "real4"
      enddo

      call read_pointer_file()          ! Read meta-data from pointer file

      return
      end

c---------------------------------------------------------------------72
      subroutine get_time(current_time)
      implicit none
      real current_time
      include 'read_hdf5_e3d.h'
      current_time = time
      return
      end

c---------------------------------------------------------------------72
      subroutine get_ndvar(n)
      implicit none
      integer n
      include 'read_hdf5_e3d.h'
      n = ndvar
      return
      end

c---------------------------------------------------------------------72
      subroutine get_dxdydz(dx, dy, dz)
      implicit none
      real dx, dy, dz
      include 'read_hdf5_e3d.h'
      dx = (xbndyR-xbndyL) / float(nx)
      dy = (ybndyR-ybndyL) / float(ny)
      dz = (zbndyR-zbndyL) / float(nz)
      return
      end

c---------------------------------------------------------------------72
      subroutine get_full_XYZrange(x0,x1,y0,y1,z0,z1)
      implicit none
      real x0,x1,y0,y1,z0,z1
      include 'read_hdf5_e3d.h'

      if(nx .lt. 0) then
        write (6,*) 'Problem: header info not set'
        stop
      endif

      x0 = xbndyL
      x1 = xbndyR
      y0 = ybndyL
      y1 = ybndyR
      z0 = zbndyL
      z1 = zbndyR

      return
      end

c---------------------------------------------------------------------72
      subroutine get_full_dimensions(nx1,ny1,nz1)
      implicit none
      integer nx1, ny1, nz1
      include 'read_hdf5_e3d.h'

      if(nx .lt. 0) then
        write (6,*) 'Problem: header info not set'
        stop
      endif

      nx1 = nx
      ny1 = ny
      nz1 = nz

      return
      end

c---------------------------------------------------------------------72
      subroutine get_vars_nams(MAX_TOF,ndvar1,cFld_tof,cNam_tof)
      implicit none
      include 'read_hdf5_e3d.h'
      integer MAX_TOF, ndvar1, i
      character*32 cFld_tof(MAX_TOF)     ! Symbol for field
      character*32 cNam_tof(MAX_TOF)     ! Name of field
      ndvar1 = ndvar
      do i = 1, ndvar
         cFld_tof(i) = cVarSymb(i)(1:32)
         cNam_tof(i) = cVarName(i)(1:32)
      enddo

      return
      end

c---------------------------------------------------------------------72
      subroutine report_header_info()
      implicit none
      include 'read_hdf5_e3d.h'
      integer ilen, i,j

      if(nx .lt. 0) then
        write (6,*) 'Problem: header info not set'
        return
      endif

      ilen = len(trim(cdumpfile))
      write (6,*) 'Data File: ', cdumpfile(1:ilen)
      write (6,9901) time, dtime, istep
9901  format("Time  = ", f10.5,/,"Dtime = ", f10.5,/,"Step  = ", i5)
      if(nx.gt.1) write (6,9902) trim(cxcoord),xbndyL,xbndyR,nx
      if(ny.gt.1) write (6,9902) trim(cycoord),ybndyL,ybndyR,ny
      if(nz.gt.1) write (6,9902) trim(czcoord),zbndyL,zbndyR,nz
9902  format(a,":  [", f9.4,",",f8.4,"]     ", "Mesh =",i4)

      write (6,*)  "  "
      write (6,*) "Symbol      Name"
      do i=1,ndvar
        write (6,9904) cVarSymb(i)(1:12), cVarName(i)(1:20)
      enddo
9904  format(" ",a, a, 1p2e15.6)

      return
      end

c---------------------------------------------------------------------72
      subroutine xrange2ixos(xlo,xhi,ixoff,ixsize)
      implicit none
      include 'read_hdf5_e3d.h'
      integer ixlo,ixhi,ixoff,ixsize, nxm
      real    xlo,xhi,xfrac,dx,eps

      nxm = nx
      eps = 0.01 * (xbndyR-xbndyL) / float(nxm)
      xfrac = (xlo-xbndyL) / (xbndyR-xbndyL)
      ixlo  = max(1, min(nxm, 1+int(float(nxm)*xfrac)))
      xfrac = max(xfrac, (xhi-xbndyL-eps) / (xbndyR-xbndyL))
      ixhi  = max(1, min(nxm, 1+int(float(nxm)*xfrac)))

      dx = (xbndyR-xbndyL) / float(nxm)
      xlo = xbndyL + dx * (ixlo - 1)
      xhi = xbndyL + dx * ixhi
      ixoff  = ixlo - 1
      ixsize = 1 + ixhi - ixlo

      return
      end

c---------------------------------------------------------------------72
      subroutine yrange2iyos(ylo,yhi,iyoff,iysize)
      implicit none
      include 'read_hdf5_e3d.h'
      integer iylo,iyhi,iyoff,iysize, nym
      real    ylo,yhi,yfrac,dy,eps

      nym = ny
      eps = 0.01 * (ybndyR-ybndyL) / float(nym)
      yfrac = (ylo-ybndyL) / (ybndyR-ybndyL)
      iylo  = max(1, min(nym, 1+int(float(nym)*yfrac)))
      yfrac = max(yfrac, (yhi-ybndyL-eps) / (ybndyR-ybndyL))
      iyhi  = max(1, min(nym, 1+int(float(nym)*yfrac)))

      dy = (ybndyR-ybndyL) / float(nym)
      ylo = ybndyL + dy * (iylo - 1)
      yhi = ybndyL + dy * iyhi
      iyoff  = iylo - 1
      iysize = 1 + iyhi - iylo

      return
      end

c---------------------------------------------------------------------72
      subroutine zrange2izos(zlo,zhi,izoff,izsize)
      implicit none
      include 'read_hdf5_e3d.h'
      integer izlo,izhi,izoff,izsize, nzm
      real    zlo,zhi,zfrac,dz, eps

      nzm = nz
      eps = 0.01 * (zbndyR-zbndyL) / float(nzm)
      zfrac = (zlo-zbndyL) / (zbndyR-zbndyL)
      izlo  = max(1, min(nzm, 1+int(float(nzm)*zfrac)))
      zfrac = max(zfrac, (zhi-zbndyL-eps) / (zbndyR-zbndyL))
      izhi  = max(1, min(nzm, 1+int(float(nzm)*zfrac)))

      dz = (zbndyR-zbndyL) / float(nzm)
      zlo = zbndyL + dz * (izlo - 1)
      zhi = zbndyL + dz * izhi
      izoff  = izlo - 1
      izsize = 1 + izhi - izlo

      return
      end

c---------------------------------------------------------------------72
      subroutine read_fld_from_one_box(ibox,cVar,iz,nrdx,nrdy,nrdz,fld)
      integer ibox,iz,nrdx,nrdy,nrdz
      real fld(nrdx,nrdy,nrdz)
      character*256 cVar
      write (6,*) "read_fld_from_one_box not implmented for hdf5"
      stop
      end

c---------------------------------------------------------------------72
      subroutine read_fields(cVars,nvars,nxsub,nysub,nzsub,
     1                       ixoff,iyoff,izoff, flds)
      implicit none
      character*256 cVars, cHost, cDir, cFile
      integer nvars,nxsub,nysub,nzsub,ixoff,iyoff,izoff
      real flds(nxsub,nysub,nzsub)         ! read only one field at a time
      include 'read_hdf5_e3d.h'
      integer ifld, nvar, ivar(MAXDVAR), ix, iy, iz

      !  Make sure header info is read in
      !  Then parse input string "cVars" for fields to be read
      !  Then allocate scratch arrays
      if(nx .lt. 0) then
        write (6,*) 'Problem: header info not set'
        stop
      endif

       ! Might eventually support more complicated boundaries, but:
      ! For now, only support reading subregions interior to full volume
      if(ixoff.lt.0 .or. ixoff+nxsub.gt.nx .or.
     1   iyoff.lt.0 .or. iyoff+nysub.gt.ny .or.
     1   izoff.lt.0 .or. izoff+nzsub.gt.nz     ) then
          write (6,*) "PROBLEM in read_hdf5_e3d.f read_fields:"
          write (6,*) "Subregion outside of full volume"
          write (6,*) "This is not currently supported"
          stop
      endif

       ! Only support reading one field at a time (no interleaved  output)
       if(nvars.gt.1) then
          write (6,*) "PROBLEM in read_hdf5_e3d.f read_fields:"
         write (6,*) "reading multiple interleaved fields not supported"
          stop
       endif

      call parse_vars(MAXDVAR,ndvar,cVarSymb,cVars,nvar,ivar)
      ifld = ivar(1)

      do iz=1,nzsub
      do iy=1,nysub
      do ix=1,nxsub
        flds(ix,iy,iz) = 2.0
      enddo
      enddo
      enddo

      call read_hdf5_sub(nxsub,nysub,nzsub,  ixoff,iyoff,izoff,
     1                   cDumpFiles(ifld),cVarName(ifld),flds)

      return
      end

C=====================================================================72
      subroutine get_HDF(iFile, cHostDirName, cHost, cDir, cName)
      character*256 cHostDirName, cHost, cDir, cName

      ! Find the index in cHostDirName berween the directory and the file name 
      ! I should be the last character in the string
      islash = 0
      ilast  = 256
      do i=1, 256
        cHost(i:i) = ' '
        cDir(i:i)  = ' '
        cName(i:i) = ' '
        if(cHostDirName(i:i)         .eq. '/') islash = i
        if(cHostDirName(257-i:257-i) .eq. ' ') ilast  = 256-i
      enddo

      iPathLen = islash - 1
      iNameLen = ilast-islash

      ! Default: the file is in the local directory
      cHost(1:9) = 'localdisk'

      if(islash .gt. 0) then
        ! cHostDirName contains path info to the file
        cDir(1:iPathLen)  = cHostDirName(1:iPathLen)
        cName(1:iNameLen) = cHostDirName(islash+1:ilast)
      else
        ! cHostDirName is just the name of the file, assume it is in cwd
        cDir(1:1)  = '.'
        cName      = cHostDirName
      endif

      ! Reformat iFile into last 3 characters of file name
      i1 =     iFile/100
      i2 = mod(iFile/10 , 10)
      i3 = mod(iFile    , 10)
      cName(iNameLen-2:iNameLen-2) = char(ichar('0') + i1)
      cName(iNameLen-1:iNameLen-1) = char(ichar('0') + i2)
      cName(iNameLen  :iNameLen  ) = char(ichar('0') + i3)

      return
      end

C=====================================================================72
      subroutine parse_vars(MAXDVAR,ndvar,cVarSymb,cVars,nvar,ivar)
      character*256 cVarSymb(MAXDVAR)
      character*256 cVars
      integer ndvar, nvar, ivar(MAXDVAR)

      integer i, j, istr, nstr
      character*256  cstr(1000)

      ! Initialize strings
      do j = 1, 1000
        do i = 1, 256
          cstr(j)(i:i) = ' '
        enddo
      enddo

      ! Break cVars into an array of strings delimited by spaces
      nstr = 0
      istr = 0
      do i = 1, 256
        if(cVars(i:i) .eq. ' ') then
          istr = 0
        else
          if(istr .eq. 0) nstr = nstr + 1
          istr = istr + 1
          cstr(nstr)(istr:istr) = cVars(i:i)
        endif
      enddo
      
      ! Get data fields to be read in
      nvar = 0
      do i = 1, nstr
        do j=1,ndvar
          if(cstr(i) .eq. cVarSymb(j)) then
            nvar = nvar + 1
            ivar(nvar) = j
c            write (6,*) "nvar,ivar(nvar)=",nvar,ivar(nvar)
           endif
        enddo
      enddo

      return
      end
c---------------------------------------------------------------------72
      subroutine set_read_header_defaults()
      implicit none
      include 'read_hdf5_e3d.h'
      integer i, j

      ! Initialize & set defaults
      !
      do i = 1, 256
        cLoXbndry(i:i) = ' '
        cLoYbndry(i:i) = ' '
        cLoZbndry(i:i) = ' '
        cHiXbndry(i:i) = ' '
        cHiYbndry(i:i) = ' '
        cHiZbndry(i:i) = ' '
      enddo
      cxcoord = "X               "
      cycoord = "Y               "
      czcoord = "Z               "

      do j=1,MAXDVAR
        do i = 1, 256
          cVarSymb(j)(i:i)  = ' '
          cVarName(j)(i:i)  = ' '
          cDumpMap(j)(i:i)  = ' '
        enddo
      enddo
      nheadsize = -1
      ncube     = 0
      ndvar     = 0
      time      = 0.0
      dtime     = 0.0
      istep     = 0.0
      nbx =1
      nby = 1
      nbz = 1
      nbricksperfile = 1

      xBndyL = -0.5
      xBndyR =  0.5
      yBndyL = -0.5
      yBndyR =  0.5
      zBndyL = -0.5
      zBndyR =  0.5

      cLoXbndry(1:12) = "continuation"
      cLoYbndry(1:12) = "continuation"
      cLoZbndry(1:12) = "continuation"
      cHiXbndry(1:12) = "continuation"
      cHiYbndry(1:12) = "continuation"
      cHiZbndry(1:12) = "continuation"

      return
      end

C=====================================================================72
      subroutine get_vec_geom(n,vec4)
      integer n, i
      real vec4(4,10000)
      n = 0
      return
      end

C=====================================================================72
c end of external interface
C=====================================================================72

      subroutine read_pointer_file()
      include "read_hdf5_e3d.h"
      character*256 line
      integer nstr, iparse_line
      character*(256) str(100)
      character*(256) ctype, cspace
      character*(256) cfield, cfile

      ! Defaults
      nthdump = 0
      dtdump  = 1.0
      nstepperdump  = 1

      ! invalid values
      time = -999.0
      dtime = -999.0
      istep = -999

      ndvar = 0
      open(unit=11,file=cDumpFile, form='formatted',err=9000)
10    call set_blank(line)
      read (11,999,end=100) line
999   format(a)
      nstr = iparse_line(line, str)
      if(nstr .gt. 1) then
99      format(a)
        if(str(1)(1: 8) .eq. "nthdump "  ) read(str(2),*) nthdump
        if(str(1)(1: 5) .eq. "file "     ) read(str(2),99) cfile
        if(str(1)(1: 9) .eq. "DATATYPE " ) read(str(2),99) ctype
        if(str(1)(1:10) .eq. "DATASPACE ") then
          read(str(2),99) cspace
          read(str(5),*) nx
          read(str(4),*) ny
          read(str(3),*) nz
        endif
        if(str(1)(1:8) .eq. "DATASET ") then
          ndvar = ndvar + 1
          read(str(2),99) cVarName(ndvar)
          read(cfile,99) cDumpFiles(ndvar)
          call gen_symb_from_name(cVarName(ndvar), cVarSymb(ndvar))
        endif
      endif
      go to 10
100   continue
      close(11)

      if(time  .eq. -999.0) time = dtdump * float(nthdump)
      if(dtime .eq. -999.0) dtime = dtdump / nstepperdump
      if(istep .eq. -999  ) istep = nstepperdump * nthdump
      return 

9000  write (6,*) "PROBLEM: could not open file"
      stop

      end

C=====================================================================72
      integer function iparse_line(line, strings)
      implicit none
      character*(*) line
      character*(256) strings(100)
      integer n, nstr, i, i0, i1, isrc, idst, itype
      logical instr   ! true if in a string
      nstr = 0
      n = len(line)
      instr = .false.
      do i = 1, n
         itype = 1                         ! NOT whte space
        if(line(i:i) .eq. ' ') itype = 0   ! white space
        if(line(i:i) .eq. '=') itype = 0   ! white space
        if(line(i:i) .eq. '}') itype = 0   ! white space
        if(line(i:i) .eq. '{') itype = 0   ! white space
        if(line(i:i) .eq. '"') itype = 0   ! white space
        if(line(i:i) .eq. ',') itype = 0   ! white space
        if(line(i:i) .eq. ')') itype = 0   ! white space
        if(line(i:i) .eq. '(') itype = 0   ! white space
        if(.not. instr) then
          if(itype .eq. 1) then
            instr = .true.
            i0 = i
          endif
        endif
        if(instr) then
          if(itype .eq. 1) i1 = i
          if(itype .ne. 1) then
            instr = .false.
            nstr = nstr + 1
            call set_blank(strings(nstr))
            strings(nstr)(1:1+i1-i0) = line(i0:i1)
          endif
        endif
      enddo
      iparse_line = nstr
      return
      end

C=====================================================================72
      subroutine gen_symb_from_name(cname, csymb)
      implicit none
      character*256 cname, csymb
      integer n, i, ic, idash
      character*1 cvec
      n = len(trim(cname))
      ! construct symble from the capical letters
      call set_blank(csymb)
      ic = 0
      idash = -1
      cvec = ' '
      do i = 1, n
        if(cname(i:i) .eq. '-') idash = i
        if(cname(i:i+2) .eq. 'vel') cvec = 'V'
        if(cname(i:i+2) .eq. 'Vel') cvec = 'V'
        if(cname(i:i) .ge. 'A' .and. cname(i:i) .le. 'Z') then
          ic = ic + 1
          csymb(ic:ic) = cname(i:i)
        endif
      enddo
      if(cvec .ne. ' ') then
        csymb(1:1) = cvec
        if(idash .eq. 2) csymb(2:2) = cname(1:1)
        if(idash .gt. 2) csymb(2:2) = cname(idash+1:idash+1)
      endif

      return
      end
C=====================================================================72
      subroutine set_blank(str)
      implicit none
      character*(*) str
      integer i
      do i = 1, len(str)
        str(i:i) = ' '
      enddo
      return
      end
      

