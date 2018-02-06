c  Read adump routines
c  David H. Porter
c  Versions:
c  Date    2010 Feb   1    Fixed offset of izs in read_fields
c  Date    2010 Apr   9    Added get_dxdydz back
c  Date    2012 May  15    addapted to e3d
c  Date    2014 Oct  22    extend to handle arm output
c
      subroutine read_header_info(cFile)
      implicit none
      include 'read_adump_e3d.h'
      character*256 cFile, cline, ckey, cLastFile
      integer i, j, nboxhead
      data nx/-1/

      if(cDumpFile .eq. cfile) return

      call set_read_header_defaults()
      cDumpFile = cFile
      ndvar = 0
      nbdryx    = 0
      nbdryy    = 0
      nbdryz    = 0
      nheadsize = 0
      nboxhead  = 0
      nboxs     = 0
      call setblank(cLastFile)

      ! Parse dump file for available variables and data format & size
      !
      open(unit=17, file=cDumpFile, form='formatted', err=180,
     1     status='old')
      rewind(17)
100   continue
      do i = 1, 256
        cline(i:i) = ' '
      enddo
      read (17,991,end=200) cline
991   format(a256)

      if(cline(1:10).eq.'time      ') read (cline,*) ckey,time
      if(cline(1:10).eq.'dtime     ') read (cline,*) ckey,dtime
      if(cline(1:10).eq.'step      ') read (cline,*) ckey,istep
      if(cline(1:10).eq.'ncube     ') read (cline,*) ckey,ncube
      if(cline(1:10).eq.'i2max     ') read (cline,*) ckey,i2max
      if(cline(1:10).eq.'i2off     ') read (cline,*) ckey,i2off
      if(cline(1:10).eq.'ncubesxyz ') read (cline,*) ckey,ncx,ncy,ncz
      if(cline(1:10).eq.'nbrickxyz ') read (cline,*) ckey,nbx,nby,nbz
      if(cline(1:10).eq.'nbpf      ') read (cline,*) ckey,nbricksperfile
      if(cline(1:10).eq.'xlimits   ') read (cline,*) ckey,xbndyL,xbndyR
      if(cline(1:10).eq.'ylimits   ') read (cline,*) ckey,ybndyL,ybndyR
      if(cline(1:10).eq.'zlimits   ') read (cline,*) ckey,zbndyL,zbndyR
      if(cline(1:10).eq.'boundary  ') then
        read (cline,*) ckey,nbdryx
        nbdryy = nbdryx
        nbdryz = nbdryx
      endif
      if(cline(1:10).eq.'nbdry_xyz ') 
     1     read (cline,*) ckey,nbdryx,nbdryy,nbdryz
      if(cline(1:10).eq.'lobndry   ') then
         read (cline,*) ckey,cLoXbndry,cLoYbndry,cLoZbndry
      endif
      if(cline(1:10).eq.'hibndry   ') then
         read (cline,*) ckey,cHiXbndry,cHiYbndry,cHiZbndry
      endif
      if(cline(1:10).eq.'field3d   ') then
        ndvar = ndvar + 1
        read (cline,*) ckey, cVarSymb(ndvar), cVarName(ndvar),
     1                 cDumpMap(ndvar),dvarmin(ndvar), dvarmax(ndvar)
        if(cVarName(ndvar)(1:2) .eq. "Z-"  .and.
     1    trim(cVarSymb(ndvar)) .eq. "Bx") cVarName(ndvar)(1:1) = 'X'
        if(cVarName(ndvar)(1:2) .eq. "Z-"  .and.
     1    trim(cVarSymb(ndvar)) .eq. "By") cVarName(ndvar)(1:1) = 'Y'
      endif
      if(cline(1:10).eq.'headsize  ') read (cline,*) ckey,nheadsize

      ! For AMR:
      if(cline(1:8).eq.'boxhead ') read (cline,*) ckey, nboxhead
      if(cline(1:5).eq.'boxb ') then
        if(nboxs .eq. 0) call calc_box_size()
        nboxs = nboxs + 1
        call setblank(cboxfile(nboxs))
        read (cline,*) ckey,iboxrefine(nboxs),(iboxoff(j,nboxs),j=1,3),
     1                       cboxfile(nboxs)
        iboxoff(1,nboxs) = nx * iboxoff(1,nboxs)
        iboxoff(2,nboxs) = ny * iboxoff(2,nboxs)
        iboxoff(3,nboxs) = nz * iboxoff(3,nboxs)
        if(trim(cLastFile) .ne. trim(cboxfile(nboxs))) then
          cLastFile = cboxfile(nboxs)
           i64boxoff(nboxs) = nboxhead
        else
           i64boxoff(nboxs) = i64boxoff(nboxs-1) + i64boxsize
        endif
      endif

      if(nheadsize .gt. 0) goto 200
      if(nboxs .gt. 0 .and. cline(1:3).ne.'box') go to 200
      goto 100

180   continue
      write (6,*) "Problem: open failed on dumpfile", cDumpFile
      stop

200   continue
      close(17)

      if(ncube .le. 0) then
        nx = ncx
        ny = ncy
        nz = ncz
      else
        nx = ncx*ncube
        ny = ncy*ncube
        nz = ncz*ncube
      endif

      boxsize(1) = (xbndyR - xbndyL) / float(nbx)
      boxsize(2) = (ybndyR - ybndyL) / float(nby)
      boxsize(3) = (zbndyR - zbndyL) / float(nbz)
      boxoffset(1) = xbndyL
      boxoffset(2) = ybndyL
      boxoffset(3) = zbndyL
      nboxmesh(1) = nx
      nboxmesh(2) = ny
      nboxmesh(3) = nz

      call find_levels()   ! For AMR

      return
      end


c---------------------------------------------------------------------72
      subroutine get_time(current_time)
      implicit none
      real current_time
      include 'read_adump_e3d.h'
      current_time = time
      return
      end

c---------------------------------------------------------------------72
      subroutine get_ndvar(n)
      implicit none
      integer n
      include 'read_adump_e3d.h'
      n = ndvar
      return
      end

c---------------------------------------------------------------------72
      subroutine get_dxdydz(dx, dy, dz)
      implicit none
      real dx, dy, dz
      include 'read_adump_e3d.h'
      dx = (xbndyR-xbndyL) / float(nx*nbx)
      dy = (ybndyR-ybndyL) / float(ny*nby)
      dz = (zbndyR-zbndyL) / float(nz*nbz)
      return
      end

c---------------------------------------------------------------------72
      subroutine get_full_XYZrange(x0,x1,y0,y1,z0,z1)
      implicit none
      real x0,x1,y0,y1,z0,z1
      include 'read_adump_e3d.h'

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
      include 'read_adump_e3d.h'

      if(nx .lt. 0) then
        write (6,*) 'Problem: header info not set'
        stop
      endif

      nx1 = nx*nbx
      ny1 = ny*nby
      nz1 = nz*nbz

      return
      end

c---------------------------------------------------------------------72
      subroutine get_vars_nams(MAX_TOF,ndvar1,cFld_tof,cNam_tof)
      implicit none
      include 'read_adump_e3d.h'
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
      include 'read_adump_e3d.h'
      integer ilen, i

      if(nx .lt. 0) then
        write (6,*) 'Problem: header info not set'
        return
      endif

      ilen = len(trim(cdumpfile))
      write (6,*) 'Data File: ', cdumpfile(1:ilen)
      if(dtime .gt. 0.0001) then
      write (6,9901) time, dtime, istep
9901  format("Time  = ", f10.5,/,"Dtime = ", f10.5,/,"Step  = ", i5)
      else
        write (6,9905) "Time  = ", time
        write (6,9905) "Dtime = ", dtime
        write (6,9906) "Step  = ", istep
9905    format(a, 1pe15.6)
9906    format(a,i10)
      endif
      write (6,9902) 'X',xbndyL,xbndyR, nx, nbx
      write (6,9902) 'Y',ybndyL,ybndyR, ny, nby
      write (6,9902) 'Z',zbndyL,zbndyR, nz, nbz
9902  format(a1,":  [", f9.4,",",f8.4,"]   ",
     1       "Brick Size =",i4,   "    # of Bricks =",i3)
      write (6,9903) nx*nbx, ny*nby, nz*nbz
9903  format("Full Dimensions: ", 3i7)
      if(nbdryx .gt. 0) write (6,9907) nbdryx,nbdryy,nbdryz
9907  format("XYZ Boundary depth brick data:", 3i4)
      if(nboxs .gt. 0) then
        write (6,*) "nboxs = ", nboxs
        write (6,*) "i64boxsize = ", i64boxsize
        call find_levels()

        write (6,9908)
9908    format(/," Level  Refinement   Number of Bricks")
        do i = 1, nlevels
          write (6,9909) i, lev_refine(i), lev_number(i)
9909      format(2i6, 2i12)
        enddo
      endif

      write (6,*)  "  "
      do i=1,ndvar
        write (6,9904) cVarSymb(i)(1:8), cVarName(i)(1:20),
     1                 cDumpMap(i)(1:10),dvarmin(i), dvarmax(i)
      enddo
9904  format(a8, a20, a10, 1p2e15.6)

      return
      end

c---------------------------------------------------------------------72
      subroutine find_levels()
      implicit none
      include 'read_adump_e3d.h'
      integer i, j, ilev
      
      nlevels = 0
      do i = 1, nboxs
        ilev = -1
        do j = 1, nlevels
          if(iboxrefine(i) .eq. lev_refine(j)) ilev = j
        enddo
        if(ilev .lt. 0) then
          nlevels = nlevels + 1
          lev_number(nlevels) = 1
          lev_refine(nlevels) = iboxrefine(i)
        else
          lev_number(ilev) = lev_number(ilev) + 1
        endif
      enddo

      return
      end
c---------------------------------------------------------------------72
      subroutine xrange2ixos(xlo,xhi,ixoff,ixsize)
      implicit none
      include 'read_adump_e3d.h'
      integer ixlo,ixhi,ixoff,ixsize
      real    xlo,xhi,xfrac,dx,eps

      eps = 0.01 * (xbndyR-xbndyL) / float(nx*nbx)
      xfrac = (xlo-xbndyL) / (xbndyR-xbndyL)
      ixlo  = max(1, min(nx*nbx, 1+int(float(nx*nbx)*xfrac)))
      xfrac = max(xfrac, (xhi-xbndyL-eps) / (xbndyR-xbndyL))
      ixhi  = max(1, min(nx*nbx, 1+int(float(nx*nbx)*xfrac)))

      dx = (xbndyR-xbndyL) / float(nx*nbx)
      xlo = xbndyL + dx * (ixlo - 1)
      xhi = xbndyL + dx * ixhi
      ixoff  = ixlo - 1
      ixsize = 1 + ixhi - ixlo

      return
      end

c---------------------------------------------------------------------72
      subroutine yrange2iyos(ylo,yhi,iyoff,iysize)
      implicit none
      include 'read_adump_e3d.h'
      integer iylo,iyhi,iyoff,iysize
      real    ylo,yhi,yfrac,dy,eps

      eps = 0.01 * (ybndyR-ybndyL) / float(ny*nby)
      yfrac = (ylo-ybndyL) / (ybndyR-ybndyL)
      iylo  = max(1, min(ny*nby, 1+int(float(ny*nby)*yfrac)))
      yfrac = max(yfrac, (yhi-ybndyL-eps) / (ybndyR-ybndyL))
      iyhi  = max(1, min(ny*nby, 1+int(float(ny*nby)*yfrac)))

      dy = (ybndyR-ybndyL) / float(ny*nby)
      ylo = ybndyL + dy * (iylo - 1)
      yhi = ybndyL + dy * iyhi
      iyoff  = iylo - 1
      iysize = 1 + iyhi - iylo

      return
      end

c---------------------------------------------------------------------72
      subroutine zrange2izos(zlo,zhi,izoff,izsize)
      implicit none
      include 'read_adump_e3d.h'
      integer izlo,izhi,izoff,izsize
      real    zlo,zhi,zfrac,dz, eps

      eps = 0.01 * (zbndyR-zbndyL) / float(nz*nbz)
      zfrac = (zlo-zbndyL) / (zbndyR-zbndyL)
      izlo  = max(1, min(nz*nbz, 1+int(float(nz*nbz)*zfrac)))
      zfrac = max(zfrac, (zhi-zbndyL-eps) / (zbndyR-zbndyL))
      izhi  = max(1, min(nz*nbz, 1+int(float(nz*nbz)*zfrac)))

      dz = (zbndyR-zbndyL) / float(nz*nbz)
      zlo = zbndyL + dz * (izlo - 1)
      zhi = zbndyL + dz * izhi
      izoff  = izlo - 1
      izsize = 1 + izhi - izlo

      return
      end

c---------------------------------------------------------------------72
      subroutine get_vec_geom(n,vec4)
      integer n, i
      real vec4(4,10000)
      n = 0
      return
      end

c---------------------------------------------------------------------72
      subroutine calc_box_size()
      implicit none
      include 'read_adump_e3d.h'
      integer(kind=8) ivaroff(MAXDVAR), i64size
      integer iv, ibpv

      if(ncube .le. 0) then
        nx = ncx
        ny = ncy
        nz = ncz
      else
        nx = ncx*ncube
        ny = ncy*ncube
        nz = ncz*ncube
      endif

      ivaroff(1) = nheadsize
      do iv = 1, ndvar
        ! Bytes Per Variable, and offset relative to start of box
        ibpv = 2
        if(cDumpMap(iv)(1:10) .eq. "real4     ") ibpv = 4

        i64size = ibpv * (nx+2*nbdryx) * (ny+2*nbdryy) * (nz+2*nbdryz)
        ivaroff(iv+1) = ivaroff(iv) + i64size
      enddo
      i64boxsize = ivaroff(ndvar+1)

      return
      end

c---------------------------------------------------------------------72
      subroutine read_fld_from_one_box(ibox,cVar,iz,nrdx,nrdy,nrdz,fld)
      implicit none
      include 'read_adump_e3d.h'
      character*256 cVar, cHost, cDir, cName
      integer ibox,iz,nrdx,nrdy,nrdz, nvar, ivar(MAXDVAR)
      real fld(nrdx,nrdy,nrdz)
      integer(kind=8) ivaroff(MAXDVAR), i64size, i64off
      real val, fac,off,fldscl(MAXDVAR),fldmin(MAXDVAR),fldmax(MAXDVAR)
      integer i2mid, ilog(MAXDVAR), ibpv(MAXDVAR),ibytezoff
      integer i, n, nbytes, iNameLen, sw_read, ix

      ! File containing this box: host, directory, & file name
      call get_HD(cDumpFile, cHost, cDir, cName, iNameLen)
      call setblank(cName)
      n = len(trim(cBoxFile(ibox)))
      cName(1:n) = trim(cBoxFile(ibox))

      ! Field array to be read in: byte offset in rile and size
      call parse_vars(MAXDVAR,ndvar,cVarsymb,cVar,nvar,ivar)
      call gen_fld_offsets(ivaroff,fldscl,fldmin,fldmax,ilog,ibpv)
      i = ivar(1)
      i64off = i64boxoff(ibox) + ivaroff(i) + ibpv(i)*nrdx*nrdy*(iz-1)
      nbytes = ibpv(i) * nrdx*nrdy*nrdz

      if(ibpv(i) .eq. 4) then
        n = sw_read(cHost,cDir,cName,i64off,nbytes,fld)
c        write (6,*) "i,i64off:",i,i64off
      else
        write (6,*) "read_fld_from_one_box: ipb.ne.4 NOT supported yet"
        stop
      endif

      return
      end

c        write (6,*) "nrdx,nrdy,nrdz:", nrdx,nrdy,nrdz
c        write (6,999) ibox, (fld(ix,1,1),ix=1,5)
c999     format("ibox=",i4, "    fld:", 1p5e12.4)

c---------------------------------------------------------------------72
      subroutine gen_fld_offsets(ivaroff,fldscl,fldmin,fldmax,ilog,ibpv)
      implicit none
      include 'read_adump_e3d.h'
      integer i2mid, ilog(MAXDVAR), ibpv(MAXDVAR),ibytezoff
      integer(kind=8) nbricksize, ibrickinfile, iv
      integer(kind=8) ivaroff(MAXDVAR), i64size, i64off, ibrickoff
      real val, fac,off,fldscl(MAXDVAR),fldmin(MAXDVAR),fldmax(MAXDVAR)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        Derive field scales and offsets, and offsets of data files    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ivaroff(1) = nheadsize
      do iv = 1, ndvar
        fldmin(iv) = dvarmin(iv)
        fldmax(iv) = dvarmax(iv)
        ilog(iv)   = 0
        if(cDumpMap(iv)(1:10) .eq. "Log_e     ") then
          fldmin(iv) = alog(dvarmin(iv))
          fldmax(iv) = alog(dvarmax(iv))
          ilog(iv)   = 1
        endif
        fldscl(iv) = (fldmax(iv)-fldmin(iv))/float(i2max)

        ! Bytes Per Variable, and offset relative to start of brick, of
        ! data in file
        ibpv(iv) = 2
        if(cDumpMap(iv)(1:10) .eq. "real4     ") ibpv(iv) = 4
        i64size = ibpv(iv)*(nx+2*nbdryx)*(ny+2*nbdryy)*(nz+2*nbdryz)  !  For AMR
        ivaroff(iv+1) = ivaroff(iv) + i64size
      enddo

      return
      end

c---------------------------------------------------------------------72
      subroutine read_fields(cVars,nvars,nxsub,nysub,nzsub,
     1                       ixoff,iyoff,izoff, flds)
      implicit none
      character*256 cVars, cHost, cDir, cFile
      integer nvars,nxsub,nysub,nzsub,ixoff,iyoff,izoff
      real flds(nvars,nxsub,nysub,nzsub)
      include 'read_adump_e3d.h'
      integer ncoefs, iilog, i, sw_read
      integer ibrick,ifile,nplanesize
      integer n,nbytes,iv,ix,iy,iz,ixs,iys,izs,ixd,iyd,izd
      integer ibx,ibxoff,ibxmin,ibxmax,ixmin,ixmax
      integer iby,ibyoff,ibymin,ibymax,iymin,iymax
      integer ibz,ibzoff,ibzmin,ibzmax,izmin,izmax
      integer nvar, ivar(MAXDVAR)
      real coefs(10)
      integer*2 slab2i
      real      slab4f
      allocatable  slab2i(:,:,:) ! (ix , iy , iz)
      allocatable  slab4f(:,:,:) ! (ix , iy , iz)

      integer i2mid, ilog(MAXDVAR), ibpv(MAXDVAR),ibytezoff
      integer(kind=8) nbricksize, ibrickinfile
      integer(kind=8) ivaroff(MAXDVAR), i64size, i64off, ibrickoff
      real val, fac,off,fldscl(MAXDVAR),fldmin(MAXDVAR),fldmax(MAXDVAR)
      integer itype
c      real valmax

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        Make sure header info is read in                              c
c        Then parse input string "cVars" for fields to be read         c
c        Then allocate scratch arrays                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(nx .lt. 0) then
        write (6,*) 'Problem: header info not set'
        stop
      endif

      call parse_vars(MAXDVAR,ndvar,cVarSymb,cVars,nvar,ivar)
      allocate (slab2i(1-nbdryx:nx+nbdryx,1-nbdryy:ny+nbdryy,nzsub))
      allocate (slab4f(1-nbdryx:nx+nbdryx,1-nbdryy:ny+nbdryy,nzsub))

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        Derive field scales and offsets, and offsets of data files    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      i2mid = i2max/2 - i2off
      call gen_fld_offsets(ivaroff,fldscl,fldmin,fldmax,ilog,ibpv)
      nbricksize = ivaroff(ndvar+1)


cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Fill fields array in subrange requeested by                      c
c     sweeping through all bricks that overlap with specified region   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ixmin = ixoff + 1            !
      iymin = iyoff + 1            !
      izmin = izoff + 1            ! XYZ ranges in full volume
      ixmax = ixoff + nxsub        !
      iymax = iyoff + nysub        !
      izmax = izoff + nzsub        !

      do ibz=1,nbz
      ibzoff = nz*(ibz-1)
      ibzmin = max( 1+ibzoff, izmin)
      ibzmax = min(nz+ibzoff, izmax)
      if(ibzmin .le. ibzmax) then
         do iby=1,nby
         ibyoff = ny*(iby-1)
         ibymin = max( 1+ibyoff, iymin)
         ibymax = min(ny+ibyoff, iymax)
         if(ibymin .le. ibymax) then
            do ibx=1,nbx
            ibxoff = nx*(ibx-1)
            ibxmin = max( 1+ibxoff, ixmin)
            ibxmax = min(nx+ibxoff, ixmax)
            if(ibxmin .le. ibxmax) then

               ! There is an overlap of spcified range with this brick

               if(nboxs .eq. 0) then
                 ibrick = ibx-1 + nbx*(iby-1 + nby*(ibz-1))   ! starts at 0
                 ifile  = ibrick / nbricksperfile             ! starts at 0
                 ibrickinfile = mod(ibrick, nbricksperfile)   ! starts at 0
                 call get_HDF(ifile, cDumpFile, cHost, cDir, cFile)
                 ibrickoff    = nbricksize * ibrickinfile
               else
                 call lev1_HDF(ibx-1,iby-1,ibz-1,cHost,cDir,cFile,
     1                         ibrickoff)
c                  write (6,*) "-------------------------------------"
c                  write (6,*) "cHost     = ", trim(cHost)
c                  write (6,*) "cDir     = ", trim(cDir)
c                  write (6,*) "cFile     = ", trim(cFile)
c                  write (6,*) "ibrickoff = ", ibrickoff
c                  write (6,*) "ibx,iby,ibz ", ibx,iby,ibz
c                  write (6,*) "-------------------------------------"
               endif
               !
               ! for amr: derive ibrickinfile and cHost,cDir,cFIle from  ibrick
               ! given offsets & file names read it from metadata file

               do i=1,nvar
                 iv = ivar(i)
                 nplanesize   = ibpv(iv) * (nx+2*nbdryx) * (ny+2*nbdryy)
                 ibytezoff    = nplanesize * (nbdryz+mod(ibzmin-1, nz))
                 i64off       = ibrickoff + ivaroff(iv) + ibytezoff
                 nbytes       = (1 + ibzmax - ibzmin) * nplanesize

                 if(ibpv(iv) .eq. 2) then
c                   write (6,*) "i64off,nbytes",i64off,nbytes
                   n = sw_read(cHost,cDir,cFile,i64off,nbytes,slab2i)

                   fac = fldscl(iv)
                   off = 0.5 * (fldmin(iv) + fldmax(iv))
                   iilog = ilog(iv)

c                   valmax = -999.999
                   do iz=ibzmin,ibzmax
c                  izs = 1 + iz - ibzmin
                   izs = iz - ibzmin + 1
                   izd = iz - izoff
                   do iy=ibymin,ibymax
                   iys = iy - ibyoff
                   iyd = iy - iyoff
                   do ix=ibxmin,ibxmax
                   ixs = ix - ibxoff
                   ixd = ix - ixoff
                   val = off + fac*float(int(slab2i(ixs,iys,izs))+i2mid)
                   if(iilog .eq. 1) val = exp(val)
                   flds(i,ixd,iyd,izd) = val
c                   valmax = max(val,valmax)
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c                   if(ixd.eq.2 .and. iyd.eq.1 .and. izd.eq.1) then
c                     write (6,*) 'val,flds = ', val,flds(1,2,1,1)
c                     write (6,*) 'off,fac  = ', off,fac
c                     write (6,*) 'i2mid    = ', i2mid
c                     write (6,*) 'slab2i   = ',int(slab2i(ixs,iys,izs))
c                     write (6,*) 'i2max,i2off=', i2max, i2off
c                     write (6,*) 'i64off   = ', i64off
c                     write (6,*) 'iv   =',iv
c                     write (6,*) 'iXYZs=',ixs,iys,izs
c                     write (6,*) 'ibzmin,ibzoff=',ibzmin,ibzoff
c                   endif
cccccccccccccccccccccccccccccccccccccccccccccccccccc
                   enddo
                   enddo
                   enddo
c                   write (6,*) 'valmax = ', valmax
c >add
                 else
                   n = sw_read(cHost,cDir,cFile,i64off,nbytes,slab4f)
c                   itype = 1
c                   valmax = -999.999
                   do iz=ibzmin,ibzmax
c                  izs = iz - ibzoff
                   izs = iz - ibzmin + 1
                   izd = iz - izoff
                   do iy=ibymin,ibymax
                   iys = iy - ibyoff
                   iyd = iy - iyoff
                   do ix=ibxmin,ibxmax
                   ixs = ix - ibxoff
                   ixd = ix - ixoff
                   flds(i,ixd,iyd,izd) = slab4f(ixs,iys,izs)
c                   valmax = max(slab4f(ixs,iys,izs),valmax)
c                   if(itype .eq. 1) then
c                     write (6,*) "i64off,nbytes",i64off,nbytes
c                     write (6,*) "iXYZ,val=", ixs,iys,izs,
c     1                                 flds(i,ixd,iyd,izd)
c                     itype = 0
c                   endif
                   enddo
                   enddo
                   enddo
c                   write (6,*) 'slab3f: valmax = ', valmax
                 endif
               enddo
            endif
            enddo
         endif
         enddo
      endif
      enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      deallocate (slab2i)
      deallocate (slab4f)

      return
      end

C=====================================================================72
      subroutine lev1_HDF(ibx,iby,ibz,cHost,cDir,cName,ibrickoff)
      implicit none
      include 'read_adump_e3d.h'
      integer ibx, iby, ibz, i, ibox, n, iNameLen, ix,iy,iz
      integer(kind=8) ibrickoff
      character*256 cHost, cDir, cName

      ! cHost & cDir will be that of the meta-data file
      ! So real data file names should be relative paths to meta-file directory
      call get_HD(cDumpFile, cHost, cDir, cName, iNameLen)

      ix = nx * ibx
      iy = ny * iby
      iz = nz * ibz
      do i = 1, nboxs
        if(iboxrefine(i) .eq. 1) then
          if( ix.eq.iboxoff(1,i) .and.
     1        iy.eq.iboxoff(2,i) .and.
     1        iz.eq.iboxoff(3,i)       ) ibox = i
        endif
      enddo

      if(ibox .lt. 0) then
        write (6,*) "PROBLEM in lev1_HDFF:"
        write (6,*) "No box at offset: ", ibx, iby, ibz
        stop
      endif

      call setblank(cName)
      n = len(trim(cBoxFile(ibox)))
      cName(1:n) = trim(cBoxFile(ibox))
      ibrickoff = i64boxoff(ibox)

      return
      end

C=====================================================================72
      subroutine get_HD(cHostDirName, cHost, cDir, cName, iNameLen)
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

      return
      end

C=====================================================================72
      subroutine get_HDF(iFile, cHostDirName, cHost, cDir, cName)
      implicit none
      integer iFile, iNameLen, i1, i2, i3
      character*256 cHostDirName, cHost, cDir, cName

      call get_HD(cHostDirName, cHost, cDir, cName, iNameLen)

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
           endif
        enddo
      enddo

      return
      end
c---------------------------------------------------------------------72
      subroutine set_read_header_defaults()
      implicit none
      include 'read_adump_e3d.h'
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

      return
      end
C=====================================================================72
