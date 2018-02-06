c  Read Tecplot ASXII (pdt) routines
c  David H. Porter
c  Versions:
c  Date    2012 Dec. 25    Adapte to teckplot
c          2013 Mar. 22    multfac & face centered fields
c
c-------------------------------------------------------------
c      call read_header_info(cDumpFile0)  ! Get header info from dump file
c      if(cOutput .eq. "show") then
c        call report_header_info()
c        stop
c      endif
c      if(cOutput .eq. "time") call report_time()
c      call gen_output_par()    ! Generate parameters/dimensions for output structures
c      call gen_ops()           ! Generate seqence ofoperations for requested field
c      call gen_read_bricks() ! Generate sequence of bricks to read in
c       call do_work()
c-------------------------------------------------------------
      subroutine read_header_info(cFile)
      implicit none
      include 'read_pdt_e3d.h'
      character*256 cFile, cline, ckey, cPDTfile, cVarDir, cVarRoots
      real dx, dy
      integer i, j, iField, maxbinvar
      integer iget_cval_from_key_file
      integer iget_fval_from_key_file
      integer iget_ival_from_key_file
      data nx/-1/

      if(cDumpFile .eq. cfile) return

      call setblank(cVarDir)
      call setblank(cVarRoots)
      maxbinvar = 0
      call set_read_header_defaults()
      cDumpFile = cFile
      i = iget_cval_from_key_file(cFile,"pdtfile",cPDTfile)
      if(i .le. 0) stop
      i = iget_cval_from_key_file(cFile,"icp.pdt",cVarDir)
      if(i .gt. 0) maxbinvar = 8000
      i = iget_cval_from_key_file(cFile,"varroots",cVarRoots)
      i = iget_ival_from_key_file(cFile,"maxbinvar",maxbinvar)
      i = iget_fval_from_key_file(cFile,"time",time)
      i = iget_fval_from_key_file(cFile,"dtime",dtime)
      i = iget_ival_from_key_file(cFile,"istep",istep)

      ! Read full ASCII Tecplot data (pdt file)
      call pdt_read(cPDTfile,cVarDir,cVarRoots,maxbinvar,istep)

      i = iget_cval_from_key_file(cFile,"xcoord",cxcoord)
      i = iget_cval_from_key_file(cFile,"ycoord",cycoord)
      i = iget_cval_from_key_file(cFile,"zcoord",czcoord)
      if(i .gt. 0) then
        write (6,*) '3D PDT not supported yet'
        stop
      endif

      call pdt_get_dims(nx,ny,nz)

      ! nx & ny include the extra zone for face centered fields
      ! nz should be 1
      ! calc multiplication factor for low res meshes
      if(nz .eq. 1) then
        multfac = max(1, 256 / min(nx-1, ny-1))
      else
        write (6,*) "3D PDT case not supported yet"
        write (6,*) "nx,ny,nz = ", nx,ny,nz
        stop
      endif

      call pdt_get_range(cxcoord, xBndyL, xBndyR, iField)
      if(iField .gt. 0) then
        dx = (xBndyR-xBndyL) / float((nx-1)*multfac)
      else
        write (6,*) "Coordinate ", trim(cxcoord)," not found"
        stop
      endif

      call pdt_get_range(cycoord, yBndyL, yBndyR, iField)
      if(iField .gt. 0) then
        dy = (yBndyR-yBndyL) / float((ny-1)*multfac)
      else
        write (6,*) "Coordinate ", trim(cxcoord)," not found"
        stop
      endif

      zbndyL = -0.5*dx
      zbndyR =  0.5*dx

      cLoXbndry(1:12) = "continuation"
      cLoYbndry(1:12) = "continuation"
      cLoZbndry(1:12) = "continuation"
      cHiXbndry(1:12) = "continuation"
      cHiYbndry(1:12) = "continuation"
      cHiZbndry(1:12) = "continuation"

      call pdt_get_field_list(cVarName, ndvar)
      do i = 1, ndvar
        cDumpMap(i)(1:5) = "ascii"
        call gen_symb_from_name(cVarName(i),cVarSymb(i))
      enddo

      ! Set times
      if(istep .lt.   0) istep = 0
      if(dtime .lt. 0.0) dtime = 1.0
      if(time  .lt. 0.0) time  = dtime * float(istep) 

      return
      end


c---------------------------------------------------------------------72
      subroutine get_time(current_time)
      implicit none
      real current_time
      include 'read_pdt_e3d.h'
      current_time = time
      return
      end

c---------------------------------------------------------------------72
      subroutine get_ndvar(n)
      implicit none
      integer n
      include 'read_pdt_e3d.h'
      n = ndvar
      return
      end

c---------------------------------------------------------------------72
      subroutine get_dxdydz(dx, dy, dz)
      implicit none
      real dx, dy, dz
      include 'read_pdt_e3d.h'
      dx = (xbndyR-xbndyL) / float((nx-1)*multfac)
      dy = (ybndyR-ybndyL) / float((ny-1)*multfac)
      dz = (zbndyR-zbndyL) / float( nz           )
      return
      end

c---------------------------------------------------------------------72
      subroutine get_full_XYZrange(x0,x1,y0,y1,z0,z1)
      implicit none
      real x0,x1,y0,y1,z0,z1
      include 'read_pdt_e3d.h'

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
      include 'read_pdt_e3d.h'

      if(nx .lt. 0) then
        write (6,*) 'Problem: header info not set'
        stop
      endif

      nx1 = (nx-1)*multfac
      ny1 = (ny-1)*multfac
      nz1 =  nz

      return
      end

c---------------------------------------------------------------------72
      subroutine get_vars_nams(MAX_TOF,ndvar1,cFld_tof,cNam_tof)
      implicit none
      include 'read_pdt_e3d.h'
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
      include 'read_pdt_e3d.h'
      integer ilen, i,j

      if(nx .lt. 0) then
        write (6,*) 'Problem: header info not set'
        return
      endif

      ilen = len(trim(cdumpfile))
      write (6,*) 'Data File: ', cdumpfile(1:ilen)
      write (6,9901) time, dtime, istep
9901  format("Time  = ", f10.5,/,"Dtime = ", f10.5,/,"Step  = ", i5)
      if(nx.gt.1) write (6,9902) trim(cxcoord),xbndyL,xbndyR,(nx-1)*nbx
      if(ny.gt.1) write (6,9902) trim(cycoord),ybndyL,ybndyR,(ny-1)*nby
      if(nz.gt.1) write (6,9902) trim(czcoord),zbndyL,zbndyR,nz*nbz
9902  format(a,":  [", f9.4,",",f8.4,"]     ", "Mesh =",i4)

      do i = 1, ndvar
        call pdt_get_range(trim(cVarName(i)),dvarmin(i),dvarmax(i),j)
        if(j .lt. 0) stop
      enddo

      if(ndvar .gt. 0) then
      write (6,*)  "  "
      write (6,*) "   Symbol      Name",
     1            "                   Min            Max"
      do i=1,ndvar
        write (6,9904) cVarSymb(i)(1:12), cVarName(i)(1:20),
     1                 dvarmin(i), dvarmax(i)
      enddo
9904  format(" ",a15, a20, 1p2e15.6)

      else
        call check_slab()
      endif

      return
      end

c---------------------------------------------------------------------72
      subroutine check_slab()
      implicit none
      include 'read_pdt_e3d.h'
      real rad, con, flx
      allocatable  rad(:,:,:) ! (ix , iy , iz)
      allocatable  con(:,:,:) ! (ix , iy , iz)
      allocatable  flx(:,:,:) ! (ix , iy , iz)
      integer icon, iflx, ix,iy,iz, irad

      allocate (rad(nx,ny,nz))
      allocate (con(nx,ny,nz))
      allocate (flx(nx,ny,nz))
      write (6,*)  "  "
c      call pdt_get_field(   "SI2H6", con, icon)
c      call pdt_get_field("FR-SI2H6", flx, iflx)
      call pdt_get_field( "R (cm)", rad, irad)
      call pdt_get_field(   "AR3S", con, icon)
      call pdt_get_field("FR-AR3S", flx, iflx)
      write (6,*) "irad = ", irad
      write (6,*) "icon = ", icon
      write (6,*) "iflx = ", iflx

      iz = 1
      iy = 14
      iy = ny/2

999   format(3i4, 1p8e12.4)
998   format(i4, 1p9e12.4)
c      do iy=ny,ny-7,-1    ! ny
      do ix=1,nx
        write (6,999) ix,iy,iz,rad(ix,iy,iz),con(ix,iy,iz),flx(ix,iy,iz)
c        write (6,998) ix,(flx(ix,iy,iz),iy=10,ny,30)
      enddo

c      write (6,*) "  "
c      do iy=8,1,-1    ! ny
c        write (6,999) ix,iy,iz,con(ix,iy,iz), flx(ix,iy,iz)
c      enddo

      deallocate (flx)
      deallocate (con)

      return
      end

c---------------------------------------------------------------------72
      subroutine xrange2ixos(xlo,xhi,ixoff,ixsize)
      implicit none
      include 'read_pdt_e3d.h'
      integer ixlo,ixhi,ixoff,ixsize, nxm
      real    xlo,xhi,xfrac,dx,eps

      nxm = (nx-1)*multfac
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
      include 'read_pdt_e3d.h'
      integer iylo,iyhi,iyoff,iysize, nym
      real    ylo,yhi,yfrac,dy,eps

      nym = (ny-1)*multfac
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
      include 'read_pdt_e3d.h'
      integer izlo,izhi,izoff,izsize, nzm
      real    zlo,zhi,zfrac,dz, eps

      nzm = 1
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
      write (6,*) "read_fld_from_one_box not implmented for techplot"
      stop
      end

c---------------------------------------------------------------------72
      subroutine read_fields(cVars,nvars,nxsub,nysub,nzsub,
     1                       ixoff,iyoff,izoff, flds)
      implicit none
      character*256 cVars
      integer nvars,nxsub,nysub,nzsub,ixoff,iyoff,izoff,j
      real flds(nvars,nxsub,nysub,nzsub)
      include 'read_pdt_e3d.h'
      integer iv,ix,iy,iz
      integer nvar, ivar(MAXDVAR), iField
      integer iface_centered, ixface, iyface
      real      slab4f, slab4fb, xoff, yoff
      integer ixs0,ixs1,iys0,iys1,izs
      integer iyd,iyo,ixd,ixo, ixmax, iymax
      real fx0,fx1,fy0,fy1, fmid, dmult, triangle_interp, bicubic_eval
      real fx, fy, fxy, a, bilinear_interp
      allocatable  slab4f(:,:,:) ! (ix , iy , iz)
      allocatable  slab4fb(:,:,:) ! (ix , iy , iz)
      allocatable  fx(:,:,:) ! (ix , iy , iz)
      allocatable  fy(:,:,:) ! (ix , iy , iz)
      allocatable  fxy(:,:,:) ! (ix , iy , iz)
      allocatable    a(:,:,:) ! (icoef, ix , iy)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        Make sure header info is read in                              c
c        Then parse input string "cVars" for fields to be read         c
c        Then allocate scratch array, read data, copy into subrange    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      if(nx .lt. 0) then
        write (6,*) 'Problem: header info not set'
        stop
      endif

      call parse_vars(MAXDVAR,ndvar,cVarSymb,cVars,nvar,ivar)
      if(nvar .ne. nvars) then
        write (6,*) "PROBLEM: could not find variables"
        write (6,*) "   ", trim(cVars)
        write (6,*) "in tecplot variable list"
        write (6,*) "  Number of variables rquested: nvars = ", nvars
        write (6,*) "  Number of variables found   : nvar  = ", nvar
        stop
      endif

      allocate (slab4f(nx,ny,nz))
      allocate (slab4fb(-1:nx+2,-1:ny+1,2))
      allocate ( fx(0:nx+1,0:ny+1,1))
      allocate ( fy(0:nx+1,0:ny+1,1))
      allocate (fxy(0:nx+1,0:ny+1,1))
      allocate (a(16,0:nx+1,0:ny+1))

      do iv = 1, nvar
         ! Get full field
c parse_var
c          write (6,*) "iv,ivar(iv)=",iv,ivar(iv)
         call pdt_get_field(cVarName(ivar(iv)), slab4f, iField)
         if(iField .lt. 1) then
          write (6,*) "PROBLEM: could not find variable"
          write (6,*) "   ", trim(cVarName(ivar(iv)))
          write (6,*) "in tecplot file"
          stop
         endif

      call fill_bndy(nx,ny,slab4f,slab4fb,cVarName(ivar(iv)))
      call calc_slopes(nx,ny,slab4fb,fx,fy,fxy)
      do iy=0,ny
      do ix=0,nx
        call bicubic_coefs(nx+1,slab4fb(ix,iy,1),fx(ix,iy,1),
     1                     fy(ix,iy,1),fxy(ix,iy,1),a(1,ix,iy))
      enddo
      enddo

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      open(unit=25,file='xy.dat',form='formatted')
c      do iys0 = 40, 70
c        write (25,*) iys0,(slab4fb(ixs0,iys0,1),ixs0=1,3)
c      enddo
c      close(25)
c                                 DEBUG
c      ixs0 = 1
c      iys0 = 50
c      write (6,*) " "
c      write (6,*) " "
c      write (6,*) "f  : a(1) ", a(1,ixs0,iys0)
c      write (6,*) "fx : a(2) ", a(2,ixs0,iys0)
c      write (6,*) "fy : a(5) ", a(5,ixs0,iys0)
c      write (6,*) "fxy: a(6) ", a(6,ixs0,iys0)
c      write (6,*) " "
c      ixs0 = 1
c      do iys0=49,52
c        write (6,9001) float(iys0-50),slab4fb(ixs0,iys0,1)
c9001    format(1p2e12.4)
c      enddo
c 
c      write (6,*) " "
c      iys0 = 50
c      fx1 = 0.0
c      do iy = 0, 40
c        fy1 = 0.025 * float(iy)
c        write (6,9001) fy1, bicubic_eval(fx1,fy1,a(1,ixs0,iys0))
c      enddo
c      write (6,*) " "
c      write (6,*) " "
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      if(iField .lt. 100000) then
c        write (6,*) "cVarName = ", trim(cVarName(iVar(iv)))
c        do ix=1,10
c          write (6,*) "ix,slab4f(ix,1,1)", ix,slab4f(ix,1,1)
c        enddo
c        stop
c      endif

          ! Determine if face or zone centered
          j = ivar(iv)
c          ixface = iface_centered(cVarName(j),cxcoord(1:1),cycoord(1:1))
c          iyface = iface_centered(cVarName(j),cycoord(1:1),cxcoord(1:1))
          ixface = 0   ! was 1
          iyface = 0   ! was 1
          if(cVarName(j)(1:2) .eq. "R ") ixface = 1
          if(cVarName(j)(1:2) .eq. "FZ") iyface = 1
c          if(cVarName(j)(1:2) .eq. "FR") ixface = 1
           ! subroutine

         ! copy subdomain into fields array
         dmult = 1.0 / float(multfac)

         if(iyface .eq. 0) then
ccc           iymax = ny-1
ccc           iyo   = multfac/2
ccc           yoff  = 1.0

           iymax = ny
           iyo   = multfac + multfac/2
           yoff  = 0.0
         else
           iymax = ny
           iyo   = multfac
           yoff  = 0.5
         endif

         if(ixface .eq. 0) then
           ixmax = nx-1
           ixo   = multfac/2
           xoff = 1.0
         else
           ixmax = nx
           ixo   = multfac
           xoff = 0.5
         endif

c         open(unit=25,file='xyf.dat',form='formatted')
         do iz = 1,nzsub
         izs = iz+izoff

c         write (6,*) "nysub = ", nysub
c         write (6,*) "iyoff = ", iyoff
         do iy = 1,nysub
         iyd = iyo + iy + iyoff - 1
         iys0 = iyd / multfac
         fy1 = dmult * (0.5+float(iyd-multfac*iys0))
ccccc         fy0 = 1.0 - fy1
         iys1 = min(iys0 + 1, iymax)
c         iys0 = max(iys0, 1)

         do ix = 1,nxsub
         ixd = ixo + ix + ixoff - 1
         ixs0 = ixd / multfac
         fx1 = dmult * (xoff+float(ixd-multfac*ixs0))
ccccc         fx0 = 1.0 - fx1
         ixs1 = min(ixs0 + 1, ixmax)
c         ixs0 = max(ixs0, 1)

c           flds(iv,ix,iy,iz) = triangle_interp(slab4fb(ixs0,iys0,izs),
c     1                                         slab4fb(ixs1,iys0,izs),
c     1                                         slab4fb(ixs0,iys1,izs),
c     1                                         slab4fb(ixs1,iys1,izs),
c     1                                         fx1, fy1)
           flds(iv,ix,iy,iz) = bilinear_interp(slab4fb(ixs0,iys0,izs),
     1                                         slab4fb(ixs1,iys0,izs),
     1                                         slab4fb(ixs0,iys1,izs),
     1                                         slab4fb(ixs1,iys1,izs),
     1                                         fx1, fy1)

c        flds(iv,ix,iy,iz) = bicubic_eval(fx1, fy1, a(1,ixs0,iys0))

c          if(ix.eq.1 .and. iy.lt.30) then
c             write (6,*) 'iy,flds:', iy,flds(iv,ix,iy,iz)
c          endif

         ! if(ix .eq. 14) this is right on ixs0=2
         if(ix .eq. 10) then
c         if(iy .eq. 121) then
c         if(iy .eq. 500) then
         if(iys0 .ge. 50 .and. iys0 .le. 70) then
c         if(ix .ge. 183  .and.  ix .le. 195) then
c           write (6,982) iys0,fx1,fy1, flds(iv,ix,iy,iz)
c982        format("iys0,fx1,fy1,f",i6, 2f10.5, 1pe12.4)
cccccccccccccccccccccc
c            write (25,982) fx1+float(ixs0),fy1+float(iys0),
c     1                      flds(iv,ix,iy,iz)
c982        format(1p3e12.4)
cccccccccccccccccccccc
c           write (6,983) slab4f(ixs0,iys0,izs),
c     1                   slab4f(ixs1,iys0,izs),
c     1                   slab4f(ixs0,iys1,izs),
c     1                   slab4f(ixs1,iys1,izs)
c
c983        format("f00,f10,f01,f11",1p4e12.4)
c           write (6,*) "  "
         endif
         endif

         enddo
         enddo
         enddo
c         close(25)
      enddo

      deallocate (slab4f)
      deallocate (slab4fb)
      deallocate (fx)
      deallocate (fy)
      deallocate (fxy)


      return
      end

c---------------------------------------------------------------------72
      real function bicubic_eval(x, y, a)
      implicit none
      integer n, m
      real x, y, a(0:3,0:3), sum
      sum = 0.0
      do m= 0, 3
      do n= 0, 3
        sum = sum + a(n,m) * (x**n) * (y**m)
      enddo
      enddo
      bicubic_eval = sum
      return
      end
c---------------------------------------------------------------------72
      subroutine bicubic_coefs(n,f,fx,fy,fxy,a)
      implicit none
      integer n
      real f(0:n+2,0:1), fx(0:n,0:1), fy(0:n,0:1),fxy(0:n,0:1)
      real a(0:3,0:3)
      real sa, sax, say, saxy, fs, fsx, fsy, fsxy

      ! f(x,y) = sum_n=0,3 { sum_m=0,3 { a(n,m)* x**n * y**m }}
      a(0,0) = f(0,0)      ! a00
      a(1,0) = fx(0,0)     ! a10
      a(0,1) = fy(0,0)     ! a01
      a(1,1) = fxy(0,0)    ! a11

      ! f10 = a00 + a10 +   a20 +   a30
      ! fx10 =      a10 + 2*a20 + 3*a30
      ! (fx10-a10) - 2*(f10-a00-a10) = a30
      a(3,0) = fx(1,0)-a(1,0) - 2.0*(f(1,0)-a(0,0)-a(1,0))
      a(2,0) = f(1,0) - (a(0,0)+a(1,0)+a(3,0))

      a(0,3) = fy(0,1)-a(0,1) - 2.0*(f(0,1)-a(0,0)-a(0,1))
      a(0,2) = f(0,1) - (a(0,0)+a(0,1)+a(0,3))

      ! fy10  = a01 + a11 +   a21 +   a31 
      ! fxy10 =       a11 + 2*a21 + 3*a31 
      a(3,1) = fxy(1,0)-a(1,1) - 2.0*(fy(1,0)-a(0,1)-a(1,1))
      a(2,1) = fy(1,0) - (a(0,1)+a(1,1)+a(3,1))

      ! fx01  = a10 + a11 +   a12 +   a13 
      ! fxy01 =       a11 + 2*a12 + 3*a13 
      a(1,3) = fxy(0,1)-a(1,1) - 2.0*(fx(0,1)-a(1,0)-a(1,1))
      a(1,2) = fx(0,1) - (a(1,0)+a(1,1)+a(1,3))


      sa  = a(0,0) + a(0,1) + a(0,2) + a(0,3)
     1    + a(1,0) + a(1,1) + a(1,2) + a(1,3)
     1    + a(2,0) + a(2,1)
     1    + a(3,0) + a(3,1)
     
      sax =     a(1,0) +     a(1,1) + a(1,2) + a(1,3)
     1    + 2.0*a(2,0) + 2.0*a(2,1)
     1    + 3.0*a(3,0) + 3.0*a(3,1)

      say = a(0,1) + 2.0*a(0,2) + 3.0*a(0,3)
     1    + a(1,1) + 2.0*a(1,2) + 3.0*a(1,3)
     1    + a(2,1)
     1    + a(3,1)

      saxy =     a(1,1) + 2.0*a(1,2) + 3.0*a(1,3)
     1     + 2.0*a(2,1)
     1     + 3.0*a(3,1)

      fs   = f(1,1) - sa
      fsx  = fx(1,1) - sax
      fsy  = fy(1,1) - say
      fsxy = fxy(1,1) - saxy

      ! f(1,1)   = sa   +   a(2,2) +   a(2,3) +   a(3,2) +   a(3,3)
      ! fx(1,1)  = sax  + 2*a(2,2) + 2*a(2,3) + 3*a(3,2) + 3*a(3,3)
      ! fy(1,1)  = say  + 2*a(2,2) + 3*a(2,3) + 2*a(3,2) + 3*a(3,3)
      ! fxy(1,1) = saxy + 4*a(2,2) + 6*a(2,3) + 6*a(3,2) + 9*a(3,3)
      !
      ! fsx-2*fs   =   a(3,2) +   a(3,3)
      ! fsxy-2*fsy = 2*a(3,2) + 3*a(3,3)
      ! fsxy-2*fsy - 2*(fsx-2*fs) = a(3,3)
      a(3,3) = fsxy-2.0*fsy - 2.0*(fsx-2*fs)
      a(3,2) = fsx-2.0*fs - a(3,3)

      ! fsy-2*fs = a(2,3) + a(3,3)
      a(2,3) = fsy-2.0*fs - a(3,3)
      a(2,2) = fs - (a(2,3) + a(3,2) + a(3,3))

      return
      end
c---------------------------------------------------------------------72
      real function bilinear_interp(v00,v10,v01,v11,fx1,fy1)
      implicit none
      real v00,v01,v10,v11,fx1,fy1, fx0,fy0
      fx0 = 1.0 - fx1
      fy0 = 1.0 - fy1
      bilinear_interp = fx0 * fy0 * v00 + fx0 * fy1 * v01
     1                + fx1 * fy0 * v10 + fx1 * fy1 * v11
      return
      end

c---------------------------------------------------------------------72
      real function triangle_interp(v00,v10,v01,v11,fx1,fy1)
      implicit none
      real v00,v01,v10,v11,fx1,fy1
      real a0011, a0110, d0011, d0110, vmm, fx0, fy0, v, f0110, f0011

      ! Set vmm to minimize variation
      a0011 = 0.5 * (v00 + v11)
      a0110 = 0.5 * (v01 + v10)
      d0011 = abs(v00-v11) + min(abs(v01-a0011),abs(v10-a0011))
      d0110 = abs(v01-v10) + min(abs(v00-a0110),abs(v11-a0110))

      if(d0011+d0110 .eq. 0.0) then
        triangle_interp = v00
        return
      endif
      f0110 = d0110 / (d0011+d0110)
      f0011 = d0011 / (d0011+d0110)
      vmm = 0.5*(f0110*(v00+v11) + f0011*(v01+v10))

      fx0 = 1.0 - fx1
      fy0 = 1.0 - fy1
      if(fx1 .gt. fy1) then
        if(fx1 .gt. fy0) then
          ! use v10, v11, vmm
          v = fy0*v10+fy1*v11 + 2.0*fx0*(vmm-0.5*(v10+v11))
        else
          ! use v00, v10, vmm
          v = fx0*v00+fx1*v10 + 2.0*fy1*(vmm-0.5*(v00+v10))
        endif
      else
        if(fx1 .gt. fy0) then
          ! use v01, v11, vmm
          v = fx0*v01+fx1*v11 + 2.0*fy0*(vmm-0.5*(v01+v11))
        else
          ! use v00, v01, vmm
          v = fy0*v00+fy1*v01 + 2.0*fx1*(vmm-0.5*(v00+v01))
        endif
      endif
      triangle_interp = v
      return
      end
c---------------------------------------------------------------------72
      integer function iface_centered(cVar1, c1, c2)
      implicit none
      include 'read_pdt_e3d.h'
      character*(*) cVar1
      character*(20) cVar2
      character*1 c1, c2
      integer i, j, n, iface
      do i=20,1,-1
        if(cVar1(i:i) .eq. ' ') n = i   ! length of string including space
      enddo

      iface = 0    ! zone centered if no coordinate match is found
      do j = 1, n
       if(cVar1(j:j) .eq. c1) then
         cVar2(1:n) = cVar1(1:n)
         cVar2(j:j) = c2
         do i = 1, ndvar
           if(cVar2(1:n) .eq. cVarName(i)(1:n)) iface = 1 ! is face cenered
         enddo
       endif
      enddo
      iface_centered = iface

      return
      end

C=====================================================================72
      subroutine calc_slopes(nx,ny,f,fx,fy,fxy)
      implicit none
      integer nx,ny, ix,iy
      real f(-1:nx+2,-1:ny+2)
      real fx(0:nx+1,0:ny+1),fy(0:nx+1,0:ny+1),fxy(0:nx+1,0:ny+1)
      logical bzero

      do iy=0,ny+1
      bzero = .true.
      do ix=nx+1,0,-1
         fx(ix,iy) = 0.5  * (f(ix+1,iy  ) - f(ix-1,iy  ))
         fy(ix,iy) = 0.5  * (f(ix  ,iy+1) - f(ix  ,iy-1))
        fxy(ix,iy) = 0.25 * ((f(ix+1,iy+1) - f(ix-1,iy+1)) -
     1                       (f(ix+1,iy-1) - f(ix-1,iy-1))  )

         ! if everything at and to right of this pt. is 0, then slopes
         ! must be zero (to prevent overshoots)
         if(f(ix,iy) .ne. 0.0) bzero = .false.
         if(bzero) then
             fx(ix,iy) = 0.0
             fy(ix,iy) = 0.0
            fxy(ix,iy) = 0.0
         endif
      enddo
      enddo


      return
      end
C=====================================================================72
      subroutine fill_bndy(nr,nz,slab4f,slab4fb,cVar)
      implicit none
      integer nr,nz, izd,izs,ird,irs
      real slab4f(nr,nz),slab4fb(-1:nr+2,-1:nz+2), fac
      character*(*) cVar

      ! Continuation BCs
      do izd=-1,nz+2
        izs = max(1,min(nz,izd))
        do ird=-1,nr+2
          irs = max(1,min(nr,ird))
          slab4fb(ird,izd) = slab4f(irs,izs)
        enddo
      enddo

      ! Mirror BCs at R=0
      fac = 1.0
      if(cVar(1:2) .eq. "FR" .or. cVar(1:2) .eq. "VR") fac = -1.0
       do izd=0,nz+1
         slab4fb( 0,izd) = fac * slab4fb(1,izd)
         slab4fb(-1,izd) = fac * slab4fb(2,izd)
       enddo

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
      subroutine gen_symb_from_name(cName,cSymb)
      character*(*) cName, cSymb
      character*256 cTmp
      integer n

      call setblank(cSymb)
      n = len(trim(cName))
      cSymb(1:n) = trim(cName)
c10    cSymb(1:n) = trim(cName)
c      if(cName(1:1) .eq. ' ') then
c        cName(1:n-1) = cSymb(2:n)
c        cName(n:n)   = ' '
c        n = n - 1
c        go to 10
c      endif

      do i = 1,n
        if(cSymb(i:i) .eq. " ") cSymb(i:i) = '_'
        if(cSymb(i:i) .eq. "(") cSymb(i:i) = '_'
        if(cSymb(i:i) .eq. ")") cSymb(i:i) = '_'
        if(cSymb(i:i) .eq. "-") cSymb(i:i) = '_'
        if(cSymb(i:i) .eq. "*") cSymb(i:i) = 's'
        if(cSymb(i:i) .eq. "+") cSymb(i:i) = 'p'
        if(cSymb(i:i) .eq. "^") cSymb(i:i) = 'p'
      enddo
      return
      end

C=====================================================================72
      subroutine parse_vars(MAXDVAR,ndvar,cVarSymb,cVars,nvar,ivar)
      character*256 cVarSymb(MAXDVAR)
      character*256 cVars
      integer ndvar, nvar, ivar(MAXDVAR)

      integer i, j, istr, nstr
      character*256  cstr(8000)

      ! Initialize strings
      do j = 1, 8000
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
      include 'read_pdt_e3d.h'
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
      time      = -999.0
      dtime     = -999.0
      istep     = -999
      nbx =1
      nby = 1
      nbz = 1
      nbricksperfile = 1

      return
      end
C=====================================================================72


C=====================================================================72
C          TECPLOT PDT ASCII FILES
C=====================================================================72
      module tecdat
      integer       :: iverbose
      integer       ::  nx,ny,nz, nvar, izoned(3), ndims, isallocated=0
      integer       :: igotdata, nstrokes
      real          :: stroke(4,10000)
      character*256 :: cFileName, cTitle, cVars(8000), cForm
      character*256 :: czonef, czonet
      real, allocatable :: data(:,:,:,:)
      end module

C=====================================================================72
      ! part of external interface
      subroutine get_vec_geom(n,vec4)
      use tecdat
      integer n, i
      real vec4(4,10000)
      n = nstrokes
      do i = 1, n
      do j = 1, 4
        vec4(j,i) = stroke(j,i)
      enddo
      enddo
      return
      end

C=====================================================================72
      subroutine pdt_read(cFile,cVarDir,cVarRoots,maxbinvar,iteration)
      use tecdat
      character*(*) cFile, cVarDir, cVarRoots
      integer maxbinvar, iteration

       ! scratch
      integer nlines, nlen, igetline, nstr
      integer i, igetstrings, istr1(8000), istr2(8000), itype, itypecur
      integer i1, i2, j1,j2, nr, igettype, ix,iy,iz,iv
      character*256  cline

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      cFileName(1:len(cFile)) = cFile
      ! Initialize
      iverbose    = 0
      !  isallocated = 0
      nvar = 0
      itypecur = 0    ! 0=none; 1=title; 2=variables; 3=zone; 4=geometry
      nstrokes = 0    ! number of "stokes" in geometry vector data

      nlines=0
      open(1,file=cFileName,status='old')
10    nlen = igetline(cline)
c      if(nlen .eq. 0) go to 10
c      if(nlen .lt. 0) go to 100
      nlines=nlines+1
      nstr = igetstrings(cline, istr1, istr2)
      i = 1                                         !  do i = 1, nstr
20    continue
c
c TITLE = "ITER= 200 : CDI_04 Ar/SiH4
c
      i1 = istr1(i)
      i2 = istr2(i)
      itype = igettype(cline(i1:i2))
      if(itype .ne. 0) then
c          write (6,*) "type=",itype,  cline(i1:i2)
        itypecur = itype
        if(itype .eq. 2) nvar = 0
        if(itype .eq. 3) then    ! reset zone properties
          nx = 1
          ny = 1
          nz = 1
          ndims = 0
          igotdata = 0
          call setblank(czonet)
          call setblank(czonef)
        endif
      else
        if(itypecur .eq. 1) then
          call strcpy(cTitle,cline(i1:i2))
          if(cTitle(1:5) .eq. "ITER=") read (cTitle(6:10),*) iteration
        endif
        if(itypecur .eq. 2) then
          nvar = nvar + 1
          if(nvar .gt. 8000) then
            write (6,*) "Number of variables exceeds limit of 8000"
            stop
          endif
          call strcpy(cVars(nvar), cline(i1:i2))
        endif
        if(itypecur .eq. 3) then
          j1 = istr1(i+1)
          j2 = istr2(i+1)
          nr = 0
          if(cline(i1:i2) .eq. "I") call geti(cline(j1:j2),nx,nr)
          if(cline(i1:i2) .eq. "J") call geti(cline(j1:j2),ny,nr)
          if(cline(i1:i2) .eq. "K") call geti(cline(j1:j2),nz,nr)
          if(cline(i1:i1) .eq. "F") call getc(cline(j1:j2),czonef,nr)
          if(cline(i1:i1) .eq. "T") call getc(cline(j1:j2),czonet,nr)
          if(cline(i1:i1) .eq. "D") call getd(cline,i,nstr,istr1,istr2,
     1                                        izoned,ndims,nr)
          i = i + nr
          if(nr .eq. 0 .and. igotdata .eq. 0) then
            if(isallocated .eq. 1) deallocate(data)
            allocate(data(nx,ny,nz,nvar+maxbinvar))
            isallocated = 1
            call getdata(nx,ny,nz,nvar,data, i, nstr, cline)
            igotdata = 1
c            write (6,*) 'i,nstr = ', i,nstr
          endif
        endif
cccccccccccccccccccccccccc geom
        if(itypecur .eq. 4) then
          call get_geometry(nstrokes,stroke)
          itypecur = 0
        endif
cccccccccccccccccccccccccc geom
      endif
      i = i + 1
      if(i .le. nstr) go to 20                      ! enddo
      if(nlen .ge. 0) go to 10
100   continue
      close(1)

      if(maxbinvar .gt. 0)
     1  call pdt_read_icp(cVarDir,cVarRoots,maxbinvar)

      return
      end

C=====================================================================72
      subroutine pdt_read_icp(cVarDir,cVarRoots,maxbinvar)
      use tecdat
      character*(*) cVarDir, cVarRoots
      integer maxbinvar, n
      character*256 cInFile, cOutFile, cline
      integer i1(8000), i2(8000), j1(8000), j2(8000)

      nnonspace = 0

      n = len(trim(cVarDir))
      call setblank(cInFile)
      cInFile(1:n+8) = trim(cVarDir) // "/icp.pdt"
      cOutFile(1:n+15) = trim(cVarDir) // "/binary_icp.bin"
      nVarRoots = igetstrings(cVarRoots, i1, i2)
      do i = 1, nVarRoots
      enddo
      open(unit=1,file=cInFile,form="formatted",status="old",err=9000)
10    if(igetline(cline) .lt. 0) go to 100
        if(cline(1:1) .ne. ' ') then
          nnonspace = nnonspace + 1
          do i=1,nVarRoots
            n = 1+i2(i)-i1(i)
            if(cline(1:n) .eq. cVarRoots(i1(i):i2(i))) then
              nvar = nvar + 1
              nstr = igetstrings(cline, j1, j2)
              call strcpy(cVars(nvar), cline(j1(1):j2(1)))
              call read_icp_values()
            endif
          enddo
        endif
      go to 10
100   continue
      close(1)
      return

9000  write (6,*) "Can not open icp.dat file: ",trim(cInFile)
      stop
      end
C=====================================================================72
      subroutine read_icp_values()
      use tecdat
      integer nxy, ic, ix, iy, iz, n, nstr, is, i1(8000), i2(8000)
      real fmin, fmax
      character*256 cline

      iz = 1
      ic = 0
      nxy = nx * (ny-1)
10    n = igetline(cline)
      nstr = igetstrings(cline, i1, i2)
      is = 0
20    is = is + 1

      ic = ic + 1
      iy = 1 + (ic-1) / nx
      ix = ic - nx*(iy-1)
      if(i1(is) .lt. 0) then
        write (6,*) "is,nstr = ", is,nstr
        write (6,*) "i1(is),i2(is) = ", i1(is),i2(is)
        write (6,*) "ic = ", ic
        write (6,*) "ix,iy,iz = ", ix,iy,iz
        stop
      endif
      read (cline(i1(is):i2(is)),*) data(ix,iy,iz,nvar)
      if(ic .lt. nxy  .and.  is .lt. nstr) go to 20
      if(ic .lt. nxy  .and.  is .ge. nstr) go to 10

c module  symbol   

      return
      end
C=====================================================================72
      subroutine pdt_header_info()
      use tecdat
      ! Debug
      print *
      write (6,*) 'TITLE = -->', trim(ctitle), "<--" 
      print *
      write (6,*) "Variables:"
      do i1=1,nvar,5
        i2 = min(nvar, i1+4)
        write (6,999) (trim(cVars(i)),i=i1,i2)
999     format(5a15)
      enddo
      print *
      write (6,*) "nx,ny,nz = ", nx, ny, nz
      write (6,*) "ZONE title  = ", trim(czonet)
      write (6,*) "ZONE format = ", trim(czonef)
      write (6,*) "ZONE dims = ", (izoned(i),i=1,ndims)
      return
      end

C=====================================================================72
      subroutine pdt_get_dims(nx1,ny1,nz1)
      use tecdat
      integer nx1, ny1, nz1
      nx1 = nx
      ny1 = ny
      nz1 = nz
      return
      end
C=====================================================================72
      subroutine pdt_get_range(cField,fmin, fmax, iField)
      use tecdat
      character*(*) cField
      real fmin, fmax
      integer i, iField
      iField = -1
      do i=1,nvar
        if(trim(cfield) .eq. trim(cVars(i))) iField = i
      enddo
      if(iField .lt. 0) then
        write (6,*) "This field is not available: ", trim(cField)
        return
      endif

      if(iverbose .gt. 0) then
        write (6,*) 'Getting range for: ', trim(cField)
        write (6,*) 'Array number : ', iField
      endif

      fmin = data(1,1,1,iField)
      fmax = data(1,1,1,iField)
      do iz=1,nz
      do iy=1,ny
      do ix=1,nx
        fmin = min(fmin, data(ix,iy,iz,ifield))
        fmax = max(fmax, data(ix,iy,iz,ifield))
      enddo
      enddo
      enddo

      return
      end
C=====================================================================72
      subroutine pdt_get_field(cField, field, iField)
      use tecdat
      character*(*) cField
      real field(nx,ny,nz), fac
      integer i, iField
      iField = -1
      do i=1,nvar
        if(trim(cfield) .eq. trim(cVars(i))) iField = i
c        if(i .lt. 20) write (6,*) "i,CVar:",i,trim(cVars(i))
      enddo
      if(iField .lt. 0) then
        write (6,*) "This field is not available: ", trim(cField)
        return
      endif

      if(iverbose .gt. 0) then
        write (6,*) 'Getting Field: ', trim(cField)
        write (6,*) 'Array number : ', iField
      endif

      do iz=1,nz
      do iy=1,ny
      do ix=1,nx
        field(ix,iy,iz) = data(ix,iy,iz,ifield)
      enddo
      enddo
      enddo

      ! Fill continuation lower boundaries when needed
      do iz=1,nz
      do iy=3,1,-1
      do ix=1,nx
        if(field(ix,iy,iz) .eq. 0.0) field(ix,iy,iz) = field(ix,iy+1,iz)
      enddo
      enddo
      enddo

      ! Fill continuation upper boundaries when needed
      do iz=1,nz
      do iy=ny-1,ny
      do ix=1,nx
        if(field(ix,iy,iz) .eq. 0.0) field(ix,iy,iz) = field(ix,iy-1,iz)
      enddo
      enddo
      enddo


c      fac = 0.0
c      do iy=1,50
c      do ix=1,9
c        fac = max(fac, field(ix,iy,1))
c      enddo
c      enddo
c      fac = 100.0 / fac
   
c      do iy=ny,ny-2,-1
c        write (6,999) iy, (field(ix,iy,1),ix=1,7)
c      enddo
c      write (6,*) "..."
c      do iy=4,1,-1
c        write (6,999) iy, (field(ix,iy,1),ix=1,7)
c      enddo
c999   format("iy="i4,"    data=" 1p7e12.4)

c        write (6,999) (fac*field(ix,iy,1),ix=1,9)
c999     format(10f9.3)

c      write (6,*) "..."
c      do iy=ny-20,ny
c        write (6,*) "iy,field(2,iy,1)=",iy,field(2,iy,1)
c      enddo

      return
      end

C=====================================================================72
      subroutine pdt_fix_cvar_string(cvar)
      character*(*) cVar
      character*256 cTmp

      n = len(trim(cVar))
10    cTmp(1:n) = trim(cVar)
      if(cVar(1:1) .eq. ' ') then
        cVar(1:n-1) = cTmp(2:n)
        cVar(n:n)   = ' '
        n = n - 1
        go to 10
      endif

      return
      end

C=====================================================================72
      subroutine pdt_get_field_list(cFields, nFields)
      use tecdat
      character*(*) cFields(*)
      integer i, n, nFields

      nfields = nvar
      do i=1,nvar
        call pdt_fix_cvar_string(cVars(i))
        n = len(trim(cVars(i)))
        cfields(i)(1:n) = trim(cVars(i))
      enddo

      return
      end

C=====================================================================72
      subroutine pdt_list_5_fields()
      use tecdat
      print *
      write (6,999) (trim(cVars(iv)),iv=1,5)
999     format(5a15)
      iy = 1
      iz = 1
      do ix = 1, nx
      write (6,998) (data(ix,iy,iz,iv),iv=1,5)
998   format(1p5e15.6)
      enddo

      return
      end

c=====================================================================72
      subroutine get_geometry(nstrokes,stroke)
      implicit none
      integer istr, nstr, nstrokes, istroke, j
      real x, y, stroke(4,10000)
      integer igetline, igetstrings, ngroups, igroup, npoints, ipoint
      integer i1(8000), i2(8000), nlen
      character*256 cline

999   format("nstr=",i3,"   cline: ",a)
c      write (6,*) "Getting a GEOMETRY"
      nlen = igetline(cline)
      nstr = igetstrings(cline, i1, i2)
c      write (6,999) nstr, trim(cline)
      read (cline(i1(1):i2(1)),*) ngroups
c      write (6,*) "ngroups = ", ngroups
      do igroup = 1, ngroups
        nlen = igetline(cline)
        nstr = igetstrings(cline, i1, i2)
        read (cline(i1(1):i2(1)),*) npoints
c        write (6,999) nstr, trim(cline)
c        write (6,*) "npoints = ", npoints
        do ipoint = 1, npoints
          nlen = igetline(cline)
          nstr = igetstrings(cline, i1, i2)
          read (cline(i1(1):i2(1)),*) x
          read (cline(i1(2):i2(2)),*) y
          nstrokes = nstrokes + 1
          if(ipoint .lt. npoints) then
            stroke(1,nstrokes) = x
            stroke(2,nstrokes) = y
          endif
          if(ipoint .gt. 1) then
            stroke(3,nstrokes-1) = x - stroke(1,nstrokes-1)
            stroke(4,nstrokes-1) = y - stroke(2,nstrokes-1)
          endif
        enddo
        nstrokes = nstrokes - 1
      enddo

c      do istroke = 1, nstrokes
c        write (6,998) (stroke(j,istroke), j=1,4)
c998     format(4f10.5)
c      enddo

      return
      end
c=====================================================================72
      subroutine getdata(nx,ny,nz,nvar,data, istr, nstr, cline)
      implicit none
      integer nx,ny,nz,nvar, istr, nstr
      real data(nx,ny,nz,nvar)
      integer igetline, igetstrings
      integer i1(8000), i2(8000), nlen, ic, ix,iy,iz,iv
      character*256 cline

c      istr = 0
c      nstr = 0
      nstr = igetstrings(cline, i1, i2)
      do iv = 1,nvar
      do iz = 1,nz
      do iy = 1,ny
      do ix = 1,nx
        if(istr .gt. nstr) then
          ic = 0
10        nlen = igetline(cline)
          ic = ic + 1
          if(ic .gt. 100) then
             write (6,*) "100 empty reads: may have run out of data"
             stop
          endif
          nstr = igetstrings(cline, i1, i2)
c          write (6,998) ix,iy,iz,iv,nstr,trim(cline)
c998       format("iXYZV=", 4i4, "     nstr=",i3,"     cline=",a)
          if(nstr .eq. 0) go to 10
          istr = 1
        endif
        read (cline(i1(istr):i2(istr)),*,err=99,end=98)
     1        data(ix,iy,iz,iv)
        istr = istr + 1
      enddo
      enddo
      enddo
      enddo
c      write (6,*) "got all of the data"
      return

98    write (6,*) "Eit end of file  reading data:"
      write (6,997) nx,ny,nz,nvar
      write (6,998) ix,iy,iz,iv,nstr,trim(cline)
      write (6,*)"nstr,istr=", nstr, istr
      stop

99    write (6,*) "Error reading data:"
      write (6,997) nx,ny,nz,nvar
997   format("iXYZV=", 4i4)
      write (6,998) ix,iy,iz,iv,nstr,trim(cline)
998   format("iXYZV=", 4i4, "     nstr=",i3,"     cline=",a)
      write (6,*)"nstr,istr=", nstr, istr
      write (6,*)"i1,i2=", i1(istr),i2(istr)
      write (6,*)"cline(i1,i2)=", cline(i1(istr):i2(istr))
      stop
      end
c=====================================================================72
      subroutine getd(cline,i0,nstr,istr1,istr2,izoned,nd,nreads)
      implicit none
      character*256 cline
      integer i0,nstr,istr1(8000),istr2(8000),izoned(3),nd,nreads
      integer i1, i2, islast
      nd = 0
10    nd = nd + 1
      i1 = istr1(i0+nd)
      i2 = istr2(i0+nd)
      islast = 0
      if(cline(i1:i1) .eq. '(') i1 = i1 + 1
      if(cline(i2:i2) .eq. ')') islast = 1
      if(cline(i2:i2) .eq. ')') i2 = i2 - 1
      read (cline(i1:i2),*) izoned(nd)

      if(nd    .ge.    3) islast = 1
      if(i0+nd .ge. nstr) islast = 1
      if(islast .ne. 1) go to 10

      nreads = nd
      return
      end
c=====================================================================72
      subroutine getc(cstr, czone, nreads)
      character*(*) cstr, czone
      call strcpy(czone, cstr)
      nreads = 1
      return
      end
c=====================================================================72
      subroutine geti(cstr, n, nreads)
      character*(*) cstr
      integer n, nreads
      read (cstr,*) n
      nreads = 1
      return
      end
c=====================================================================72
      subroutine strcpy(s1,s2)
      character*(*) s1, s2
      integer n1, n2, i
      n1 = len(s1)
      n2 = len(s2)
      do i = 1, min(n1,n2)
        s1(i:i) = s2(i:i)
      enddo
      do i = min(n1,n2)+1, n1
        s1(i:i) = ' '
      enddo
      return
      end
c=====================================================================72
      integer function igettype(str)
      character*(*) str
      igettype = 0
      if(str .eq. "TITLE"    ) igettype = 1
      if(str .eq. "VARIABLES") igettype = 2
      if(str .eq. "ZONE"     ) igettype = 3
      if(str .eq. "GEOMETRY" ) igettype = 4
      return
      end
c=====================================================================72
      integer function igetstrings(cline, istr1, istr2)
      character*256  cline
      integer istr1(8000), istr2(8000)
      integer n, i, nstr, instr, inquote
      logical bwhite
      n = len(trim(cline))
c      write (6,*) "len(trim(cline))=", n
      nstr = 0
      instr1qt2 = 0
      do i = 1, n
        bwhite = .false.
        if(cline(i:i).eq.'  ') bwhite = .true.
        if(cline(i:i) .eq. ',') bwhite = .true.
        if(cline(i:i) .eq. '=') bwhite = .true.
        if(instr1qt2  .eq. 0) then
          if(.not. bwhite) then
            nstr = nstr+1
            istr2(nstr) = -1
            if(cline(i:i).eq.'"') then
              istr1(nstr) = i+1
              instr1qt2 = 2
            else
              istr1(nstr) = i
              instr1qt2 = 1
            endif
          endif
        else if(instr1qt2 .eq. 2) then
          if(cline(i:i).eq.'"') then
            istr2(nstr) = i-1
            instr1qt2 = 0
          endif
        else  ! instr1qt2 has to be 1 (= in a string but not in a quoate) 
          if(bwhite) then
            istr2(nstr) = i-1
            instr1qt2 = 0
          endif
          if(cline(i:i).eq.'"') then
            istr2(nstr) = i-1
            instr1qt2 = 0
            istr1(nstr) = i+1
            instr1qt2 = 2
          endif
        endif
      enddo
      if(nstr .gt. 0) then
        if(istr2(nstr) .lt. 0) istr2(nstr) = n
      endif

      igetstrings = nstr
      return
      end

c=====================================================================72
      integer function igetline(cline)
      character*256 cline

       read (1,999,end=20) cline
999   format(256a) 
      igetline = len(trim(cline))
      return

20    igetline = -1
      return
      end
c=====================================================================72
c=====================================================================72
c           PARSE
c=====================================================================72
      integer function iget_cval_from_key_file(cFile, cKey, cVal)
      character*(*) cFile, cKey, cVal
      character*64 ckey1
      character*256 cline

      n = len(cKey)
      nfill = 0
      open(unit=17, file=cFile, form="formatted", err=900)
10    continue
      do i = 1, 256
        cline(i:i) = ' '
      enddo
      read (17,990,end=1000) cline
990   format(a256)
      if(cline(1:n) .eq. ckey) then
        read (cline,*) ckey1, cVal
        nfill = 1
      endif
      go to 10

900   continue
      write (6,*) 'error opening file:', trim(cFile)
      iget_cval_from_key_file = -1
      return

1000  continue
      close(17)
      iget_cval_from_key_file = nfill

      return
      end


      integer function iget_ival_from_key_file(cFile, cKey, iVal)
      character*(*) cFile, cKey
      nfill = iget_fval_from_key_file(cFile, cKey, fVal)
      if(nfill .gt. 0) iVal = int(fVal + 0.5)
      iget_ival_from_key_file = nfill
      return
      end
      

      integer function iget_fval_from_key_file(cFile, cKey, fVal)
      character*(*) cFile, cKey
      character*64 ckey1
      character*256 cline

      n = len(cKey)
      nfill = 0
      open(unit=17, file=cFile, form="formatted", err=900)
10    continue
      do i = 1, 256
        cline(i:i) = ' '
      enddo
      read (17,990,end=1000) cline
990   format(a256)
      if(cline(1:n) .eq. ckey) then
        read (cline,*) ckey1, fVal
        nfill = 1
      endif
      go to 10

900   continue
      write (6,*) 'error opening file:', trim(cFile)
      iget_fval_from_key_file = -1
      return

1000  continue
      close(17)
      iget_fval_from_key_file = nfill
      return
      end
