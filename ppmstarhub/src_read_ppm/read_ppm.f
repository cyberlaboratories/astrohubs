
      implicit none
      include "read_ppm.h"
      integer i, iv,im,ir
      real flds, moms
      allocatable flds(:,:,:,:)
      allocatable moms(:,:,:,:)

      call get_inputs()
      call gen_geometry()

      allocate(flds(nxfb,nyfb,nzfb,nReadVars))
      call fill_fields(ix0f,iy0f,iz0f,nxf,nyf,nzf,flds)

      if(cOutFile(1:5) .eq. "bofs:") then
        call dump_bofs(flds)
      else
        call dump_fields(flds)
      endif

      stop
      end
      
C=====================================================================72
      subroutine dump_bofs(flds)
      implicit none
      include "read_ppm.h"
      real flds(nxfb,nyfb,nzfb,nReadVars)
      character*256 cfile
      integer n, i, iv, m
c      real xyzranges(6)

      write (6,*) " "
      write (6,901)
901   format(/, "Writing fields in real*4 Brick-Of-Float format.")

      write (6,*) "      Min       Max      Mesh"
c      call get_xyzranges(xyzranges)
      write (6,902) "X: ", xyzranges(1), xyzranges(4), nxfb
      write (6,902) "Y: ", xyzranges(2), xyzranges(5), nxfb
      write (6,902) "Z: ", xyzranges(3), xyzranges(6), nxfb
902   format(a3,2f10.5,i7)

      write (6,*) " "
      n = len(trim(cOutFile))
      do i = 1,nReadVars
        iv = iReadVars(i)
        m  = len(trim(cppmvars(iv)))
        call set_str_blank(cfile)
        cfile(1:n-5) = cOutFile(6:n)
        cfile(n-4:n+m+1) = "-" // trim(cppmvars(iv)) // ".bof"
        write (6,*) "Writing BOF file: ", trim(cfile)
        call write_real_array(nxfb*nyfb*nzfb, flds(1,1,1,i), cfile)
      enddo

      return
      end

C=====================================================================72
      subroutine write_real_array(n, a, cfile)
      implicit none
      integer n
      real a(n)
      character*(*) cfile

#ifndef ISGFORTRAN
      open(unit=12,file=cfile,form="binary")
#else
      open(unit=12,file=cfile,form="unformatted",access="stream")
#endif
      write (12) a
      close(12)

      return
      end

C=====================================================================72
c      subroutine get_xyzranges(xyzranges)
c      implicit none
c      include "read_ppm.h"
c      real xyzranges(6)
c      real xsize,ysize,zsize, dx, dy, dz, x0f,y0f,z0f
c
c      dx = bsizex / float(ncpbx)
c      dy = bsizey / float(ncpby)
c      dz = bsizez / float(ncpbz)
c      xsize =  bsizex * nfilex * nbrkx
c      ysize =  bsizey * nfiley * nbrky
c      zsize =  bsizez * nfilez * nbrkz
c      x0f = -0.5 * xsize
c      y0f = -0.5 * ysize
c      z0f = -0.5 * zsize
c      xyzranges(1) = x0f + dx * float(ix0f-1)
c      xyzranges(2) = y0f + dy * float(iy0f-1)
c      xyzranges(3) = z0f + dz * float(iz0f-1)
c      xyzranges(4) = xyzranges(1) + dx * float(nxf)
c      xyzranges(5) = xyzranges(2) + dy * float(nyf)
c      xyzranges(6) = xyzranges(3) + dz * float(nzf)
c
c      return
c      end
c
C=====================================================================72
      subroutine gen_geometry()
      implicit none
      include "read_ppm.h"

      dx = bsizex / float(ncpbx)
      dy = bsizey / float(ncpby)
      dz = bsizez / float(ncpbz)
      xsize =  bsizex * nfilex * nbrkx
      ysize =  bsizey * nfiley * nbrky
      zsize =  bsizez * nfilez * nbrkz
      x0f = dx * float(ix0f-1) - 0.5 * xsize
      y0f = dy * float(iy0f-1) - 0.5 * ysize
      z0f = dz * float(iz0f-1) - 0.5 * zsize
      xyzranges(1) = x0f
      xyzranges(2) = y0f
      xyzranges(3) = z0f
      xyzranges(4) = xyzranges(1) + dx * float(nxf)
      xyzranges(5) = xyzranges(2) + dy * float(nyf)
      xyzranges(6) = xyzranges(3) + dz * float(nzf)

      return
      end

C=====================================================================72
      subroutine dump_fields(flds)
      implicit none
      include "read_ppm.h"
      real flds(nxfb,nyfb,nzfb,nReadVars)

      integer NHEAD
      parameter (NHEAD = 16*1024)
      integer*2 abuf(NHEAD)
      integer ndump, istep, itx,ity,itz, ntx,nty,ntz, i, iv, nbytes
      real dtime

      call remove_nan(99.0, nxfb,nyfb,nzfb,nReadVars, flds)

      ndump = 0
      dtime = time / float(ncycle)
      istep = ncycle
      itx   = 1
      ity   = 1
      itz   = 1
      ntx   = 1
      nty   = 1
      ntz   = 1

      ! Set dump parameters, including a default list variable to dump & what they are called
      call set_dump_parameters(ntx,nty,ntz)

      do i=1,nReadVars
        iv = iReadVars(i)
        call set_dump_var(i,cppmvars(iv),cppmvars(iv),"real4",0.0,0.0)
      enddo


      ! Format meta-data (header info) into a buffer (abuf).
      ! to customize for your 3D data array format.

      nbytes = 2*NHEAD
      call format_adump_header(abuf,ndump,time,dtime,istep,xyzranges,
     1             nbytes,nReadVars,nxfb,nyfb,nzfb,
     1             itx,ity,itz,ntx,nty,ntz)


      ! Write header buffer to a file
#ifndef ISGFORTRAN
      open(unit=12,file=cOutFile,form="binary")
#else
      open(unit=12,file=cOutFile,form="unformatted",access="stream")
#endif
      write (12) abuf
      write (12) flds
      close(12)

      return
      end
      
C=====================================================================72
      real function var_face(n1,n2,face)
      integer n1, n2, i,j
      real face(100,100), sum

      sum = 0.0
      do j = 1,n2-1
      do i = 1,n1-1
        sum = sum + (face(i+1,j)-face(i,j))**2
        sum = sum + (face(i,j+1)-face(i,j))**2
      enddo
      enddo

      var_face = sqrt(sum / float(2*(n1-1)*(n2-1)))
  
      return
      end

      
C=====================================================================72
      real function diff_face(n1,n2,face1,face2)
      integer n1, n2, i,j
      real face1(100,100)
      real face2(100,100), sum

      sum = 0.0
      do j = 1,n2
      do i = 1,n1
        sum = sum + (face1(i,j)-face2(i,j))**2
      enddo
      enddo

      diff_face = sqrt(sum / float(n1*n2))
  
      return
      end

C=====================================================================72
      subroutine fill_fields(ix0,iy0,iz0,nx,ny,nz,flds)
      implicit none
      include "read_ppm.h"

      ! Arguments
      integer ix0,iy0,iz0,nx,ny,nz,ioff(3)
      real flds(nxfb,nyfb,nzfb,nReadVars)

      ! scratch
      integer i, iv, ifi, ifi_last, ib, ib1, limits(6), iloc, idel
      character*256 cfile
      real fld, samp(100)
      character*1 cflds
      allocatable cflds(:,:,:,:,:)
      allocatable fld(:,:,:)

      allocate (cflds(ncpbx,ncpby,ncpbz,nbpb,nvpc))
      allocate (fld(ncpbx,ncpby,ncpbz))

      ifi_last = -1
      do ifi  = 1, nfilex * nfiley * nfilez
      do ib   = 1, nbrkx  * nbrky  * nbrkz

        call map_brick(ib, ib1)

        call overlap(ifi,ib1,ix0,iy0,iz0,nx,ny,nz,cfile,limits,ioff)
        if(cFile(1:4) .ne. "none") then
          if(ifi .ne. ifi_last) then
            write (6,*) trim(cFile)
            ifi_last = ifi
          endif
          call read_a_brick(ifi, cfile, ib, cflds)
          do iv=1,nReadVars
c >>>
            i = iReadVars(iv)
            iloc = ilocppmvars(i)
            idel = idelppmvars(i)
            if(nresppmvars(iloc) .eq. 1) then
              call ppm_var_map(iloc,cflds, fld)
            else
              call ppm_var_map8(iloc,cflds(1,1,1,1,iloc), fld)
            endif
c            write (6,*) 'iv,fld(1,1,1) = ', iv,fld(1,1,1)
            call copy_fld_subreg(ioff,limits,idel,fld,nx,ny,nz,iv,flds)
          enddo
        endif

       enddo
       enddo

      call read_a_brick(-1, cfile, ib, cflds)  ! close the last file
      deallocate(cflds)
      deallocate(fld)

      return
      end

C=====================================================================72
      subroutine map_brick(ib, ib1)
      implicit none
      include "read_ppm.h"
      integer ib, ib1,  nnbx,nnby,nnbz,nnb, in0, ini
      integer nndbx,nndby,nndbz, inx, iny, inz, ib0
      integer inix, iniy, iniz, ibx, iby, ibz

      nnbx = nbrkx / nodsx
      nnby = nbrky / nodsy
      nnbz = nbrkz / nodsz
      nnb  = nnbx * nnby * nnbz   ! number of bricks in a node

      ib1 = ib
      if(nnb .eq. 1) return

      ! count from 0
      ib0 = ib - 1
      in0 = ib0 / nnb

      ! XYZ of the node
      inx = mod(in0, nodsx)
      iny = mod(in0/nodsx, nodsy )
      inz = in0/(nodsx*nodsy)

      ! XYZ inside the node
      ini  = ib0 - in0 * nnb
      inix = mod(ini            , nnbx)
      iniy = mod(ini/ nnbx      , nnby)
      iniz =     ini/(nnbx*nnby)

      ! XYZ of brick
      ibx = inix + inx*nnbx
      iby = iniy + iny*nnby
      ibz = iniz + inz*nnbz

      ib1 = 1 + ibx + nbrkx*(iby + nbrky*ibz)


      return
      end
C=====================================================================72
      subroutine copy_fld_subreg(ioff,limits,idel,fld,nx,ny,nz,iv,flds)
      implicit none
      include "read_ppm.h"
      integer ioff(3), limits(6), idel, nx,ny,nz,iv
      real fld(ncpbx,ncpby,ncpbz)
      real flds(nxfb,nyfb,nzfb,nReadVars)
      integer ix,iy,iz, ix1,iy1,iz1, ixd,iyd,izd
      real factor, sum, delta, radav

      factor = 1.0 / float(nblend**3)
      do iz=limits(5),limits(6), nblend
      do iy=limits(3),limits(4), nblend
      do ix=limits(1),limits(2), nblend
        sum = 0.0
        do iz1 = iz, iz+nblend-1
        do iy1 = iy, iy+nblend-1
        do ix1 = ix, ix+nblend-1
          sum = sum + fld(ix,iy,iz)
        enddo
        enddo
        enddo

        ixd = 1 + (ix+ioff(1)-1)/nblend
        iyd = 1 + (iy+ioff(2)-1)/nblend
        izd = 1 + (iz+ioff(3)-1)/nblend
c >>> if  use radial profile (rrho or rprs) goes here
        if(idel .le. 0) then
          flds(ixd,iyd,izd,iv) = factor * sum
        else
          delta = factor * sum
          if(idel .eq. 1) call gen_radav(ixd,iyd,izd,nrrho,rrho,radav)
          if(idel .eq. 2) call gen_radav(ixd,iyd,izd,nrprs,rprs,radav)
          flds(ixd,iyd,izd,iv) = radav * (1.0 + delta)
        endif
      enddo
      enddo
      enddo
      return
      end
C=====================================================================72
      subroutine gen_radav(ixd,iyd,izd,n,rval,radav)
      implicit none
      include "read_ppm.h"
      integer ixd,iyd,izd,n,      ir0,ir1,ir2, nmomav, ix,iy,iz
      real rval(maxrad,2), radav,   dr,drinv,fnmax,x,y,z,ri,f1,f2
      real fac, dx1,dy1,dz1,off, x0,y0,z0, asum

      dr = (rval(n,1) - rval(1,1)) / float(n-1)
      drinv = 1.0 / dr
      fnmax = float(n-2)

      nmomav = 4
      fac = 1.0 / float(nmomav)
      dx1 = dx * fac
      dy1 = dy * fac
      dz1 = dz * fac
      off = 0.5 * float(nmomav+1)

      x0  = x0f + dx * (float(ixd) - 0.5)
      y0  = y0f + dy * (float(iyd) - 0.5)
      z0  = z0f + dz * (float(izd) - 0.5)

      asum = 0.0
      do iz = 1, nmomav
      z = z0 + dz1 * (float(iz) - off) 
      do iy = 1, nmomav
      y = y0 + dy1 * (float(iy) - off) 
      do ix = 1, nmomav
      x = x0 + dx1 * (float(ix) - off) 

      ri  = min(fnmax, drinv*sqrt(x**2 + y**2 + z**2))
      ir0 = int(ri)
      ir1 = ir0 + 1
      ir2 = ir0 + 2
      f2  = ri - float(ir0)
      f1  = 1.0 - f2
      asum = asum +  f1*rval(ir1,2) + f2*rval(ir2,2)
      enddo
      enddo
      enddo
      radav = asum *fac*fac*fac

      return
      end

C=====================================================================72
      subroutine overlap(ifile,ib,ix0,iy0,iz0,nx,ny,nz,cfile,
     1                   limits,ioff)
      implicit none
      include "read_ppm.h"
      integer ifile,ib,ix0,iy0,iz0,nx,ny,nz,limits(6)
      character*256 cfile
      integer ifx,ify,ifz, ibx,iby,ibz

      integer mefirst, i, ilen, ioff(3), idebug
      common /dgbover/ mefirst
      data mefirst/1/
      idebug = 0

      ! ib = 1 + (ibx-1)*ndbx + (iby-1)*ndby + (ibz-1)*ndbz

      ifx = 1 + mod(ifile-1,nfilex)
      ibx = 1 + mod(ib   -1, nbrkx)
cccc      ibx = 1 + mod((ib-1)/ndbx, nbrkx)

      call overlap1d(ifx,ibx, nbrkx,ncpbx, ix0,nx, limits(1),ioff(1))

      ify = 1 + mod((ifile-1)/nfilex, nfiley)
      iby = 1 + mod((ib   -1)/ nbrkx,  nbrky)
cccc      iby = 1 + mod(ib-1, nbrky)
      call overlap1d(ify,iby, nbrky,ncpby, iy0,ny, limits(3),ioff(2))

      ifz = 1 + (ifile-1) / (nfilex*nfiley)
      ibz = 1 + (ib-1) / ndbz
      call overlap1d(ifz,ibz, nbrkz,ncpbz, iz0,nz, limits(5),ioff(3))

      if(limits(1).le.limits(2) .and. limits(3).le.limits(4) .and.
     1  limits(5).le.limits(6)) then
          ! there is an overlap
          cfile = cRootFile
          ilen = len(trim(cRootFile))
          cfile(ilen-2:ilen-2) = char(ichar('a')+ifx-1)
          cfile(ilen-1:ilen-1) = char(ichar('a')+ify-1)
          cfile(ilen  :ilen  ) = char(ichar('a')+ifz-1)

         if(idebug .gt. 999) then
           if(mefirst .eq. 1) then
             mefirst = 0
             write (6,*)
     1         "File          ibx iby ibz    XYZ limits  OFFSET"
           endif
           write (6,999) cfile(ilen-11:ilen),
     1                   ibx,iby,ibz,(limits(i),i=1,6), (ioff(i),i=1,3)
999        format(a12, "   "3i4, 3("   ", 2i3), "   ", 3i5)
         endif
      else
          ! there is no overlap
          cFile(1:4) = "none"
      endif


      return
      end

C=====================================================================72
      subroutine overlap1d(ifd,ibd, nbrk,ncpb, irange0,nrange, lim,ioff)
      implicit none
      integer ifd,ibd, nfiled,ncpbd, id0,nd, lim(2)
      integer nbrk,ncpb, nrange, irange0, irange1, ibrick0, ibrick1
      integer ibrk, ioff
      integer ioverlap0, ioverlap1

c      integer mefirst
c      common /dgbover/ mefirst
c      data mefirst/1/

      irange1 = irange0 + nrange-1

      ibrk    = ibd + nbrk*(ifd-1)
      ibrick0 = 1 + ncpb * (ibrk - 1)
      ibrick1 =     ncpb *  ibrk
      ioverlap0 = max(irange0, ibrick0)
      ioverlap1 = min(irange1, ibrick1)
      lim(1) = 1 + ioverlap0 - ibrick0
      lim(2) = 1 + ioverlap1 - ibrick0
      ioff = ibrick0 - irange0

c      if(mefirst .eq. 1) then
c        write (6,*) "ibrk    brick            range           ",
c     1            "  overlap          limits"
c        mefirst = 0
c      endif
c      write (6,888) ibrk, ibrick0,ibrick1, irange0,irange1,
c     1                   ioverlap0,ioverlap1, lim(1),lim(2)
c888   format(i5, 2i7,"    ",2i7,"    ",2i7, "    ", 2i7)

      return
      end

C=====================================================================72
      subroutine read_a_brick(ifi, cfile, ib, cflds)
      implicit none
      include "read_ppm.h"

      ! Arguments
      character*256 cFile
      integer ifi, ib
      character*1 cflds(ncpbx,ncpby,ncpbz,nbpb,nvpc)
c      real        fld(ncpbx,ncpby,ncpbz)

      ! scratch
      INTEGER SEEK_SET, SEEK_CUR, SEEK_END, ierr
      integer ibrick, nb_brick
      integer(8) i8a, i8b, nb_offset, ifi_last
      common /rabdat/ ifi_last
      data ifi_last/-1/

      if(ifi .lt. 0) then
        close(11)
        ifi_last = -1
        return
      endif

      SEEK_SET = 0                                     ! Other options: SEEK_CUR = 1 SEEK_END = 2
      nb_brick = ncpbx * ncpby * ncpbz * nbpb * nvpc   ! bytes per brick
      i8a = ib - 1
      i8b = nb_brick
      nb_offset =  i8a * i8b           ! Offset into file to read the reques ted brick

      if(ifi .ne .ifi_last) then
        if(ifi_last .gt. 0) close(11)
#ifndef ISGFORTRAN
        open(11,file=cfile,form='binary',status='old', err=200,
     1       access="sequential", position='rewind')
#else
        open(11,file=cfile,form='unformatted',status='old', err=200,
     1       access="stream", position='rewind')
#endif
        ifi_last = ifi
      endif
      CALL FSEEK(11, nb_offset, SEEK_SET, ierr)
      read (11,err=111) cflds

      return

111   write (6,*) 'hit end of cfile'
      close(11)
      stop

200   write (6,*) 'open failed'
      stop
      end

c      call myopen(cfile)
c      call myread(nb_offset, cflds, nb_brick)
c      call myclose()
c      write (6,*) "read: ", trim(cfile)

C=====================================================================72
c      real axy(100,100)
c      character*256 cVNMfile
c      integer i,iv, ix, iy, iz, ib, ibrick,ival,ic,nny
c      integer axy(100,100)
c      real s1, s2, s10
c      integer is1, is2

c        write (6,*) "br,ix,iy,iz,ib,iv,cflds"
c        do iv = 1,24
c        do ib = 1,2
c
c        ic = 0
c        do iz = 1,ncpbz,ncpbz/4
c        do iy = 1,ncpby,ncpby/4
c        do ix = 1,ncpbx,ncpbx/4
c          ival = iachar(cflds(ix,iy,iz,ib,iv))
c          if(ival .ne. 222 .and. ic .lt. 100) then
c            ic = ic + 1
c            write (6,999) ibrick,ix,iy,iz,ib,iv,ival
c          endif
c999       format(6i3, i5)
c       enddo
c       enddo
c       enddo
c       if(ic .gt. 0) write (6,*) '---------'
c
c       enddo
c       enddo


C=====================================================================72
      subroutine ppm_var_map(iv,cflds,fld)

      include "read_ppm.h"
      integer iv
      character*1 cflds(ncpbx*ncpby*ncpbz,nbpb,nvpc)
      real          fld(ncpbx*ncpby*ncpbz)
      real offset, scale, s1, s2, y, p0, p1inv
      integer ididmap, i1, i2, ixyz, ncxyz
      real y1,y2,y3,y4
      real*8  pp0, pp1inv, yy, yy1, yy2, e, c, t, a, d, g, xx
      real fvmin, fvmax, ffvimg
      real faclog10, fvscale, fvlgfc

c      write (6,*) "Field:  ", trim(cppmvars(iv))

      ncxyz   = ncpbx*ncpby*ncpbz
      ididmap = 0

      if(cppmmap(iv)(1:12) .eq. "ppm-posdef  ") then
        ! map for positive deffinite fields
        offset = pvar000(iv)
        scale  = 1.0 / pvar001inv(iv)
        do ixyz = 1, ncxyz
           s1 =       float(iachar(cflds(ixyz,1,iv)))
           s2 = 0.5 + float(iachar(cflds(ixyz,2,iv)))
           y = (s1 + s2/250.499)/250.499
           y1 = y
           if(offset .gt. 0) y = 2.0*y - 1.0
           y2 = y
           y = sqrt((y+1.0)/(1.0-y))
           y3 = y
           y = (y*y - 1.0)/(2.0*y)
           y4 = y
           fld(ixyz) = offset + scale*y
        enddo
        ididmap = 1
      endif

      if(cppmmap(iv)(1:12) .eq. "ppm-signed  ") then
        ! map for singed fields
        scale  = 1.0 / pvar001inv(iv)
        do ixyz = 1, ncxyz
           s1 =       float(iachar(cflds(ixyz,1,iv)))
           s2 = 0.5 + float(iachar(cflds(ixyz,2,iv)))
           y = (s1 + s2/250.499)/250.499
           y = sqrt(y/(1-y))
           y = (y*y - 1.0)/(2.0*y)  ! check
           fld(ixyz) = scale*y
        enddo
        ididmap = 1
      endif

      if(cppmmap(iv)(1:12) .eq. "sak-posdef  ") then
        ! map for positive deffinite fields for Sakurai's star model (4908)
        pp0    = pvar000(iv)
        pp1inv = pvar001inv(iv)
        if (pp0 .le. 0.d0) then
          do ixyz = 1, ncxyz
            i1 = iachar(cflds(ixyz,1,iv))
            i2 = iachar(cflds(ixyz,2,iv))
c           call invposmapstar(p0, p1inv, i1,i2, fld(ixyz))
            yy1  =         dble(i1 - 2)
            yy2  = 0.5d0 + dble(i2 - 2)
            yy  = yy1 + yy2 / 250.499d0
            e = yy / 250.499d0
            c  = (1.0d0 + e) / (1.0d0 - e)
            t  = sqrt(c)
            a  = 0.5d0 * (t*t - 1.0d0) / t
            fld(ixyz) = (pp0 + a / pp1inv)
          enddo
        else
          do ixyz = 1, ncxyz
            i1 = iachar(cflds(ixyz,1,iv))
            i2 = iachar(cflds(ixyz,2,iv))
c           call invposmapstar(p0, p1inv, i1,i2, fld(ixyz))
            yy1  =         dble(i1 - 2)
            yy2  = 0.5d0 + dble(i2 - 2)
            yy  = yy1 + yy2 / 250.499d0
            e = (2.0d0/250.499d0) * yy - 1.0d0
            c  = (1.0d0 + e) / (1.0d0 - e)
            t  = sqrt(c)
            a  = 0.5d0 * (t*t - 1.0d0) / t
            fld(ixyz) = (pp0 + a / pp1inv)
          enddo
        endif


        ididmap = 1
      endif

      if(cppmmap(iv)(1:12) .eq. "sak-signed  ") then
        ! map for positive deffinite fields for Sakurai's star model (4908)
        pp1inv = pvar001inv(iv)
        do ixyz = 1, ncxyz
           i1 = iachar(cflds(ixyz,1,iv))
           i2 = iachar(cflds(ixyz,2,iv))
c           call invsgnmapbb(p0, p1inv, i1,i2, fld(ixyz))
           yy1 =         dble(i1 - 2)
           yy2 = 0.5d0 + dble(i2 - 2)
           yy  = yy1 + yy2 / 250.499d0
           a = yy / 250.499d0
           c = a / (1.0d0 - a)
           d = sqrt(c)
           g = 0.5d0*(d*d - 1)/d
           fld(ixyz) = g / pp1inv
        enddo

        ididmap = 1
      endif

      if(cppmmap(iv)(1:12) .eq. "sak-logmap  ") then
        ! map for positive deffinite fields for Sakurai's star model (4908)
        fvmin = pvar000(iv)
        fvmax = pvar001inv(iv)
        faclog10 = 1. / log(10.)
        fvscale = faclog10 * log(fvmax/fvmin)
        fvlgfc = 250.499 / fvscale
        do ixyz = 1, ncxyz
           i1 = iachar(cflds(ixyz,1,iv))
           i2 = iachar(cflds(ixyz,2,iv))
           y1 =       float(i1 - 2)
           y2 = 0.5 + float(i2 - 2)
           ffvimg  = y1 + y2 / 250.499
           fld(ixyz) = exp((ffvimg/fvlgfc - fvscale) / faclog10)
        enddo

        ididmap = 1
      endif

      if(ididmap .ne. 1) then
        write (6,*) "This map is not supported:"
        write (6,*) "-->", cppmmap(iv)(1:12), "<---   iv=",iv
        stop
      endif

      return
      end

C=====================================================================72
      subroutine ppm_var_map8(iv,cflds,fld)

      include "read_ppm.h"
      integer iv
      character*1 cflds(8*ncpbx*ncpby*ncpbz,nbpb)
      real          fld(  ncpbx*ncpby*ncpbz)
      real offset, scale, s1, s2, y, p0, p1inv
      integer ididmap, i1, i2, ixyz, ncxyz
      real y1,y2,y3,y4
      real*8  pp0, pp1inv, yy, yy1, yy2, e, c, t, a, d, g, xx
      real fvmin, fvmax, ffvimg
      real faclog10, fvscale, fvlgfc
      real fld8
      allocatable fld8(:)
      allocate(fld8(8*ncpbx*ncpby*ncpbz))

c      write (6,*) "Field:  ", trim(cppmvars(iv))

      ncxyz   = ncpbx*ncpby*ncpbz
      ididmap = 0

      if(cppmmap(iv)(1:12) .eq. "sak-logmap  ") then
        ! map for positive deffinite fields for Sakurai's star model (4908)
        fvmin = pvar000(iv)
        fvmax = pvar001inv(iv)
        faclog10 = 1. / log(10.)
        fvscale = faclog10 * log(fvmax/fvmin)
        fvlgfc = 250.499 / fvscale
        do ixyz = 1, 8*ncxyz
           i1 = iachar(cflds(ixyz,1))
           i2 = iachar(cflds(ixyz,2))
           y1 =       float(i1 - 2)
           y2 = 0.5 + float(i2 - 2)
           ffvimg  = y1 + y2 / 250.499
           if(i1+i2 .gt. 4) then
             fld8(ixyz) = exp((ffvimg/fvlgfc - fvscale) / faclog10)
           else
             fld8(ixyz) = 0.0
           endif
        enddo

        ididmap = 1
      endif


      if(ididmap .ne. 1) then
        write (6,*) "This map is not supported:"
        write (6,*) "-->", cppmmap(iv)(1:12), "<---   iv=",iv
        stop
      endif

      call blend8(ncpbx,ncpby,ncpbz, fld8, fld)
      deallocate(fld8)

      return
      end

C=====================================================================72
      subroutine blend8(nx,ny,nz,f8,f)
      implicit none
      integer nx, ny, nz, ix, iy, iz, ix1,iy1,iz1,ix2,iy2,iz2
      real f8(2*nx,2*ny,2*nz), f(nx,ny,nz)

      do iz = 1, nz
        iz2 = 2*iz
        iz1 = iz2 - 1
      do iy = 1, ny
        iy2 = 2*iy
        iy1 = iy2 - 1
      do ix = 1, nx
        ix2 = 2*ix
        ix1 = ix2 - 1

        f(ix,iy,iz) = 0.125 * (f8(ix1,iy1,iz1) + f8(ix2,iy1,iz1) +
     1                         f8(ix1,iy2,iz1) + f8(ix2,iy2,iz1) +
     1                         f8(ix1,iy1,iz2) + f8(ix2,iy1,iz2) +
     1                         f8(ix1,iy2,iz2) + f8(ix2,iy2,iz2)  )

      enddo
      enddo
      enddo

      return
      end

C=====================================================================72
      subroutine inv_fv_map(fvmin, fvmax, fv, ffvimg, imgffv1, imgffv2)
      implicit none
      real fvmin, fvmax, fv, ffvimg
      real faclog10, fvscale, fvlgfc, y1, y2
      integer imgffv1, imgffv2, ffvimg1, ffvimg2

      y1 =       float(imgffv1 - 2)
      y2 = 0.5 + float(imgffv2 - 2)
      ffvimg  = y1 + y2 / 250.499

      faclog10 = 1. / log(10.)
      fvscale = faclog10 * log(fvmax/fvmin)
      fvlgfc = 250.499 / fvscale

      fv = exp((ffvimg/fvlgfc - fvscale) / faclog10)

      return
      end

C=====================================================================72
      subroutine invsgnmapbb (s0, s1inv, i1,i2, x)

         y1 =       float(i1 - 2)
         y2 = 0.5 + float(i2 - 2)
         y  = y1 + y2 / 250.499

         a = y / 250.499
         c = a / (1.0 - a)
         d = sqrt(c)
         g = 0.5*(d*d - 1)/d
         x = g / s1inv

      return
      end

C=====================================================================72
      subroutine invposmapstar(p0, p1inv, i1,i2, x)
      implicit none
      real    p0, p1inv, x
      integer i1,i2
      real*8  pp0, pp1inv, yy, p1, p2, e, c, t, a, xx

      pp0    = p0
      pp1inv = p1inv

      !  im1 = y
      !  p1 = im1
      !  p2 = (y - p1) * 250.499
      !  im2 = p2
      !  i1 = im1 + 2
      !  i2 = im2 + 2

      p1  =         dble(i1 - 2)
      p2  = 0.5d0 + dble(i2 - 2)
      yy  = p1 + p2 / 250.499d0

      if (p0 .le. 0.) e = yy / 250.499d0
      if (p0 .gt. 0.) e = (2.0d0/250.499d0) * yy - 1.0d0
      c  = (1.0d0 + e) / (1.0d0 - e)
      t  = sqrt(c)
      a  = 0.5d0 * (t*t - 1.0d0) / t
      xx = p0 + a / p1inv
      x  = xx

      return
      end

C=====================================================================72
      subroutine get_lines(nlines, clines)
      implicit none
      integer i, nlines
      character*256 clines(256)

      do nlines = 1, 256
      do i = 1, 256
        clines(nlines)(i:i) = ' '
      enddo
      enddo

      open(17,file="read_ppm.in",form='formatted',status='old',err=230)

      nlines = 1
100   read (17,999,end=200) clines(nlines)
999   format(a256)
      nlines = nlines + 1
      if(nlines .gt. 256) then
        write (6,*) 'Subroutine setoutput:'
        write (6,*) 'Number of line in input exceeds the limit'
        write (6,*) 'Check input or increase the dims of clines'
        stop
      endif
      go to 100
200   continue
      close(17)
      return

230   write (6,*) 'failed to open read_ppm.in'
      stop
      end

C=====================================================================72
      subroutine set_defaults()
      implicit none
      include "read_ppm.h"
      integer i
      do i = 1, 256
        cRootFile(i:i) = " "
        cOutFile(i:i)  = " "
        cReadVars(i:i) = " "
      enddo

      ! invalid values
      nzf  = -1
      ix0f = -1
      ncpbx = -1
      ncpby = -1
      ncpbz = -1
      nvpc  = 0
      nblend = 1
      nfmoms = 1
      nfffvars = 0
      nodsx = 1
      nodsy = 1
      nodsz = 1
      time = 1.0
      ncycle = 1
      nrrho = -1
      nrprs = -1

      return
      end

C=====================================================================72
      subroutine gen_iloc_idel()
      implicit none
      include "read_ppm.h"
      character*16 cVar
      integer i, idel

      nppmvars = nvpc
      do i = 1, nvpc 
        ilocppmvars(i) = i
        idelppmvars(i) = 0
        cvar = cppmvars(i)
        idel = 0
        if(trim(cvar).eq."dRho" .and. nrrho.gt.0) idel = 1
        if(trim(cvar).eq."dPrs" .and. nrprs.gt.0) idel = 2
        if(idel .gt. 0) then
          nppmvars = nppmvars + 1
          nresppmvars(nppmvars) = nresppmvars(i)
          ilocppmvars(nppmvars) = i
          idelppmvars(nppmvars) = idel
          cppmvars(nppmvars) = "                "
          cppmvars(nppmvars)(1:3) = cvar(2:4)
        endif
      enddo
c >>>
      return
      end
      
c      write (6,*) "nppmvars = ", nppmvars
c      do i = 1, nppmvars 
c         write (6,999) i,ilocppmvars(i),idelppmvars(i)
c999      format("i,iloc,idel = ", 3i5)
c      enddo


C=====================================================================72
      subroutine get_read_var_list()
      implicit none
      include "read_ppm.h"
      character*16 cVar
      integer i0, i, irv

      call gen_iloc_idel()

      nReadVars = 0
      i0 = 1

10    continue
      cVar = "                "
      read (cReadVars(i0:256),*, end=100) cVar
      irv = -1
      do i = 1, nppmvars    ! nvpc
        if(trim(cVar) .eq. trim(cppmvars(i))) irv = i
      enddo
      if(irv .lt. 0) then
        write (6,*) "Requested variable: ==>", trim(cvar), "<=="
        write (6,*) "is not available in:"
        do i = 1, nppmvars    ! nvpc
          write (6,*) i, "==>", trim(cppmvars(i)), "<=="
        enddo
        stop
      endif

      nReadVars = nReadVars + 1
      iReadVars(nReadVars) = irv
c >>>
      i0 = i0 + len(trim(cVar)) + 1
      if(nReadVars .lt. 100) go to 10

100   continue
      write (6,*) "cReadVars = ", trim(cReadVars)
      write (6,999) (iReadVars(i),i=1,nReadVars)
999   format(" Read Vars: ", 100i4)


      return
      end

C=====================================================================72
      subroutine get_inputs()
      implicit none
      include "read_ppm.h"
      integer i, j, nlines, ix1f, iy1f, iz1f, ispace
      character*256 clines(256), str, ckey
      character*16 cvar

      call get_lines(nlines, clines)
      call set_defaults()
      do i = 1, nlines
        str = clines(i)
        if(str(1:8).eq."readvars") read (str,*) ckey,cReadVars
        if(str(1:8).eq."ixyz0   ") read (str,*) ckey,ix0f,iy0f,iz0f
        if(str(1:8).eq."nxyz    ") read (str,*) ckey,nxf ,nyf ,nzf
        if(str(1:8).eq."nblend  ") read (str,*) ckey,nblend
        if(str(1:8).eq."nfmoms  ") read (str,*) ckey,nfmoms
        if(str(1:8).eq."file    ") read (str,*) ckey,cRootFile
        if(str(1:8).eq."outfile ") read (str,*) ckey,cOutFile
        if(str(1:8).eq."time    ") read (str,*) ckey,time
        if(str(1:8).eq."ncycle  ") read (str,*) ckey,ncycle
        if(str(1:8).eq."rho_rad ") call read_rad(str,maxrad,nrrho,rrho)
        if(str(1:8).eq."prs_rad ") call read_rad(str,maxrad,nrprs,rprs)

        if(str(1:8).eq."nbrkx   ") read (str,*) ckey,nbrkx
        if(str(1:8).eq."nbrky   ") read (str,*) ckey,nbrky
        if(str(1:8).eq."nbrkz   ") read (str,*) ckey,nbrkz
        if(str(1:8).eq."nodsx   ") read (str,*) ckey,nodsx
        if(str(1:8).eq."nodsy   ") read (str,*) ckey,nodsy
        if(str(1:8).eq."nodsz   ") read (str,*) ckey,nodsz
        if(str(1:8).eq."bsizex  ") read (str,*) ckey,bsizex
        if(str(1:8).eq."bsizey  ") read (str,*) ckey,bsizey
        if(str(1:8).eq."bsizez  ") read (str,*) ckey,bsizez
        if(str(1:8).eq."ndbx    ") read (str,*) ckey,ndbx
        if(str(1:8).eq."ndby    ") read (str,*) ckey,ndby
        if(str(1:8).eq."ndbz    ") read (str,*) ckey,ndbz
        if(str(1:8).eq."nfilex  ") read (str,*) ckey,nfilex
        if(str(1:8).eq."nfiley  ") read (str,*) ckey,nfiley
        if(str(1:8).eq."nfilez  ") read (str,*) ckey,nfilez
        if(str(1:8).eq."ncpbx   ") read (str,*) ckey,ncpbx
        if(str(1:8).eq."ncpby   ") read (str,*) ckey,ncpby
        if(str(1:8).eq."ncpbz   ") read (str,*) ckey,ncpbz
        if(str(1:8).eq."nbpb    ") read (str,*) ckey,nbpb
        if(str(1:8).eq."field   ") then
          nvpc = nvpc + 1
          read (str,*) ckey,cppmvars(nvpc),cppmmap(nvpc),
     1                 pvar000(nvpc),pvar001inv(nvpc)
          nresppmvars(nvpc) = 1
        endif
        if(str(1:8).eq."field8  ") then
          nvpc = nvpc + 1
          read (str,*) ckey,cppmvars(nvpc),cppmmap(nvpc),
     1                 pvar000(nvpc),pvar001inv(nvpc)
          nresppmvars(nvpc) = 2
          do ispace = 2, 8
            nvpc = nvpc + 1
            cppmvars(nvpc)(1:8) = "#       "
          enddo
        endif
      enddo
      call get_read_var_list()
      call validate_inputs()
      nxfb = nxf / nblend
      nyfb = nyf / nblend
      nzfb = nzf / nblend
      return
      end

C=====================================================================72
      subroutine read_rad(str, maxrad, nrad, rval)
      implicit none
      character*(*) str
      character*256 ckey, cfile, cline
      integer       maxrad, nrad
      real          rval(maxrad,2)

      read (str,*) ckey, cfile
      open(unit=11,file=cfile,form="formatted",err=1001)

      nrad = 0
100   read (11,999,end=200) cline
999   format(a256)
      if(cline(1:1) .ne. '#') then
        nrad = nrad + 1
        read (cline,*,err=1002) rval(nrad, 1), rval(nrad,2)
        go to 100
      endif
200   continue

      close(11)
      write (6,*) "read_rad: ", trim(cfile), nrad
      return

1001  write (6,*) "Problem in read_rad: can not open", trim(cfile)
      stop

1002  write (6,*) "Problem in read_rad: can not read line", trim(cfile)
      stop

      end

C=====================================================================72
      subroutine validate_inputs()
      implicit none
      include "read_ppm.h"
      integer i

      if(nzf.le.0               ) write (6,*) "need xyzoff"
      if(              ix0f.le.0) write (6,*) "need yxzsize"
      if(nzf.le.0 .or. ix0f.le.0) stop
      if(cOutFile(1:4) .eq. "    ") then
        write (6,*) 'outfile not set'
        stop
      endif
      if(mod(nxf,nblend).ne.0 .or. mod(nyf,nblend).ne.0 .or.
     1   mod(nyf,nblend).ne.0 .or. mod(ncpbx,nblend).ne.0 .or.
     1   mod(ncpby,nblend).ne.0 .or. mod(ncpbz,nblend).ne.0) then
       write (6,*) 'n?f and ncpb? must all be multiples of nblend'
       stop
      endif

      ! Debug
      write (6,*) " "
      write (6,*) "RootFile   = ", trim(cRootFile)
      write (6,*) "xyzoffset  = ", ix0f,iy0f,iz0f
      write (6,*) "xyzsize    = ", nxf,nyf,nzf
      write (6,*) "nfile[xyz] = ", nfilex,nfiley,nfilez
      write (6,*) "nbrk[xyz]  = ", nbrkx,nbrky,nbrkz
      write (6,*) "nods[xyz]  = ", nodsx,nodsy,nodsz
      write (6,*) "ncpb[xyz]  = ", ncpbx,ncpby,ncpbz
      write (6,*) "nbpb,nvpc  = ", nbpb,nvpc
      write (6,*) " "
      write (6,*) "  i   Name           Map             ",
     1                  "       pvar000     pvar001inv"
      do i = 1, nvpc
         if(cppmvars(i)(1:1) .ne. "#")
     1   write (6,901) i,cppmvars(i),cppmmap(i),pvar000(i),pvar001inv(i)
901      format(i4, "   ", 2a15, 1p2e15.6)
      enddo
      write (6,*) " "

      return
      end

C=====================================================================72
      subroutine set_ppm_scalings(cVNMfile)
      implicit none
      character*256 cVNMfile
      character*5 cstr5
      include "read_ppm.h"
      integer i

      write (6,*) trim(cVNMfile)
      ! Read in argsimg info
      open(11,file=cVNMfile,form='formatted',status='old',err=230)

      ! Debug
      read (11,801) cstr5, nvpc
801   format(a5,i5)

      write (6,*) "nvpc = ", nvpc
      do i = 1, nvpc
      read (11,802) cppmvars(i),cppmmap(i),pvar000(i),pvar001inv(i)
802   format(2a16,1p2e16.6)
      enddo
      close(11)

      write (6,*) " "
      write (6,*) "  i   Name           Map             ",
     1                  "       pvar000     pvar001inv"
      do i = 1, nvpc
         write (6,901) i,cppmvars(i),cppmmap(i),pvar000(i),pvar001inv(i)
901      format(i4, "   ", 2a15, 1p2e15.6)
      enddo
      write (6,*) " "

      return

230   write (6,*) 'failed to open ', trim(cVNMfile)
      stop

      end

c=====================================================================72
      subroutine remove_nan(val, nx,ny,nz,nv, flds)
      implicit none
      integer nx,ny,nz,nv
      real val, flds(nx,ny,nz,nv)
      integer ix,iy,iz,iv
      logical is_nan

      do iv = 1, nv
      do iz = 1, nz
      do iy = 1, ny
      do ix = 1, nx
        if(is_nan(flds(ix,iy,iz,iv))) flds(ix,iy,iz,iv) = val
      enddo
      enddo
      enddo
      enddo

      return
      end

c=====================================================================72
      logical function is_nan(v)
      real v
      logical is_not_real
      is_not_real = .true.
      if(v .gt. 0.0) is_not_real = .false.
      if(v .lt. 1.0) is_not_real = .false.
      is_nan = is_not_real
      return
      end

C=====================================================================72
      subroutine set_str_blank(c)
      implicit none
      character*(*) c
      integer i
      do i=1,len(c)
        c(i:i) = ' '
      enddo
      return
      end

C=====================================================================72

