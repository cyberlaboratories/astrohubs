
C=====================================================================72

      subroutine sofplot(d,nx,ny,nDim,xmin,xmax,ymin,ymax,coutfile,
     1                   cx,cy, vmin,vmax,ix1,ix2,iy1,iy2,
     1                   ncbarsteps, cbarString, srin)
      implicit none
      character*(*) coutfile, cx, cy, cbarstring
      character*64  cpltfile
      integer nx,ny,nDim, ix1,ix2,iy1,iy2
      real d(nx,ny),xmin,xmax,ymin,ymax, vmin, vmax, v1, v2,scale,ex
      integer ix,iy, i,j,n,m, imap, nx1, ny1, ixx, iyy, nbar
      real  bof, xtic, auto_tic, x1, x2, y1, y2, fmap, cbar(4,100)
      real srused, srin
      integer nstrokes, ncbarsteps
      real stroke(4,10000)

      allocatable bof(:,:)

      v1 = vmin
      v2 = vmax
      call autoscale(d,nx,ny,v1,v2,imap,scale)

      allocate (bof(0:1+ix2-ix1,0:1+iy2-iy1))
      do iyy=iy1,iy2
      iy = iyy
      if(iyy .lt. 1) iy = 1 - iyy
      do ixx=ix1,ix2
      ix = ixx
      if(ixx .lt. 1) ix = 1 - ixx
        if(imap .ne. 2) then
          bof(ixx+1-ix1,iyy+1-iy1) = max(v1,d(ix,iy))
        else
          bof(ixx+1-ix1,iyy+1-iy1) = fmap(d(ix,iy),scale,v1,v2)
        endif
      enddo
      enddo
      bof(0,0) = float(1+ix2-ix1)

      nx1 = ix2+1-ix1
      ny1 = iy2+1-iy1
      do ix=1,nx1
        bof(ix,0) = xmin+(xmax-xmin)*(float(ix-1)+0.5)/float(nx1)
      enddo
      do iy=1,ny1
        bof(0,iy) = ymin+(ymax-ymin)*(float(iy-1)+0.5)/float(ny1)
      enddo
      x1 = bof(  1,0)
      x2 = bof(nx1,0)
c      write (6,*) 'x1,x2=', x1,x2
c      write (6,*) 'ix1,ix2=', ix1,ix2
c      write (6,*) 'nx1=', nx1
      y1 = bof(0,  1)
      y2 = bof(0,ny1)

#ifndef ISGFORTRAN
      open(unit=21,file=coutfile,form="binary")
#else
      open(unit=21,file=coutfile,form="unformatted",access="stream")
#endif
      write (21) bof
      close(21)
      deallocate (bof)

      ! Generate gnuplot script
      cpltfile( 1:32) = "                                "
      cpltfile(33:64) = "                                "
      n = len(trim(coutfile))
      cpltfile(1:n) = trim(coutfile)
      cpltfile(n+1:n+6) = ".plt  "
      open(unit=23,file=cpltfile,form="formatted")
999   format(a, a, a, a, a, a, a)
998   format(a, " [", e15.6,":",e15.6,"]")
c      write (23,999) "#unset title"
      srused = (y2-y1)/(x2-x1)
      if(srin .gt. 0.0) srused = srin
      write (23,996) "set size ratio", srused
#ifndef ISGFORTRAN
996   format(a, f, f)
#else
996   format(a, f6.2, f6.2)
#endif
c       ((y2-y1)/(x2-x1) .gt. 2.0)
      if(srused .gt. 2.0) then
        xtic = auto_tic(x1, x2)
        write (23,996) "set xtic", xtic, xtic
      endif

      write (23,999) "set pm3d map"
      if(imap .ne. 2) write (23,998) "set cbrange", v1, v2
      if(imap .eq. 2) write (23,998) "set cbrange", -1.0, 1.0
      if(imap .eq. 1) write (23,999) "set logscale zcb"

      if(cbarString(1:1) .ne. " ") then
        call string_to_cbar(cbarString, cbar, nbar)
       else
        call set_cbar(imap,v1,v2,nbar,cbar)
      endif
      if(ncbarsteps.gt.1) call cbar_substeps(nbar,ncbarsteps,cbar)
      call write_cbar(23, nbar, ncbarsteps, cbar)

      write (23,999) "set xlabel '", trim(cx), "'"
      if(imap .eq. 2) then
        write (23,999) "set cbtic 1 1"
c        write (23,999) "set cbtics add ('  '  1.0)"
c        write (23,999) "set cbtics add ('  ' -1.0)"
        call add_fmap_tics(23,scale,v1,v2)
c        write (23,997) "set xlabel '",trim(cx), scale
c997     format(a, a,"   [htan map, scale=",1pe8.1,"]'")
      endif
      write (23,999) "set ylabel '", trim(cy), "'"
      write (23,998) "set xrange", x1, x2
      write (23,998) "set yrange", y1, y2

      call get_vec_geom(nstrokes,stroke)
      if(nstrokes .le. 0) then
        if(ndim .eq. 1) then
          write (23,999) "splot '", trim(coutfile), "' binary notitle"
        else
c >>>
          write (23,999) "set style arrow 1 head filled ",
     1                    "size screen 0.009,10,15 ",
     1                    " linecolor rgb 'black' lw 2"
          write (23,999) "splot '", trim(coutfile),
     1      "' binary notitle with pm3d,",
     1      " 'vecfield.vec' using 1:2:(1):3:4:(1)",
     1      " notitle with vectors arrowstyle 1"
c <<<
        endif
      else
        if(ix1 .lt. 1) call mirror_strokes(nstrokes,stroke,-1.0,+1.0)
        if(iy1 .lt. 1) call mirror_strokes(nstrokes,stroke,+1.0,-1.0)
        call myclip(nstrokes,stroke,x1,x2,y1,y2)
        open(unit=37,file="geometry.4vec",form="formatted")
        do i=1,nstrokes
          write (37,*) (stroke(j,i),j=1,4)
        enddo
        close(37)

        write (23,999) "set clip two"
        write (23,999) "#"
        write (23,999) "# ack 0.018,20,30"
        write (23,999) "#"
        write (23,999) "set style arrow 1 head filled ",
     1                "size screen 0.009,10,15 ",
     1                " linecolor rgb 'black' lw 2"
        write (23,999) "set style arrow 7 nohead",
     1                " linecolor rgb 'red' lw 2"

        if(nDim .eq. 1) then
          write (23,999) "splot '", trim(coutfile), "' binary notitle",
     1        " with pm3d, 'geometry.4vec' using 1:2:(1):3:4:(1)",
     1        " notitle with vectors arrowstyle 7"
        else
          write (23,999) "splot '", trim(coutfile), "' binary notitle",
     1      " with pm3d, 'geometry.4vec' using 1:2:(1):3:4:(1)",
     1      " notitle with vectors arrowstyle 7, ",
     1      "'vecfield.vec' using 1:2:(1):3:4:(1) notitle with vectors",
     1      " arrowstyle 1"
        endif

c splot 'cdi00-0000-VZ_AE1_C.RZsof' binary notitle with pm3d, "strokes.dat" using 1:2:(10):3:4:(10) notitle with vectors arrowstyle 7

      endif
      close(23)
      return
      end

C=====================================================================72
      subroutine cbar_substeps(n, ncbarsteps, cbar)
      implicit none
      integer iunit, n,m, i0, j,k, ncbarsteps
      real cbar(4,100), rebar(4,100), f0, f1

      m = 0
      do i0=1,n-1
        do j = 0,ncbarsteps-1
          f1 = float(j) / float(ncbarsteps)
          f0 = 1.0 - f1
          m = m + 1
          do k = 1, 4
            rebar(k,m) = f0*cbar(k,i0) + f1*cbar(k,i0+1)
          enddo
        enddo
      enddo
      m = m + 1
      do k = 1, 4
        rebar(k,m) = cbar(k,n)
      enddo

      do i0 = 1, m
        do k = 1, 4
          cbar(k,i0) = rebar(k,i0)
        enddo
      enddo
      n = m
      return
      end

C=====================================================================72
      subroutine write_cbar(iunit, n, ncbarsteps, cbar)
      implicit none
      integer iunit, n, i, j,k, ncbarsteps
      real cbar(4,n)    !  (index,r,g,b)
      real fmid

      write (iunit,999) "set pm3d explicit at b"
      write (iunit,999) "set palette defined ( \"
      if(ncbarsteps .ge. 1) then
        write (iunit,998) (cbar(j,1), j=1,4), ", \"
        do i = 1, n-1
          fmid = 0.5 * (cbar(1,i) + cbar(1,i+1))
          write (iunit,998) fmid,(cbar(j,i  ), j=2,4), ", \"
          write (iunit,998) fmid,(cbar(j,i+1), j=2,4), ", \"
        enddo
        write (iunit,998) (cbar(j,n), j=1,4), " )"
      else
        do i = 1, n
         if(i.lt.n) write (iunit,998) (cbar(j,i), j=1,4), ", \"
         if(i.ge.n) write (iunit,998) (cbar(j,i), j=1,4), " )"
        enddo
      endif
999   format(a, a, a, a, a, a, a)
998   format(4f10.5, a)
      
      return
      end

C=====================================================================72
      subroutine set_cbar_row(fred,fgrn,fblu, n, cbar)
      implicit none
      integer n
      real find,fred,fgrn,fblu, cbar(4,100)
      n = n + 1
      cbar(1,n) = float(n)
      cbar(2,n) = fred
      cbar(3,n) = fgrn
      cbar(4,n) = fblu
      return
      end

C=====================================================================72
      subroutine set_cbar(imap,v1,v2,n,cbar)
      implicit none
      integer imap,n
      real v1,v2,cbar(4,100)

      n = 0

      if(imap .eq. 0) then
        if(v1 .ge. 0) then
          call set_cbar_row(0.3,0.3,0.3, n,cbar)
          call set_cbar_row(0.0,1.0,0.0, n,cbar)
          call set_cbar_row(1.0,0.0,0.0, n,cbar)
          call set_cbar_row(1.0,1.0,0.0, n,cbar)
          call set_cbar_row(1.0,1.0,1.0, n,cbar)
        else
          call set_cbar_row(0.0,1.0,1.0, n,cbar)
          call set_cbar_row(0.0,0.0,1.0, n,cbar)
          call set_cbar_row(0.0,0.5,0.0, n,cbar)
          call set_cbar_row(1.0,0.0,0.0, n,cbar)
          call set_cbar_row(1.0,1.0,0.0, n,cbar)
        endif
      endif

      if(imap .eq. 1) then
          call set_cbar_row(0.3,0.3,0.3, n,cbar)
          call set_cbar_row(0.0,0.0,1.0, n,cbar)
          call set_cbar_row(0.0,0.8,0.0, n,cbar)
          call set_cbar_row(1.0,0.0,0.0, n,cbar)
          call set_cbar_row(0.0,1.0,1.0, n,cbar)
          call set_cbar_row(1.0,0.0,1.0, n,cbar)
          call set_cbar_row(1.0,1.0,0.0, n,cbar)
      endif

      if(imap .eq. 2) then
        write (6,*) "doing imap=2"
        call set_cbar_row(0.0,1.0,1.0, n,cbar)
        call set_cbar_row(0.0,0.0,1.0, n,cbar)
        call set_cbar_row(0.0,0.5,0.0, n,cbar)
        call set_cbar_row(1.0,0.0,0.0, n,cbar)
        call set_cbar_row(1.0,1.0,0.0, n,cbar)

c        call set_cbar_row(0.0,1.0,1.0, n,cbar)
c        call set_cbar_row(0.0,0.0,1.0, n,cbar)
c        call set_cbar_row(0.3,0.3,0.3, n,cbar)
c        call set_cbar_row(1.0,0.0,0.0, n,cbar)
c        call set_cbar_row(1.0,1.0,0.0, n,cbar)
      endif

      return
      end

C=====================================================================72
      subroutine string_to_cbar(string, cbar, n)
      implicit none
      character*(*) string
      real cbar(4,100)
      integer n, nlen, i
      nlen = len(trim(string))
      n = 0
      do i = 1, nlen
        if(string(i:i).eq.'0') call set_cbar_row(0.0,0.0,0.0, n,cbar)
        if(string(i:i).eq.'1') call set_cbar_row(0.1,0.1,0.1, n,cbar)
        if(string(i:i).eq.'2') call set_cbar_row(0.2,0.2,0.2, n,cbar)
        if(string(i:i).eq.'3') call set_cbar_row(0.3,0.3,0.3, n,cbar)
        if(string(i:i).eq.'4') call set_cbar_row(0.4,0.4,0.4, n,cbar)
        if(string(i:i).eq.'5') call set_cbar_row(0.5,0.5,0.5, n,cbar)
        if(string(i:i).eq.'6') call set_cbar_row(0.6,0.6,0.6, n,cbar)
        if(string(i:i).eq.'7') call set_cbar_row(0.7,0.7,0.7, n,cbar)
        if(string(i:i).eq.'8') call set_cbar_row(0.8,0.8,0.8, n,cbar)
        if(string(i:i).eq.'9') call set_cbar_row(0.9,0.9,0.9, n,cbar)
        if(string(i:i).eq.'w') call set_cbar_row(1.0,1.0,1.0, n,cbar)
        if(string(i:i).eq.'r') call set_cbar_row(1.0,0.0,0.0, n,cbar)
        if(string(i:i).eq.'g') call set_cbar_row(0.0,1.0,0.0, n,cbar)
        if(string(i:i).eq.'b') call set_cbar_row(0.0,0.0,1.0, n,cbar)
        if(string(i:i).eq.'a') call set_cbar_row(0.0,1.0,1.0, n,cbar)
        if(string(i:i).eq.'p') call set_cbar_row(1.0,0.0,1.0, n,cbar)
        if(string(i:i).eq.'y') call set_cbar_row(1.0,1.0,0.0, n,cbar)
      enddo
      return
      end

C=====================================================================72
      subroutine myclip(n,v4,x1,x2,y1,y2)
      implicit none
      integer n, i
      real v4(4,10000), x1,x2,y1,y2
      do i=1,n
        call clip_to_limits(v4(1,i),x1,x2,y1,y2)
      enddo
      return
      end
C=====================================================================72
      subroutine clip_to_limits(v4,x1e,x2e,y1e,y2e)
      implicit none
      integer n, i
      real x1e,x2e,y1e,y2e
      real v4(4), x1,x2,y1,y2, smin, smax, s
      real x, y, xnew0,ynew0,xnew1,ynew1
      x1 = x1e + 0.001*(x2e-x1e)
      x2 = x2e - 0.001*(x2e-x1e)
      y1 = y1e + 0.001*(y2e-y1e)
      y2 = y2e - 0.001*(y2e-y1e)
      smax = -1.0
      smin = +2.0
      do i=-1,10001
        s = 0.0001*float(i)
        x = v4(1) + s * v4(3)
        y = v4(2) + s * v4(4)
        if(x.gt.x1 .and. x.lt.x2 .and. y.gt.y1 .and. y.lt.y2) then
          smax = max(smax, s)
          smin = min(smin, s)
        endif
      enddo
      if(smin .le. 0.0 .and. smax .ge. 1.0) return   ! segment is in plot
      if(smin .gt. smax) return   ! segment enirely outside of plot
      xnew0 = v4(1)
      ynew0 = v4(2)
      xnew1 = v4(1) + v4(3)
      ynew1 = v4(2) + v4(4)
      if(smin .gt. 0.0) then
        xnew0 = v4(1) + smin * v4(3)
        ynew0 = v4(2) + smin * v4(4)
      endif
      if(smin .lt. 1.0) then
        xnew1 = v4(1) + smax * v4(3)
        ynew1 = v4(2) + smax * v4(4)
      endif
      v4(1) = xnew0
      v4(2) = ynew0
      v4(3) = xnew1 - xnew0
      v4(4) = ynew1 - ynew0

      return
      end

C=====================================================================72
      subroutine mirror_strokes(nstrokes,stroke,xfac,yfac)
      implicit none
      integer nstrokes, i
      real stroke(4,10000), xfac, yfac
      if(xfac .ne. 0) then
      do i=1,nstrokes
        if(stroke(1,i) .eq. 0.0 .and. stroke(3,i) .eq. 0.0) then
          stroke(2,i) = stroke(4,i)   ! remove if on symmetry axis
        endif
      enddo
      endif
      do i=1,nstrokes
        stroke(1,nstrokes+i) = xfac * stroke(1,i)
        stroke(2,nstrokes+i) = yfac * stroke(2,i)
        stroke(3,nstrokes+i) = xfac * stroke(3,i)
        stroke(4,nstrokes+i) = yfac * stroke(4,i)
      enddo
      nstrokes = 2 * nstrokes
      return
      end

C=====================================================================72
      subroutine add_fmap_tics(iunit,scale,v1,v2)
      implicit none
      integer iunit
      real scale,v1,v2,x0, x,y,fmap, ylast
      integer is, ix, ic
      character*7 cx

      ylast = -1.0
      ic = 0
      x = 10.0**(int(100.1+alog(scale)/alog(10.0)) - 101)
10    x = 10.0 * x
      y = fmap(x,scale,v1,v2)
c      write (6,*) "x,y = ", x, y
      if(y .le. 0.9 .and. y-ylast.ge. 0.1) then
         call add_tic(iunit, x, y)
         call add_tic(iunit,-x,-y)
         ylast = y
      endif
      ic = ic + 1
      if(y .lt. 0.99 .and. ic .lt. 100) go to 10

      y = fmap(v2,scale,v1,v2)
      call add_tic(iunit, v2, y)
      call add_tic(iunit,-v2,-y)

      return
      end
C=====================================================================72
      subroutine add_tic(iunit,x,y)
      implicit none
      integer iunit
      real x,y
      character*7 cx

      cx = "     "
      write (cx,998) x
      write (iunit,999) cx, y

998   format(1pe7.0)
#ifndef ISGFORTRAN
999   format("set cbtics add ('",a,"' ",f,")")
#else
999   format("set cbtics add ('",a,"' ",f6.2,")")
#endif
      return
      end

C=====================================================================72
      real function auto_tic(xmin, xmax)
      implicit none
      real xmin, xmax, del, del10, tic
      del = xmax - xmin
      if(abs(xmax+xmin) .lt. 0.1*del) del = xmax
      del10 = 10.0**(float(int(1000.0+alog(del)/alog(10.0))-1000))
c      write (6,*) "del = ", del
c      write (6,*) "del10 = ", del10
      tic = del10
      if(del/del10 .gt. 2.0) tic = 2.0*del10
      if(del/del10 .gt. 3.0) tic = 3.0*del10
      if(del/del10 .gt. 4.0) tic = 4.0*del10
      if(del/del10 .gt. 5.0) tic = 5.0*del10
      if(del/del10 .gt. 6.0) tic = 6.0*del10
      if(del/del10 .gt. 7.0) tic = 7.0*del10
      if(del/del10 .gt. 8.0) tic = 8.0*del10
      if(del/del10 .gt. 9.0) tic = 9.0*del10
c      write (6,*) "tic = ", tic
      auto_tic = tic
      return
      end
C=====================================================================72
      real function fmap(x,scale,v1,v2)
      implicit none
      real x,scale,v1,v2, y

      if(abs(x) .le. scale) then
        y = 0.1 * x / scale
      else
            y = (alog(abs(x))-alog(scale))/(alog(v2)-alog(scale))
            y = min(1.0, 0.1 + 0.9 * y)
            if(x .lt. 0.0) y = -y
      endif
      fmap = y

      return
      end

c      ! 0.54930614433405484569
c      ex = exp(0.5493061*x/scale)
c      fmap = (ex - 1.0/ex) / (ex + 1.0/ex)

C=====================================================================72
      subroutine autoscale(d,nx,ny,v1,v2,imap, scale)
      implicit none
      integer nx,ny, imap, nsmall, ix,iy, iround, nonzero, i10pc
      real d(nx,ny),v1, v2, round_to_nice_value,av20,rms,scale,avmax
      real*8 s2
      real eps, dbin, vmedian, val, scl, vabsmin, v1a, v10pc
      integer mbin, ibin, imedian, isum, imin
      parameter(dbin=0.1)
      parameter(mbin=380)  ! single precision limit
      integer idist(-mbin:mbin)

      if(v1 .lt. v2) then      ! limits have been externally set
        if(v1 .eq. 0.0) then
          imap = 0   ! linear
          return
        endif

        if(0.0.lt.v1  .and.  v1.lt.0.1*v2) then
          imap = 1   ! log
          return
        endif 

        if(v1.lt.0.0  .and.  abs(v1).lt.0.1*v2) then
          imap = 2   !  double log
          scale = abs(v1)
          return
        endif 

        imap = 0   ! linear
        return
      endif

      iround = 0
      if(v1 .gt. v2) then
         ! Find min & max values in d if v1 & v2 have not been set
         do ibin=-mbin,mbin
           idist(ibin) = 0
         enddo
         nonzero = 0
         v1 = d(1,1)
         v2 = d(1,1)
         do iy=1,ny
         do ix=1,nx
           val = d(ix,iy)
           v1 = min(v1,val)
           v2 = max(v2,val)
           if(val .ne. 0) then
             ibin = int(float(mbin)+alog(abs(val))/dbin) - mbin
             ibin = max(-mbin,min(mbin,ibin))
             idist(ibin) = idist(ibin) + 1
             nonzero = nonzero + 1
           endif
         enddo
         enddo
         if(v1 .eq. v2) then
           if(v1 .eq. 0) then
             v1 = -1.0
             v2 =  1.0
           else
             v1 = v1 - 0.1 * abs(v1)
             v2 = v2 + 0.1 * abs(v2)
           endif
         endif
         iround = 1
      endif
      imin = mbin + 1
      imedian = 0
      i10pc   = 0
      isum    = 0
      do ibin=-mbin,mbin
        isum = isum + idist(ibin)
        if(2*isum .lt. nonzero) imedian = ibin
        if(10*isum .lt. nonzero) i10pc  = ibin
        if(imin.gt.mbin .and. idist(ibin).gt.0) imin = ibin
      enddo
      v10pc   = exp(dbin*float(  i10pc))
      vmedian = exp(dbin*float(imedian))
      vabsmin = exp(dbin*float(imin))

c      ! gard agains roundoff
c      if(abs(v1) .lt. 1.0e-36) v1 = 0.0
c      if(abs(vabsmin) .lt. 1.0e-36) vabsmin = 0.0

      avmax = max(abs(v1),abs(v2))
      av20 = 0.05 * avmax
      eps  = 0.00001 * avmax
      nsmall = 0
      s2 = 0.0
      nonzero = 0
      do iy=1,ny
      do ix=1,nx
        if(abs(d(ix,iy)) .ge. eps) then
          if(abs(d(ix,iy)) .lt. av20) nsmall = nsmall + 1
          s2 = s2 + dble(d(ix,iy))**2
          nonzero = nonzero + 1
        endif
      enddo
      enddo
      rms = 0.0
      if(s2 .gt. 0.0) rms = sqrt(s2 / dble(nonzero))


      imap = 0                     ! safe default

      if(nonzero .gt. 2*nsmall .or. rms .eq. 0.0) then
        imap = 0                   ! linear map (dynamic range NOT large)
      else
        if(v1 .eq. 0.0) v1 = vabsmin
        if(v1 .gt. 0.0) imap = 1   ! log map  (dynamic range is large + values only)
        if(v1 .le. 0.0) imap = 2   ! htan map (dynamic range is large +&- values)
      endif
c      write (6,*) "avmax = ", avmax
c      write (6,*) "av20  = ", av20
c      write (6,*) "eps   = ", eps
c      write (6,*) "nonzero = ", nonzero
c      write (6,*) "nsmall  = ", nsmall
c      write (6,*) "rms     = ", rms
c      write (6,*) "median value = ", vmedian
c      write (6,*) "v1           = ", v1
c      write (6,*) "imap         = ", imap

      if(imap .eq. 2) then         ! double log scale
        scale = 10.0**(float(int(1000.0+alog(v10pc)/alog(10.0))-1000))
        v2     = max(abs(v1), abs(v2))
        v2    = 10.0**(float(int(1001.0+alog(v2   )/alog(10.0))-1000))
c        write (6,*) "v10pc   = ", v10pc
c        write (6,*) "scale   = ", scale
c        write (6,*) "v1      = ", v1
c        write (6,*) "v2      = ", v2

cccccccccccccccccccccccccccccccccccccccccccccccccc
c        scl = min(av20, rms)
c        scale = round_to_nice_value(vmedian,0.5)
cccc        if(v1 .lt. -eps) then
cccc          v1 = -1.0
cccc        else
cccc          v1 = 0.0
cccc        endif
cccc        v2 =  1.0
cccccccccccccccccccccccccccccccccccccccccccccccccc
        return
      endif

      if(iround .eq. 0) return   ! min,max values were externally set
      ! round auto-scale scale range to "nice" values

      if(imap .eq. 0) then
        if(v2.gt.0.0 .and. v1.lt.0.0) then
          v2 = round_to_nice_value(avmax,0.99999)
          v1 = -v2
          return
        endif
        if(v2 .gt. 0.0 .and. v1 .lt. 0.5*v1) then
          v2 = round_to_nice_value(avmax,0.99999)
          v1 = 0.0
          return
        endif
        v1 = round_to_nice_value(v1,0.0)
        v2 = round_to_nice_value(v2,0.99999)
      endif
      if(imap .eq. 1) then
        v1 = 10.0**(float(int(1000.0+alog(v1)/alog(10.0))-1000))
        v2 = 10.0**(float(int(1001.0+alog(v2)/alog(10.0))-1000))
      endif

      return
      end

C=====================================================================72
      real function round_to_nice_value(v,r)
      implicit none
      real v,r, av, sv, av10, avn
      if(v .eq. 0.0) then
        round_to_nice_value = 0.0
        return
      endif
      sv = 1.0
      if(v .lt. 0) sv = -1.0
      av = sv * v

      av10 = 10.0**float(int(1000.5 + alog(av)/alog(10.))-1000)
      avn = float(int(r + 10.0*av/av10)) / 10.0

      round_to_nice_value = sv * avn * av10

      return
      end
C=====================================================================72                             

