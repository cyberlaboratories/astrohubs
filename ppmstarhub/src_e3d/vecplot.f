
C=====================================================================72
      subroutine gen_vecfield_file(d,nx,ny,xmin,xmax,ymin,ymax,ix1,iy1)
      implicit none
      integer nx,ny, ix1,iy1, ix,iy, ic, i, j, nstrokes, ndel,ixoff,iix
      real d(nx,ny,2),xmin,xmax,ymin,ymax,scale, stroke(4,10000)
      real eps, dx, dy, x0, y0, x, y, signx
      real percentile_value   ! (d,n,p)
      real v2, vv, v80
      allocatable v2(:)

      x0 = xmin
      y0 = ymin
      if(ix1 .lt. 1) x0 = 0.0
      if(iy1 .lt. 1) y0 = 0.0

      dx = (xmax-x0)/float(nx)
      dy = (ymax-y0)/float(ny)
      ndel = max(2, int(max(nx,ny)/30))

      ! Calculate scales & sizes
      allocate (v2(nx*ny))
      ic = 0
      do iy = 1, ny
      do ix = 1, nx
        vv = d(ix,iy,1)**2 + d(ix,iy,2)**2
        if(vv .gt. 0.0) then
          ic = ic + 1
          v2(ic) = vv
        endif
      enddo
      enddo
      v80 = sqrt(percentile_value(v2,ic,80.0))
      deallocate (v2)

      eps = v80 / 1000000.0
      scale = 2.0 * float(ndel) * min(dx,dy) / v80

      ! Generate & mirror vectors as needed
      nstrokes = 0
      ixoff = 1
      do iy=ndel/2,ny,ndel
      ixoff = ixoff + max(1, ndel/5) 
      if(ixoff .gt. ndel) ixoff = 1

      do iix=ix1+ixoff,nx,ndel
        ix = iix
        signx = scale
        if(iix .lt. 1) then
          ix = 1-iix
          signx  = -1.0 * scale
        endif
        if(sqrt(d(ix,iy,1)**2+d(ix,iy,2)**2) .ge. eps) then
          nstrokes = nstrokes + 1
          stroke(1,nstrokes) = x0 + dx*(0.5+float(iix-1))
          stroke(2,nstrokes) = y0 + dy*(0.5+float(iy -1))
          stroke(3,nstrokes) = signx * d(ix,iy,1)
          stroke(4,nstrokes) = scale * d(ix,iy,2)
        endif
      enddo
      enddo
c      if(ix1 .lt. 1) call mirror_strokes(nstrokes,stroke,-1.0,+1.0)
c      if(iy1 .lt. 1) call mirror_strokes(nstrokes,stroke,+1.0,-1.0)

      ! Write vectors to file
      open(unit=23,file="vecfield.vec",form="formatted")
      do i = 1, nstrokes
        write (23,999) (stroke(j,i), j=1,4)
999     format(1p4e15.6)
      enddo
      close(23)

      return
      end

C=====================================================================72
      ! Returns the value of the p-th percentile in d(1:n)
      real function percentile_value(d,n,p)
      implicit none
      real d(n), p, dbin, val, aval, pval
      integer n, i, mbin, ibin, icp, ic, ibinp
      parameter(mbin=10000)
      integer idist(-mbin:mbin)

      dbin = 0.01
      do ibin=-mbin,mbin
        idist(ibin) = 0
      enddo
      ic = 0

      do i = 1, n
        val = d(i)
        if(val .eq. 0) then
          idist(0) = idist(0) + 1
        else
          aval = abs(val)
          ibin = int(float(mbin)+alog(abs(aval))/dbin) - mbin
          ibin = max(-mbin/2,min(mbin/2,ibin))
          if(val .gt. 0.0) ibin = ibin + mbin/2
          if(val .lt. 0.0) ibin = ibin - mbin/2
          idist(ibin) = idist(ibin) + 1
        endif
        ic = ic + 1
      enddo

      icp = int(0.01 * p * float(ic))
      ibinp = -mbin
      ic = 0
      do ibin=-mbin,mbin
       if(ic .lt. icp) ibinp = ibin
       ic = ic + idist(ibin)
      enddo

      if(ibinp .eq. 0) pval = 0.0
      if(ibinp .gt. 0) pval = exp( dbin*float(ibinp - mbin/2))
      if(ibinp .lt. 0) pval = exp(-dbin*float(ibinp + mbin/2))

      percentile_value = pval
      return
      end

C=====================================================================72
      subroutine vecplot(d,nx,ny,xmin,xmax,ymin,ymax,coutfile,
     1                   cx,cy, vmin,vmax,ix1,ix2,iy1,iy2)
      implicit none
      character*(*) coutfile, cx, cy
      character*64  cpltfile
      integer nx,ny, ix1,ix2,iy1,iy2
      real d(nx,ny,2),xmin,xmax,ymin,ymax, vmin, vmax, v1, v2,scale,ex
      integer ix,iy, i,j,n,m, imap, nx1, ny1, ixx, iyy
      real  xtic, auto_tic, x1, x2, y1, y2, fmap, vx1,vx2,vy1,vy2
      integer nstrokes
      real stroke(4,10000)

      vx1 = vmin
      vx2 = vmax
      call autoscale(d(1,1,1),nx,ny,vx1,vx2,imap,scale)
      vy1 = vmin
      vy2 = vmax
      call autoscale(d(2,1,1),nx,ny,vy1,vy2,imap,scale)
      v1 = 0.0
      v2 = max(vx2,vy2)


      ! Generate gnuplot script
      cpltfile( 1:32) = "                                "
      cpltfile(33:64) = "                                "
      n = len(trim(coutfile))
      cpltfile(1:n) = trim(coutfile)
      cpltfile(n+1:n+6) = ".plt  "
      open(unit=23,file=cpltfile,form="formatted")
999   format(a, a, a, a, a, a)
998   format(a, " [", e15.6,":",e15.6,"]")
      write (23,996) "set size ratio", (y2-y1)/(x2-x1)
#ifndef ISGFORTRAN
996   format(a, f, f)
#else
996   format(a, f6.2, f6.2)
#endif
      if((y2-y1)/(x2-x1) .gt. 2.0) then
        xtic = auto_tic(x1, x2)
        write (23,996) "set xtic", xtic, xtic
      endif

      write (23,999) "set xlabel '", trim(cx), "'"
      write (23,999) "set ylabel '", trim(cy), "'"
      write (23,998) "set xrange", x1, x2
      write (23,998) "set yrange", y1, y2

      call get_vec_geom(nstrokes,stroke)
      if(nstrokes .le. 0) then
        write (23,999) "splot 'vecfield.vec' using 1:2:3:4"
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
        write (23,999) "set style arrow 7 nohead ls 1"
#ifndef ISGFORTRAN
        write (23,999) "splot 'vecfield.vec using' 1:2:3:4"
     1                 ", 'geometry.4vec' using 1:2:3:4",
     1                 " notitle with vectors arrowstyle 7"
#else
        write (23,999) "splot 'vecfield.vec using' 1:2:3:4",
     1                 ", 'geometry.4vec' using 1:2:3:4",
     1                 " notitle with vectors arrowstyle 7"
#endif

      endif
      close(23)
      return
      end

C=====================================================================72

