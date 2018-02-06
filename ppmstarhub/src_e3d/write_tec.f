c=====================================================================72
      subroutine write_tec(cfile,ct,cx,cy,cv, nx,ny,sof)
      implicit none
      character*(*) cfile,ct,cx,cy,cv
      integer nx, ny
      real sof(0:nx, 0:ny)
      integer i1, i2, i, j

      write (6,*) "at top of write_tec"

c      cv = "PSD     "
c      cx = "x       "
c      cy = "y       "
c      ct = "Is MINE "
c      cfile = "a.pdt   "

      write (6,*) trim(ct)
      write (6,*) trim(cx)
      write (6,*) trim(cy)
      write (6,*) trim(cv)
      write (6,*) trim(cfile)

      open(22,file=cfile,form="formatted")
900   format(a,a,a,a,a)
      write (22,900) " TITLE = ", '"', trim(ct), '"'
      write (22,900) "     VARIABLES=", '"', trim(cx), '"', ','
      write (22,900)                    '"', trim(cy), '"', ','
      write (22,900)                    '"', trim(cv), '"'

901   format("     ZONE  I=",i4,", J=",i4,", F=BLOCK,")
      write (22,901) nx, ny
      write (22,900) "           T=",'"', "2-D DATA", '"'

902   format(1p7e11.3)
      do j  = 1, ny
      do i1 = 1, nx, 7
        i2 = min(i1+6, nx)
        write (22,902) (sof(i,0),i=i1,i2)
      enddo
      enddo

      do j  = 1, ny
      do i1 = 1, nx, 7
        i2 = min(i1+6, nx)
        write (22,902) (sof(0,j),i=i1,i2)
      enddo
      enddo

      do j  = 1, ny
      do i1 = 1, nx, 7
        i2 = min(i1+6, nx)
        write (22,902) (sof(i,j),i=i1,i2)
      enddo
      enddo
      close(22)

      return
      end
c=====================================================================72

