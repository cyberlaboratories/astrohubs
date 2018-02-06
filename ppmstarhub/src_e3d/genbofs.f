
      parameter(n=128)
      real x(n,n,n), y(n,n,n), z(n,n,n)

      dx = 10.0 / float(n)
      do k = 1, n
      do j = 1, n
      do i = 1, n
        x(i,j,k) = dx * float(i) - 0.5 * dx - 5.0
        y(i,j,k) = dx * float(j) - 0.5 * dx - 5.0
        z(i,j,k) = dx * float(k) - 0.5 * dx - 5.0
      enddo
      enddo
      enddo

      open(unit=11,file="xbof",form="binary")
      write (11) x
      close(11)

      open(unit=11,file="ybof",form="binary")
      write (11) y
      close(11)

      open(unit=11,file="zbof",form="binary")
      write (11) z
      close(11)

      stop
      end
