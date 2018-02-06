C=====================================================================72
      subroutine strfnout(n1,n2,ms,ns,k,nd,idisp,dvari,cOut,cChan,
     1                   ixdisp, iydisp, delx, time)
      integer idisp(n1:n2,ms,ns)
      character*256 cOut, cOutScl
      character*1   cChan
      parameter (MAXDISPSETS=10)
      parameter (MAXNDISPINSET=100)
      integer ixdisp(MAXNDISPINSET,MAXDISPSETS)
      integer iydisp(MAXNDISPINSET,MAXDISPSETS)
      real*8 sums(8), avs(8), var, dsum, var_to_the_mom

      ! Find limits of bins for output
      write (6,*) 'ms,ns,k =', ms,ns,k
      write (6,*) 'n1,n2,nd=', n1, n2,nd
      m1 = n2
      m2 = n1
      do j=1,nd
      do i=n1,n2
        if(idisp(i,j,k) .gt. 0) then
          m1 = min(i, m1)
          m2 = max(i, m2)
        endif
      enddo
      enddo
      write (6,*) 'm1,m2 = ', m1, m2

      ! Write output
      i = len(trim(cOut))
      cOut(i:i) = cChan
      write (6,*) "strfn3vout: ", cOut(1:80)
      open(unit=11, file=cOut, form="formatted")
      write (11,997) "#time     ", time
      write (11,997) "#delx     ", delx
997   format(a10, f12.6)
      write (11,998) "#         "
      write (11,998) "#xdisp    ", (ixdisp(j,k), j=1,nd)
      write (11,998) "#ydisp    ", (iydisp(j,k), j=1,nd)
      write (11,998) "#         "
998   format(a10, 20i12)
      do i=m1,m2
        var = (0.5 + float(i)) / dvari
        write (11,999) var,(idisp(i,j,k), j=1,nd)
999     format(f10.6, 20i12)
      enddo
      close(11)


      ! Calculate & output structure functions

      do i = 1, 256
        cOutScl(i:i) = ' '
      enddo
      ilen = len(trim(cOut))
      cOutscl(1:ilen+4) = trim(cOut)//"-scl"
      open(unit=11,file=cOutScl,form="formatted")
      write (11,997) "#time     ", time
      write (11,997) "#delx     ", delx
      write (11,994) "#           "
      write (11,994) "#Moment     ", (mom, mom=1,8)
      write (11,994) "#           "
994   format(a12, 8i12)
      do j=1,nd
        dsum = 0.0d0
        do mom=1,8
          sums(mom) = 0.0d0
        enddo
        do i=m1,m2
          dsum = dsum + dble(idisp(i,j,k))
          var = (0.5 + float(i)) / dvari
          var = dabs(var)
          var_to_the_mom = 1.0d0
          do mom=1,8
            var_to_the_mom = var_to_the_mom * var
            sums(mom) = sums(mom) + var_to_the_mom*dble(idisp(i,j,k))
          enddo
        enddo
        do mom=1,8
          avs(mom) = sums(mom) / dsum
        enddo
        delta = delx * sqrt(float(ixdisp(j,k)**2 + iydisp(j,k)**2))
        write (11,995) delta, (avs(mom),mom=1,8)
995     format(1p9e12.4)
      enddo
      close(11)

      if(n1 .ge. 0) return


      ! If signed scalar (n1<0): odd moments can cancel.
      !  This leads to a different set of structure functions.
      cOutscl(1:ilen+5) = trim(cOut)//"-scls"
      open(unit=11,file=cOutScl,form="formatted")
      write (11,997) "#time     ", time
      write (11,997) "#delx     ", delx
      write (11,994) "#           "
      write (11,994) "#Moment     ", (mom, mom=1,8)
      write (11,994) "#           "
      do j=1,nd
        dsum = 0.0d0
        do mom=1,8
          sums(mom) = 0.0d0
        enddo
        do i=m1,m2
          dsum = dsum + dble(idisp(i,j,k))
          var = (0.5 + float(i)) / dvari
          var_to_the_mom = 1.0d0
          do mom=1,8
            var_to_the_mom = var_to_the_mom * var
            sums(mom) = sums(mom) + var_to_the_mom*dble(idisp(i,j,k))
          enddo
        enddo
        do mom=1,8
          avs(mom) = sums(mom) / dsum
        enddo
        delta = delx * sqrt(float(ixdisp(j,k)**2 + iydisp(j,k)**2))
        write (11,995) delta, (avs(mom),mom=1,8)
      enddo
      close(11)

      return
      end
C=====================================================================72
