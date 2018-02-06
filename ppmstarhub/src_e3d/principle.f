
c=====================================================================72
      subroutine princip(a,ye,ee, ierr)
      real*8 a(3,3), ye(3), ee(3,3)
      real*8 P,Q,R, aa, bb, pi, modulus, argument, ae(3,3), y
      real*8 fNorm, temp, x1, x2, x3, theta, ddot
      logical debug
      debug = .false.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    Solve for eigenvalues
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! Calculate P, Q, R invarients of eigenvalue equations
      P = a(1,1) + a(2,2) + a(3,3)
      Q = a(1,1)*a(2,2) - a(1,2)*a(2,1) +
     1    a(1,1)*a(3,3) - a(1,3)*a(3,1) +
     1    a(2,2)*a(3,3) - a(2,3)*a(3,2)
      R = a(1,1)*a(2,2)*a(3,3) - a(1,1)*a(2,3)*a(3,2) +
     1    a(1,2)*a(2,3)*a(3,1) - a(1,2)*a(2,1)*a(3,3) +
     1    a(1,3)*a(2,1)*a(3,2) - a(1,3)*a(2,2)*a(3,1)

      if(debug) write (6,*) 'PQR =', P, Q, R

      ! Solve eigenvalue equation:  det(Sd-yI) = -(y**3) + P*(y**2) - Q*y + R = 0
      ! => y**3 + p*(y**2) + q*y + r = 0
      ! where p = -P   q = Q   r = -R
      ! Let y = x - p/3 = x + P/3
      ! solve:  x**3 + a*x + b = 0
      ! where a = (3*q - p**2)/3  b = (2*p**3 - 9*p*q + 27*r)/27
      aa = (3.0d0 * Q - P*P) / 3.0d0
      bb = (9.0d0*P*Q  - 2.0d0*P*P*P - 27.0d0 * R) / 27.0d0


      if(.not. (aa .lt. 0.0d0)) then
c      write (6,*) "=============================================="
c        write (6,*) 'aa,bb =', aa, bb
c        do i=1,3
c          write (6,*) "a(i,j): ",(a(i,j),j=1,3)
c        enddo
c        write (6,*) "P = ", P
c        write (6,*) "Q = ", Q
c        write (6,*) "aa = (3.0d0 * Q - P*P) / 3.0d0"
c        write (6,*) "aa = ", aa
c        write (6,*) "aa > 0 in a symmetric matix"
c      write (6,*) "=============================================="
c        ierr = -1

        ierr = 0
        do i = 1, 3
          do j = 1, 3
            ee(i,j) = 0.0d+00
          enddo
          ye(i) = 0.0d+00
          ee(i,i) = 1.0d+00
        enddo
        return
      endif

      pi = 2.0d0 * acos(0.0d0)
      modulus = 2.0d0 * sqrt(-aa / 3.0d0)
      argument = 3.0d0 * bb / (aa * modulus)
      if(argument .ge. 1.0d0) then
        theta = 0.0d0
      else if(argument .le. -1.0d0) then
        theta = pi
      else
        theta = acos(argument) / 3.0d0
      endif

      if(debug) write (6,*) 'Modulus,theta =', modulus,theta
      if(debug) write (6,*) '3.0d0*bb/(aa*modulus) =', argument

       x1 = modulus * cos(theta)
       x2 = modulus * cos(theta + 2.0d0 * pi / 3.0d0)
       x3 = modulus * cos(theta + 4.0d0 * pi / 3.0d0)

      if(debug) write (6,*) 'x123 =', x1,x2,x3

       ye(1) = x1 + P / 3.0d0
       ye(2) = x2 + P / 3.0d0
       ye(3) = x3 + P / 3.0d0

      if(debug) write (6,*) 'ye123 =', (ye(i),i=1,3)

       ! Sort: largest to smallest
       do i = 1, 2
         imax = i
         do j = i+1, 3
           if(ye(j) .gt. ye(imax)) imax = j
         enddo
         if(imax .ne. i) then
           temp = ye(i)
           ye(i) = ye(imax)
           ye(imax) = temp
         endif
       enddo

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     Now solve for eigenvectors
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! eigenvalues 1 & 3 are farthest appart in value
      do ign = 1, 3, 2
        y = ye(ign)
        do j=1,3
        do i=1,3
          ae(i,j) = a(i,j)
          if(i .eq. j) ae(i,j) = a(i,j) - y
        enddo
        enddo

        if(debug) then
        write (6,*) 'ign=', ign
        do i=1,3
          write (6,*) (ae(i,j),j=1,3)
        enddo
        write (6,*) ' '
        endif

        ! Find 1st pivot
        a1max = -1.0d0
        do j=1,3
        do i=1,3
          if(abs(ae(i,j)) .gt. a1max) then
            a1max = abs(ae(i,j))
            i1    = i
            j1    = j
          endif
        enddo
        enddo

        ! reduce column of pivot (j1)
        do i=1,3
          if(i .ne. i1) then
            factor = ae(i,j1) / ae(i1,j1)
            do j=1,3
              ae(i,j) = ae(i,j) - factor * ae(i1,j)
            enddo
          endif
        enddo

        if(debug) then
          write (6,*) 'i1,j1=', i1,j1
          do i=1,3
            write (6,*) (ae(i,j),j=1,3)
          enddo
          write (6,*) ' '
         endif

        ! Find 2nd pivot
        a2max = -1.0d0
        do j=1,3
        do i=1,3
          if(abs(ae(i,j)).gt.a2max .and. i.ne.i1 .and. j.ne.j1) then
            a2max = abs(ae(i,j))
            i2    = i
            j2    = j
          endif
        enddo
        enddo

        ! Find remaining column
        if(j1.ne.1 .and. j2.ne.1) j3 = 1
        if(j1.ne.2 .and. j2.ne.2) j3 = 2
        if(j1.ne.3 .and. j2.ne.3) j3 = 3

c        write (6,*) 'j123=', j1, j2, j3

        ! Can now solve for eigenvector (ae(i2,j1) was row reduced to 0)
        ! 0 =                        ae(i2,j2)*ee(j2,ign) + ae(i2,j3)*ee(j3,ign)
        ! 0 = ae(i1,j1)*ee(j1,ign) + ae(i1,j2)*ee(j2,ign) + ae(i1,j3)*ee(j3,ign)
        !
        ee(j3,ign) = 1.0
        ee(j2,ign) = -ae(i2,j3)/ae(i2,j2)
        ee(j1,ign) = -(ae(i1,j2)*ee(j2,ign) + ae(i1,j3)*ee(j3,ign))
     1             /  ae(i1,j1)

        if(debug) then
          write (6,998) "ee       =", (ee(j,ign), j=1,3)
          write (6,998) "ae(i1,:) =", (ae(i1, j), j=1,3)
998       format(a10, 3f10.5)
          write (6,*) ' '
          write (6,*) ' '
        endif

c        ireduce = 1
c        do i=ireduce+1,3
c          factor = ae(i,ireduce) / ae(ireduce,ireduce)
c          do j=1,3
c            ae(i,j) = ae(i,j) - factor * a(ireduce,j)
c          enddo
c        enddo
c        ee(2,ign) = 1.0d0
c        ee(3,ign) = -ae(2,2) /ae(2,3)
c        ee(1,ign) = -(ae(1,2)*ee(2,ign) + ae(1,3)*ee(3,ign)) / ae(1,1)

         call dnormalize(ee(1,ign))
      enddo

      if(debug) then
        write (6,*) "before graham-shmidtt for ee3:  ee1, ee3"
        write (6,*) (ee(i,1),i=1,3)
        write (6,*) (ee(i,3),i=1,3)
        write (6,*) " "
      endif

      ! Make sure e3 is normal to e1
      d13 = ddot(ee(1,1),ee(1,3))
      do i=1,3
        ee(i,3) = ee(i,3) - d13*ee(i,1)
      enddo

      if(debug) then
        write (6,*) "after ee3 = ee3 - d13*ee1:  ee1 then ee3"
        write (6,*) "d13 = ", d13
        write (6,*) (ee(i,1),i=1,3)
        write (6,*) (ee(i,3),i=1,3)
        write (6,*) " "
      endif

      call dnormalize(ee(1,3))

      if(debug) then
        write (6,*) "after normailize ee3:  ee3"
        write (6,*) (ee(i,3),i=1,3)
        write (6,*) " "
      endif

      ! Construct e2 normal to e1 and e3
      ee(1,2) = ee(2,3)*ee(3,1) - ee(3,3)*ee(2,1)
      ee(2,2) = ee(3,3)*ee(1,1) - ee(1,3)*ee(3,1)
      ee(3,2) = ee(1,3)*ee(2,1) - ee(2,3)*ee(1,1)

      if(debug) then
        write (6,*) "before  normailization of ee2:  ee2"
        write (6,*) (ee(i,2),i=1,3)
        write (6,*) " "
      endif

      call dnormalize(ee(1,2))
      if(debug) write (6,*) "anorm ee2:", (ee(i,2),i=1,3)
      if(debug) then
        write (6,*) "after  normailization of ee2:  ee2"
        write (6,*) (ee(i,2),i=1,3)
        write (6,*) " "
      endif

      ierr = 0
      return
      end

c---------------------------------------------------------------------72
      subroutine dnormalize(v)
      real*8 v(3), fNorm
      fNorm = 1.0d0 / sqrt(v(1)**2+v(2)**2+v(3)**2)
      v(1) = fNorm * v(1)
      v(2) = fNorm * v(2)
      v(3) = fNorm * v(3)
      return
      end

c---------------------------------------------------------------------72
      real*8 function ddot(a,b)
      real*8 a(3),b(3)
      ddot = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
      return
      end

c---------------------------------------------------------------------72
