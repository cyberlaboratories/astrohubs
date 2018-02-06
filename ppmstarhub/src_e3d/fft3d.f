c---------------------------------------------------------------------72

c           call fftxy(nxbof,nybof,ntslab,nzbof,nv,itslab_zoff,ftnv)
      subroutine fftxy(nx,ny,nz,nzfull,nv,izoff,ftnv)
      implicit none
      integer nx,ny,nz,nzfull,izoff, ix,iy,iz,iv,nv
      real ftnv(nx,ny,nz,nv)
      integer nxmax,nymax,nzmax, nx2,ntx,nty,nxt,nyt
      common /com_fftxy/ nxmax,nymax,nzmax, nx2,ntx,nty,nxt,nyt
      real*8 vri
      allocatable vri(:,:,:,:,:)   ! (ix,iri,iv,iy,iz)

      nx2 = nx / 2
      allocate (vri(nx2,2,nv,ny,nz))

      ! This routines assumes nx & ny span XY

      if(izoff .eq. 0) then
        call check_dim(nx,nxmax)  ! check that dimentions are supported
        call check_dim(ny,nymax)
        call check_dim(nzfull,nzmax)
        call gen_tiles(nx,ny,nzfull,nv,ntx,nty,nxt,nyt)
        call open_fft_files()
      endif

      call do_fft_x(nx,nx2,ny,nz,nxmax,nv,ftnv,vri)
      call do_fft_y(nx,nx2,ny,nz,nymax,nv,vri)
      call write_fft_data(nx,ny,nz,nv,vri)

      deallocate (vri)

      if(izoff+nz .ge. nzfull) call close_fft_files()

      return
      end

c      if(izoff .eq. 0) open(31,file="ft3v_vs_z.dat",form="formatted")
c      do iz = 1, nz
c        write (31,999) iz+izoff,ft3v(1,1,iz,1),ft3v(nx,ny,iz,1),
c     1                          ft3v(1,1,iz,2),ft3v(nx,ny,iz,2),
c     1                          ft3v(1,1,iz,3),ft3v(nx,ny,iz,3)
c999     format(i4, 1p6e12.4)
c      enddo
c      if(izoff+nz .ge. nzfull) close(31)
c---------------------------------------------------------------------72
      subroutine do_fft_z_spc(nzfull,nv,icounts,spc)
      implicit none
      integer nzfull, nv
      integer nxmax,nymax,nzmax, nx2,ntx,nty,nxt,nyt
      common /com_fftxy/ nxmax,nymax,nzmax, nx2,ntx,nty,nxt,nyt
      integer itx,ity,iu, iunit, ix,iy,iz,iv
      real*8 vri
      allocatable vri(:,:,:,:,:)
      character*12 cfile
      real*8 ar3(4096,3), ai3(4096,3)
      real*8 fkx,fky,fkz, fnfft22, fkmag, fkmag2
      real*8 spccmpr,spccmpi,spccmp,spctot
      real*8 spc(0:4096,3)
      integer nfft, ixf, iyf, izf, k, icounts(0:4096)

      allocate (vri(nxt/2,2,nv,nyt,nzfull))

      call fftset(nzmax)
      nfft = nzfull
      fnfft22 = float((nfft/2)**2)

      do k = 0, nfft
        spc(k,1)   = 0.0d0
        spc(k,3)   = 0.0d0
        icounts(k) = 0
      enddo

      do ity = 0,nty-1
      do itx = 0,ntx-1
        iu    = itx + ntx*ity
        call gen_name(iu,cfile)
        write (6,*) "reading file: ", cfile
#ifndef ISGFORTRAN
        open(30,file=cfile,form="binary")
#else
        open(30,file=cfile,form="unformatted", access="stream")
#endif
        read (30) vri
        close(30)

        do iy = 1,nyt
        do ix = 1,nxt/2
          do iv = 1,nv
            do iz = 1, nzfull
              ar3(iz,iv) = vri(ix,1,iv,iy,iz)
              ai3(iz,iv) = vri(ix,2,iv,iy,iz)
            enddo
            call fft(ar3(1,iv), ai3(1,iv), -nzmax)
          enddo

          ixf = ix + nxt*itx
          iyf = iy + nyt*ity
          do izf = 1, nfft
            fkx = dble(ixf-1)
            fky = dble(iyf-1)
            fkz = dble(izf-1)
            if(ixf .gt. nfft/2) fkx = dble(ixf-nfft-1)
            if(iyf .gt. nfft/2) fky = dble(iyf-nfft-1)
            if(izf .gt. nfft/2) fkz = dble(izf-nfft-1)

              fkmag2 = fkx**2 + fky**2 + fkz**2
            if(fkmag2 .lt. fnfft22) then

              fkmag = sqrt(fkmag2)
              k     = int(0.5d0 + fkmag)
              icounts(k) = icounts(k) + 1

              if(nv .eq. 1) then
                spctot = ar3(izf,1)**2 + ai3(izf,1)**2
                spc(k,3) = spc(k,3) + spctot  ! total L2 "energy" in spectrum
              else
                spccmpr = fkx*ar3(izf,1)+fky*ar3(izf,2)+fkz*ar3(izf,3)
                spccmpi = fkx*ai3(izf,1)+fky*ai3(izf,2)+fkz*ai3(izf,3)
                spccmp  = (spccmpr**2 + spccmpi**2) / max(0.01, fkmag2)
                spctot = ar3(izf,1)**2 + ar3(izf,2)**2 + ar3(izf,3)**2
     1                 + ai3(izf,1)**2 + ai3(izf,2)**2 + ai3(izf,3)**2
                spc(k,1) = spc(k,1) + spccmp  ! comprissional component
                spc(k,3) = spc(k,3) + spctot  ! total L2 "energy" in spectrum
            endif

            endif
          enddo  ! izf loop
        enddo    ! ix loop
        enddo    ! iy loop
      enddo      ! itx loop
      enddo      ! ity loop

      deallocate(vri)

      return
      end
c---------------------------------------------------------------------72
      subroutine close_fft_files()
      implicit none
      integer nxmax,nymax,nzmax, nx2,ntx,nty,nxt,nyt
      common /com_fftxy/ nxmax,nymax,nzmax, nx2,ntx,nty,nxt,nyt
      integer itx,ity,iu, iunit

      do ity = 0,nty-1
      do itx = 0,ntx-1
        iu    = itx + ntx*ity
        iunit = 30 + iu
        close(iunit)
      enddo
      enddo
      return
      end
c---------------------------------------------------------------------72
      subroutine gen_name(iu,cfile)
      implicit none
      integer iu
      character*12 cfile
      cfile = "fftdata_0000"
      cfile( 9: 9) = char(ichar('0') +     iu/1000   )
      cfile(10:10) = char(ichar('0') + mod(iu/100,10))
      cfile(11:11) = char(ichar('0') + mod(iu/10 ,10))
      cfile(12:12) = char(ichar('0') + mod(iu    ,10))
      return
      end
c---------------------------------------------------------------------72
      subroutine open_fft_files()
      implicit none
      integer nxmax,nymax,nzmax, nx2,ntx,nty,nxt,nyt
      common /com_fftxy/ nxmax,nymax,nzmax, nx2,ntx,nty,nxt,nyt
      integer itx,ity,iu, iunit
      character*12 cfile
      do ity = 0,nty-1
      do itx = 0,ntx-1
        iu    = itx + ntx*ity
        iunit = 30 + iu
        call gen_name(iu,cfile)
#ifndef ISGFORTRAN
        open(unit=iunit,file=cfile,form='binary')
#else
        open(unit=iunit,file=cfile,form='unformatted',access="stream")
#endif
      enddo
      enddo

      return
      end

c---------------------------------------------------------------------72
      subroutine do_fft_x(nx,nx2,ny,nz,nxmax,nv,ftnv,vri)
      implicit none
      integer nx,nx2,ny,nz,nxmax,  ix,iy,iz,iv,nv
      real ftnv(nx,ny,nz,nv)
      real*8 vri(nx2,2,nv,ny,nz)
      real*8 ar(4097), ai(4097)

      call fftset(nxmax)

      do iv = 1, nv
      do iz = 1, nz
      do iy = 1, ny
        do ix = 1, nx
          ar(ix) = ftnv(ix,iy,iz,iv)
          ai(ix) = 0.0d0
        enddo
        call fft(ar,ai,-nxmax)
        do ix = 1, nx2
          vri(ix,1,iv,iy,iz) = ar(ix)
          vri(ix,2,iv,iy,iz) = ai(ix)
        enddo
      enddo
      enddo
      enddo

      return
      end

c---------------------------------------------------------------------72
      subroutine do_fft_y(nx,nx2,ny,nz,nymax,nv,vri)
      implicit none
      integer nx,nx2,ny,nz,nymax,  ix,iy,iz,iv,nv
      real*8 vri(nx2,2,nv,ny,nz)
      real*8 ar(4097), ai(4097)

      call fftset(nymax)

      do iv = 1, nv
      do iz = 1, nz
      do ix = 1, nx2
        do iy = 1, ny
          ar(iy) = vri(ix,1,iv,iy,iz)
          ai(iy) = vri(ix,2,iv,iy,iz)
        enddo
        call fft(ar,ai,-nymax)
        do iy = 1, ny
          vri(ix,1,iv,iy,iz) = ar(iy)
          vri(ix,2,iv,iy,iz) = ai(iy)
        enddo
      enddo
      enddo
      enddo

      return
      end
c---------------------------------------------------------------------72
      subroutine write_fft_data(nx,ny,nz,nv,vri)
      implicit none
      integer nx,ny,nz,  ix,iy,iz,iv,nv
      integer nxmax,nymax,nzmax, nx2,ntx,nty,nxt,nyt
      common /com_fftxy/ nxmax,nymax,nzmax, nx2,ntx,nty,nxt,nyt
      real*8 vri(nx2,2,nv,ny,nz)
      integer itx,ity,iu, iunit

      if(ntx*nty .eq. 1) then
        write (30) vri
        return
      endif

c      do ity = 0,nty-1
c      do itx = 0,ntx-1
c        iu    = itx + ntx*ity
c        iunit = 30 + iu
c        write (iunit) ...
c      enddo
c      enddo

      write (6,*) "write_fft_data ntx*nty>1 not supported yet"
      stop
      end

c---------------------------------------------------------------------72
      subroutine gen_tiles(nx,ny,nzfull,nv,ntx,nty,nxt,nyt)
      implicit none
      integer nx,ny,nzfull,ntx,nty,nxt,nyt, nv
      integer mb_full, mb_max

c      mb_max = 1024   ! Max size (in MB) ofa3ri array
      mb_max = 3100   ! Max size (in MB) ofa3ri array
      mb_full = int(0.1+(8d0*dble(nx*ny)*dble(nzfull*nv)/1024d0**2))
      ntx = 1
      nty = 1

10    if(mb_full .gt. mb_max*nty) nty = nty * 2
      if(mb_full .gt. mb_max*nty .and. nty.lt.ny) go to 10

20    if(mb_full .gt. mb_max*ntx*nty) ntx = ntx * 2
      if(mb_full .gt. mb_max*ntx*nty .and. ntx.lt.nx) go to 20

      nxt = nx / ntx
      nyt = ny / nty

      if(ntx*nty .gt. 1) then
       write (6,*) "mb_max = ", mb_max
       write (6,*) "mb_full = ", mb_full
       write (6,*) "nx,ny,nzfull = ", nx,ny,nzfull
       write (6,*) "nv = ", nv
       write (6,*) "ntx,nty = ", ntx,nty
      endif

      return
      end
c---------------------------------------------------------------------72
      subroutine check_dim(n,nmax)
      implicit none
      integer n,nmax, i,j,ic, ifa, nfacs, ifac(100), nfac(100)
      common /supported_factors/ nfacs, ifac, nfac

      call set_supported_factors()
      ic = 0
      j = n
10    i = j
      do ifa = 1, nfacs
        if(ifac(ifa)*(j/ifac(ifa)) .eq. j) then
          j = j/ifac(ifa)
          nfac(ifa) = nfac(ifa) + 1
        endif
      enddo
      ic = ic + 1
      if(i .gt. j  .and.  j .gt. 1  .and.  ic.lt.100) go to 10

      if(j .ne. 1) then
        write (6,*) "PROBLEM: dimension n is not supported"
        write (6,*) "   n = ", n
        write (6,*) "   Supported facgtors:"
        do ifa = 1, nfacs
          write (6,*) ifac(ifa)
        enddo
        stop
      endif
      nmax = nfac(1)   ! n = 2**nmax
   
      return
      end

c---------------------------------------------------------------------72
      subroutine set_factor(ifac_in)
      implicit none
      integer ifac_in, nfacs, ifac(100), nfac(100)
      common /supported_factors/ nfacs, ifac, nfac
      nfacs = nfacs + 1
      ifac(nfacs) = ifac_in
      nfac(nfacs) = 0
      return
      end

c---------------------------------------------------------------------72
      subroutine set_supported_factors()
      implicit none
      integer nfacs, ifac(100), nfac(100)
      common /supported_factors/ nfacs, ifac, nfac
      nfacs = 0
      call set_factor(2)     ! Currently, only powers of 2 are currently supported
      return
      end

c---------------------------------------------------------------------72
c        A very short, cheap, and dirty implementation of a
c        FAST FOURIER TRANSFORM
c
c        Note: w0, and w1 need to be set previously in
c              subroutine fftset(mmax)
c
      subroutine fft(ar,ai,nmax)
c     ******************************************************************
c     * THIS CODE WAS WRITTEN BY:  DAVID PORTER  ALL RIGHTS RESERVED   *
c     ******************************************************************
      parameter (mp2 = 4097)
      double precision ar(1), ai(1), br(mp2), bi(mp2), fninv
      double precision w0r(mp2,20),w1r(mp2,20),w0i(mp2,20),w1i(mp2,20)
      common /fftbin/ w0r, w1r, w0i, w1i, lastmmax

      mmax = iabs(nmax)
      n = 2 ** mmax
      n2 = n / 2

      do 1000 m = 1,mmax
      ifund = 2 ** (m-1)
      do 500 ii = 1,n
      i = ii - 1
      i1 = i / ifund
      i2 = i - ((i1+1) / 2) * ifund   + 1
      br(ii) = w0r(ii,m) * ar(i2)    - w0i(ii,m) * ai(i2)
     1       + w1r(ii,m) * ar(i2+n2) - w1i(ii,m) * ai(i2+n2)
      bi(ii) = w0r(ii,m) * ai(i2)    + w0i(ii,m) * ar(i2)
     1       + w1r(ii,m) * ai(i2+n2) + w1i(ii,m) * ar(i2+n2)
500   continue

      do 600 ii = 1,n
      ar(ii) = br(ii)
600   ai(ii) = bi(ii)

1000  continue

      if(nmax .gt. 0) return
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c            do below for Fourier analysis                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      fninv = 1. / float(n)
      do 1200 ii = 1,n
      ar(ii) = ar(ii) * fninv
1200  ai(ii) = ai(ii) * fninv

      return
      end
c---------------------------------------------------------------------72
c
      subroutine fftset(mmax)
c     ******************************************************************
c     * THIS CODE WAS WRITTEN BY:  DAVID PORTER  ALL RIGHTS RESERVED   *
c     ******************************************************************
      parameter (mp2 = 4097)
      double precision w0r(mp2,20),w1r(mp2,20),w0i(mp2,20),w1i(mp2,20)
      double precision factor, arg, pi, fninv
      common /fftbin/ w0r, w1r, w0i, w1i, lastmmax

      data lastmmax/0/
      if(mmax .eq. lastmmax) return
      lastmmax = mmax

      pi = acos(-1.)

      n = 2 ** mmax
      factor = 2. * pi / float(n)
      do 1000 m = 1,mmax
      ifund = 2 ** (m-1)
      do 500 ii = 1,n
      i = ii - 1
      i1 = i / ifund
      ind = mod(i1, 2) * ifund * (i1 / 2)
c      arg = factor * float(ind)
      arg = factor * float(ind)
      w0r(ii,m) = dcos(arg)
      w1r(ii,m) = w0r(ii,m) * sign(1., float(-mod(i1, 2)))
      w0i(ii,m) = dsin(arg)
      w1i(ii,m) = w0i(ii,m) * sign(1., float(-mod(i1, 2)))
500   continue
1000  continue

      return
      end
c---------------------------------------------------------------------72
