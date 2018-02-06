c---------------------------------------------------------------------72
c        Program: vert2stats
c
c      x            yav         ydisp        yrms        ymin         ymax
c  -7.4219e-01  -5.4644e+00   1.1451e+00   5.5831e+00  -5.9999e+00  -1.2335e+00

      character*256 cLine
      nc = 0
      rmssum = 0.
      avrsum = 0.
10    continue
      read (5,998,end=100) cLine
998   format(a256)
      if(cLine(1:1) .eq. '#') go to 10
      read (cLine,*) x, yav, ydisp, yrms, ymin, ymax
      nc = nc + 1
      avrsum = avrsum + yav
      fmssum = fmssum + yrms**2
      if(nc .eq. 1) then
        fmin = ymin
        fmax = ymax
        xmax = x
      else
        if(ymax .gt. fmax) xmax = x
        fmin = amin1(fmin, ymin)
        fmax = amax1(fmax, ymax)
      endif
      go to 10

100   continue

      avr   = avrsum / float(nc)
      rms2  = fmssum / float(nc)
      rms   = sqrt(rms2)
      disp2 = rms2 - avr**2
      disp  = 0.0
      if(disp2 .gt. 0) disp = sqrt(disp2)
      if(fmin .eq. 0.0  .and.  fmax .eq. 0.0) then
        disp = 0.0
        rms  = 0.0
      endif

      write (6,999) avr, disp, rms, fmin, fmax, xmax
999   format("  ", 1p6e14.6)

      stop
      end
