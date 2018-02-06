
      subroutine get_dump_header_info()
      implicit none
      include 'header_info_e3d.h'
      character*256 cline, ckey
      integer i

      ! Parse dump file for available variables and data format & size
      !
      open(unit=17, file=cDumpFile, form='formatted', err=180,
     1     status='old')
100   continue
      do i = 1, 256
        cline(i:i) = ' '
      enddo
      read (17,991,end=200) cline
991   format(a256)

      if(cline(1:10).eq.'time      ') read (cline,*) ckey,time
      if(cline(1:10).eq.'dtime     ') read (cline,*) ckey,dtime
      if(cline(1:10).eq.'step      ') read (cline,*) ckey,istep
      if(cline(1:10).eq.'ncube     ') read (cline,*) ckey,ncube
      if(cline(1:10).eq.'i2max     ') read (cline,*) ckey,i2max
      if(cline(1:10).eq.'i2off     ') read (cline,*) ckey,i2off
      if(cline(1:10).eq.'ncubesxyz ') read (cline,*) ckey,ncx,ncy,ncz
      if(cline(1:10).eq.'nbrickxyz ') read (cline,*) ckey,nbx,nby,nbz
      if(cline(1:10).eq.'nbpf      ') read (cline,*) ckey,nbricksperfile
      if(cline(1:10).eq.'xlimits   ') read (cline,*) ckey,xbndyL,xbndyR
      if(cline(1:10).eq.'ylimits   ') read (cline,*) ckey,ybndyL,ybndyR
      if(cline(1:10).eq.'zlimits   ') read (cline,*) ckey,zbndyL,zbndyR
      if(cline(1:10).eq.'lobndry   ') then
         read (cline,*) ckey,cLoXbndry,cLoYbndry,cLoZbndry
      endif
      if(cline(1:10).eq.'hibndry   ') then
         read (cline,*) ckey,cHiXbndry,cHiYbndry,cHiZbndry
      endif
      if(cline(1:10).eq.'field3d   ') then
        ndvar = ndvar + 1
        read (cline,*) ckey, cVarSymb(ndvar), cVarName(ndvar),
     1                 cDumpMap(ndvar),dvarmin(ndvar), dvarmax(ndvar)
      endif
      if(cline(1:10).eq.'headsize  ') read (cline,*) ckey,nheadsize

      if(nheadsize .gt. 0) goto 200
      goto 100

180   continue
      write (6,*) "Problem: open failed on dumpfile: ", trim(cDumpFile)
      stop

200   continue
      close(17)

      if(ncube .le. 0) then
        nx = ncx
        ny = ncy
        nz = ncz
      else
        nx = ncx*ncube
        ny = ncy*ncube
        nz = ncz*ncube
      endif

      return
      end
