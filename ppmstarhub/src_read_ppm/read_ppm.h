
      integer    maxppmvar, maxrad
      parameter (maxppmvar = 100)
      parameter (maxrad    = 20001)
      real argsimg(maxppmvar), pvar000(maxppmvar)
      real pvar001inv(maxppmvar), svar001inv(maxppmvar)
      real time, rrho(maxrad,2), rprs(maxrad,2)
      real xsize,ysize,zsize, dx, dy, dz, x0f,y0f,z0f, xyzranges(6)

      integer ncpbx, ncpby, ncpbz, nbpb, nvpc
      integer nbrkx, nbrky, nbrkz, nfilex, nfiley, nfilez
      integer nodsx, nodsy, nodsz, nppmvars
      integer nReadVars, iReadVars(maxppmvar)
      integer ilocppmvars(maxppmvar), idelppmvars(maxppmvar)
      integer ix0f,iy0f,iz0f,nxf,nyf,nzf, nxff,nyff,nzff
      integer ndbx, ndby, ndbz
      integer nblend, nxfb,nyfb,nzfb, nfmoms, nfffvars
      integer ncycle, nrrho, nrprs
      real bsizex, bsizey, bsizez
      integer jRho, jRhoUx, jRhoUy,jRhoUz,jFv,nfilter,nfvars
      integer nresppmvars(maxppmvar)

      character*16 cppmvars(maxppmvar)
      character*16 cfffvars(2*maxppmvar)
      character*16 cppmmap(maxppmvar)
      character*256 cRootFile, cReadVars, cOutFile
      common /ppm_scalings/ argsimg, pvar000, pvar001inv, svar001inv
      common /ppm_scalings/ time, rrho, rprs
      common /ppm_scalings/ xsize,ysize,zsize, dx, dy, dz
      common /ppm_scalings/ x0f,y0f,z0f, xyzranges
      common /ppm_scalings/ ncpbx, ncpby, ncpbz, nbpb, nvpc
      common /ppm_scalings/ nbrkx, nbrky, nbrkz, nfilex, nfiley, nfilez
      common /ppm_scalings/ nodsx, nodsy, nodsz
      common /ppm_scalings/ ix0f,iy0f,iz0f,nxf,nyf,nzf, ndbx,ndby,ndbz
      common /ppm_scalings/ nxff,nyff,nzff,nfvars, nppmvars
      common /ppm_scalings/ nReadVars, iReadVars,nblend, nxfb,nyfb,nzfb
      common /ppm_scalings/ ilocppmvars, idelppmvars
      common /ppm_scalings/ nresppmvars, ncycle, nrrho, nrprs
      common /ppm_scalings/ jRho, jRhoUx, jRhoUy,jRhoUz,jFv,nfilter
      common /ppm_scalings/ bsizex, bsizey, bsizez,nfmoms,nfffvars
      common /ppm_scalings/ cppmvars, cppmmap, cRootFile, cReadVars
      common /ppm_scalings/ cfffvars, cOutFile
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

