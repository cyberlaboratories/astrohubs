
      integer MAXBOXS, MAXLEVELS
      parameter(MAXBOXS=100000)            ! For AMR: Maximum number of bricks
      parameter(MAXLEVELS=1000)            ! For AMR: Maximum number of brick levels
      real    boxsize(3)                   ! For AMR: xyz size of each box on level 1
      real    boxoffset(3)                 ! For AMR: xyz offset of full volume
      integer(kind=8) i64boxoff(MAXBOXS)   ! For AMR: box offset into file
      integer(kind=8) i64boxsize           ! For AMR: size of one box
      integer iboxoff(3,MAXBOXS)           ! For AMR: 3D box offset
      integer iboxrefine(MAXBOXS)          ! For AMR: box refinment factor
      integer nboxs                        ! For AMR: Number of boxs in this dump
      integer icurrent_box                 ! For AMR: current box number
      integer nboxmesh(3)                  ! For AMR: xyz mesh interior of each box
      integer nbdryx,nbdryy,nbdryz         ! For AMR: Boundary depth in data (all sides)
      integer nlevels                      ! For AMR: number of box levels
      integer lev_number(MAXLEVELS)        ! For AMR: number of bricks on level
      integer lev_refine(MAXLEVELS)        ! For AMR: refinment factor of level
      character*256 cboxfile(MAXBOXS)      ! For AMR: box file

      common /amrinfo/ i64boxsize, i64boxoff, boxsize, boxoffset
      common /amrinfo/ iboxoff, iboxrefine, nbdryx,nbdryy,nbdryz 
      common /amrinfo/ nboxs, nboxmesh
      common /amrinfo/ nlevels, lev_number, lev_refine, icurrent_box
      common /amrinfo/ cboxfile

