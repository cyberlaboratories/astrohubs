
      integer MAX_FORMULAS, MAXDVAR
      parameter (MAX_FORMULAS = 4000)
      character*256 cFormulas(MAX_FORMULAS) ! List of all formulas
      character*256 cFunction(MAX_FORMULAS) ! Function performed on input varialbles
      character*256 cFld1,cFld2            ! Fields for correlate
      character*256 cField                 ! Field to generate
      character*256 cFieldLabel            ! Label for output Field
      character*256 cOutput                ! Output format
      character*32  cPart                  ! Part (subset or moment) of output
      character*256 cOutSymb               ! Name of output field(s)
      character*256 cOutDir                ! Output directory
      character*256 cView                  ! View input file
      character*256 cRegion                ! Region name for labeling
      integer nFunction                    ! Number of functions to perform
      integer nFormulas                    ! Number of formulas
      integer nslab                        ! Z interior thickness of slab
      integer ntslab                       ! slab thickness for building block-trees
      integer nzbslab                      ! blended slab thickness
      integer nxbof,nybof,nzbof            ! BOF dimensions
      integer(kind=8) n64slabsize          ! 64-bit size of bof write
      integer itslab_zoff                  ! slab offset in z-direction
      integer idebug                       ! debug level
      integer ilogscale                    ! (1/0) ==> (use/don't use) log scale in plots
      integer isizelogscale                ! (1/0) ==> (use/don't use) log scale on sectional size axis
      integer isRZ                         ! is 2D radial (R), axial (Z) coordinates
      real    vview(3,3)                   ! Horizontal, Vertical, to_eye unit vectors of view
      real    HVMIN, HVMAX, hvscale        ! limits of HV range & scale
      real    xmin,ymin,zmin               ! Min of XYZ range to evaluate in
      real    xmax,ymax,zmax               ! Max of XYZ range to evaluate in
      real    endpt0(3), endpt1(3)         ! end points of line segment
      real    extent                       ! +/- extent of line segment
      real    endpt_vec(3)                 ! Unit vector in direction from endpt0 to endpt1
      real    THRESHOLD                    ! a specified threashold
      real    X00, Y00, Z00                ! X, Y, and/or Z coordinates to evaluate at
      integer iX00, iY00, iZ00             ! X, Y, and/or Z indexes to evaluate at
      real finvalid                        ! not a valid value for anything
      integer NBLEND                       ! Blend facgtor for output
      integer nsmooth                      ! smooth over [ix-nsmooth:ix+nsmooth]
      integer ncbarsteps                   ! number of cbar steps
      character*256 cbarString             ! Color bar for plots
      character*256 cSectCbar              ! Color bar for sectional concentration plots
      character*256 COUTHOST               ! AIO target host
      character*256 COUTFILE               ! AIO output file
      character*256 CINVAR                 ! Input variable name
      integer MAXRES                       ! Max resolution before auto-blending
      integer iflip1                       ! 0: no flip; 1: flip ordinate axis
      integer(kind=8) I64OUTPUT            ! Size of output in bytes
      real SMALL                           ! a general purpose lower bound
      real QCMIXFAC, QCDIVLIM0, QCDIVLIM   ! factors for quality curl calc
      character*256 cCurlCalc              ! face centered calc or cell centered calc
      character*256 cDumpFile0             ! Path to dump file
      character*256 cVersion               ! e3d Version string
      character*16 cxcoord,cycoord,czcoord ! XYZ coordinate labels
      character*6  cxgrad ,cygrad ,czgrad  ! XYZ derivatives: grad_?
      character*256 c2label, c2map         ! Alternative label & map of 2nd coordinate

      common /postinfo/ I64OUTPUT,n64slabsize
      common /postinfo/ SMALL, QCMIXFAC, QCDIVLIM0, QCDIVLIM, MAXRES
      common /postinfo/ vview
      common /postinfo/ x00, y00, z00, iX00, iY00, iZ00, iflip1
      common /postinfo/ cFormulas,cOutput,cPart,cOutSymb,cOutDir,cView
      common /postinfo/ cRegion
      common /postinfo/ COUTHOST, COUTFILE, CINVAR, cField, cSectCbar
      common /postinfo/ cFieldLabel, cbarString, cFld1,cFld2
      common /postinfo/ nFormulas, finvalid, isRZ, isizelogscale
      common /postinfo/ nFunction,cFunction,nzbslab,ncbarsteps
      common /postinfo/ nslab, ntslab, itslab_zoff, idebug, ilogscale
      common /postinfo/ nxbof,nybof,nzbof,xmin,xmax,zmin,zmax
      common /postinfo/ endpt0, endpt1, extent, endpt_vec
      common /postinfo/ hvmin,hvmax, hvscale, THRESHOLD, ymin,ymax
      common /postinfo/ nblend,nsmooth,cDumpFile0,cVersion,cCurlCalc
      common /postinfo/ cxcoord,cycoord,czcoord
      common /postinfo/ cxgrad ,cygrad ,czgrad, c2label, c2map

      integer iop_noop, iop_read, iop_curl, iop_norm,iop_radc,iop_getc
      integer iop_cross, iop_div, iop_dot, iop_plus, iop_minus
      integer iop_times, iop_divide, iop_scale, iop_add, iop_pullb
      integer iop_agradb, iop_grad, iop_cminus, iop_curld, iop_smooth
      integer iop_grad_1,iop_grad_2,iop_grad_3, iop_const, iop_trsub
      integer iop_log, iop_log10, iop_sqrt, iop_sin, iop_cos, iop_exp
      integer iop_trace, iop_ballflt, iop_comp, iop_secsum, iop_det
      integer iop_step, iop_abs, iop_igenv_s, iop_setval
      integer iop_coord1, iop_coord2,  iop_coord3
      parameter (iop_noop   =  0)
      parameter (iop_read   =  1)
      parameter (iop_curl   =  3)
      parameter (iop_norm   =  4)
      parameter (iop_cross  =  5)
      parameter (iop_div    =  6)
      parameter (iop_dot    =  7)
      parameter (iop_pullb  =  8)
      parameter (iop_agradb =  9)
      parameter (iop_grad   = 10)
      parameter (iop_grad_1 = 11)
      parameter (iop_grad_2 = 12)
      parameter (iop_grad_3 = 13)
      parameter (iop_curld  = 14)
      parameter (iop_smooth = 15)
      parameter (iop_trsub  = 16)
      parameter (iop_trace  = 17)
      parameter (iop_ballflt= 18)
      parameter (iop_comp   = 19)
      parameter (iop_secsum = 20)
      parameter (iop_det    = 21)
      parameter (iop_radc   = 22)
      parameter (iop_igenv_s = 23)
      parameter (iop_getc = 24)

      parameter (iop_log    = 101)
      parameter (iop_log10  = 102)
      parameter (iop_sqrt   = 103)
      parameter (iop_sin    = 104)
      parameter (iop_cos    = 105)
      parameter (iop_exp    = 106)
      parameter (iop_const  = 107)
      parameter (iop_step   = 108)
      parameter (iop_abs    = 109)

      parameter (iop_plus   = 201)
      parameter (iop_minus  = 202)
      parameter (iop_times  = 203)
      parameter (iop_divide = 204)
      parameter (iop_scale  = 205)
      parameter (iop_add    = 206)
      parameter (iop_cminus = 207)
      parameter (iop_setval = 208)
      parameter (iop_coord1 = 211)
      parameter (iop_coord2 = 212)
      parameter (iop_coord3 = 213)

      ! Table Of Formulas
      integer MAX_TOF, MAX_DO_TOF
      parameter(MAX_TOF = 4000)
      parameter(MAX_DO_TOF = 4000)
      character*32  cFld_tof(MAX_TOF)     ! Symbol for field
      character*32  cNam_tof(MAX_TOF)     ! Name of field
      integer       nDim_tof(MAX_TOF)     ! Dimension of field (number of work arrays)
      integer       nDim_wrk(MAX_TOF)     ! Dimension of work field (number of work arrays that go together)
      character*32  cOpr_tof(MAX_TOF)     ! Name of operaion
      integer       iOpr_tof(MAX_TOF)     ! Operation number (for computed goto)
      integer      nArg_tof(MAX_TOF)     ! Number of argument to operation
      integer       nBdy_tof(3,MAX_TOF)   ! Number of argument to operation
      character*32  cArg_tof(MAX_TOF,4000) ! Names of argument variables
      character*256 cFile_tof(MAX_TOF)    ! Names of dump file to read from
      character*256 cPar_tof(MAX_TOF)     ! TOF parameter string
      integer       icomp_tof(MAX_TOF)    ! component offset into tensor
      integer       iWrk_tof(MAX_TOF)     ! Work array assigned to field
      integer       iArg_tof(MAX_TOF,4000)  ! Work array assigned to argument variable
      integer       iArg_dim(MAX_TOF,4000)  ! Dimenstion of argument variable
      real          fVal_tof(MAX_TOF)     ! scalar value associated with this field
      integer       iRda_tof(MAX_TOF)     ! Read array assigned to this field (to copy from)
      integer       itof_Rda(MAX_TOF)     ! itof associated with read array (to find file, field, ...)
      integer       ntof                  ! number of fields in table
      integer       nrda                  ! number of read fields in table
      integer       nrdx, nrdy, nrdz      ! Current Dimensions of read array
      integer       nIndexValues          ! Number of index values
      character*2   cIndexValues(64)      ! Array of index values
      integer       nIndexVars            ! Number of index variables
      character*2   cIndexVars(64)        ! Array of index variables
      character*256 c_add_tof_problem    ! Potential reason for a read error
      common /tof/ cFld_tof, nDim_tof, cOpr_tof, iOpr_tof, nDim_wrk
      common /tof/ cNam_tof
      common /tof/ nArg_tof,cArg_tof,iWrk_tof,iArg_tof,nrda,icomp_tof
      common /tof/ nBdy_tof, ntof,iRda_tof,itof_rda,iArg_dim
      common /tof/ cFile_tof,cPar_tof
      common /tof/ nIndexValues, nIndexVars, cIndexValues, cIndexVars
      common /tof/ nrdx, nrdy, nrdz, fVal_tof, c_add_tof_problem

      ! List of operations
      integer MAX_OPS
      parameter (MAX_OPS = 4000)
      integer nops, itof_ops(MAX_OPS)
      common /operations/ nops, itof_ops


      ! Sequence of bricks
      integer MAXBRICKS
      parameter (MAXBRICKS = 400000)
      integer n_rb                         ! Number of bricks in read brick list (rb)
      integer n_wb                         ! Number of bricks in work brick list (wb)
      integer ioff_rb(3,MAXBRICKS)         ! Offset of read brick into full volume
      integer ioff_wb(3,MAXBRICKS)         ! Offset of work brick into full volume
      integer nf_rb(3,MAXBRICKS)           ! Interior size of read brick
      integer n0_wb(3,MAXBRICKS)           ! Interior size of work brick
      integer nboxsused                    ! For AMR: number of boxs to be read in 
      integer  ibox_list(MAXBRICKS)        ! For AMR: list of boxs (w/oundaries) to be read in
      integer  maxrefine                   ! For AMR: maximum amr level used
      integer  maxboxs_db                  ! For AMR: max number of boxs used
      integer nbdy_max(3)                  ! Max boundary depth in each direction
      integer nxfull, nyfull, nzfull       ! Full dimension of problem
      integer izgbllo, izgblhi             ! Global lo & high Z-range needed for output
      integer nfx,nfy,nfz,nfv              ! Full XYZ dimensions of flds array
      integer NFX0, NFY0, NFZ0             ! Internal dims of work array
      integer nfbndx,nfbndy,nfbndz         ! Boundary depth in XYZ directions
      integer nadvect                      ! Maximum number of zones can advect
      common /brickdim/ nbdy_max, nxfull, nyfull, nzfull
      common /brickdim/ izgbllo, izgblhi
      common /brickdim/ nfx,nfy,nfz,nfv, NFX0, NFY0, NFZ0
      common /brickdim/ nfbndx,nfbndy,nfbndz,n_wb,ioff_wb,n0_wb
      common /brickdim/ nboxsused, ibox_list, maxrefine, maxboxs_db
      common /brickdim/ n_rb,ioff_rb,nf_rb, nadvect

      ! read brick to work brick copy info (for a given pair or read and work bricks)
      integer n_cp             ! Number of 3D copy blocks
      integer iwb0_cp(3,27)    ! 3D starting point of of copy block in work array
      integer iwb1_cp(3,27)    ! 3D ending   point of of copy block in work array
      integer irb0_cp(3,27)    ! 3D starting point of of copy block in read array
      integer idel_cp(3,27)    ! 3D increment of copy in read array (0 for continuation)
      integer iboundary(3)     ! Boundary types in XYZ directions (0=periodic, 1=continuation)
      common /rb_to_wb/ n_cp, irb0_cp, iwb0_cp, iwb1_cp, idel_cp
      common /rb_to_wb/ iboundary

      ! output stucture dims
      integer iStr0              ! 1st character of output type in output file
      integer maxcols            ! Max number of columns in profile output
      integer maxprof            ! Max number of rows    in profile output
      integer nbdy_output        ! boundary depth needed in output structure
      real    delbin1            ! bin size of profile or 1st bin
      real    profcoordmin       ! Lower bound of output coordinates
      real binmin, binscale      ! Hist/Dist bin minimum & scale factor
      real bin10,bin11,bin20,bin21 ! Bins for correlate output
      integer nbin1, nbin2         ! bin res for correlate output
      integer iout_rprof
      integer iout_xprof
      integer iout_yprof
      integer iout_zprof
      integer iout_xsect
      integer iout_ysect
      integer iout_hv
      integer iout_xyzavs
      integer iout_vrz, iout_vry, iout_vrx
      integer iout_bof
      integer iout_strip
      integer iout_xstrip
      integer iout_ystrip
      integer iout_xyslice
      integer iout_dist
      integer iout_corr
      integer iout_spc3v
      integer iout_view
      integer iout_spcnv, nvec
      real    dx,dy,dz           ! Mesh spacing in XYZ
      real    dx2i,dy2i,dz2i     ! For centered derivatives
      real    fLimits(6)         ! XYZ min and max of full domain
      integer ixout1, iyout1, izout1  ! XYZ lower bounds of output
      integer ixout2, iyout2, izout2  ! XYZ upper bounds of output
      integer nxout, nyout, nzout     ! XYZ size of output

      common /outstruct/ iStr0, maxcols, maxprof, nbdy_output, delbin1
      common /outstruct/ binmin, binscale, bin10,bin11,bin20,bin21
      common /outstruct/ nbin1, nbin2
      common /outstruct/ profcoordmin,iout_zprof,iout_yprof,iout_xprof
      common /outstruct/ iout_bof, iout_xstrip, iout_xyslice,iout_hv
      common /outstruct/ iout_strip, iout_rprof
      common /outstruct/ iout_ystrip, iout_dist, iout_spc3v, iout_spcnv
      common /outstruct/ iout_xsect, iout_ysect, iout_corr
      common /outstruct/ iout_vrz, iout_vry, iout_vrx, iout_view
      common /outstruct/ dx,dy,dz, dx2i,dy2i,dz2i, fLimits,iout_xyzavs
      common /outstruct/ ixout1,iyout1,izout1, ixout2,iyout2,izout2
      common /outstruct/ nxout, nyout, nzout, nvec


      ! Structure functions
      integer MAXDISPSETS, MAXNDISPINSET
      parameter (MAXDISPSETS=10)
      parameter (MAXNDISPINSET=100)
      integer ndisp_in_set(MAXDISPSETS)
      integer ixdisp(MAXNDISPINSET,MAXDISPSETS)
      integer iydisp(MAXNDISPINSET,MAXDISPSETS)
      integer ndisp_sets, max_disp_in_set, ndiff_sampl
      integer iout_strfn, nstructfn
      real del_vari,del_vari2
      character*256 cdisp_names(MAXDISPSETS)
      common /structfn/ ndisp_in_set,ixdisp,iydisp,ndisp_sets
      common /structfn/ max_disp_in_set, ndiff_sampl,del_vari,del_vari2
      common /structfn/ iout_strfn, nstructfn
      common /structfn/ cdisp_names


      ! Work array
      integer nfx0_max, nfy0_max, nfz0_max, nfb0_max
      parameter (nfx0_max = 2048)
      parameter (nfy0_max = 32)
      parameter (nfz0_max = 32)
      parameter (nfb0_max =  2)
c      parameter (nfv0_max =  40000)
      integer           iwxlo,iwxhi,iwylo,iwyhi,iwzlo,iwzhi,nfv0_max
      common /wrkarray/ iwxlo,iwxhi,iwylo,iwyhi,iwzlo,iwzhi,nfv0_max

      ! Constant values for offsets & coefficients 
      integer    MAX_VALUES
      parameter (MAX_VALUES =  100)
      integer        nValues
      character*256 cValues(MAX_VALUES)
      real          fValues(MAX_VALUES)
      common /const_vals/ nValues, fValues, cValues

      integer    MAX_DUMPS
      parameter (MAX_DUMPS =  10000)
      integer        nDumps
      character*256 cDumps(MAX_DUMPS)
      real          fTimes(MAX_DUMPS)
      common /const_dumps/ nDumps, cDumps, fTimes

      ! Sectional information
      integer nSecFields, MAXSECDIMS
      parameter (MAXSECDIMS=10)
      real fSecF0(MAXSECDIMS), fSecF1(MAXSECDIMS)
      integer isect_min, isect_max, isect_bin(4000)
      character*32 cSecFilter, cPsigns, cMsigns
      character*32 cSecFormat(MAXSECDIMS),cSecName(MAXSECDIMS)
      character*32 cSecUnit(MAXSECDIMS),cSecMap(MAXSECDIMS)
      character*256 cSecTitle, cSecLoc(MAXSECDIMS)
      common /sections_com/ fSecF0, fSecF1, nSecFields
      common /sections_com/ isect_min, isect_max, isect_bin
      common /sections_com/ cSecFilter, cPsigns, cMsigns
      common /sections_com/ cSecFormat, cSecName, cSecUnit, cSecMap
      common /sections_com/ cSecTitle, cSecLoc

      integer nrep_used
      common /amr_com/ nrep_used

