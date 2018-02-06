c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c         E3D - Extensible analysis tool for 3-D data                  c
c----------------------------------------------------------------------c
c  COPYRIGHT David H. Porter      2012 May 12 -  2011 June 21          c
c----------------------------------------------------------------------c
c Versions                                                             c
c    Date    2012 May   12    Initial version                          c
c    Date    2012 June  21    grad operator implemented                c
c    Date    2012 June  22    fixed bug in do_agradb                   c
c    Date    2012 June  28    xyslice implemented                      c
c    Date    2012 Dec   18    xyzavs  implemented                      c
c    Date    2012 Dec   21    data reader routines separated out       c
c    Date    2012 Dec   27    Tecplot data reader implemented          c
c    Date    2013 Apr   27    div for 2D & radial symmetry             c
c    Date    2013 May   30    Face centered curl calc.                 c
c    Date    2013 June  10    HDF5_simple reader                       c
c    Date    2013 June  11    fixed 2* problem in face centered curl   c
c    Date    2013 Oct    2    test/develop face centered curl "curlf"  c
c    Date    2013 Oct    6    added tesnor math in  do_math & gen_ops  c
c    Date    2013 Oct    9    implemented sectional output             c
c    Date    2013 Oct   10    implemented constant field input         c
c    Date    2013 Oct   14    Sectional color bar: cSectCbar           c
c    Date    2013 Dec   30    Implemented spc3v                        c
c    Date    2014 Jan    4    Implemented sectional Histograms (Hsect) c
c    Date    2014 Jan    9    Implemented ball smoothing filter        c
c    Date    2014 Jan   10    Extended dot to do tensor contraction    c
c    Date    2014 Jan   10    Implemented matrix trace                 c
c    Date    2014 Jan   13    implemented gaussian smooth              c
c    Date    2014 Jan   18    Add component noop reference             c
c    Date    2014 Jan   22    Add index reference to arguments         c
c    Date    2014 Jan   23    sect_sum/sect_av function implemented    c
c    Date    2014 Jan   24    index variables implemneted              c
c    Date    2014 Jan   25    nbdy_tof dependencies for components     c
c    Date    2014 Jan   26    Implemented cPar_tof & do_sectional_args c
c    Date    2014 Mar   23    Fixed bug in add_tofs                    c
c    Date    2014 Apr   24    Fixed allocation of flds                 c
c    Date    2014 May   15    Increased range for structure fns        c
c    Date    2014 May   20    Try-cubick interp in pullback (is good)  c
c    Date    2014 May   20    Try-cubick interp in pullback (is good)  c
c    Date    2014 July  19    view fzx ==> fxz so that pairs will work c
c    Date    2014 Aug    5    grad_?L (Left edge gradient) implemented c
c    Date    2014 Sept  16    flip1 implemented for profiles           c
c    Date    2014 Sept  18    flip1 implemented for stip plots         c
c    Date    2014 Sept  22    fied sofplot low value edge case         c
c    Date    2014 Oct.  23    support mlti-level AMR for view          c
c    Date    2015 Jan.  10    support mlti-level AMR for strip         c
c    Date    2015 May   25    Implemented radial profiles (radprof)    c
c    Date    2015 June  25    Implemented vrz volume rendering         c
c    Date    2015 July  11    Implemented vry volume rendering         c
c    Date    2015 July  20    Implemented vrx volume rendering         c
c    Date    2015 Oct.  24    cbar= option to set cbarstring           c
c    Date    2015 Nov.  11    eigenv_of_s: Eigenvalues of sym. matrix  c
c    Date    2015 Nov.  17    added getc (get component as function)   c
c    Date    2015 Nov.  21   added  correlate & vixed dist             c
c    Date    2015 Nov.  30   fixed read brick z range                  c
c    Date    2015 Nov.  30   increased work brick X range to max       c
c    Date    2015 Nov.  30   added progres by read brick & work brick  c
c    Date    2015 Dec.   3   fixed Z ranges in read and work bricks    c
c    Date    2015 Dec.  10   show generates region & vars list (ipnb)  c
c    Date    2016 Jan.  21   idebug >= 30 zranges & reads              c
c    Date    2016 Jan.  21   idebug >= 31 report work bricks           c
c    Date    2016 Jan.  21   idebug == 21 just formula analysis        c
c    Date    2016 Jan.  21   idebug >= 10 TOF values reported          c
c    Date    2016 Mar.   5   implmented coordinate fields iop_coord?  c
c    Date    2016 Mar.  11   implmented smooth_x01y02z03(f)           c
c    Date    2016 Mar.  19   implmented itime--                       c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      implicit none
      include 'postproc_info_e3d.h'

      cVersion = "e3d version: 2013 Oct 14"
      idebug = 0

      call set_postproc_defaults()

      call read_formulas_file()   ! Get formulas for all relevant fields

      call read_input_file()   ! Get data file to read, field to generate, what to output, ..

      call read_header_info(cDumpFile0)  ! Get header info from dump file

      if(cOutput .eq. "show") then
        call report_header_info()
        call gen_regions_vars()
        stop
      endif
      if(cOutput .eq. "time") call report_time()

      call gen_output_par()    ! Generate parameters/dimensions for output structures

      call gen_ops()           ! Generate seqence ofoperations for requested field

      call gen_read_bricks() ! Generate sequence of bricks to read in

       call do_work()

      ! Cleanup: Finalize analysis, close files, ...

      stop
      end

C ====================================================================72
      subroutine do_work()
      implicit none
      include 'postproc_info_e3d.h'
      integer i, j, k, iwrk0, iwrk1, ix, iy, iz, iout_vr
      real*8 tr2w0, t_r2w

      integer*2 slab2i
      real      slab4f
      real      rdas, fout, fout3, bofslab, fxy,fyz,fxz,ft3v, ftnv,fnv
      real      frgb
      real*8    sym33(3,3)
      real*8    profile, average,rms2,disp2,rms,disp,dispmax,dout
      integer   ivbrick, icorrelate
      integer   idispsumL, idispsumT, idispsum2, idispsumS
      character*1 bout
      allocatable ivbrick(:,:,:)   ! (ibx, iby, iv)
      allocatable  slab2i(:,:,:,:) ! (ix , iy , iz , ivbrick)
      allocatable  slab4f(:,:,:,:) ! (ix , iy , iz , ivbrick)
      allocatable    fout(:,:,:)   ! (ix , iy , iz)
      allocatable  fxy(:,:,:)     ! (ix, iy,ivar)
      allocatable  fyz(:,:,:)     ! (iy, iz,ivar)
      allocatable  fxz(:,:,:)     ! (ix, iz,ivar)
      allocatable   fout3(:,:,:,:) ! (ix , iy , iz, 3)
      allocatable    bout(:,:,:)   ! (ix , iy , iz)
      allocatable bofslab(:,:,:)   ! (ix , iy , iz)
      allocatable icorrelate(:,:)  ! (field1, field2)
      allocatable    ft3v(:,:,:,:) ! (ix,iy,iz,iv)
      allocatable    ftnv(:,:,:,:) ! (ix,iy,iz,iv)
      allocatable     fnv(:,:,:,:) ! (ix,iy,iz,iv)
      allocatable idispsumS(:,:,:) ! (idiff,idisp,iset)
      allocatable idispsum2(:,:,:) ! (idiff,idisp,iset)
      allocatable idispsumL(:,:,:) ! (idiff,idisp,iset)
      allocatable idispsumT(:,:,:) ! (idiff,idisp,iset)
c count
      real           flds
      allocatable    flds(:,:,:,:) ! (ix , iy , iz , iv)
      allocatable    rdas(:,:,:,:) ! (ix , iy , iz , iv)
      allocatable profile(:,:)     ! (coordinate , channel)

      ! scratch 
      integer i_rb, iop, iv, i_wb
      integer ido, itof
      real*8 OMP_GET_WTIME, t_start, t_end, t0, t_read, t_work, t_out
      real*8 t_full

      t_start = OMP_GET_WTIME()
      t_r2w  = 0.0d0
      t_read = 0.0d0
      t_work = 0.0d0
      t_out  = 0.0d0

      ! ***************************************************************
      ! *               INTTIALIZE OUTPUT STURCTURES                  *
      ! ***************************************************************
      if(maxprof  .gt. 0) then
        allocate (profile(maxprof,maxcols))   ! (coordinate , channel)

        do i=1,maxprof
           profile(i,1) = profcoordmin + delbin1 * (0.5 + float(i-1))
        enddo
        if(maxcols .eq. 2) then
          do i=1,maxprof
            profile(i,2) = 0.0
          enddo
        endif
        if(maxcols .eq. 7) then
          do i=1,maxprof
            profile(i,2) = 0.0
            profile(i,3) = 0.0
            profile(i,5) =  1.0d+33
            profile(i,6) = -1.0d+33
            profile(i,7) = 0.0d0
          enddo
        endif
      endif

      if(iout_corr .eq. 1) then
        allocate (icorrelate(nbin1,nbin2))
        do i = 1, nbin1
        do j = 1, nbin2
          icorrelate(i,j) = 0
        enddo
        enddo
      endif
      if(iout_view .eq. 1) then
        do i = 1, 3
        do j = 1, 3
          sym33(i,j) = 0.0d0
        enddo
        enddo
      endif

      if(iout_hv  .eq. 1) allocate (bout(nxfull,nyfull,ntslab)) ! (ix , iy , iz)
      if(iout_xyslice  .eq. 1) allocate (bout(nxfull,nyfull,1)) ! (ix , iy)
      iout_vr = max(iout_vrz, iout_vry, iout_vrx)
      if(iout_vr .ge. 1) then
        allocate (fxy(nxout,nyout,3))
        allocate (fxz(nxout,nzout,3))
        allocate (fyz(nyout,nzout,3))
      endif
      if(iout_xyzavs+iout_xstrip+iout_ystrip  .ge. 1) then
        allocate (fxy(nxout,nyout,3))
        allocate (fyz(nyout,nzout,3))
        allocate (fxz(nxout,nzout,3))
      endif

      if(iout_bof .eq. 1) then
        allocate (fout(nxbof,nybof,nzbslab))    ! (ix , iy , iz)
        do iz = 1, nzbslab
        do iy = 1, nybof
        do ix = 1, nxbof
          fout(ix,iy,iz) = 0.0
        enddo
        enddo
        enddo
      endif

      if(iout_spc3v .eq. 1) then
        allocate (ft3v(nxbof,nybof,ntslab,3))    ! (ix , iy , iz, iv)
      endif

      if(iout_spcnv .eq. 1) then
        allocate (ftnv(nxbof,nybof,ntslab,nvec))    ! (ix , iy , iz)
      endif

      if(iout_strfn .eq. 1) then
        allocate (fnv(nxbof,nybof,ntslab,nvec))    ! (ix , iy , iz)
        if(nvec .eq. 1) then
        allocate (idispsumS(-ndiff_sampl:ndiff_sampl,
     1                      max_disp_in_set, ndisp_sets))     ! (idiff,idisp,iset)
        do iz=1,ndisp_sets
        do iy=1,max_disp_in_set
        do ix=-ndiff_sampl,ndiff_sampl
          idispsumS(ix,iy,iz) = 0
        enddo
        enddo
        enddo

        endif
      endif

      ! Allocate flds work array
      iwxlo = 1-nbdy_max(1)
      iwylo = 1-nbdy_max(2)
      iwzlo = 1-nbdy_max(3)
      iwxhi = min(nxfull, nfx0_max)+nbdy_max(1)
      iwyhi = nfy0_max+nbdy_max(2)
      iwzhi = min(nfz0_max,1+izgblhi-izgbllo) + 2*nbdy_max(3)
      nfv0_max = nfv
      if(idebug .ge. 30) then
        write (6,*) "iwxlo,iwxhi = ", iwxlo,iwxhi
        write (6,*) "iwylo,iwyhi = ", iwylo,iwyhi
        write (6,*) "iwzlo,iwzhi = ", iwzlo,iwzhi
        write (6,*) "nfv0_max    = ", nfv0_max
        write (6,*) "flds [MB],= ", float(1+iwxhi-iwxlo) *
     1                              float(1+iwyhi-iwylo) *
     1                              float(1+iwzhi-iwzlo) *
     1                              float(nfv0_max) / (1024.0**2)
      endif
      allocate (flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max))

      ! ***************************************************************
      ! *             READ & WORK LOOPS                               *
      ! ***************************************************************
      nrdx = 0
      nrdy = 0
      nrdz = 0

      ! Loop of read-brick (_rb) regions
      do i_rb = 1, n_rb
        t0 = OMP_GET_WTIME()
        if(idebug .ge. 30) then
          if(idebug .ge. 31) write (6,*) " "
          write (6,999) ioff_rb(3,i_rb)+1,ioff_rb(3,i_rb)+nf_rb(3,i_rb)
999       format("Read full Z-range:     ", 2i5)
        endif
        ! allocate/reallocate read arrays
        if(nrda .gt. 0) then
          if(nrdx .ne. nf_rb(1,i_rb) .or.
     1       nrdy .ne. nf_rb(2,i_rb) .or.
     1       nrdz .ne. nf_rb(3,i_rb)     ) then
            if(nrdx .gt. 0) deallocate(rdas)
            nrdx = nf_rb(1,i_rb)
            nrdy = nf_rb(2,i_rb)
            nrdz = nf_rb(3,i_rb)
            allocate (rdas(nrdx,nrdy,nrdz,nrda))   ! (ix,iy,iz,iv)
          endif
         endif

        ! Read in all fields in this read brick sub-region
        do iv = 1, nrda
          call read_fld(iv, i_rb, rdas(1,1,1,iv))
        enddo
        t_read = t_read + (OMP_GET_WTIME() - t0)

        ! Generate & loop over work bricks in this read brick sub-region
        call gen_work_bricks(i_rb)

        do i_wb = 1, n_wb
          ! Generate list of 3D block copies for current pair of read & work bricks
          call gen_cp_list(i_rb,i_wb)

          if(idebug .ge. 31) call print_wbxyz(i_wb)

          ! Loop over operations to genrate field
          t0 = OMP_GET_WTIME()
          do iop = 1, nops
            itof = itof_ops(iop)
            ido = iOpr_tof(itof)
            if(ido .eq. iop_read ) then
               tr2w0 = OMP_GET_WTIME()
               call r2w(itof,i_rb,i_wb,rdas,flds)
               t_r2w = t_r2w + (OMP_GET_WTIME() - tr2w0)
            endif
            if(ido .eq. iop_norm ) call do_norm(itof, i_wb, flds)
            if(ido .eq. iop_getc ) call do_getc(itof, i_wb, flds)
            if(ido .eq. iop_radc ) call do_radc(itof, i_wb, flds)
            if(ido .eq. iop_curl ) call do_curl(itof, i_wb, flds)
            if(ido .eq. iop_curld) call do_curld3(itof, i_wb, flds)
            if(ido .eq. iop_igenv_s ) call do_eigenv_s(itof, i_wb, flds)
            if(ido .eq. iop_div  ) call do_div(itof, i_wb, flds)
            if(ido .eq. iop_grad ) call do_grad(itof, i_wb, flds)
            if(ido .eq. iop_grad_1) call do_grad_1(itof, i_wb, flds)
            if(ido .eq. iop_grad_2) call do_grad_2(itof, i_wb, flds)
            if(ido .eq. iop_grad_3) call do_grad_3(itof, i_wb, flds)
            if(ido .eq. iop_agradb) call do_agradb(itof, i_wb, flds)
            if(ido .eq. iop_cross) call do_cross(itof, i_wb, flds)
            if(ido .eq. iop_dot  ) call do_dot(itof, i_wb, flds)
            if(ido .eq. iop_pullb) call do_pullb(itof, i_wb, flds)
            if(ido .eq. iop_smooth) call do_gauss_filt(itof,i_wb,flds)
            if(ido .eq. iop_ballflt) call do_ball_filter(itof,i_wb,flds)
            if(ido .eq. iop_trsub) call do_subtrace(itof, i_wb, flds)
            if(ido .eq. iop_trace) call do_trace(itof, i_wb, flds)
            if(ido .eq. iop_det  ) call do_determinant(itof, i_wb, flds)
            if(ido .eq. iop_secsum) call do_secsum(itof, i_wb, flds)
            if(ido.ge.100 .and. ido.le.199)
     1              call do_function(ido, itof, i_wb, flds)      ! uniary fns
            if(ido.ge.200 .and. ido.le.299)
     1             call do_math(itof, ido, i_wb, flds)       ! binary fns

             if(idebug .gt. 10 .and. i_wb .eq. 1 .and. i_rb .eq. 1) then
               iwrk0 = iWrk_tof(itof)
               iwrk1 = iWrk_tof(itof) + nDim_tof(itof) - 1
               write (6,998) trim(cFld_tof(itof)),
     1                      (flds(1,1,1,j),j=iwrk0,iwrk1)
998            format(a20,"  Value:", 1p300e15.6)
             endif
          enddo  ! end iop loop
          if(idebug .eq. 100) stop
          t_work = t_work + (OMP_GET_WTIME() - t0)
 
          ! Feed output structure from results in work array
          t0 = OMP_GET_WTIME()
          if(nrep_used .eq. 1) then
            if(iout_xstrip .eq. 1) call fill_xstrip(profile,i_wb,flds)
            if(iout_ystrip .eq. 1) call fill_ystrip(profile,i_wb,flds)
          else 
            if(iout_xstrip+iout_ystrip .eq. 1) then
              call fill_avs(fxy,fyz,fxz, i_wb, flds)
            endif
          endif
          if(iout_strip  .eq. 1) call fill_strip(profile,i_wb,flds)
          if(iout_rprof  .eq. 1) call fill_rprof(profile,i_wb, flds)
          if(iout_xprof  .eq. 1) call fill_xprof(profile,i_wb, flds)
          if(iout_yprof  .eq. 1) call fill_yprof(profile,i_wb, flds)
          if(iout_zprof  .eq. 1) call fill_zprof(profile,i_wb, flds)
          if(iout_xsect  .eq. 1) call fill_xsect(profile,i_wb, flds)
          if(iout_ysect  .eq. 1) call fill_ysect(profile,i_wb, flds)
          if(iout_dist   .eq. 1) call fill_dist(profile,i_wb, flds)
          if(iout_corr   .eq. 1) call fill_corr(icorrelate,i_wb, flds)
          if(iout_view   .eq. 1) call fill_view(sym33, i_wb, flds)
          if(iout_xyslice .eq. 1) call fill_xyslice(bout,i_wb,flds)
          if(iout_hv    .eq. 1) call fill_hv(bout,i_wb, flds)
          if(iout_bof   .eq. 1) call fill_bof(fout,i_wb, flds)
          if(iout_spc3v .eq. 1) call fill_spc3v(ft3v,i_wb, flds)
          if(iout_spcnv .eq. 1) call fill_spcnv(ftnv,i_wb, flds)
          if(iout_strfn .eq. 1) call fill_strfn(fnv,i_wb,flds,idispsumS)
          if(iout_xyzavs.eq. 1) call fill_avs(fxy,fyz,fxz, i_wb, flds)
          if(iout_vr    .eq.  1) call fill_vr(fxy,fxz,fyz,i_wb,flds)
          t_out = t_out + (OMP_GET_WTIME() - t0)
        enddo   ! end I_wb loop

      enddo     ! end i_rb loop

c      write (6,*) ' '
c      do i = 1, nrda
c        write (6,*) 'rd field::', rdas(4,4,4,i)
c      enddo

      ! Finalize output
      t0 = OMP_GET_WTIME()
      if(iout_xstrip+iout_ystrip .ge. 1 .and. nrep_used .gt. 1) then
        call fxy_to_prof(fxy,profile)
      endif
      if(iout_xstrip .eq. 1) call write_strip(profile)
      if(iout_ystrip .eq. 1) call write_strip(profile)
      if(iout_strip  .eq. 1) call write_strip(profile)
      if(iout_xyslice .eq. 1) call write_xyslice(bout)
      if(iout_rprof .eq. 1) call write_prof(profile)
      if(iout_xprof .eq. 1) call write_prof(profile)
      if(iout_yprof .eq. 1) call write_prof(profile)
      if(iout_zprof .eq. 1) call write_prof(profile)
      if(iout_xyzavs.eq. 1) call write_avs(fxy,fyz,fxz)
      if(iout_vrz .eq. 1) call write_vr(fxy,nxout,nyout)
      if(iout_vry .eq. 1) call write_vr(fxz,nxout,nzout)
      if(iout_vrx .eq. 1) call write_vr(fyz,nyout,nzout)
      if(iout_ysect+iout_xsect .eq. 1) call write_sect(profile)
      if(iout_dist  .eq. 1) call write_dist(profile)
      if(iout_corr  .eq. 1) call write_corr(icorrelate)
      if(iout_view  .eq. 1) call write_view_from_sym33(sym33)
      t_out = t_out + (OMP_GET_WTIME() - t0)

      t_end = OMP_GET_WTIME()
      t_full = t_end - t_start
      if(idebug .gt. 0 .or. t_full .gt. 20.0) then
        write (6,997) t_read, t_work, t_out, t_full
997     format("#Wall Clock: read,work,out,full:", 4f10.3)
c        write (6,*) "t_r2w = ", t_r2w
      endif

      return
      end
 
C ====================================================================72
      subroutine gen_regions_vars()
      implicit none
      include 'postproc_info_e3d.h'
      real x0,x1,y0,y1,z0,z1, zmid, del
      integer i, ndvar
      character*32 cFile

      ! Generate list of variabies in compressed cump
      call get_vars_nams(MAX_TOF, ndvar, cFld_tof, cNam_tof)
      open(unit=21,file='cdvars.txt',form='formatted')
      do i = 1, ndvar
        write (21,980) trim(cFld_tof(i)), trim(cNam_tof(i))
980     format(a, " = ", a)
      enddo
      close(21)

      ! Generate base set of rgions of iPython Notebook GUI (ipnb)
      call get_full_dimensions(nxfull, nyfull, nzfull)
      call get_full_XYZrange(x0,x1,y0,y1,z0,z1)
      open(unit=21,file='full.rgn',form='formatted')
      write (21,990) cxcoord(1:1)
      write (21,990) cycoord(1:1)
      write (21,990) czcoord(1:1)
      close(21)
990   format("ee ", a1, "range")
991   format("ee ", a1, "range", 2f10.5)
      zmid = (z0 + z1) / 2.0
      del = (z1-z0)/float(nzfull)
      if(zmid .lt. 0.5*del) zmid = 0.0
      if(nzfull .gt. 1) then
        if(zmid .lt. 0.5*del) then
          cFile = "z=0.rgn                         "
        else
          cFile = "z_mid_plane.rng                 "
        endif
        cFile(1:1) = czcoord(1:1)
        open(unit=21,file=cFile,form='formatted')
        write (21,990) cxcoord(1:1)
        write (21,990) cycoord(1:1)
        write (21,991) czcoord(1:1), zmid, zmid
        close(21)
      else
        z0 = zmid
        z1 = zmid
      endif

      open(unit=21,file="full_dims.txt",form='formatted')
      write (21,991) cxcoord(1:1), x0, x1
      write (21,991) cycoord(1:1), y0, y1
      write (21,991) czcoord(1:1), z0, z1
      close(21)

      return
      end

C ====================================================================72
      subroutine print_wbxyz(i_wb)
      implicit none
      include 'postproc_info_e3d.h'
      integer i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1

      ix0 = ioff_wb(1,i_wb) + 1
      ix1 = ioff_wb(1,i_wb) + n0_wb(1,i_wb)
      iy0 = ioff_wb(2,i_wb) + 1
      iy1 = ioff_wb(2,i_wb) + n0_wb(2,i_wb)
      iz0 = ioff_wb(3,i_wb) + 1
      iz1 = ioff_wb(3,i_wb) + n0_wb(3,i_wb)

      write (6,999) ix0,ix1, iy0,iy1, iz0,iz1
999   format("Work Brick  x:", 2i5, "     y:", 2i5,"    z:", 2i5)

      return
      end

C ====================================================================72
      subroutine read_fld_from_a_box(irda, i_rb, fld, cVar)
      implicit none
      include 'postproc_info_e3d.h'
      include 'amr_e3d.h'
      real fld(nrdx,nrdy,nrdz)
      integer irda, i_rb, ibox
      character*256 cVar
      integer iz, nz

      ! Since this will read ONLY from one box, which will alwyas have
      ! boundaries deep enough, read full extent of box in X and Y, and
      ! specify offset (in Z only) relative to lower bound of the box.
      ! As alwyas, ALL mesh sizes & offsets are in units of the box.
      ! Note: nrdx and nrdy should be full box widths, but nrdz is JUST
      ! the Z mesh size of what will be read in.

      ibox = ibox_list(i_rb)
      iz = ioff_rb(3,i_rb) + nbdryz - iboxoff(3,ibox)
      if(ibox .ne. icurrent_box) then
        dx = boxsize(1) / nboxmesh(1)
        dy = boxsize(2) / nboxmesh(2)
        dz = boxsize(3) / nboxmesh(3)
        dx2i = 1.0 / (2.0 * dx)
        dy2i = 1.0 / (2.0 * dy)
        dz2i = 1.0 / (2.0 * dz)
      endif
      icurrent_box = ibox

      if(idebug .eq. 21 .and. irda .eq. 1) then
        write (6,999) ioff_rb(3,i_rb)+1, ioff_rb(3,i_rb)+nrdz,
     1                ibox, iz-nbdryz, iz+nrdz-(1+nbdryz)
999     format("   Reading Z-subrange: ", 2i6,
     1         "      from ibox,box_iz[12]:", 3i6)
      endif

c      write (6,*) "read_fld_from_one_box yet to be implemented"
      call read_fld_from_one_box(ibox,cVar,iz,nrdx,nrdy,nrdz,fld)

      return
      end
C ====================================================================72
      subroutine read_fld(irda, i_rb, fld)
      implicit none
      include 'postproc_info_e3d.h'
      real fld(nrdx,nrdy,nrdz)
      integer irda, i_rb
      character*256 cVar
      integer itof, i, ixoff, iyoff, izoff, nz, iz, izoffs(4),iat

      itof = itof_rda(irda)

      if(trim(cFile_tof(itof)) .eq. "no_file") then
        write (6,*)  " "
        write (6,*) trim(c_add_tof_problem)
        write (6,*)  " "
        write (6,*) "Check directory for needed data files."
        write (6,*) "You will probably need to generate"
        write (6,*) "or re-generate dump_times by running"
        write (6,*) "      ee times"
        stop
      endif

       if(cdumpfile0 .ne. cFile_tof(itof)) then   ! reset current dump file
          if(idebug .gt. 100) then
            write (6,*) "Switching dumps:"
            write (6,*) "  old dump:", trim(cdumpfile0)
            write (6,*) "  New dump:", trim(cFile_tof(itof))
          endif
          cdumpfile0 = cFile_tof(itof)
          call read_header_info(cDumpFile0)
        endif

      ! Get field string & trim off any @...
      do i = 1, 256
        cVar(i:i) = ' '
      enddo
      cVar(1:32) = cFld_tof(itof)
      iat = -1
      do i=32,1,-1
        if(cVar(i:i) .eq. '@') iat = i
      enddo
      if(iat .gt. 0) then
        do i=iat,32
          cVar(i:i) = ' '
        enddo
      endif

      if(nboxsused .gt. 0) then
        call read_fld_from_a_box(irda, i_rb, fld, cVar)
        return
      endif

      ixoff = ioff_rb(1,i_rb)
      iyoff = ioff_rb(2,i_rb)
      izoffs(1) =        ioff_rb(3,i_rb)
      izoffs(2) = max(0, ioff_rb(3,i_rb))
      izoffs(3) = min(nzfull, izoffs(1) + nf_rb(3,i_rb))
      izoffs(4) =             izoffs(1) + nf_rb(3,i_rb)

      if(nrdx .ne. nxfull .or. nrdy .ne. nyfull .or.
     1  ixoff .ne. 0      .or. iyoff .ne. 0          ) then
         write (6,*) "Read Field case not currently supported"
         write (6,*) "  Reads do not span full X-Y plane in dump"
         stop
      endif

      do i = 1, 3
        nz    = izoffs(i+1) - izoffs(i)
        izoff = izoffs(i)
        iz    = 1 + izoff - ioff_rb(3,i_rb)
        if(nz .gt. 0) then
          if(iboundary(3) .eq. 0) then
            ! Periodic in Z for reading data in
            if(izoff .lt.      0) izoff = izoff + nzfull
            if(izoff .ge. nzfull) izoff = izoff - nzfull
          else
            ! Continuation in Z for reading data in
            if((izoff.lt.0 .or. izoff.ge.nzfull) .and. idebug.eq.31)
     1         write (6,*) "continuation & izoff = ", izoff
            if(izoff .lt.      0) izoff = 0
            if(izoff .ge. nzfull) izoff = nzfull-nz
          endif

          if(idebug .ge. 31 .and. irda .eq. 1) then
            write (6,*) "------------------------------------"
            write (6,999) "X", ixoff+1,ixoff+nrdx
            write (6,999) "Y", iyoff+1,iyoff+nrdy
            write (6,999) "Z", izoff+1,izoff+nz
999         format("   Reading ",a1,"-subrange: ", 2i5)
            write (6,*) "nrdx,nrdy,nrdz = ", nrdx,nrdy,nrdz
            write (6,*) "calling read_fields"
          endif

          call read_fields(cVar,1,nrdx,nrdy,nz,
     1                     ixoff,iyoff,izoff, fld(1,1,iz))

          if(idebug .ge. 31 .and. irda .eq. 1) then
            write (6,*) "returned from read_fields"
            write (6,*) "------------------------------------"
          endif
        endif
      enddo

      if(idebug .ge. 30) then
         write (6,*) "Read: ",trim(cVar),"  1st value: ",fld(1,1,1)
      endif

      return
      end

C ====================================================================72
      subroutine gen_cp_list(i_rb,i_wb)
      implicit none
      include 'postproc_info_e3d.h'
      integer itof, i_rb, i_wb

      integer irda, iwrk, j, icp_rb0(3), icp_rb1(3)
      integer ioff_wr(3)  ! add this to convert work index to read index

      integer i, imid
      integer i0, nfull(3), ir(3), irx, iry, irz
      integer irange0(3,4)

      integer IRB0, IRB1, IAVAIL0_WB(3), IAVAIL1_WB(3)
      integer ICP0_FULL_WB(3), ICP1_FULL_WB(3), ICP0_FULL_RB(3)
      integer ICP1_FULL_RB(3)
   
      ! Generate list of 3D block copies to fill work array from read array
      ! 3D blocks will be in work array index space

      ! 3D Index OFFset from Work to Read indexs
      ! add  ioff_wr(j) to convert Work array index to corresponding Read array index
      do j = 1, 3
        ioff_wr(j) = ioff_wb(j,i_wb) - ioff_rb(j,i_rb)
      enddo

      ! Available range in read brick data (converted to work brick indexes)
      do j = 1, 3
       irb0 = 1
       irb1 = nf_rb(j,n_rb)
       iavail0_wb(j) = irb0 - ioff_wr(j)
       iavail1_wb(j) = irb1 - ioff_wr(j)
      enddo

      ! Full range of this copy in work brick and in read brick indexes
      do j = 1, 3
        icp0_full_wb(j) = 1             - nbdy_max(j)
        icp1_full_wb(j) = n0_wb(j,i_wb) + nbdy_max(j)

        icp0_full_rb(j) = icp0_full_wb(j) + ioff_wr(j)
        icp1_full_rb(j) = icp1_full_wb(j) + ioff_wr(j)
      enddo

      ! 1st 3D block copy is always the overlap region (ubiquitous & simple)
      ! copy do loops will be indexed in work brick space
      n_cp = 1
      do j = 1, 3
        iwb0_cp(j,1) = max(iavail0_wb(j), icp0_full_wb(j))
        iwb1_cp(j,1) = min(iavail1_wb(j), icp1_full_wb(j))

        irb0_cp(j,1) = iwb0_cp(j,n_cp) + ioff_wr(j)
        idel_cp(j,1) = 1
      enddo

      ! Find 3 ranges in each of x y and z directions
      do j = 1, 3
        irange0(j,1) = icp0_full_wb(j)
        irange0(j,2) = iwb0_cp(j,1)
        irange0(j,3) = iwb1_cp(j,1)    + 1
        irange0(j,4) = icp1_full_wb(j) + 1
      enddo 

      ! Impose periodic for now 
      do j = 1, 3
        if(isRZ .eq. 1) iboundary(j) = 1     ! 1=continuation
      enddo
      nfull(1) = nxfull
      nfull(2) = nyfull
      nfull(3) = nzfull

      ! enumerate every combination of ranges in X Y and Z
      do irz = 1, 3
      do iry = 1, 3
      do irx = 1, 3
        imid = 0
        if(irx .eq. 2) imid = imid + 1
        if(iry .eq. 2) imid = imid + 1
        if(irz .eq. 2) imid = imid + 1
c        write (6,998) irx,iry,irz,irange0(1,irx+1)-irange0(1,irx),
c     1                            irange0(2,iry+1)-irange0(2,iry),
c     1                            irange0(2,irz+1)-irange0(2,irz), imid
998     format("irxyz=", 3i2,"   nxyz=",3i4,"   imid=",i2)
        if(irange0(1,irx) .lt. irange0(1,irx+1) .and.
     1     irange0(2,iry) .lt. irange0(2,iry+1) .and.
     1     irange0(3,irz) .lt. irange0(3,irz+1) .and.
     1     imid .lt. 3                               ) then
          ! New non-empty region to fill
          ir(1) = irx
          ir(2) = iry
          ir(3) = irz
          n_cp = n_cp + 1
          do j = 1, 3
            iwb0_cp(j,n_cp) = irange0(j,ir(j)  ) 
            iwb1_cp(j,n_cp) = irange0(j,ir(j)+1) - 1

            ! Default (center overlap) case
            irb0_cp(j,n_cp) = iwb0_cp(j,n_cp) + ioff_wr(j)
            idel_cp(j,n_cp) = 1

            ! Edge cases
            if(ir(j) .ne. 2) then
              if(iboundary(j) .eq. 0) then
                ! In this direction: periodic AND boundareis not read in
                ! Can only assume that read brick spans this direction
                i0 = irb0_cp(j,n_cp)
                irb0_cp(j,n_cp) = 1 + mod(i0+nfull(j)-1,nfull(j))
              else
                ! Continuation
                idel_cp(j,n_cp) = 0
                if(ir(j) .eq. 1) then
                  irb0_cp(j,n_cp) = iwb1_cp(j,n_cp)+1 + ioff_wr(j)
                else
                  irb0_cp(j,n_cp) = iwb0_cp(j,n_cp)-1 + ioff_wr(j)
                endif
              endif
            endif
          enddo
        endif
      enddo
      enddo
      enddo

      if(idebug .gt. 100) then
      write (6,*) ' '
      write (6,*) 'i_wb, n_cp: ', i_wb, n_cp
      write (6,*) "Range       Xmin Xmax   Ymin Ymax   Zmin Zmax"
      write (6,999) "Interior",(ioff_wb(j,i_wb)+1,
     1                          ioff_wb(j,i_wb)+n0_wb(j,i_wb), j=1,3)
999   format(a8, "  ", 3("  ",2i5))

c      write (6,*) ' '
      write (6,999) "iavail  ",(iavail0_wb  (j),iavail1_wb(j),j=1,3)
      write (6,999) "cp_full ",(icp0_full_wb(j),icp1_full_wb(j),j=1,3)

      write (6,*) "Range       Xmin Xmax   Ymin Ymax   Zmin Zmax",
     1            "   RB0x RB0y RB0z   Delx Dely Delz"
      do i = 1, n_cp
        write (6,997) "3D-block", (iwb0_cp(j,i), iwb1_cp(j,i), j=1,3),
     1                (irb0_cp(j,i), j=1,3), (idel_cp(j,i), j=1,3)
997     format(a8, "  ", 3("  ",2i5), 2("  ", 3i5))
      enddo
      endif

      return
      end

C ====================================================================72
      subroutine r2w(itof,i_rb,i_wb,rdas,flds)
      implicit none
      include 'postproc_info_e3d.h'
      real rdas(nrdx,nrdy,nrdz,nrda)
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_rb, i_wb, irda, iwrk
      integer ixs, iys, izs, ixd, iyd, izd, i_cp
      
      irda = iRda_tof(itof)
      iwrk = iWrk_tof(itof)

      ! Loop over 3D block copies
      do i_cp = 1, n_cp

c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
      izs = irb0_cp(3,i_cp)
      do izd=iwb0_cp(3,i_cp), iwb1_cp(3,i_cp)
        iys = irb0_cp(2,i_cp)
        do iyd=iwb0_cp(2,i_cp), iwb1_cp(2,i_cp)
          ixs = irb0_cp(1,i_cp)
          do ixd=iwb0_cp(1,i_cp), iwb1_cp(1,i_cp)
            flds(ixd,iyd,izd,iwrk) = rdas(ixs,iys,izs,irda)
            ixs = ixs + idel_cp(1,i_cp)
          enddo
          iys = iys + idel_cp(2,i_cp)
        enddo
        izs = izs + idel_cp(3,i_cp)
      enddo

      enddo

      if(idebug .gt. 100) then
        i_cp = 1
        ixs = irb0_cp(1,i_cp)
        iys = irb0_cp(2,i_cp)
        ixd = iwb0_cp(1,i_cp)
        iyd = iwb0_cp(2,i_cp)
        write (6,*) ' '
        write (6,*) ' '
        write (6,*) 'irda,iwrk    = ', irda, iwrk
        write (6,*) 'ixs,ixd      = ', ixs, ixd
        write (6,*) 'iys,iyd      = ', iys, iyd
        write (6,*) 'idel_cp(3,1) = ', idel_cp(3,i_cp)
        write (6,*) ' '
        write (6,*) ' izs    rdas       izd    flds'
        izs = irb0_cp(3,i_cp)
        do izd=iwb0_cp(3,i_cp), iwb1_cp(3,i_cp)
          write (6,999) izs, rdas(ixs,iys,izs,irda),
     1                  izd, flds(ixd,iyd,izd,iwrk)
999       format(i5,f10.5,i10,f10.5)
          izs = izs + idel_cp(3,i_cp)
        enddo

        write (6,*) 'stopping'
        stop
      endif

      return
      end

C ====================================================================72
      subroutine do_function(ido, itof, i_wb,flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer ido, itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1, iWrk, iArg, ix,iy,iz
      real    fac, value, val

      iwrk = iWrk_tof(itof)
      iArg = iArg_tof(itof,1)

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      if(ido .eq. iop_log10) then
      fac = 1.0 / alog(10.0)
      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
        flds(ix,iy,iz,iwrk) = fac * alog(flds(ix,iy,iz,iArg))
      enddo
      enddo
      enddo
      endif
      if(ido .eq. iop_log) then
      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
        flds(ix,iy,iz,iwrk) = alog(flds(ix,iy,iz,iArg))
      enddo
      enddo
      enddo
      endif
      if(ido .eq. iop_sqrt) then
      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
        flds(ix,iy,iz,iwrk) = sqrt(flds(ix,iy,iz,iArg))
      enddo
      enddo
      enddo
      endif
      if(ido .eq. iop_step) then
      value = fVal_tof(itof)
      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
        if(flds(ix,iy,iz,iArg) .lt. value) then
          flds(ix,iy,iz,iwrk) = 0.0
        else
          flds(ix,iy,iz,iwrk) = 1.0
        endif
      enddo
      enddo
      enddo
      endif
      if(ido .eq. iop_abs) then
      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
          flds(ix,iy,iz,iwrk) = abs(flds(ix,iy,iz,iArg))
      enddo
      enddo
      enddo
      endif
      if(ido .eq. iop_sin) then
      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
        flds(ix,iy,iz,iwrk) = sin(flds(ix,iy,iz,iArg))
      enddo
      enddo
      enddo
      endif
      if(ido .eq. iop_cos) then
      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
        flds(ix,iy,iz,iwrk) = cos(flds(ix,iy,iz,iArg))
      enddo
      enddo
      enddo
      endif
      if(ido .eq. iop_exp) then
      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
        val = flds(ix,iy,iz,iArg)
c        val = min(88.0, max(-88.0, val))
        flds(ix,iy,iz,iwrk) = exp(val)
      enddo
      enddo
      enddo
      endif
      if(ido .eq. iop_const) then
      value = fVal_tof(itof)
      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
        flds(ix,iy,iz,iwrk) = value
      enddo
      enddo
      enddo
      endif

      return
      end

C ====================================================================72
      subroutine do_norm(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1, ndim1, iArg
      integer iWrk, iArg1,iArg2,iArg3, ix,iy,iz
      real    fac, s
      
      iwrk = iWrk_tof(itof)
      iArg1 = iArg_tof(itof,1)
      iArg2 = iArg1 + 1
      iArg3 = iArg1 + 2
      nDim1 = iArg_dim(itof,1)

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      if(ndim1 .eq. 3) then
        do iz = iz0, iz1
        do iy = iy0, iy1
        do ix = ix0, ix1
          flds(ix,iy,iz,iwrk) = sqrt(flds(ix,iy,iz,iArg1)**2 +
     1                               flds(ix,iy,iz,iArg2)**2 +
     1                               flds(ix,iy,iz,iArg3)**2  )
        enddo
        enddo
        enddo
      else
        do iz = iz0, iz1
        do iy = iy0, iy1
        do ix = ix0, ix1
          s = 0.0
          do iArg = iArg1, iArg1+ndim1-1
            s = s + flds(ix,iy,iz,iArg)**2
          enddo
          flds(ix,iy,iz,iwrk) = sqrt(s)
        enddo
        enddo
        enddo
      endif

      return
      end

C ====================================================================72
      subroutine do_getc(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1, ndim1,iarg,iArg1,iarg2,iarg3
      integer iWrk, ix,iy,iz, ivc
      real    fac, s
      integer lentof, mycomp
      
      iwrk = iWrk_tof(itof)
      nDim1 = iArg_dim(itof,1)

      ivc = -1
      lentof = len(trim(cOpr_tof(itof)))    ! getc
      if(cOpr_tof(itof)(5:6) .eq. '_v') then
        read (cOpr_tof(itof)(7:lentof),*) ivc
        iArg1 = iArg_tof(itof,1)
        iArg2 = iArg_tof(itof,1)+1
        iArg3 = iArg_tof(itof,1)+2
      else
        read (cOpr_tof(itof)(5:lentof),*) mycomp
        iArg = iArg_tof(itof,1) + mycomp
      endif

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      if(ivc .le. 0) then
        do iz = iz0, iz1
        do iy = iy0, iy1
        do ix = ix0, ix1
          flds(ix,iy,iz,iwrk) = flds(ix,iy,iz,iArg)
        enddo
        enddo
        enddo
      else
        do iz = iz0, iz1
        do iy = iy0, iy1
        do ix = ix0, ix1
          flds(ix,iy,iz,iwrk) = vview(1,ivc) * flds(ix,iy,iz,iArg1)
     1                        + vview(2,ivc) * flds(ix,iy,iz,iArg2)
     1                        + vview(3,ivc) * flds(ix,iy,iz,iArg3)
        enddo
        enddo
        enddo
      endif

      return
      end

C ====================================================================72
      subroutine do_radc(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1, ndim1, iArg
      integer iWrk, iArg1,iArg2,iArg3, ix,iy,iz
      real    xx,yy,zz,xyzi,xyzmid0
      
      iwrk = iWrk_tof(itof)
      iArg1 = iArg_tof(itof,1)
      iArg2 = iArg1 + 1
      iArg3 = iArg1 + 2
      nDim1 = iArg_dim(itof,1)

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)
      xyzmid0 = 0.5 * float(nxfull + 1)

      do iz = iz0, iz1
      zz = float(iz + ioff_wb(3,i_wb)) - xyzmid0
      do iy = iy0, iy1
      yy = float(iy + ioff_wb(2,i_wb)) - xyzmid0
      do ix = ix0, ix1
      xx = float(ix + ioff_wb(1,i_wb)) - xyzmid0
 
        xyzi = 1.0 / sqrt(xx**2 + yy**2 + zz**2)
        flds(ix,iy,iz,iwrk) = xyzi*(xx*flds(ix,iy,iz,iArg1) +
     1                              yy*flds(ix,iy,iz,iArg2) +
     1                              zz*flds(ix,iy,iz,iArg3)  )
      enddo
      enddo
      enddo

      return
      end

C ====================================================================72
      subroutine do_pullb(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1, ia20, i
      integer iWrk, ia00, ia11,ia12,ia13, ia21,ia22,ia23, ix,iy,iz
      real    fac, dt, x1(3),x2(3),v1(3),v2(3),dinv(3), q2

      ! debug
c      integer i0(3), istop
c      real frac(3,0:1)
c      common /frac_at/ i0, frac
      
      iwrk = iWrk_tof(itof)

      ia20 = iArg_tof(itof,1)

      ia11 = iArg_tof(itof,2)
      ia12 = ia11 + 1
      ia13 = ia11 + 2

      ia21 = iArg_tof(itof,3)
      ia22 = ia21 + 1
      ia23 = ia21 + 2

      dt = fVal_tof(itof)
      dinv(1) = 1.0 / dx
      dinv(2) = 1.0 / dy
      dinv(3) = 1.0 / dz

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

c      istop = 0
      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1

        ! get time averaged vel
        x1(1)  = dx * (float(ix) - 0.5)
        x1(2)  = dy * (float(iy) - 0.5)
        x1(3)  = dz * (float(iz) - 0.5)
        v1(1) = flds(ix,iy,iz,ia11)
        v1(2) = flds(ix,iy,iz,ia12)
        v1(3) = flds(ix,iy,iz,ia13)

        do i = 1, 3
          x2(i) = x1(i) + v1(i) * dt
        enddo
        call set_frac_at(x2,dinv)
        call get_val_at(flds,ia21,v2(1))
        call get_val_at(flds,ia22,v2(2))
        call get_val_at(flds,ia23,v2(3))
        do i = 1, 3
          x2(i) = x1(i) + 0.5*(v1(i)+v2(i)) * dt
        enddo
        call set_frac_at(x2,dinv)
        call get_val_at(flds,ia20,q2)
        flds(ix,iy,iz,iwrk) = q2

c      if(iy00 .eq. iy + ioff_wb(2,i_wb) .and.
c     1   iz00 .eq. iz + ioff_wb(3,i_wb)        ) then
cc          write (6,*) 'ix,value = ', ix, q2
c          write (6,999) ix, (i0(i),i=1,3), (x2(i),i=1,3), q2,
c     1                      (frac(i,0),i=1,3)
c999      format("ix=",i2,"   i0=",3i3,"  x2=",3f10.5,"  val=",f10.5,
c     1          "   frac=",3f10.5)
c         istop = 1
c      endif

      enddo
      enddo
      enddo
c      if(istop .eq. 1) stop

      return
      end
C ====================================================================72
      subroutine set_frac_at(x,dinv)
      implicit none
      real x(3), dinv(3)
      real xdinv
      integer i
      integer i0(3)
      real frac(3,0:1)
      common /frac_at/ i0, frac

c      integer idebug
c      idebug = 10

      do i = 1, 3
        xdinv = x(i) * dinv(i)
        i0(i) = int(99.5 + xdinv) - 99
        frac(i,1) = xdinv + 0.5 - float(i0(i))
        frac(i,0) = 1.0 - frac(i,1)

c         if(idebug .gt. 1000) then
c           write (6,999) i,i0(i),x(i),dinv(i),frac(i,0)
c999        format("set_frac: i,i0,x,dinv,frac(i,0):",2i4,3f10.5)
c        endif
      enddo
c       if(idebug .gt. 1000) then
c         write (6,999) (i0(i),i=1,3),(x(i),i=1,3)
c999      format("set_frac: i0=",3i3,"   x=",3f10.5)
c      endif

      return
      end

C ====================================================================72
      real function cubic_interp(v,x)
      implicit none
      real v(-1:2), x, s0, s1, d0, c0, dv

      s0 = 0.5 * (v(1) - v(-1)) 
      s1 = 0.5 * (v(2) - v( 0)) 
      ! v(x) = v(0) + s0*x + c0*(x**2) + d0*x**3)
      ! v(1) = v(0) + s0 +   c0 +   d0
      ! s1   =        s0 + 2*c0 + 3*d0
      ! s1-2*v(1) = -2*v(0) - s0 + d0
      dv = v(1) - v(0)
      d0 = s0 + s1 - 2.0*dv
      c0 = dv - s0 - d0
      cubic_interp = v(0) + x*(s0 + x*(c0 + x*d0))

      return
      end

C ====================================================================72
c tri-cubic version
      subroutine get_val_at(flds,ia,val)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer ia
      real val
      integer i0(3)
      real frac(3,0:1)
      common /frac_at/ i0, frac
      integer ix,iy,iz, ix0, iy0, iz0
      real fracx, fracy, fracz
      real vyz(-1:2, -1:2), vz(-1:2), cubic_interp

      ix0 = i0(1)
      iy0 = i0(2)
      iz0 = i0(3)
      fracx = frac(1,1)
      fracy = frac(2,1)
      fracz = frac(3,1)

      do iz = -1, 2
      do iy = -1, 2
        vyz(iy,iz) = cubic_interp(flds(ix0-1,iy0+iy,iz0+iz,ia),fracx)
      enddo
      enddo

      do iz = -1, 2
        vz(iz) = cubic_interp(vyz(-1,iz),fracy)
      enddo

      val = cubic_interp(vz(-1),fracz)

      return
      end

C ====================================================================72
c Tri-linear version:
      subroutine get_val_at_tl(flds,ia,val)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer ia
      real val
      integer i0(3)
      real frac(3,0:1)
      common /frac_at/ i0, frac
      integer ix,iy,iz, iix, iiy, iiz
      real fracx, fracy, fracz

      val = 0.0

      do iix = 0, 1
      fracx = frac(1,iix)
      ix = i0(1) + iix

      do iiy = 0, 1
      fracy = frac(2,iiy)
      iy = i0(2) + iiy

      do iiz = 0, 1
      fracz = frac(3,iiz)
      iz = i0(3) + iiz

        if(idebug .gt. 1000) then
          write (6,888) ix,iy,iz,fracx,fracy,fracz,flds(ix,iy,iz,ia)
888       format("iXYZ:",3i4,"    fracXYZ:",3f10.5,"  flds:",f10.5)
        endif
        val = val + fracx*fracy*fracz*flds(ix,iy,iz,ia)

      enddo
      enddo
      enddo

        if(idebug .gt. 1000) then
         write (6,*) "val = ", val
         stop
        endif

      return
      end

C ====================================================================72
      subroutine gen_v_perp(ix,iy,iz,vx,flds,g,v)
      implicit none
      include 'postproc_info_e3d.h'
      integer ix,iy,iz,vx
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      real g(3), v(3), gdv
      v(1) = flds(ix,iy,iz,vx  )
      v(2) = flds(ix,iy,iz,vx+1)
      v(3) = flds(ix,iy,iz,vx+2)
      gdv = g(1)*v(1) + g(2)*v(2) + g(3)*v(3)
      v(1) = v(1) - gdv * g(1)
      v(2) = v(2) - gdv * g(2)
      v(3) = v(3) - gdv * g(3)
      return
      end

C ====================================================================72
      subroutine gen_v_minus(ix,i,iy,j,iz,k,vx,flds,g,v,vofs)
      implicit none
      include 'postproc_info_e3d.h'
      integer ix,iy,iz,vx, i,j,k, is0, is1
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      real g(3), v(3), s, sc, vofs(3), w0, w1, vv
      s = float(i)*g(1) + float(j)*g(2)+float(k)*g(3)
      is0 = int(2.5 + s)
      sc = float(is0-2)
      is1 = is0 + 1
      if(s .lt. sc .or. is0 .eq. 3) is1 = is0 - 1
      w0 = 1.0 - abs(s-sc)
      w1 = 1.0 - w0
      if(is0 .lt. 1  .or.  is1 .lt. 1) then
        write (6,*) "is0,is1", is0,is1
        write (6,*) "s,sc", s,sc
        write (6,*) "i,j,k", i,j,k
      endif
      vv = w0*vofs(is0) + w1*vofs(is1)
      
      v(1) = flds(ix+i,iy+j,iz+k,vx  ) - vv * g(1)
      v(2) = flds(ix+i,iy+j,iz+k,vx+1) - vv * g(2)
      v(3) = flds(ix+i,iy+j,iz+k,vx+2) - vv * g(3)
c      v(1) = flds(ix+i,iy+j,iz+k,vx  )
c      v(2) = flds(ix+i,iy+j,iz+k,vx+1)
c      v(3) = flds(ix+i,iy+j,iz+k,vx+2)
      return
      end
C ====================================================================72
      subroutine cross(a,b,c)
      implicit none
      real a(3), b(3), c(3)
      a(1) = b(2)*c(3) - b(3)*c(2)
      a(2) = b(3)*c(1) - b(1)*c(3)
      a(3) = b(1)*c(2) - b(2)*c(1)
      return
      end
C ====================================================================72
      subroutine vecsub(a,b,c)
      implicit none
      real a(3), b(3), c(3)
      a(1) = b(1) - c(1)
      a(2) = b(2) - c(2)
      a(3) = b(3) - c(3)
      return
      end
C ====================================================================72
      subroutine ortho(g,v)
      implicit none
      real g(3), v(3), gdv
      gdv = g(1)*v(1) + g(2)*v(2) + g(3)*v(3)
      v(1) = v(1) - gdv * g(1)
      v(2) = v(2) - gdv * g(2)
      v(3) = v(3) - gdv * g(3)
      return
      end
C ====================================================================72
      real function vecnorm(v)
      implicit none
      real v(3)
      vecnorm = sqrt(v(1)**2 + v(2)**2 + v(3)**2)
      return
      end
C ====================================================================72
      subroutine normalize(v)
      implicit none
      real v(3), vinv
      vinv = 1.0 / sqrt(v(1)**2 + v(2)**2 + v(3)**2)
      v(1) = vinv * v(1)
      v(2) = vinv * v(2)
      v(3) = vinv * v(3)
      return
      end
C ====================================================================72
      subroutine gen_vofs(s,v,vofs)
      implicit none
      real s(27), v(27), vofs(3)
      integer is0, i, n(3)
      do i=1,3
        n(i) = 0
        vofs(i) = 0.0
      enddo
      do i=1,27
        is0 = int(2.5 + s(i))
        if(is0 .ge. 1 .and. is0 .le. 3) then
          n(is0) = n(is0) + 1
          vofs(is0) = vofs(is0) + v(i)
        endif
      enddo
      do i=1,3
        if(n(i) .eq. 0) vofs(i) = 0.0
        if(n(i) .gt. 0) vofs(i) = vofs(i) / float(n(i))
      enddo

      return
      end

C ====================================================================72
      subroutine do_curld3(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1
      integer iWrk1,iWrk2,iWrk3, ro, vx,vy,vz, ix,iy,iz
      real wx, wy, wz, divv, g(3)
      real vxp(3),vxm(3),vyp(3),vym(3),vzp(3),vzm(3), w(3), gdv
      integer i,j,k, ijk
      real s(27), v(27), vofs(3)
      
      iwrk1 = iWrk_tof(itof)
      iwrk2 = iWrk1 + 1
      iwrk3 = iWrk1 + 2
      ro = iArg_tof(itof,1)
      vx = ro + 1
      vy = ro + 2
      vz = ro + 3

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
        ! not in a shock
        w(1) =  dy2i*(flds(ix  ,iy+1,iz  ,vz)-flds(ix  ,iy-1,iz  ,vz))
     1        - dz2i*(flds(ix  ,iy  ,iz+1,vy)-flds(ix  ,iy  ,iz-1,vy))

        w(2) =  dz2i*(flds(ix  ,iy  ,iz+1,vx)-flds(ix  ,iy  ,iz-1,vx))
     1        - dx2i*(flds(ix+1,iy  ,iz  ,vz)-flds(ix-1,iy  ,iz  ,vz))

        w(3) =  dx2i*(flds(ix+1,iy  ,iz  ,vy)-flds(ix-1,iy  ,iz  ,vy))
     1        - dy2i*(flds(ix  ,iy+1,iz  ,vx)-flds(ix  ,iy-1,iz  ,vx))

        divv = dx2i*(flds(ix+1,iy  ,iz  ,vx)-flds(ix-1,iy  ,iz  ,vx))
     1       + dy2i*(flds(ix  ,iy+1,iz  ,vy)-flds(ix  ,iy-1,iz  ,vy))
     1       + dz2i*(flds(ix  ,iy  ,iz+1,vz)-flds(ix  ,iy  ,iz-1,vz))
        if(divv .lt. threshold) then
          w(1) = 0.0
          w(2) = 0.0
          w(3) = 0.0
        endif

        flds(ix,iy,iz,iwrk1) = w(1)
        flds(ix,iy,iz,iwrk2) = w(2)
        flds(ix,iy,iz,iwrk3) = w(3)
      enddo
      enddo
      enddo

      return
      end

C ====================================================================72
      real function vlav(ix,iy,iz,idir,vv,ro)
      implicit none
      include 'postproc_info_e3d.h'
      integer ix,iy,iz,idir, ix1,iy1,iz1
      real vv(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi)
      real ro(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi)
      real v(-1:1), m(-1:1), sm, sv

      if(idir .eq. 1) then
        do ix1 = ix-1, ix+1
c         v(ix1-ix) = 0.0
          v(ix1-ix) = 0.25*(vv(ix1,iy  ,iz  ) + vv(ix1,iy-1,iz  )
     1              +       vv(ix1,iy  ,iz-1) + vv(ix1,iy-1,iz-1))
          m(ix1-ix) = 0.25*(ro(ix1,iy  ,iz  ) + ro(ix1,iy-1,iz  )
     1              +       ro(ix1,iy  ,iz-1) + ro(ix1,iy-1,iz-1))
        enddo
      endif
      if(idir .eq. 2) then
        do iy1 = iy-1, iy+1
c          v(iy1-iy) = 0.0
          v(iy1-iy) = 0.25*(vv(ix  ,iy1 ,iz  ) + vv(ix-1,iy1 ,iz  )
     1              +       vv(ix  ,iy1 ,iz-1) + vv(ix-1,iy1 ,iz-1))
          m(iy1-iy) = 0.25*(ro(ix  ,iy1 ,iz  ) + ro(ix-1,iy1 ,iz  )
     1              +       ro(ix  ,iy1 ,iz-1) + ro(ix-1,iy1 ,iz-1))
        enddo
      endif
      if(idir .eq. 3) then
        do iz1 = iz-1, iz+1
c          v(iz1-iz) = 1.0
          v(iz1-iz) = 0.25*(vv(ix  ,iy  ,iz1) + vv(ix  ,iy-1,iz1)
     1              +       vv(ix-1,iy  ,iz1) + vv(ix-1,iy-1,iz1))
          m(iz1-iz) = 0.25*(ro(ix  ,iy  ,iz1) + ro(ix  ,iy-1,iz1)
     1              +       ro(ix-1,iy  ,iz1) + ro(ix-1,iy-1,iz1))
        enddo
      endif

      ! momentum(i) = m(i) * v(i) = integral{m(x)*v(x)}
      ! in zone0: m(x) = m(0) + sm*x    -0.5 <= x <= 0.5
      ! in zone0: v(x) = v0   + sv*x    -0.5 <= x <= 0.5
      ! m(0)*v(0) = m(0)*v0 + sm*sv/4
      !      v(0) =      v0 + sm*sv/(4.0*m(0))
c      sm = 0.5 * (m(1) - m(-1))     ! dx=1 here
c      sv = 0.5 * (v(1) - v(-1))     ! dx=1 here
c      vlav = v(0) - sm*sv/(4.0*m(0))
c      vlav = v(0)
       vlav = 0.25*(v(-1)+v(1)) + 0.5*v(0)

      return
      end

C ====================================================================72
      subroutine do_curld2(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1
      integer iWrk1,iWrk2,iWrk3, ro, vx,vy,vz, ix,iy,iz
      real wx, wy, wz, divv, g(3)
      real vxp(3),vxm(3),vyp(3),vym(3),vzp(3),vzm(3), w(3), gdv
      integer i,j,k, ijk
      real s(27), v(27), vofs(3)
      
      iwrk1 = iWrk_tof(itof)
      iwrk2 = iWrk1 + 1
      iwrk3 = iWrk1 + 2
      ro = iArg_tof(itof,1)
      vx = ro + 1
      vy = ro + 2
      vz = ro + 3

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1

        divv = dx2i*(flds(ix+1,iy  ,iz  ,vx)-flds(ix-1,iy  ,iz  ,vx))
     1       + dy2i*(flds(ix  ,iy+1,iz  ,vy)-flds(ix  ,iy-1,iz  ,vy))
     1       + dz2i*(flds(ix  ,iy  ,iz+1,vz)-flds(ix  ,iy  ,iz-1,vz))

          ! not in a shock
        w(1) =  dy2i*(flds(ix  ,iy+1,iz  ,vz)-flds(ix  ,iy-1,iz  ,vz))
     1        - dz2i*(flds(ix  ,iy  ,iz+1,vy)-flds(ix  ,iy  ,iz-1,vy))

        w(2) =  dz2i*(flds(ix  ,iy  ,iz+1,vx)-flds(ix  ,iy  ,iz-1,vx))
     1        - dx2i*(flds(ix+1,iy  ,iz  ,vz)-flds(ix-1,iy  ,iz  ,vz))

        w(3) =  dx2i*(flds(ix+1,iy  ,iz  ,vy)-flds(ix-1,iy  ,iz  ,vy))
     1        - dy2i*(flds(ix  ,iy+1,iz  ,vx)-flds(ix  ,iy-1,iz  ,vx))

        if(divv .lt. threshold) then
          ! in a shock: porject out component of v in shock direction
          g(1) = dx2i*(flds(ix+1,iy  ,iz  ,ro)-flds(ix-1,iy  ,iz  ,ro))
          g(2) = dy2i*(flds(ix  ,iy+1,iz  ,ro)-flds(ix  ,iy-1,iz  ,ro))
          g(3) = dz2i*(flds(ix  ,iy  ,iz+1,ro)-flds(ix  ,iy  ,iz-1,ro))
          call normalize(g)

           ijk = 0
           do k=-1,1
           do j=-1,1
           do i=-1,1
             ijk = ijk + 1
             s(ijk) = float(i)*g(1) + float(j)*g(2)+float(k)*g(3)
             v(ijk) = g(1)*flds(ix+i,iy+j,iz+k,vx)
     1              + g(2)*flds(ix+i,iy+j,iz+k,vy)
     1              + g(3)*flds(ix+i,iy+j,iz+k,vz)
cccccccccccccccccccccccccccccccccccccc
c             if(divv .lt. -15.0) then
c                write (6,*) flds(ix+i,iy+j,iz+k,vx),
c     1                      flds(ix+i,iy+j,iz+k,vy),
c     1                      flds(ix+i,iy+j,iz+k,vz)
c             endif
cccccccccccccccccccccccccccccccccccccc
           enddo
           enddo
           enddo
           call gen_vofs(s,v,vofs)

cccccccccccccccccccccccccccccccccccccc
c           if(divv .lt. -15.0) then
c            write (6,*) "   "
c            write (6,*) "divv:", divv
c            write (6,*) "g   :", (g(i),i=1,3)
c            write (6,*) "   "
c
c           do i=1,27
c             write (6,*) "s,v", s(i), v(i)
c            enddo
c            write (6,*) "   "
c            write (6,*) "   "
c            do i = 1,3
c              write (6,*) "i,vofs(i)", i, vofs(i)
c            enddo
c            stop
c          endif
cccccccccccccccccccccccccccccccccccccc

          call gen_v_minus(ix,+1,iy, 0,iz, 0,vx,flds,g,vxp,vofs)
          call gen_v_minus(ix,-1,iy, 0,iz, 0,vx,flds,g,vxm,vofs)
          call gen_v_minus(ix, 0,iy,+1,iz, 0,vx,flds,g,vyp,vofs)
          call gen_v_minus(ix, 0,iy,-1,iz, 0,vx,flds,g,vym,vofs)
          call gen_v_minus(ix, 0,iy, 0,iz,+1,vx,flds,g,vzp,vofs)
          call gen_v_minus(ix, 0,iy, 0,iz,-1,vx,flds,g,vzm,vofs)

c           if(divv .lt. -15.0) then
c            write (6,*) "vyp(3),flds", vyp(3),flds(ix,iy+1,iz,vz)
c            write (6,*) "vym(3),flds", vym(3),flds(ix,iy-1,iz,vz)
c            write (6,*) "vzp(2),flds", vzp(2),flds(ix,iy,iz+1,vy)
c            write (6,*) "vzm(2),flds", vzm(2),flds(ix,iy,iz-1,vy)
c
c            write (6,*) "  "
c            write (6,*) "  "
c            write (6,*) "OLD w:", (w(i),i=1,3)
c          endif

          w(1) =  dy2i*(vyp(3)-vym(3)) - dz2i*(vzp(2)-vzm(2))
          w(2) =  dz2i*(vzp(1)-vzm(1)) - dx2i*(vxp(3)-vxm(3))
          w(3) =  dx2i*(vxp(2)-vxm(2)) - dy2i*(vyp(1)-vym(1))

c           if(divv .lt. -15.0) then
c            write (6,*) "NEW w:", (w(i),i=1,3)
c            stop
c          endif

c          gdv = g(1)*w(1) + g(2)*w(2) + g(3)*w(3)
c          w(1) = gdv * g(1)
c          w(2) = gdv * g(2)
c          w(3) = gdv * g(3)
        endif

        flds(ix,iy,iz,iwrk1) = w(1)
        flds(ix,iy,iz,iwrk2) = w(2)
        flds(ix,iy,iz,iwrk3) = w(3)
      enddo
      enddo
      enddo

      return
      end

C ====================================================================72
      subroutine do_curlf(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1
      integer iWrk1,iWrk2,iWrk3, ro,vx,vy,vz, ix,iy,iz
      real vxL, vyL, vzL, wxL, wyL, wzL, fac, vlav
      allocatable vxL(:,:,:), vyL(:,:,:), vzL(:,:,:)
      allocatable wxL(:,:,:), wyL(:,:,:), wzL(:,:,:)

      iwrk1 = iWrk_tof(itof)
      iwrk2 = iWrk1 + 1
      iwrk3 = iWrk1 + 2
      ro = iArg_tof(itof,1)
      vx = ro + 1
      vy = ro + 2
      vz = ro + 3
      
      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      allocate (vxL(ix0:ix1+1, iy0:iy1+1, iz0:iz1+1))
      allocate (vyL(ix0:ix1+1, iy0:iy1+1, iz0:iz1+1))
      allocate (vzL(ix0:ix1+1, iy0:iy1+1, iz0:iz1+1))
      allocate (wxL(ix0:ix1+1, iy0:iy1+1, iz0:iz1+1))
      allocate (wyL(ix0:ix1+1, iy0:iy1+1, iz0:iz1+1))
      allocate (wzL(ix0:ix1+1, iy0:iy1+1, iz0:iz1+1))

      ! Edge centered velocities based on all sourounding cells
      do iz = iz0, iz1+1
      do iy = iy0, iy1+1
      do ix = ix0, ix1
        vxL(ix,iy,iz) = vlav(ix,iy,iz,1,flds(iwxlo,iwylo,iwzlo,vx),
     1                                  flds(iwxlo,iwylo,iwzlo,ro))
      enddo
      enddo
      enddo

      do iz = iz0, iz1+1
      do iy = iy0, iy1
      do ix = ix0, ix1+1
        vyL(ix,iy,iz) = vlav(ix,iy,iz,2,flds(iwxlo,iwylo,iwzlo,vy),
     1                                  flds(iwxlo,iwylo,iwzlo,ro))
      enddo
      enddo
      enddo

      do iz = iz0, iz1
      do iy = iy0, iy1+1
      do ix = ix0, ix1+1
        vzL(ix,iy,iz) = vlav(ix,iy,iz,3,flds(iwxlo,iwylo,iwzlo,vz),
     1                                  flds(iwxlo,iwylo,iwzlo,ro))
      enddo
      enddo
      enddo
      ! Face centered circulation (*1/dx)
      ! Assumes dx=dy=dz=const
      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1+1
        wxL(ix,iy,iz) = vyL(ix  ,iy  ,iz  ) + vzL(ix  ,iy+1,iz  )
     1                - vyL(ix  ,iy  ,iz+1) - vzL(ix  ,iy  ,iz  )
      enddo
      enddo
      enddo
      do iz = iz0, iz1
      do iy = iy0, iy1+1
      do ix = ix0, ix1
        wyL(ix,iy,iz) = vxL(ix  ,iy  ,iz  ) + vzL(ix+1,iy  ,iz  )
     1                - vxL(ix  ,iy  ,iz+1) - vzL(ix  ,iy  ,iz  )
      enddo
      enddo
      enddo
      do iz = iz0, iz1+1
      do iy = iy0, iy1
      do ix = ix0, ix1
        wzL(ix,iy,iz) = vxL(ix  ,iy  ,iz  ) + vyL(ix+1,iy  ,iz  )
     1                - vxL(ix  ,iy+1,iz  ) - vyL(ix  ,iy  ,iz  )
      enddo
      enddo
      enddo

      ! Face centered vorticity averaged back to cell center
      ! Sumes, above are (1/dx)*cericulation around each face
      ! Two of these face centered values are used here & need
      ! to divide by face area dx**2)
      fac = 1.0/(2.0*dx)

      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
        flds(ix,iy,iz,iwrk1) = fac * (wxL(ix,iy,iz)+wxl(ix+1,iy,iz))
        flds(ix,iy,iz,iwrk2) = fac * (wyL(ix,iy,iz)+wyl(ix,iy+1,iz))
        flds(ix,iy,iz,iwrk3) = fac * (wzL(ix,iy,iz)+wzl(ix,iy,iz+1))
      enddo
      enddo
      enddo

      deallocate (vxL)
      deallocate (vyL)
      deallocate (vzL)
      deallocate (wxL)
      deallocate (wyL)
      deallocate (wzL)

      return
      end

C ====================================================================72
      subroutine do_curl(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1
      integer iWrk1,iWrk2,iWrk3, vx,vy,vz, ix,iy,iz
      real vxL, vyL, vzL, wxL, wyL, wzL, fac
      allocatable vxL(:,:,:), vyL(:,:,:), vzL(:,:,:)
      allocatable wxL(:,:,:), wyL(:,:,:), wzL(:,:,:)
      
      iwrk1 = iWrk_tof(itof)
      iwrk2 = iWrk1 + 1
      iwrk3 = iWrk1 + 2
      vx = iArg_tof(itof,1)
      vy = vx + 1
      vz = vx + 2

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      if(cCurlCalc(1:4) .eq. "cell") then
      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
        
        flds(ix,iy,iz,iwrk1) =
     1          dy2i*(flds(ix  ,iy+1,iz  ,vz)-flds(ix  ,iy-1,iz  ,vz))
     1        - dz2i*(flds(ix  ,iy  ,iz+1,vy)-flds(ix  ,iy  ,iz-1,vy))

        flds(ix,iy,iz,iwrk2) =
     1          dz2i*(flds(ix  ,iy  ,iz+1,vx)-flds(ix  ,iy  ,iz-1,vx))
     1        - dx2i*(flds(ix+1,iy  ,iz  ,vz)-flds(ix-1,iy  ,iz  ,vz))

        flds(ix,iy,iz,iwrk3) =
     1          dx2i*(flds(ix+1,iy  ,iz  ,vy)-flds(ix-1,iy  ,iz  ,vy))
     1        - dy2i*(flds(ix  ,iy+1,iz  ,vx)-flds(ix  ,iy-1,iz  ,vx))

      enddo
      enddo
      enddo
      else

      ! Face centered curl calculation: average result back to center
      allocate (vxL(ix0:ix1+1, iy0:iy1+1, iz0:iz1+1))
      allocate (vyL(ix0:ix1+1, iy0:iy1+1, iz0:iz1+1))
      allocate (vzL(ix0:ix1+1, iy0:iy1+1, iz0:iz1+1))
      allocate (wxL(ix0:ix1+1, iy0:iy1+1, iz0:iz1+1))
      allocate (wyL(ix0:ix1+1, iy0:iy1+1, iz0:iz1+1))
      allocate (wzL(ix0:ix1+1, iy0:iy1+1, iz0:iz1+1))

      ! Edge centered velocities (*16) based on all sourounding cells
      do iz = iz0, iz1+1
      do iy = iy0, iy1+1
      do ix = ix0, ix1
      vxL(ix,iy,iz) = flds(ix-1,iy  ,iz  ,vx) + flds(ix-1,iy-1,iz  ,vx)
     1              + flds(ix-1,iy  ,iz-1,vx) + flds(ix-1,iy-1,iz-1,vx)
     1         + 2.0*(flds(ix  ,iy  ,iz  ,vx) + flds(ix  ,iy-1,iz  ,vx)
     1              + flds(ix  ,iy  ,iz-1,vx) + flds(ix  ,iy-1,iz-1,vx))
     1              + flds(ix+1,iy  ,iz  ,vx) + flds(ix+1,iy-1,iz  ,vx)
     1              + flds(ix+1,iy  ,iz-1,vx) + flds(ix+1,iy-1,iz-1,vx)
      enddo
      enddo
      enddo

      do iz = iz0, iz1+1
      do iy = iy0, iy1
      do ix = ix0, ix1+1
      vyL(ix,iy,iz) = flds(ix  ,iy-1,iz  ,vy) +f lds(ix-1,iy-1,iz  ,vy)
     1              + flds(ix  ,iy-1,iz-1,vy) +f lds(ix-1,iy-1,iz-1,vy)
     1        + 2.0 *(flds(ix  ,iy  ,iz  ,vy) +f lds(ix-1,iy  ,iz  ,vy)
     1              + flds(ix  ,iy  ,iz-1,vy) +f lds(ix-1,iy  ,iz-1,vy))
     1              + flds(ix  ,iy+1,iz  ,vy) +f lds(ix-1,iy+1,iz  ,vy)
     1              + flds(ix  ,iy+1,iz-1,vy) +f lds(ix-1,iy+1,iz-1,vy)
      enddo
      enddo
      enddo

      do iz = iz0, iz1
      do iy = iy0, iy1+1
      do ix = ix0, ix1+1
      vzL(ix,iy,iz) = flds(ix  ,iy  ,iz-1,vz) +f lds(ix-1,iy  ,iz-1,vz)
     1              + flds(ix  ,iy-1,iz-1,vz) +f lds(ix-1,iy-1,iz-1,vz)
     1         + 2.0*(flds(ix  ,iy  ,iz  ,vz) +f lds(ix-1,iy  ,iz  ,vz)
     1              + flds(ix  ,iy-1,iz  ,vz) +f lds(ix-1,iy-1,iz  ,vz))
     1              + flds(ix  ,iy  ,iz+1,vz) +f lds(ix-1,iy  ,iz+1,vz)
     1              + flds(ix  ,iy-1,iz+1,vz) +f lds(ix-1,iy-1,iz+1,vz)
      enddo
      enddo
      enddo
      ! Face centered circulation (*16/dx)
      ! Assumes dx=dy=dz=const
      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1+1
        wxL(ix,iy,iz) = vyL(ix  ,iy  ,iz  ) + vzL(ix  ,iy+1,iz  )
     1                - vyL(ix  ,iy  ,iz+1) - vzL(ix  ,iy  ,iz  )
      enddo
      enddo
      enddo
      do iz = iz0, iz1
      do iy = iy0, iy1+1
      do ix = ix0, ix1
        wyL(ix,iy,iz) = vxL(ix  ,iy  ,iz  ) + vzL(ix+1,iy  ,iz  )
     1                - vxL(ix  ,iy  ,iz+1) - vzL(ix  ,iy  ,iz  )
      enddo
      enddo
      enddo
      do iz = iz0, iz1+1
      do iy = iy0, iy1
      do ix = ix0, ix1
        wzL(ix,iy,iz) = vxL(ix  ,iy  ,iz  ) + vyL(ix+1,iy  ,iz  )
     1                - vxL(ix  ,iy+1,iz  ) - vyL(ix  ,iy  ,iz  )
      enddo
      enddo
      enddo

      ! Face centered vorticity averaged back to cell center
      ! Sumes, above are (16/dx)*cericulation around each face
      ! Two of these face centered values are used here & need
      ! to divide by face area dx**2)
      fac = 1.0/(32.0*dx)

      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
        flds(ix,iy,iz,iwrk1) = fac * (wxL(ix,iy,iz)+wxl(ix+1,iy,iz))
        flds(ix,iy,iz,iwrk2) = fac * (wyL(ix,iy,iz)+wyl(ix,iy+1,iz))
        flds(ix,iy,iz,iwrk3) = fac * (wzL(ix,iy,iz)+wzl(ix,iy,iz+1))
      enddo
      enddo
      enddo

      deallocate (vxL)
      deallocate (vyL)
      deallocate (vzL)
      deallocate (wxL)
      deallocate (wyL)
      deallocate (wzL)

      endif

      return
      end

C ====================================================================72
      subroutine do_eigenv_s(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1
      integer iWrk1,iWrk2,iWrk3, ia00, ix,iy,iz,iv, i,j, ierr
      real*8 a(3,3), ye(3), ee(3,3)
      
      iwrk1 = iWrk_tof(itof)
      iwrk2 = iWrk1 + 1
      iwrk3 = iWrk1 + 2
      ia00 = iArg_tof(itof,1)

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
        iv = ia00
        do j = 1,3
        do i = 1,3
          a(i,j) = flds(ix,iy,iz,iv)
          iv = iv + 1
        enddo
        enddo
        call princip(a,ye,ee, ierr)
        if(ierr .lt. 0) then
          write (6,*) "Error in principle called from do_eigenv_s"
          stop
        endif
        flds(ix,iy,iz,iwrk1) = ye(1)
        flds(ix,iy,iz,iwrk2) = ye(2)
        flds(ix,iy,iz,iwrk3) = ye(3)
      enddo
      enddo
      enddo

      return
      end

C ====================================================================72
      subroutine do_grad(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1
      integer iWrk1,iWrk2,iWrk3, s0,s00, ix,iy,iz, ndim1, iv
      

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)
      nDim1 = iArg_dim(itof,1)

      if(nDim1 .eq. 1) then
        s0 = iArg_tof(itof,1)
        iwrk1 = iWrk_tof(itof)
        iwrk2 = iWrk1 + 1
        iwrk3 = iWrk1 + 2

        do iz = iz0, iz1
        do iy = iy0, iy1
        do ix = ix0, ix1
        
          flds(ix,iy,iz,iwrk1) =
     1          dx2i*(flds(ix+1,iy  ,iz  ,s0)-flds(ix-1,iy  ,iz  ,s0))
          flds(ix,iy,iz,iwrk2) =
     1          dy2i*(flds(ix  ,iy+1,iz  ,s0)-flds(ix  ,iy-1,iz  ,s0))
          flds(ix,iy,iz,iwrk3) =
     1          dz2i*(flds(ix  ,iy  ,iz+1,s0)-flds(ix  ,iy  ,iz-1,s0))

        enddo
        enddo
        enddo
      else
        do iv = 0, nDim1-1
          s0 = iArg_tof(itof,1) + iv
          iwrk1 = iWrk_tof(itof) + iv
          iwrk2 = iwrk1 +   nDim1
          iwrk3 = iwrk1 + 2*nDim1

          do iz = iz0, iz1
          do iy = iy0, iy1
          do ix = ix0, ix1
        
            flds(ix,iy,iz,iwrk1) =
     1          dx2i*(flds(ix+1,iy  ,iz  ,s0)-flds(ix-1,iy  ,iz  ,s0))
            flds(ix,iy,iz,iwrk2) =
     1          dy2i*(flds(ix  ,iy+1,iz  ,s0)-flds(ix  ,iy-1,iz  ,s0))
            flds(ix,iy,iz,iwrk3) =
     1          dz2i*(flds(ix  ,iy  ,iz+1,s0)-flds(ix  ,iy  ,iz-1,s0))

          enddo
          enddo
          enddo
        enddo
      endif

      return
      end

C ====================================================================72



C ====================================================================72
      subroutine do_grad_1(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max), dxi
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1
      integer iWrk1,iWrk2,iWrk3, s0, ix,iy,iz
      
      iwrk1 = iWrk_tof(itof)
      s0 = iArg_tof(itof,1)
      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      if(cOpr_tof(itof)(7:7) .ne. "L") then
        do iz = iz0, iz1
        do iy = iy0, iy1
        do ix = ix0, ix1
          flds(ix,iy,iz,iwrk1) =
     1            dx2i*(flds(ix+1,iy,iz,s0)-flds(ix-1,iy,iz,s0))
        enddo
        enddo
        enddo
      else
        dxi = 1.0 / dx
        do iz = iz0, iz1
        do iy = iy0, iy1
        do ix = ix0, ix1
          flds(ix,iy,iz,iwrk1) =
     1          dxi*(flds(ix,iy,iz,s0)-flds(ix-1,iy,iz,s0))
        enddo
        enddo
        enddo
      endif

      return
      end
C ====================================================================72
      subroutine do_grad_2(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max), dyi
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1
      integer iWrk1,iWrk2,iWrk3, s0, ix,iy,iz
      
      iwrk1 = iWrk_tof(itof)
      s0 = iArg_tof(itof,1)
      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      if(cOpr_tof(itof)(7:7) .ne. "L") then
        do iz = iz0, iz1
        do iy = iy0, iy1
        do ix = ix0, ix1
          flds(ix,iy,iz,iwrk1) =
     1          dy2i*(flds(ix  ,iy+1,iz  ,s0)-flds(ix  ,iy-1,iz  ,s0))
        enddo
        enddo
        enddo
      else
        dyi = 1.0 / dy
        do iz = iz0, iz1
        do iy = iy0, iy1
        do ix = ix0, ix1
          flds(ix,iy,iz,iwrk1) =
     1          dyi*(flds(ix,iy,iz,s0) - flds(ix,iy-1,iz,s0))
        enddo
        enddo
        enddo
      endif

      return
      end
C ====================================================================72
      subroutine do_grad_3(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max), dzi
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1
      integer iWrk1,iWrk2,iWrk3, s0, ix,iy,iz
      
      iwrk1 = iWrk_tof(itof)
      s0 = iArg_tof(itof,1)
      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      if(cOpr_tof(itof)(7:7) .ne. "L") then
        do iz = iz0, iz1
        do iy = iy0, iy1
        do ix = ix0, ix1
          flds(ix,iy,iz,iwrk1) =
     1          dz2i*(flds(ix  ,iy  ,iz+1,s0)-flds(ix  ,iy  ,iz-1,s0))
        enddo
        enddo
        enddo
      else
        dzi = 1.0 / dz
        do iz = iz0, iz1
        do iy = iy0, iy1
        do ix = ix0, ix1
          flds(ix,iy,iz,iwrk1) =
     1          dzi*(flds(ix  ,iy  ,iz,s0)-flds(ix  ,iy  ,iz-1,s0))
        enddo
        enddo
        enddo
      endif

      return
      end

C ====================================================================72
      subroutine do_div(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1
      integer iWrk1, vx,vy,vz, ix,iy,iz

      if(isRZ .eq. 1) then
        call do_div_RZ(itof, i_wb, flds)
        return
      endif
      
      iwrk1 = iWrk_tof(itof)
      vx = iArg_tof(itof,1)
      vy = vx + 1
      vz = vx + 2

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
        
        flds(ix,iy,iz,iwrk1) =
     1          dx2i*(flds(ix+1,iy  ,iz  ,vx)-flds(ix-1,iy  ,iz  ,vx))
     1        + dy2i*(flds(ix  ,iy+1,iz  ,vy)-flds(ix  ,iy-1,iz  ,vy))
     1        + dz2i*(flds(ix  ,iy  ,iz+1,vz)-flds(ix  ,iy  ,iz-1,vz))

      enddo
      enddo
      enddo

      return
      end

C ====================================================================72
      subroutine do_div_RZ(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1
      integer iWrk1, vx,vy, ix,iy,iz
      real r0,r1,areaZface,volinv,vrL,vrR,vzL,vzR
      real zfac(-9:3000), rfacL(-9:3000), rfacR(-9:2000)

      iwrk1 = iWrk_tof(itof)
      vx = iArg_tof(itof,1)
      vy = vx + 1

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)

      do ix=ix0,ix1+1
        r0 = dx * float(ix-1 + ioff_wb(1,i_wb))
        r1 = dx * float(ix   + ioff_wb(1,i_wb))
        areaZface = r1**2 - r0**2
        volinv    = 1.0 / (dy * areaZface)
        zfac(ix) = areaZface * volinv
        rfacL(ix) = 2.0 * r0 * dy * volinv
        rfacR(ix) = 2.0 * r1 * dy * volinv
      enddo

      iz = 1
      do iy = iy0, iy1
      do ix = ix0, ix1
        vrL = 0.5*(flds(ix  ,iy  ,iz,vx)+flds(ix-1,iy  ,iz,vx))
        vrR = 0.5*(flds(ix+1,iy  ,iz,vx)+flds(ix  ,iy  ,iz,vx))
        vzL = 0.5*(flds(ix  ,iy  ,iz,vy)+flds(ix  ,iy-1,iz,vy))
        vzR = 0.5*(flds(ix  ,iy+1,iz,vy)+flds(ix  ,iy  ,iz,vy))
        flds(ix,iy,iz,iwrk1) = (rfacR(ix)*vrR - rfacL(ix)*vrL)
     1                       + ( zfac(ix)*vzR -  zfac(ix)*vzL)
      enddo
      enddo

      return
      end
cccccccccccccccccc
cc      yyy = dy*(iy+ioff_wb(2,i_wb))
cc      if(yyy .gt. 0.7 .and.  ix.lt.10) then
ccc      if(yyy .gt. 0.7 .and.  ix.lt.2) then
cc ranges
c      if(iy+ioff_wb(2,i_wb).eq.iyout1 .and.  ix.lt.5) then
c        write (6,*) "========================================"
c        write (6,*) "div,Fr ", flds(ix,iy,iz,iwrk1),
c     1                         flds(ix,iy,iz,vx)
cc      if(iy+ioff_wb(2,i_wb).ge.iyout1 .and.  ix.eq.1 .and.
cc         iy+ioff_wb(2,i_wb).le.iyout2) then
cc        write (6,*) "========================================"
cc        write (6,*) "div,Fz ", flds(ix,iy,iz,iwrk1),
cc     1                         flds(ix,iy,iz,vy)
c        write (6,*) "div_r  ", rfacR(ix)*vrR - rfacL(ix)*vrL
c        write (6,*) "div_z  ", zfac(ix)*vzR -  zfac(ix)*vzL
c        write (6,*) "rfacRL ", rfacR(ix), rfacL(ix)
c        write (6,*) "vr[RL] ", vrR, vrL
c        write (6,*) "zfac    ", zfac(ix)
c        write (6,*) "vz[RL] ", vzR, vzL
cc        write (6,998) zfac(ix),vzR, vzL,
cc     1                flds(ix,iy,iz,vy),
cc     1                flds(ix,iy,iz,iwrk1)
cc998     format("zfac,vzR,vzL,f,f:",1p5e12.4)
cc
cc        write (6,999) rfac(ix+1),vrR, rfac(ix),vrL,
cc     1                flds(ix,iy,iz,vx)
cc999     format("rfacR,vrR,rfacL,vrL,f:",1p5e12.4)
cc
cc   1/dZ = 1 / 0.00775 = 129.03225806451612903225
c      endif
cccccccccccccccccc


C ====================================================================72
      subroutine do_cross(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1, ix,iy,iz
      integer iax,iay,iaz, ibx,iby,ibz, icx,icy,icz

      iax = iWrk_tof(itof)
      iay = iax + 1
      iaz = iax + 2

      ibx = iArg_tof(itof,1)
      iby = ibx + 1
      ibz = ibx + 2

      icx = iArg_tof(itof,2)
      icy = icx + 1
      icz = icx + 2

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
        
        flds(ix,iy,iz,iax) = flds(ix,iy,iz,iby)*flds(ix,iy,iz,icz)
     1                     - flds(ix,iy,iz,ibz)*flds(ix,iy,iz,icy)

        flds(ix,iy,iz,iay) = flds(ix,iy,iz,ibz)*flds(ix,iy,iz,icx)
     1                     - flds(ix,iy,iz,ibx)*flds(ix,iy,iz,icz)

        flds(ix,iy,iz,iaz) = flds(ix,iy,iz,ibx)*flds(ix,iy,iz,icy)
     1                     - flds(ix,iy,iz,iby)*flds(ix,iy,iz,icx)

      enddo
      enddo
      enddo

      return
      end


C ====================================================================72
      subroutine do_agradb(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1, ix,iy,iz
      integer iax,iay,iaz, ibx,iby,ibz, icx,icy,icz
      real adx, ady, adz, ax, ay, az

      icx = iWrk_tof(itof)
      icy = icx + 1
      icz = icx + 2

      iax = iArg_tof(itof,1)
      iay = iax + 1
      iaz = iax + 2

      ibx = iArg_tof(itof,2)
      iby = ibx + 1
      ibz = ibx + 2

c      write (6,*) 'icx,icy,icz=',icx,icy,icz
c      write (6,*) 'iax,iay,iaz=',iax,iay,iaz
c      write (6,*) 'ibx,iby,ibz=',ibx,iby,ibz
c      write (6,*) 'dx2i=',dx2i
c      write (6,*) 'dy2i=',dy2i
c      write (6,*) 'dz2i=',dz2i

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
         
        adx = dx2i * flds(ix,iy,iz,iax)
        ady = dy2i * flds(ix,iy,iz,iay)
        adz = dz2i * flds(ix,iy,iz,iaz)

        ax = adx*(flds(ix+1,iy  ,iz  ,ibx) - flds(ix-1,iy  ,iz  ,ibx))
     1     + ady*(flds(ix  ,iy+1,iz  ,ibx) - flds(ix  ,iy-1,iz  ,ibx))
     1     + adz*(flds(ix  ,iy  ,iz+1,ibx) - flds(ix  ,iy  ,iz-1,ibx))

        ay = adx*(flds(ix+1,iy  ,iz  ,iby) - flds(ix-1,iy  ,iz  ,iby))
     1     + ady*(flds(ix  ,iy+1,iz  ,iby) - flds(ix  ,iy-1,iz  ,iby))
     1     + adz*(flds(ix  ,iy  ,iz+1,iby) - flds(ix  ,iy  ,iz-1,iby))

        az = adx*(flds(ix+1,iy  ,iz  ,ibz) - flds(ix-1,iy  ,iz  ,ibz))
     1     + ady*(flds(ix  ,iy+1,iz  ,ibz) - flds(ix  ,iy-1,iz  ,ibz))
     1     + adz*(flds(ix  ,iy  ,iz+1,ibz) - flds(ix  ,iy  ,iz-1,ibz))

        flds(ix,iy,iz,icx) = ax
        flds(ix,iy,iz,icy) = ay
        flds(ix,iy,iz,icz) = az

c         write (6,999)
c     1      adx*(flds(ix+1,iy  ,iz  ,ibx) - flds(ix-1,iy  ,iz ,ibx)),
c     1      flds(ix,iy,iz,6) * flds(ix,iy,iz,9)
c999      format("(VgradB3)_xx=",f10.5,"      VgradBx_x=",f10.5)
c         write (6,999)
c     1      ady*(flds(ix  ,iy+1,iz  ,ibx) - flds(ix  ,iy-1,iz  ,ibx)),
c     1      flds(ix,iy,iz,7) * flds(ix,iy,iz,10)
c999      format("(VgradB3)_xy=",f10.5,"      VgradBx_y=",f10.5)
c         write (6,999)
c     1      adz*(flds(ix  ,iy  ,iz+1,ibx) - flds(ix  ,iy  ,iz-1,ibx)),
c     1      flds(ix,iy,iz,8) * flds(ix,iy,iz,11)
c999      format("(VgradB3)_xz=",f10.5,"      VgradBx_z=",f10.5)
c         write (6,999)
c     1      dz2i*(flds(ix  ,iy  ,iz+1,ibx) - flds(ix  ,iy  ,iz-1,ibx)),
c     1      flds(ix,iy,iz,11)
c999      format("agradb:Bx,xz=",f10.5,"      Bx,=",f10.5)
c         write (6,999) adz, dz2i*flds(ix,iy,iz,8)
c999      format("agradb:adz",f10.5,"    dz2i*Vz,=",f10.5)

      enddo
c      if(idebug .gt. 0) stop
      enddo
      enddo


      return
      end

C ====================================================================72
      subroutine do_dot(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1, ix,iy,iz
      integer iax, ibx,iby,ibz, icx,icy,icz, nDim1, nDim2
      integer i1s,i2s,j1,j2,ib0,ic0,ia,ib,ic, i
      real ss

      iax = iWrk_tof(itof)
      ibx = iArg_tof(itof,1)
      icx = iArg_tof(itof,2)
      nDim1 = iArg_dim(itof,1)
      nDim2 = iArg_dim(itof,2)

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      if(ndim1 .eq. 3  .and.  nDim2 .eq. 3) then
        iby = ibx + 1
        ibz = ibx + 2
        icy = icx + 1
        icz = icx + 2
        do iz = iz0, iz1
        do iy = iy0, iy1
        do ix = ix0, ix1
          flds(ix,iy,iz,iax) = flds(ix,iy,iz,ibx)*flds(ix,iy,iz,icx)
     1                       + flds(ix,iy,iz,iby)*flds(ix,iy,iz,icy)
     1                       + flds(ix,iy,iz,ibz)*flds(ix,iy,iz,icz)
        enddo
        enddo
        enddo
      else
        i1s = nDim1 / 3
        i2s = nDim2 / 3

        do iz = iz0, iz1
        do iy = iy0, iy1
        do ix = ix0, ix1
          do j1 = 0, i1s - 1
          do j2 = 0, i2s - 1
            ia  = iax + j1 + i1s * j2
            ib0 = ibx + j1
            ic0 = icx + j2 * i2s
            ss = 0.0
            do i = 0, 2
             ib = ib0 + i1s*i   ! always contract on last index of 1st operand
             ic = ic0 +     i   ! always contract on 1st  index of 2nd operand
             ss = ss + flds(ix,iy,iz,ib) * flds(ix,iy,iz,ic)
            enddo
            flds(ix,iy,iz,ia) = ss
          enddo
          enddo
        enddo
        enddo
        enddo
      endif

      if(idebug .gt. 1000) then
      write (6,*) 'in dot:'
      write (6,*) 'ix0,ix1: ', ix0, ix1
      write (6,*) 'iy0,iy1: ', iy0, iy1
      write (6,*) 'iz0,iz1: ', iz0, iz1
      write (6,*) 'iax,flds:',iax,flds(1,1,1,iax)
      write (6,*) 'ibx,flds:',ibx,flds(1,1,1,ibx)
      write (6,*) 'iby,flds:',iby,flds(1,1,1,iby)
      write (6,*) 'ibz,flds:',ibz,flds(1,1,1,ibz)
      write (6,*) 'icx,flds:',icx,flds(1,1,1,icx)
      write (6,*) 'icy,flds:',icy,flds(1,1,1,icy)
      write (6,*) 'icz,flds:',icz,flds(1,1,1,icz)
      if(iax .ge. 0) stop
      endif

      return
      end

C ====================================================================72
      subroutine get_sect_vals_from_string(str,rad,crg)
      implicit none
      include 'postproc_info_e3d.h'
      character*(*) str
      real rad, crg, value_sect, value, fsect
      integer isect, ic , isec_value

      rad = -1000
      crg = -1000
      do ic = 1, 2
       isect = isec_value(str,cSecFormat(ic),cmsigns)
       fsect = float(isect)
       if(cSecLoc(ic) .eq. "low")  fsect = float(isect) - 0.5
       if(cSecLoc(ic) .eq. "high") fsect = float(isect) + 0.5
       value = value_sect(fsect,ic)
        if(cSecName(ic) .eq. "Diameter") rad = 0.5 * value
        if(cSecName(ic) .eq. "Radius") rad = value
        if(cSecName(ic) .eq. "Charge") crg = value
      enddo
c gen_sect_bins

      return
      end

C ====================================================================72
      subroutine do_secsum(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1, ix,iy,iz
      integer ia, ib, nDim1, iv, j, istof, imatch_tof, n
      real rad, crg, pi, rfac, rpow, cfac
      real*8 wt, s0, ss, val(1000)
      logical bAverage


      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      ia = iWrk_tof(itof)
      ib = iArg_tof(itof,1)
      istof = imatch_tof(cArg_tof(itof,1),ntof,cFld_tof)

      ! Default is just sum of concentrations
      rfac = 1.0
      rpow = 0.0
      cfac = 0.0
      bAverage = .false.
      pi = 4.0 * atan2(1.0,1.0)
      n = len(trim(cOpr_tof(itof)))
      if(cOpr_tof(itof)(n-5:n) .eq. "radius" ) then
         rfac = 1.0
         rpow = 1.0
         bAverage = .true.
       endif
      if(cOpr_tof(itof)(n-6:n) .eq. "surface") then
         rfac = 4.0 * pi
         rpow = 2.0
         bAverage = .true.
       endif
      if(cOpr_tof(itof)(n-5:n) .eq. "volume" ) then
         rfac = (4.0/3.0) * pi
         rpow = 3.0
         bAverage = .true.
       endif
      if(cOpr_tof(itof)(n-5:n) .eq. "charge" ) then
         rfac = 0.0
         cfac = 1.0
         bAverage = .true.
      endif

      ! Generate array with sectional values
      nDim1 = iArg_dim(itof,1)
      do iv = 1, nDim1
        call get_sect_vals_from_string(cArg_tof(istof,iv),rad,crg)
        val(iv) = rfac * rad**rpow + cfac*crg
      enddo

      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
        s0 = 0.0
        ss = 0.0
        do iv = 1, nDim1
          wt = flds(ix,iy,iz,ib+iv)
          s0 = s0 + wt
          ss = ss + wt * val(iv)
        enddo
        flds(ix,iy,iz,ia) = ss
        if(s0 .gt. 0.0 .and. bAverage) flds(ix,iy,iz,ia) = ss / s0
      enddo
      enddo
      enddo

      return
      end

C ====================================================================72
      subroutine do_determinant(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      real a11,a12,a13,a21,a22,a23,a31,a32,a33
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1, ix,iy,iz
      integer ia, ib

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)
      ia  = iWrk_tof(itof)
      ib = iArg_tof(itof,1)

      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
        a11 = flds(ix,iy,iz,ib  )
        a21 = flds(ix,iy,iz,ib+1)
        a31 = flds(ix,iy,iz,ib+2)
        a12 = flds(ix,iy,iz,ib+3)
        a22 = flds(ix,iy,iz,ib+4)
        a32 = flds(ix,iy,iz,ib+5)
        a13 = flds(ix,iy,iz,ib+6)
        a23 = flds(ix,iy,iz,ib+7)
        a33 = flds(ix,iy,iz,ib+8)
        flds(ix,iy,iz,ia) = a11*a22*a33 + a12*a23*a31 + a13*a21*a32
     1                    - a11*a23*a32 - a12*a21*a33 - a13*a22*a31
      enddo
      enddo
      enddo

      return
      end

C ====================================================================72
      subroutine do_trace(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1, ix,iy,iz
      integer ia, ib0,ib1,ib2, nDim1, iv

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      nDim1 = iArg_dim(itof,1)
      do iv = 0, (nDim1-1)/9
        ia  = iWrk_tof(itof) + iv
        ib0 = iArg_tof(itof,1) + 3*iv
        ib1 = ib0 + 1+nDim1/3
        ib2 = ib1 + 1+nDim1/3
        do iz = iz0, iz1
        do iy = iy0, iy1
        do ix = ix0, ix1
          flds(ix,iy,iz,ia) = flds(ix,iy,iz,ib0)
     1                      + flds(ix,iy,iz,ib1)
     1                      + flds(ix,iy,iz,ib2)
        enddo
        enddo
        enddo
      enddo

      return
      end

C ====================================================================72
      subroutine do_subtrace(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1, ix,iy,iz
      integer ia, ib, nDim1
      real    tr3

      ia = iWrk_tof(itof)
      ib = iArg_tof(itof,1)
      nDim1 = iArg_dim(itof,1)

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
        tr3 = (flds(ix,iy,iz,ib) + flds(ix,iy,iz,ib+4)
     1                            + flds(ix,iy,iz,ib+8)) / 3.0
        flds(ix,iy,iz,ia  ) = flds(ix,iy,iz,ib  ) - tr3
        flds(ix,iy,iz,ia+1) = flds(ix,iy,iz,ib+1)
        flds(ix,iy,iz,ia+2) = flds(ix,iy,iz,ib+2)
        flds(ix,iy,iz,ia+3) = flds(ix,iy,iz,ib+3)
        flds(ix,iy,iz,ia+4) = flds(ix,iy,iz,ib+4) - tr3
        flds(ix,iy,iz,ia+5) = flds(ix,iy,iz,ib+5)
        flds(ix,iy,iz,ia+6) = flds(ix,iy,iz,ib+6)
        flds(ix,iy,iz,ia+7) = flds(ix,iy,iz,ib+7)
        flds(ix,iy,iz,ia+8) = flds(ix,iy,iz,ib+8) - tr3
      enddo
      enddo
      enddo

      return
      end

C ====================================================================72
      subroutine do_ball_filter(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1, ix,iy,iz
      integer ia, ib, iac, ibc, jx,jy,jz, n, nDim1, idim
      real    fac, s, fsmooth

      ia = iWrk_tof(itof)
      ib = iArg_tof(itof,1)
      nDim1 = iArg_dim(itof,1)

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      fsmooth = 0.5 + float(nsmooth)
      n = 0
      do jz = -nsmooth, nsmooth
      do jy = -nsmooth, nsmooth
      do jx = -nsmooth, nsmooth
        if(float(jx**2 + jy**2 + jz**2) .lt. fsmooth) n = n + 1
      enddo
      enddo
      enddo
      fac = 1.0 / float(n)

      do idim = 1, ndim1
      iac = ia + idim - 1
      ibc = ib + idim - 1
      do iz = iz0, iz1
      do iy = iy0, iy1
      do ix = ix0, ix1
        s = 0.0
        do jz = -nsmooth, nsmooth
        do jy = -nsmooth, nsmooth
        do jx = -nsmooth, nsmooth
          if(float(jx**2 + jy**2 + jz**2) .lt. fsmooth) then
            s = s + flds(ix+jx,iy+jy,iz+jz,ibc)
          endif
        enddo
        enddo
        enddo
        flds(ix,iy,iz,iac) = fac * s
      enddo
      enddo
      enddo
      enddo

      return
      end


C ====================================================================72
      subroutine do_gauss_filt(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1, ix,iy,iz, m, isrc, j, n1,n2,n3
      integer ia, ib, iac, ibc, jx,jy,jz, n, nDim1, idim, nchar,m1,m2,m3
      real    fac, s, fsmooth, a
      allocatable a(:,:,:,:)

      !                   smooth_???????
      nchar = len(trim(cOpr_tof(itof)))
      call parse_xyz_pars(cOpr_tof(itof)(8:nchar), n1,n2,n3)
      if(max(n1,n2,n3) .le. 0) then
        call do_gauss_filt_old(itof, i_wb, flds)
        return
      endif

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      n = max(n1,n2,n3)
      allocate(a(ix0-n:ix1+n,iy0-n:iy1+n,iz0-n:iz1+n,2))

      nDim1 = iArg_dim(itof,1)
      do idim = 0, nDim1-1
        ia = idim + iWrk_tof(itof)
        ib = idim + iArg_tof(itof,1)
        call cp_f2a(ib,flds,ix0,ix1,iy0,iy1,iz0,iz1,n,n1,n2,n3,1,a)
        m1 = n1
        m2 = n2
        m3 = n3
        isrc = 1
        do m = 1, n
          if(m1 .gt. 0)
     1      call filt_xdir(ix0,ix1,iy0,iy1,iz0,iz1,n,m1,m2,m3,isrc,a)
          if(m2 .gt. 0)
     1      call filt_ydir(ix0,ix1,iy0,iy1,iz0,iz1,n,m1,m2,m3,isrc,a)
          if(m3 .gt. 0)
     1      call filt_zdir(ix0,ix1,iy0,iy1,iz0,iz1,n,m1,m2,m3,isrc,a)
        enddo
        call cp_to_flds(ia,flds,ix0,ix1,iy0,iy1,iz0,iz1,n,0,isrc,a)
      enddo

      deallocate(a)

      return
      end


C ====================================================================72
      subroutine do_gauss_filt_old(itof, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1, ix,iy,iz, m, idst, j
      integer ia, ib, iac, ibc, jx,jy,jz, n, nDim1, idim
      real    fac, s, fsmooth, a
      allocatable a(:,:,:,:)

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

       n = nsmooth
      allocate(a(ix0-n:ix1+n,iy0-n:iy1+n,iz0-n:iz1+n,2))
c      write (6,*) "in do_gauss_filt ix1,n=",ix1,n
c      write (6,999) ix0,ix1,iy0,iy1,iz0,iz1,n
c999   format("do_gauss_filt  x:",2i5,"    y:",2i5,"    z:",2i5,
c     1       "    n:",i4)

      nDim1 = iArg_dim(itof,1)
      do idim = 0, nDim1-1
        ia = idim + iWrk_tof(itof)
        ib = idim + iArg_tof(itof,1)
        call cp_from_flds(ib,flds,ix0,ix1,iy0,iy1,iz0,iz1,n,n,1,a)
        idst = 1
        do m = n-1, 0, -1
           do j=1,3
             idst = 3 - idst
             call filt_1dir(ix0,ix1,iy0,iy1,iz0,iz1,n,m,j,idst,a)
           enddo
        enddo
        call cp_to_flds(ia,flds,ix0,ix1,iy0,iy1,iz0,iz1,n,0,idst,a)
      enddo

      deallocate(a)
c      write (6,*) "do_gauss_filt done"

      return
      end

C ====================================================================72
      subroutine cp_f2a(in,flds,ix0,ix1,iy0,iy1,iz0,iz1,n,m1,m2,m3,io,a)
      implicit none
      include 'postproc_info_e3d.h'
      integer in,ix0,ix1,iy0,iy1,iz0,iz1,n,m1,m2,m3,io
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      real a(ix0-n:ix1+n,iy0-n:iy1+n,iz0-n:iz1+n,2)
      integer ix,iy,iz

c      write (6,*) "in cp_from_flds: iz1,m=",iz1,m

      do iz = iz0-m3, iz1+m3
      do iy = iy0-m2, iy1+m2
      do ix = ix0-m1, ix1+m1
        a(ix,iy,iz,io) = flds(ix,iy,iz,in)
      enddo
      enddo
      enddo

      return
      end
      
C ====================================================================72
      subroutine cp_from_flds(in,flds,ix0,ix1,iy0,iy1,iz0,iz1,n,m,io,a)
      implicit none
      include 'postproc_info_e3d.h'
      integer in,ix0,ix1,iy0,iy1,iz0,iz1,n,m,io
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      real a(ix0-n:ix1+n,iy0-n:iy1+n,iz0-n:iz1+n,2)
      integer ix,iy,iz

c      write (6,*) "in cp_from_flds: iz1,m=",iz1,m

      do iz = iz0-m, iz1+m
      do iy = iy0-m, iy1+m
      do ix = ix0-m, ix1+m
        a(ix,iy,iz,io) = flds(ix,iy,iz,in)
      enddo
      enddo
      enddo

      return
      end
      
C ====================================================================72
      subroutine filt_xdir(ix0,ix1,iy0,iy1,iz0,iz1,n,m1,m2,m3,isrc,a)
      implicit none
      include 'postproc_info_e3d.h'
      integer in,ix0,ix1,iy0,iy1,iz0,iz1,n,m1,m2,m3,isrc
      real a(ix0-n:ix1+n,iy0-n:iy1+n,iz0-n:iz1+n,2), fac
      integer ix,iy,iz, idst, m1a

      m1a = m1 - 1
      fac = 1.0 / 3.0
      idst = 3 - isrc
      do iz = iz0-m3 , iz1+m3
      do iy = iy0-m2 , iy1+m2
      do ix = ix0-m1a, ix1+m1a
        a(ix,iy,iz,idst) = fac*(a(ix-1,iy  ,iz  ,isrc) +
     1                          a(ix  ,iy  ,iz  ,isrc) +
     1                          a(ix+1,iy  ,iz  ,isrc)  )
      enddo
      enddo
      enddo
      m1   = m1a
      isrc = idst

      return
      end
c-----------------------------------------------------------------------
      subroutine filt_ydir(ix0,ix1,iy0,iy1,iz0,iz1,n,m1,m2,m3,isrc,a)
      implicit none
      include 'postproc_info_e3d.h'
      integer in,ix0,ix1,iy0,iy1,iz0,iz1,n,m1,m2,m3,isrc
      real a(ix0-n:ix1+n,iy0-n:iy1+n,iz0-n:iz1+n,2), fac
      integer ix,iy,iz, idst, m2a

      m2a = m2 - 1
      fac = 1.0 / 3.0
      idst = 3 - isrc
      do iz = iz0-m3 , iz1+m3
      do iy = iy0-m2a, iy1+m2a
      do ix = ix0-m1 , ix1+m1
        a(ix,iy,iz,idst) = fac*(a(ix  ,iy-1,iz  ,isrc) +
     1                          a(ix  ,iy  ,iz  ,isrc) +
     1                          a(ix  ,iy+1,iz  ,isrc)  )
      enddo
      enddo
      enddo
      m2   = m2a
      isrc = idst

      return
      end
c-----------------------------------------------------------------------
      subroutine filt_zdir(ix0,ix1,iy0,iy1,iz0,iz1,n,m1,m2,m3,isrc,a)
      implicit none
      include 'postproc_info_e3d.h'
      integer in,ix0,ix1,iy0,iy1,iz0,iz1,n,m1,m2,m3,isrc
      real a(ix0-n:ix1+n,iy0-n:iy1+n,iz0-n:iz1+n,2), fac
      integer ix,iy,iz, idst, m3a

      m3a = m3 - 1
      fac = 1.0 / 3.0
      idst = 3 - isrc
      do iz = iz0-m3a, iz1+m3a
      do iy = iy0-m2 , iy1+m2
      do ix = ix0-m1 , ix1+m1
        a(ix,iy,iz,idst) = fac*(a(ix  ,iy  ,iz-1,isrc) +
     1                          a(ix  ,iy  ,iz  ,isrc) +
     1                          a(ix  ,iy  ,iz+1,isrc)  )
      enddo
      enddo
      enddo
      m3   = m3a
      isrc = idst

      return
      end
c-----------------------------------------------------------------------
      
C ====================================================================72
      subroutine filt_1dir(ix0,ix1,iy0,iy1,iz0,iz1,n,m,idir,idst,a)
      implicit none
      include 'postproc_info_e3d.h'
      integer in,ix0,ix1,iy0,iy1,iz0,iz1,n,m,idir,idst
      real a(ix0-n:ix1+n,iy0-n:iy1+n,iz0-n:iz1+n,2), fac
      integer ix,iy,iz, isrc, m1

      m1 = m + 1
      fac = 1.0 / 3.0
      isrc = 3 - idst
      if(idir .eq. 1) then
        do iz = iz0-m1, iz1+m1
        do iy = iy0-m1, iy1+m1
        do ix = ix0-m , ix1+m
          a(ix,iy,iz,idst) = fac*(a(ix-1,iy  ,iz  ,isrc) +
     1                            a(ix  ,iy  ,iz  ,isrc) +
     1                            a(ix+1,iy  ,iz  ,isrc)  )
        enddo
        enddo
        enddo
      else if(idir .eq. 2) then
        do iz = iz0-m1, iz1+m1
        do iy = iy0-m , iy1+m
        do ix = ix0-m , ix1+m
          a(ix,iy,iz,idst) = fac*(a(ix  ,iy-1,iz  ,isrc) +
     1                            a(ix  ,iy  ,iz  ,isrc) +
     1                            a(ix  ,iy+1,iz  ,isrc)  )
        enddo
        enddo
        enddo
      else
        do iz = iz0-m, iz1+m
        do iy = iy0-m, iy1+m
        do ix = ix0-m, ix1+m
          a(ix,iy,iz,idst) = fac*(a(ix  ,iy  ,iz-1,isrc) +
     1                            a(ix  ,iy  ,iz  ,isrc) +
     1                            a(ix  ,iy  ,iz+1,isrc)  )
        enddo
        enddo
        enddo
      endif

      return
      end
      
C ====================================================================72
      subroutine cp_to_flds(io,flds,ix0,ix1,iy0,iy1,iz0,iz1,n,m,in,a)
      implicit none
      include 'postproc_info_e3d.h'
      integer in,ix0,ix1,iy0,iy1,iz0,iz1,n,m,io
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      real a(ix0-n:ix1+n,iy0-n:iy1+n,iz0-n:iz1+n,2)
      integer ix,iy,iz

      do iz = iz0-m, iz1+m
      do iy = iy0-m, iy1+m
      do ix = ix0-m, ix1+m
        flds(ix,iy,iz,io) = a(ix,iy,iz,in)
      enddo
      enddo
      enddo

      return
      end
      
C ====================================================================72
      subroutine do_math(itof, ido, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, ido, i_wb
      integer ix0,ix1,iy0,iy1,iz0,iz1, ix,iy,iz
      integer ia, ib, ic, iv, nDim1, nDim2
      real    fac, add_const, div_const, xx

      ia = iWrk_tof(itof)
      ib = iArg_tof(itof,1)
      ic = iArg_tof(itof,2)

      ix0 = 1             - NBdy_tof(1,itof)
      ix1 = n0_wb(1,i_wb) + NBdy_tof(1,itof)
      iy0 = 1             - NBdy_tof(2,itof)
      iy1 = n0_wb(2,i_wb) + NBdy_tof(2,itof)
      iz0 = 1             - NBdy_tof(3,itof)
      iz1 = n0_wb(3,i_wb) + NBdy_tof(3,itof)

      if(ido .eq. iop_scale) then
        nDim1 = max(iArg_dim(itof,1), iArg_dim(itof,2))
        fac = fVal_tof(itof)
        if(ic .eq. 0) ic = ib
        do iv = 0, nDim1-1
        do iz = iz0, iz1
        do iy = iy0, iy1
        do ix = ix0, ix1
          flds(ix,iy,iz,ia+iv) = fac * flds(ix,iy,iz,ic+iv)
        enddo
        enddo
        enddo
        enddo
      else if(ido .eq. iop_plus) then
        nDim1 = iArg_dim(itof,1)
        ia = iWrk_tof(itof)   - 1
        ic = iArg_tof(itof,2) - 1
        do ib = iArg_tof(itof,1), iArg_tof(itof,1)+nDim1-1
          ia = ia + 1
          ic = ic + 1
          do iz = iz0, iz1
          do iy = iy0, iy1
          do ix = ix0, ix1
            flds(ix,iy,iz,ia) = flds(ix,iy,iz,ib) + flds(ix,iy,iz,ic)
          enddo
          enddo
          enddo
        enddo
      else if(ido .eq. iop_minus) then
        nDim1 = iArg_dim(itof,1)
        ia = iWrk_tof(itof)   - 1
        ic = iArg_tof(itof,2) - 1
        do ib = iArg_tof(itof,1), iArg_tof(itof,1)+nDim1-1
          ia = ia + 1
          ic = ic + 1
          do iz = iz0, iz1
          do iy = iy0, iy1
          do ix = ix0, ix1
            flds(ix,iy,iz,ia) = flds(ix,iy,iz,ib) - flds(ix,iy,iz,ic)
          enddo
          enddo
          enddo
        enddo
      else if(ido .eq. iop_add) then
        add_const = fVal_tof(itof)
        if(ic .eq. 0) ic = ib
        do iz = iz0, iz1
        do iy = iy0, iy1
        do ix = ix0, ix1
          flds(ix,iy,iz,ia) = add_const + flds(ix,iy,iz,ic)
        enddo
        enddo
        enddo
      else if(ido .eq. iop_cminus) then
        add_const = fVal_tof(itof)
        if(ic .eq. 0) ic = ib
        do iz = iz0, iz1
        do iy = iy0, iy1
        do ix = ix0, ix1
          flds(ix,iy,iz,ia) = add_const - flds(ix,iy,iz,ic)
        enddo
        enddo
        enddo
      else if(ido .eq. iop_times) then
        nDim1 = iArg_dim(itof,1)
        nDim2 = iArg_dim(itof,2)

        ia = iWrk_tof(itof) - 1
        do ic = iArg_tof(itof,2), iArg_tof(itof,2)+nDim2-1
        do ib = iArg_tof(itof,1), iArg_tof(itof,1)+nDim1-1
          ia = ia + 1
          do iz = iz0, iz1
          do iy = iy0, iy1
          do ix = ix0, ix1
            flds(ix,iy,iz,ia) = flds(ix,iy,iz,ib) * flds(ix,iy,iz,ic)
          enddo
          enddo
          enddo
        enddo
        enddo
c        write (6,*) "itof = ", itof
c        write (6,*) "iArg1,iarg2=", iArg_tof(itof,1),iArg_tof(itof,2)
c        write (6,*) "nDim1,nDim2 = ", nDIm1,nDIm2
c        write (6,*) 'ia,ib,ic=', ia, ib, ic
c        write (6,*) 'flds(1,1,1,ia) = ', flds(1,1,1,ia)
c        write (6,*) 'flds(1,1,1,ib) = ', flds(1,1,1,ib)
c        write (6,*) 'flds(1,1,1,ic) = ', flds(1,1,1,ic)
c        stop
      else if(ido .eq. iop_divide) then
        ib = iArg_tof(itof,1)
        ic = iArg_tof(itof,2)
        div_const = fVal_tof(itof)
        if(ib .eq. 0) then
          do iz = iz0, iz1
          do iy = iy0, iy1
          do ix = ix0, ix1
            flds(ix,iy,iz,ia) = div_const / flds(ix,iy,iz,ic)
          enddo
          enddo
          enddo
        else if(ic .eq. 0) then
          do iz = iz0, iz1
          do iy = iy0, iy1
          do ix = ix0, ix1
            flds(ix,iy,iz,ia) = flds(ix,iy,iz,ib) / div_const
          enddo
          enddo
          enddo
        else
          nDim1 = iArg_dim(itof,1)
          ia = iWrk_tof(itof) - 1
          ic = iArg_tof(itof,2)
          do ib = iArg_tof(itof,1), iArg_tof(itof,1)+nDim1-1
            ia = ia + 1
            do iz = iz0, iz1
            do iy = iy0, iy1
            do ix = ix0, ix1
              flds(ix,iy,iz,ia) = flds(ix,iy,iz,ib) / flds(ix,iy,iz,ic)
            enddo
            enddo
            enddo
          enddo
        endif
      else if(ido .eq. iop_setval) then
        add_const = fVal_tof(itof)
        if(ic .eq. 0) ic = ib
        do iz = iz0, iz1
        do iy = iy0, iy1
        do ix = ix0, ix1
          flds(ix,iy,iz,ia) = add_const
        enddo
        enddo
        enddo
      else if(ido .eq. iop_coord1) then
        add_const = fVal_tof(itof)
        if(ic .eq. 0) ic = ib
        do iz = iz0, iz1
        do iy = iy0, iy1
        do ix = ix0, ix1
        xx = dx * (float(ix + ioff_wb(1,i_wb)) - 0.5) + flimits(1)
          flds(ix,iy,iz,ia) = xx - add_const
        enddo
        enddo
        enddo
      else if(ido .eq. iop_coord2) then
        add_const = fVal_tof(itof)
        if(ic .eq. 0) ic = ib
        do iz = iz0, iz1
        do iy = iy0, iy1
        xx = dy * (float(iy + ioff_wb(2,i_wb)) - 0.5) + flimits(3)
        do ix = ix0, ix1
          flds(ix,iy,iz,ia) = xx - add_const
        enddo
        enddo
        enddo
      else if(ido .eq. iop_coord3) then
        add_const = fVal_tof(itof)
        if(ic .eq. 0) ic = ib
        do iz = iz0, iz1
        xx = dz * (float(iz + ioff_wb(3,i_wb)) - 0.5) + flimits(5)
        do iy = iy0, iy1
        do ix = ix0, ix1
          flds(ix,iy,iz,ia) = xx - add_const
        enddo
        enddo
        enddo

      endif

      return
      end

c--add
c     integer      nDim_tof(MAX_TOF)     ! Dimension of field (number of work arrays)
c     character*32 cOpr_tof(MAX_TOF)     ! Name of operaion
c     integer      iOpr_tof(MAX_TOF)     ! Operation number (for computed goto)
c     integer      nArg_tof(MAX_TOF)     ! Number of argument to operation
c     integer      nBdy_tof(3,MAX_TOF)   ! Number of argument to operation
c     character*32  cArg_tof(MAX_TOF,10) ! Names of argument variables
c     character*256 cFile_tof(MAX_TOF)   ! Names of dump file to read from
c     integer      iWrk_tof(MAX_TOF)     ! Work array assigned to field
c     integer      iArg_tof(MAX_TOF,10)  ! Work array assigned to argument variable
c     integer      iRda_tof(MAX_TOF)     ! Read array assigned to this field (to copy from)
c     integer      itof_Rda(MAX_TOF)     ! itof associated with read array (to find file, field, ...)


C ====================================================================72
      subroutine init_prof(profile, icounts)
      implicit none
      include 'postproc_info_e3d.h'
      integer icounts(maxprof), i
      real*8 profile(maxprof,maxcols)
       do i=1,maxprof
           profile(i,2) = 3.3
           profile(i,3) = 2.2
           profile(i,5) =  1.0d+33
           profile(i,6) = -1.0d+33
           icounts(i)   = 0
        enddo
       write (6,*) ' '
       write (6,*) 'In init: maxprof =  ', maxprof
       do i = 1, 3
         write (6,*) 'i,x,prof:', i, profile(i,2),profile(i,3)
        enddo

        return
        end

C ====================================================================72
      subroutine gen_box_limits(i,boxlimits)
      implicit none
      include 'postproc_info_e3d.h'
      include 'amr_e3d.h'
      integer i, k, n1
      real boxlimits(2,3), d1

      do k = 1, 3
        n1 = nboxmesh(k)
        d1 = boxsize(k) / float(iboxrefine(i)*n1)
        boxlimits(1,k) = boxoffset(k) + d1*float(iboxoff(k,i)   )
        boxlimits(2,k) = boxoffset(k) + d1*float(iboxoff(k,i)+n1)
      enddo

      return
      end

C ====================================================================72
      logical function boxs_overlap(i,j)
      implicit none
      include 'amr_e3d.h'
      integer i, j, k
      real xyzi(2,3), xyzj(2,3)
      logical b

      boxs_overlap = .false.
      if(i.lt.1 .or. j.lt.1) return

      call gen_box_limits(i,xyzi)
      call gen_box_limits(j,xyzj)

      b = .true.
      do k = 1, 3
       if(xyzi(1,k).ge.xyzj(2,k) .or. xyzj(1,k).ge.xyzi(2,k)) b=.false.
      enddo

      boxs_overlap = b

      return
      end

c      if(b) then
c        write (6,*) "boxs overlap: ", i,j
c        write (6,*) "xyzi_min:", (xyzi(1,k), k=1,3)
c        write (6,*) "xyzi_max:", (xyzi(2,k), k=1,3)
c        write (6,*) " "
c        write (6,*) "xyzj_min:", (xyzj(1,k), k=1,3)
c        write (6,*) "xyzj_max:", (xyzj(2,k), k=1,3)
c        write (6,*) " "
c        write (6,*) " "
c      endif

C ====================================================================72
      subroutine gen_box_list()
      implicit none
      include 'postproc_info_e3d.h'
      include 'amr_e3d.h'
      integer i, j, k, nri, nrk
      logical box_in_vol, box_is_ok, boxs_overlap
      real del, xyz(2,3)

      if(nbdryx .lt. nbdy_max(1) .or.
     1   nbdryy .lt. nbdy_max(2) .or.
     1   nbdryz .lt. nbdy_max(3)     ) then
        write (6,*) "PROBLEM: box boundarys not deep enough for formula"
        stop
      endif

      ! Generate list of boxs in bounding region closest to nrep_used
      ! asume boxs align & higher level boxs are covered by higher ones
      nboxsused = 0
      do i = 1, nboxs
        box_is_ok = .false.
        if(box_in_vol(i)) then
          box_is_ok = .true.
          j = 1
          do while(j .le. nboxsused .and. box_is_ok)
            k = ibox_list(j)
            if(boxs_overlap(i,k)) then
              nri = iboxrefine(i)
              nrk = iboxrefine(k)
              if(nrk.eq.nrep_used) box_is_ok = .false.
              if(nrk.lt.nrep_used .and. nri.lt.nrk) box_is_ok = .false.
              if(nrk.gt.nrep_used) then
                if(nri .gt. nrk      ) box_is_ok = .false.
                if(nri .lt. nrep_used) box_is_ok = .false.
              endif
              if(box_is_ok) ibox_list(j) = -1  ! take k out of list
            endif ! end of boxs overlap if-block
            j = j + 1
          enddo   ! end of while run over box list loop
        endif     ! end of box i is in volume if-block

        ! Collapse list if boxs were removed
        k = 0
        do j = 1, nboxsused
          do while(j+k.le.nboxsused .and. ibox_list(j+k) .lt. 0)
            k = k+1
          enddo
          if(k.gt.0 .and. k+j.le.nboxsused) ibox_list(j)=ibox_list(j+k)
        enddo
        nboxsused = nboxsused - k

        ! add i to end of list if it passed all tests (so far)
        if(box_is_ok) then
          nboxsused = nboxsused + 1
          ibox_list(nboxsused) = i
        endif
      enddo       ! if of i = 1, nboxs   loop

      if(maxboxs_db .gt. 0) then
         nboxsused = min(nboxsused, maxboxs_db)
         write (6,*) "nboxsused reset to:", nboxsused
      endif

      if(idebug .ge. 10) write (6,*) "nboxsused = ", nboxsused

      open(unit=32,file="box_list.dat",form="formatted")
      do i = 1, nboxsused
        j = ibox_list(i)
        if(idebug .ge. 10) write (6,*) "ibox,nrefine:", j,iboxrefine(j)
        del = dx / float(iboxrefine(j))
        call gen_box_limits(j,xyz)
        write (32,*) xyz(1,1)+del, xyz(1,2)+del
        write (32,*) xyz(2,1)-del, xyz(1,2)+del
        write (32,*) xyz(2,1)-del, xyz(2,2)-del
        write (32,*) xyz(1,1)+del, xyz(2,2)-del
        write (32,*) xyz(1,1)+del, xyz(1,2)+del
        write (32,*) " "
      enddo
      close(32)
c for amr

      return
      end
C ====================================================================72
      subroutine read_bricks_from_box_list()
      implicit none
      include 'postproc_info_e3d.h'
      include 'amr_e3d.h'
      integer i, j, k, ibox, irep, izoff, ioff
      integer nxread0, nyread0, nzread0, ix1,ix2,iy1,iy2,iz1,iz2
      integer ixyzmin(3), ixyzmax(3)

      nzread0 = min(nzfull, nfx0_max)

      if(idebug .ge. 10)
     1  write (6,998) "       ibox     irep         k",
     1                "   ixyzmin   ixyzmax",
     1                "   ioff_rb     nf_rb       1z1       iz2"
998   format(3a)
      n_rb = nboxsused
      do i= 1, nboxsused
        ibox = ibox_list(i)
        irep = iboxrefine(ibox)

        ! These are all INTIERIOR ranges excluding boundaries of any! kind
        call gen_mesh_range(irep,ix1,ix2,iy1,iy2,iz1,iz2)
        ixyzmin(1) = max(ix1, iboxoff(1,ibox)+1)
        ixyzmin(2) = max(iy1, iboxoff(2,ibox)+1)
        ixyzmin(3) = max(iz1, iboxoff(3,ibox)+1)
        ixyzmax(1) = min(ix2, iboxoff(1,ibox)+nboxmesh(1))
        ixyzmax(2) = min(iy2, iboxoff(2,ibox)+nboxmesh(2))
        ixyzmax(3) = min(iz2, iboxoff(3,ibox)+nboxmesh(3))

c        ixyzmin(1) = max(ix1, iboxoff(1,ibox)+nbdry+1)
c        ixyzmin(2) = max(iy1, iboxoff(2,ibox)+nbdry+1)
c        ixyzmin(3) = max(iz1, iboxoff(3,ibox)+nbdry+1)
c        ixyzmax(1) = min(ix2, iboxoff(1,ibox)+nbdry+nboxmesh(1))
c        ixyzmax(2) = min(iy2, iboxoff(2,ibox)+nbdry+nboxmesh(2))
c        ixyzmax(3) = min(iz2, iboxoff(3,ibox)+nbdry+nboxmesh(3))


c      write (6,*) "nbdry",nbdry
c        write (6,*) "iboxoff(3,ibox),nboxmesh(3)=",
c     1               iboxoff(3,ibox),nboxmesh(3)
c      nbdry           2
c   iboxoff(3,ibox),nboxmesh(3)=          64          32
c         ibox     irep         k   ixyzmin   ixyzmax   ioff_rb     nf_rb       1z1       iz2
c          17         2         3        67        65        66        -1        65        65

       ! Read full extent in X and Y (for efficient IO)
        do k = 1, 2
          if(k .eq. 1) then
            ioff_rb(k,i) = iboxoff(k,ibox) - nbdryx
            nf_rb(k,i)   = nboxmesh(k) + 2 * nbdryx
          else
            ioff_rb(k,i) = iboxoff(k,ibox) - nbdryy
            nf_rb(k,i)   = nboxmesh(k) + 2 * nbdryy
          endif
        enddo

        ! read trimmed range in Z
        k = 3
        ioff         = ixyzmin(k) - 1
        ioff_rb(k,i) = ioff - nbdy_max(3)
        nf_rb(k,i)   = (ixyzmax(k) - ioff) + 2 * nbdy_max(3)
      if(idebug .ge. 10)
     1    write (6,999) ibox, irep,k,ixyzmin(k),ixyzmax(k),
     1                  ioff_rb(k,i),nf_rb(k,i),iz1,iz2
999     format(12i10)
      enddo

        ! read trimmed range in Z
        ! offset & range must take irep into account
c        if(irep .eq. 1) then
c          izoff = iz1-1 - nbdry
c          ioff_rb(3,i) = izoff
c          nf_rb(3,i) = min(nzread0,izgblhi-izoff)+2*nbdy_max(3)
c        else
c          write (6,*) "PROBLEM: read_bricks_from_box_list: "
c          write (6,*) "         still need to implment for irpe>1"
c          stop
c        endif

      if(idebug .ge. 10)
     1   write (6,*) "Finished read_bricks_from_box_list"

      return
      end

c        ioff_rb(1,i) = 0
c        ioff_rb(2,i) = 0
c        izoff = izgbllo-1 + nzread0*(i-1) - nbdy_max(3)
c        ioff_rb(3,i) = izoff
c        nf_rb(1,i) = nxfull
c        nf_rb(2,i) = nyfull
c        nf_rb(3,i) = min(nzread0,izgblhi-izoff)+2*nbdy_max(3)

C ====================================================================72
      subroutine gen_read_bricks()
      implicit none
      include 'postproc_info_e3d.h'
      include 'amr_e3d.h'
      integer i, j
      integer nxread0, nyread0, nzread0, izoff

      ! Generate default brick list
       call get_full_dimensions(nxfull, nyfull, nzfull)

      if(nrep_used .gt. 1) then
        call gen_box_list()
        call read_bricks_from_box_list()
        return
      endif

      ! ***************************************************************
      ! *                    GENERATE: READ BRICKS                    *
      ! ***************************************************************
      ! Read case 1: read arrays span full X & Y
      ! Internal read slab size (in z)
      nxread0 = nxfull
      nyread0 = nyfull
      nzread0 = ntslab    ! itwas32:      nzread0 = 32
c      write (6,*) "nzread0,ntslab", nzread0,ntslab

      n_rb = 0
      do i = izgbllo, izgblhi, nzread0
        n_rb = n_rb + 1
        ioff_rb(1,n_rb) = 0
        ioff_rb(2,n_rb) = 0
        izoff = izgbllo-1 + nzread0*(n_rb-1) - nbdy_max(3)
        ioff_rb(3,n_rb) = izoff

        nf_rb(1,n_rb) = nxfull
        nf_rb(2,n_rb) = nyfull
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
        nf_rb(3,n_rb) = min(nzread0,izgblhi-izoff)+2*nbdy_max(3)
      enddo
 
      if(idebug .gt. 1000) then
        write (6,*) ' '
        write (6,*) 'Read brick full ranges '
        write (6,*) "Brick     Xmin Xmax   Ymin Ymax   Zmin Zmax"
        do i = 1, n_rb
          write (6,999) i,(ioff_rb(j,i)+1,ioff_rb(j,i)+nf_rb(j,i),j=1,3)
999       format(i6, "  ", 3("  ",2i5))
        enddo

c         call gen_work_bricks(3)
      endif

      return
      end

C ====================================================================72
      subroutine gen_work_bricks(i_rb)
      implicit none
      include 'postproc_info_e3d.h'
      include 'amr_e3d.h'
      integer i, j, k, m, i_rb
      integer ixoff, iyoff, izoff, nxw0, nyw0
      integer nbx,nby,nbz, nmaxx,nmaxy,nmaxz, nwx0,nwy0,nwz0

      ! ***************************************************************
      ! *      GENERATE: WORK BRICKS FOR A GIVEN READ BRICK           *
      ! ***************************************************************
      ! work case 1: work arrays interiors are cubes
      ! Internal work size in x, y, & z matches read internal size
      ! assume read arrays span X & Y
      if(nboxsused .le. 0) then
      if(nf_rb(1,i_rb) .ne. nxfull .or. nf_rb(2,i_rb) .ne. nyfull) then
        write (6,*) "PROBLEM in gen_work_briks"
        write (6,*) "  Case not supported: read brick does not span x&y"
        stop
      endif
      endif

      if(nboxsused .le. 0) then
        ! For level 1 only (full span in XY)
        ! TODO: make sure work bricks fit in cache & long in X
c xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
c        write (6,*) "nf_rb(3,i_rb),nbdy_max(3)=",
c     1               nf_rb(3,i_rb),nbdy_max(3)
        nwz0 = nf_rb(3,i_rb) - 2*nbdy_max(3)
        nyw0 = ntslab    ! itwa32: nyw0 = 32
        nxw0 = min(nxfull, nfx0_max)
    
        n_wb = 0
        izoff = ioff_rb(3,i_rb) + nbdy_max(3)
        do iyoff = 0, nyfull-1, nyw0
        do ixoff = 0, nxfull-1, nxw0
          n_wb = n_wb + 1
          ioff_wb(1,n_wb) = ixoff
          ioff_wb(2,n_wb) = iyoff
          ioff_wb(3,n_wb) = izoff

          n0_wb(1,n_wb) = min(nxw0,nxfull-ioff_wb(1,n_wb))
          n0_wb(2,n_wb) = min(nyw0,nyfull-ioff_wb(2,n_wb))
          n0_wb(3,n_wb) = nwz0
        enddo
        enddo
      else
        ! For AMR  (see read_bricks_from_box_list)
        ! For AMR: mesh offsets & sizes are ALWAYS in "units" of the box level
        
        ! work brick interior sizes
c        nxw0 = min(nxfull, nfx0_max)
        nwx0 = nf_rb(1,i_rb) - 2*nbdryx
        nwy0 = ntslab     !  itwas32:      nwy0 = 32
        nwz0 = ntslab     !  itwas32:      nwz0 = 32

        ! READ (!) brick boundaries
        nbx = nbdryx
        nby = nbdryy
        nbz = nbdy_max(3)

        nmaxx = ioff_rb(1,i_rb) +  nf_rb(1,i_rb) - nbx
        nmaxy = ioff_rb(2,i_rb) +  nf_rb(2,i_rb) - nby
        nmaxz = ioff_rb(3,i_rb) +  nf_rb(3,i_rb) - nbz
    
        n_wb = 0
        do izoff = ioff_rb(3,i_rb)+nbz, nmaxz-1, nwz0
        do iyoff = ioff_rb(2,i_rb)+nby, nmaxy-1, nwy0
        do ixoff = ioff_rb(1,i_rb)+nbx, nmaxx-1, nwx0
          n_wb = n_wb + 1
          ioff_wb(1,n_wb) = ixoff
          ioff_wb(2,n_wb) = iyoff
          ioff_wb(3,n_wb) = izoff

          n0_wb(1,n_wb) = min(nwx0,nmaxx-ixoff)
          n0_wb(2,n_wb) = min(nwy0,nmaxy-iyoff)
          n0_wb(3,n_wb) = min(nwz0,nmaxz-izoff)
      enddo
      enddo
      enddo

      endif
 
      if(idebug .gt. 100) then
        write (6,*) ' '
        write (6,*) 'work brick interior ranges for read brick ', i_rb
        write (6,*) 'number of work bricks: n_wb = ', n_wb
        write (6,*) "WorkBirck Xmin Xmax   Ymin Ymax   Zmin Zmax"
        do i = 1, n_wb
          write (6,999) i,(ioff_wb(j,i)+1,ioff_wb(j,i)+n0_wb(j,i),j=1,3)
999       format(i6, "  ", 3("  ",2i5))
        enddo
      endif

      return
      end
c        write (6,*) "-----------------------------"
c        write (6,998) "i_rb   :", i_rb
c        write (6,998) "ioff_rb:", (ioff_rb(k,i_rb),k=1,3)
c        write (6,998) "nf_rb  :", (nf_rb(k,i_rb),k=1,3)
c        write (6,998) "nb[xyz]:", nbx, nby, nbz
c998     format(a,3i8)
c        write (6,*) "-----------------------------"

C ====================================================================72
      logical function b_has_char(cstr,char)
      implicit none
      character*(*) cstr
      character*1 char
      integer n, i
      logical isin
      isin = .false.
      n = len(cstr)
      do i = 1, n
        if(cstr(i:i) .eq. char) isin = .true.
      enddo
      b_has_char = isin
      return
      end

C ====================================================================72
      integer function ilast_numeric_char_in_str(str, ic1)
      implicit none
      character*(*) str
      integer ic1, ic2, n, i, isnum
      n = len(trim(str))

      i = ic1
      ic2 =  -1
10    isnum = 0
      if(str(i:i) .eq. ' ') isnum = 1
      if(str(i:i) .eq. '-') isnum = 1
      if(str(i:i) .eq. '0') isnum = 1
      if(str(i:i) .eq. '.') isnum = 1
      if(str(i:i) .ge. '1'  .and.  str(i:i) .le. '9') isnum = 1
      if(isnum .eq. 1) ic2 = i
      i = i + 1
      if(isnum .gt. 0  .and.  i .le. n) go to 10

      ilast_numeric_char_in_str = ic2

      return
      end

C ====================================================================72
      subroutine get_par_val(cPar,cLine,fDefault,fVal)
      implicit none
      character*(*) cPar, cLine
      integer ic, ic1, ic2
      integer ilast_numeric_char_in_str, imatch_str_in_str
      real fDefault, fVal

      fVal = fDefault

ccc      write (6,*) ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>"
ccc      write (6,*) "cPar_tof=",trim(cLine)
      ic = imatch_str_in_str(cPar,cLine)
ccc      write (6,999) cPar, ic
ccc999   format("cPar=", a, "     ic=",i)
      if(ic .gt. 0) then
        ic1 = ic + len(cPar)
        ic2 = ilast_numeric_char_in_str(cLine, ic1)
cc        write (6,*) "ic1,ic2 = ", ic1, ic2
        if(ic1 .le. ic2) read (cLine(ic1:ic2),*) fVal
cc        if(ic1 .le. ic2) write (6,*) "cNumField=", cLine(ic1:ic2)
      endif
ccc      write (6,*) "fVal = ", fVal
ccc      write (6,*) "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<"

      return
      end
      
      
C ====================================================================72
      subroutine do_sectional_args()
      implicit none
      include 'postproc_info_e3d.h'
      integer itof, iarg, i1,i2,  i_have_sectional_wildcards, i, ic
      integer imatch_filter
      logical b_has_char, isin
      character*32 ctmpargs(1000), cFilt
      real rad_low, rad_high, rad, crg

      if(i_have_sectional_wildcards() .eq. 0) return
      call get_section_info()

      do itof = 1, ntof
        isin = .false.
        do iarg = 1, nArg_tof(itof)
          if(b_has_char(cArg_tof(itof,iarg),'?')) isin = .true.
        enddo

        if(isin) then        ! redo entire list of arguments

          ! set/get limits on sections
          call get_par_val("rad_low=" ,cPar_tof(itof),    0.0, rad_low)
          call get_par_val("rad_high=",cPar_tof(itof),1.0e+12,rad_high)

          do iarg = 1, nArg_tof(itof)
           ctmpargs(iarg) = cArg_tof(itof,iarg)
          enddo

          ic = 0
          do iarg = 1, nArg_tof(itof)
            if(b_has_char(cArg_tof(itof,iarg),'?')) then
              ! Insert expanded list of sectional fields
              cFilt = cArg_tof(itof,iarg)
              ic = ic - 1     ! overwrite wildcard argument
              do i = 1, ntof
                if(imatch_filter(cFld_tof(i), cFilt) .eq. 1) then
                  call get_sect_vals_from_string(cFld_tof(i),rad,crg)
                  if(rad_low.le.rad .and. rad.lt.rad_high) then
cccc                    write (6,999) cFld_tof(i)(1:12),rad,rad_low,rad_high
cccc999                 format("cFld_tof=", a12, "  rads:", 1p3e12.4)
                    ic = ic + 1
                    cArg_tof(itof,iarg+ic) = cFld_tof(i)
cccc                 else
cccc                    write (6,999) cFld_tof(i)(1:12),rad
                 endif

                endif
              enddo
            else
              cArg_tof(itof,iarg+ic) = ctmpargs(iarg)
            endif
          enddo   ! iarg = 1, nArg_tof(itof) loop
          nArg_tof(itof) = nArg_tof(itof) + ic
          nDim_tof(itof) = nArg_tof(itof)

        endif

      enddo     ! itof = 1, ntof loop

      return
      end
      
C ====================================================================72
      integer function i_have_sectional_wildcards()
      implicit none
      include 'postproc_info_e3d.h'
      integer isin, itof, iarg, i1,i2, i
      isin = 0
      do itof = 1, ntof
        do iarg = 1, nArg_tof(itof)
          i1 = 1
          i2 = len(trim(cArg_tof(itof,iarg)))
          do i = i1, i2
            if(cArg_tof(itof,iarg)(i:i) .eq. "?") isin = 1
          enddo
        enddo
      enddo
      i_have_sectional_wildcards = isin
      return
      end
C ====================================================================72
      integer function in_tof_args(cstr)
      implicit none
      character*(*) cstr
      include 'postproc_info_e3d.h'
      integer isin, itof, iarg, i1,i2, ilenstr

      ilenstr = len(trim(cstr))
c      write (6,*) "cstr,ilenstr=",trim(cstr),ilenstr
      isin = 0
      do itof = 1, ntof
        do iarg = 1, nArg_tof(itof)
          i2 = len(trim(cArg_tof(itof,iarg)))
          i1 = max(1, 1 + i2 - ilenstr)
c           write (6,*) "cArg_tof(itof,iarg)=",trim(cArg_tof(itof,iarg))
c           write (6,*) "i1,i2=",i1,i2
          if(cArg_tof(itof,iarg)(i1:i2) .eq. trim(cstr)) isin = 1
        enddo
      enddo
      
c       write (6,*) "ntof = ", ntof
c       write (6,*) "isin = ", isin
c       write (6,*) "cstr = ", trim(cstr)
c       write (6,*) "nArg_tof(ntof) = ", nArg_tof(ntof)
c       do iarg = 1, nArg_tof(ntof)
c         i2 = len(trim(cArg_tof(itof,iarg)))
c         i1 = 1 + i2 - ilenstr
c         write (6,*) "iArg,cArg,cArg_end = ",
c     1                            iarg, trim(cArg_tof(ntof,iarg)),
c     1                            cArg_tof(ntof,iarg)(i1:i2)
c       enddo

      in_tof_args = isin

      return
      end

C ====================================================================72
      subroutine read_dump_times()
      implicit none
      include 'postproc_info_e3d.h'
      real time
      integer i, j

      if(nDumps .gt. 0) return    ! dump time have already been read in

      do i = 1, MAX_DUMPS
      do j = 1, 256
        cDumps(i)(j:j) = ' '
      enddo
      enddo

      i =  -1000
      open(unit=17, file='dump_times', form='formatted',
     1     status='old', err=210)
      i = 0
10    i = i + 1
      read (17,*,end=200) cDumps(i), fTimes(i)
      goto 10

200   close(17)
210   nDumps = i - 1

      ! Error handling
      if(nDumps .lt. 0) then
        c_add_tof_problem(1:25) = "Fail3d to open dump_times"
c      write (6,*) "PROBLEM: need to generate the dump_times file"
c      write (6,*) " RUN:  'ee times'; and try again"
c      stop
      endif
      if(nDumps .eq. 0) then
        c_add_tof_problem(1:27) = "Got 0 dumps from dump_times"
      endif
      if(nDumps .le. 0) then
        call get_time(time)
        nDumps = 1
        cDumps(1) = cdumpfile0
        fTimes(1) = time
      endif

      if(idebug .gt. 100) then
        write (6,*) 'nDumps = ', nDumps
        do i = 1, nDumps
          write (6,991) trim(cDumps(i)), fTimes(i)
991       format("   ", a, f10.5)
        enddo
      endif

      return
      end

C ====================================================================72
      subroutine print_values()
      implicit none
      include 'postproc_info_e3d.h'
      integer i

      write (6,*) " "
      write (6,*) "List of named constant/scalar values:"
      do i = 1, nValues
         write (6,888) trim(cValues(i)), fValues(i)
888      format(a20,"  ", 1pe15.6)
      enddo

      return
      end

C ====================================================================72
      subroutine get_section_info()
      implicit none
      include 'postproc_info_e3d.h'
      character*256 cline, ckey
      integer i, j, iStrLen, n, imatch_filter, ii,  isec_value
      real fval(10), f0, f1, value_sect

      n = 0
      cPsigns    = "+p                              "
      cMsigns    = "-_m                             "
      cSecFilter = "                                "
      call setblank(cSectCbar)
      cSectCbar(1:5) = "wgrb0"
      do i = 1, MAXSECDIMS
        
      enddo

      ! read sectional config file
      open(unit=17,file="sections.config",form="formatted",
     1             status="old",err=900)

10    continue
      do i = 1, 256
        cline(i:i) = ' '
      enddo
      read (17,990,end=90) cline
990   format(a256)
      if(cline(1:10).eq.'filter    ') read (cline,*) ckey,cSecFilter
      if(cline(1:10).eq.'title     ') read (cline,*) ckey,cSecTitle
      if(cline(1:10).eq.'cbar      ') read (cline,*) ckey,cSectCbar
      if(cline(1:10).eq.'field     ') then
        n = n + 1
        read (cline,*) ckey,cSecFormat(n),cSecName(n),cSecUnit(n),
     1     cSecMap(n), fSecF0(n), fSecF1(n), cSecLoc(n)
        nSecFields = n
      endif
      go to 10
90    continue
      close(17)

c --------- debug ------------------
c      write (6,*) "Section filter: ", trim(cSecFilter)
c      write (6,*) "Number or section dimensions: ", nSecFields
c      write (6,*) "ntof =  ", ntof
c      write (6,*) "maps: ", (trim(cSecMap(j)), j=1,nSecFields)
c      write (6,991) (cSecName(j)(1:10), j=1,nSecFields)
c991   format("Variable         ", 10a10)
c      do i = 1, ntof
c        if(imatch_filter(cFld_tof(i), cSecFilter) .eq. 1) then
c          do j=1,nSecFields
c            ii = isec_value(cFld_tof(i),cSecFormat(j),cmsigns)
c            if(j .eq. 1) ii = ii - 1
c             fval(j) = value_sect(float(ii),j)
c          enddo
c          write (6,992) cFld_tof(i)(1:12), (fval(j),j=1,nSecFields)
c992       format(a12, 10f10.3)
c        endif
c      enddo
c      write (6,*) " "
c      write (6,*) " "
c      write (6,*) "Title: ", trim(cSecTitle)
c --------- debug ------------------
      return

900   write (6,*) "Could not open sections.config"
      stop
      end

C ====================================================================72
      integer function isec_value(cFld,cFormat,cmsigns)
      implicit none
      character*(*) cFld,cFormat,cmsigns
      integer nmsigns, i, nlen, ival, j, isgn,jval,ifac, k
      nmsigns = len(trim(cmsigns))
      nlen = len(trim(cformat))
      isgn = 1
      ival = 0
      do j=1,nlen
        if(cFormat(j:j) .eq. "?") then
          jval = ichar(cFld(j:j)) - ichar("0")
          if(0 .le. jval  .and.  jval .le. 9) ival = 10 * ival + jval
          do k = 1,nmsigns
            if(cFld(j:j) .eq. cmsigns(k:k)) isgn = -1
          enddo
        endif
      enddo

      isec_value = isgn * ival
      return
      end

C ====================================================================72
      integer function imatch_str_in_str(s1,s2)
      implicit none
      character*(*) s1, s2
      integer n1, n2, i, imatch

      n1 = len(trim(s1))
      n2 = len(trim(s2))
      imatch_str_in_str = -1
      if(n2 .lt. n1) return

      do i = 1, 1+n2-n1
        if(s1(1:n1) .eq. s2(i:n1+i-1)) imatch_str_in_str = i
      enddo

      return
      end
      
C ====================================================================72
      integer function imatch_filter(cStr, cFilter)
      implicit none
      character*(*) cStr, cFilter
      integer nlen, i, imatch
      nlen = len(trim(cFilter))
      imatch = 1
      do i=1,nlen
        if(cFilter(i:i) .ne. '?' .and.
     1     cfilter(i:i) .ne. cStr(i:i)) imatch = 0
      enddo
        if(cstr(nlen+1:nlen+1) .ne. ' ') imatch = 0
      imatch_filter = imatch
      return
      end
C ====================================================================72
      subroutine add_tofs(cStr,ntof0)
      implicit none
      include 'postproc_info_e3d.h'
      character*(*) cStr
      character*32 cdumpfile1
      integer ntor0, itof, itof0, itime, i, ntof0, idtime
      integer im, imatch_tof, j, isin, i0, i1, i2, ilenstr
      real val
      integer itime1

      ! Read list of dumps (with thier times) that are available
      call read_dump_times()

      idtime = 0
      if(trim(cStr) .eq. "@^") idtime = 1
      if(trim(cStr) .eq. "@v") idtime = -1

      if(idtime .ne. 0) then
        ! Find next or previous dump in sequence
        ! Find currnt dump in this list
        itime = -1
        do i = 1, nDumps
          if(trim(cdumpfile0) .eq. trim(cDumps(i))) itime = i
        enddo

         ! Set next time and evaluate & set dTime constant
         if(itime .lt. 0 .and. c_add_tof_problem .eq. "none") then
           write (6,*) "Problem: current dump is not in dump_times"
           stop
         endif
         itime1 = itime + idtime
         if((itime1 .lt. 1)  .or.  (itime1 .gt. nDumps)) then
            cDumps(nDumps+1)(1:7) = "no_file"
            fTimes(nDumps+1)      = 0.0
            if(c_add_tof_problem .eq. "none") then
             c_add_tof_problem(1:30) = "No matching dump in dump_times"
           endif
         endif
         call add_dt_vals(idTime, itime, itime1)

        if(idebug .gt. 100) then
          write (6,*) ' '
          write (6,*) 'Current  Dump = ', trim(cDumpFile0)
          if(idtime .gt. 0)
     1      write (6,*) 'Next     Dump = ', trim(cDumps(itime1))
          if(idtime .lt. 0)
     1      write (6,*) 'Previous Dump = ', trim(cDumps(itime1))
          call print_values()
        endif

      else if(trim(cStr) .eq. "@o") then

        ! Find dump with different root name and same sequence number
        itime1 = -1
        i0 = len(trim(cdumpfile0))
        do i = 1, nDumps
          i1 = len(trim(cDumps(i)))
          if(trim(cdumpfile0) .ne. trim(cDumps(i)) .and.
     1       cdumpfile0(i0-8:i0) .eq. cdumps(i)(i1-8:i1) ) itime1 = i
c           write (6,*) "-------------------------------------------"
c           write (6,*) "i,itime1 = ", i,itime1
c           write (6,*) "i0,i1 = ", i0,i1
c           write (6,*) trim(cdumpfile0)," -- ",trim(cDumps(i))
c           write (6,*) cdumpfile0(i0-8:i0)," -- ", cdumps(i)(i1-8:i1)
c           write (6,*) "-------------------------------------------"
        enddo

        if(itime1 .lt. 0) then
          itime1 = nDumps+1
          cDumps(itime1)(1:7) = "no_file"
          fTimes(itime1)      = 0.0
          if(c_add_tof_problem .eq. "none") then
            c_add_tof_problem(1:30) = "No matching dump in dump_times"
          endif
c          write (6,*) " "
c          write (6,*) "Problem: '@o'  in formula.e3d requests"
c          write (6,*) "a dump with the same sequence number as"
c          write (6,*) "   ", trim(cdumpfile0)
c          write (6,*) "form a different run"
c          write (6,*) "Such a dump is missing in dump_times"
c          write (6,*) "Either collect dumps from diffeent runs in"
c          write (6,*) "the same directory and run"
c          write (6,*) "    ee times"
c          write (6,*) "or remove the field with @o from formulas.e3d"
c          write (6,*) " "
c          stop
        endif

        if(idebug .gt. 100) then
          write (6,*) ' '
          write (6,*) 'Current Dump = ', trim(cDumpFile0)
          write (6,*) 'Other   Dump = ', trim(cDumps(itime1))
        endif

      else
        write (6,*) " "
        write (6,*) "Problem: ", trim(cStr)
        write (6,*) "in formula.e3d is not a recognized option for"
        write (6,*) "specifying input from a different data file"
        write (6,*) "listed in dump times"
        write (6,*) " use ...@^ to specify the next time"
        write (6,*) " use ...@v to specify the previous time"
        write (6,*) " use ...@o to specify a different run"
        write (6,*) " "
        stop
      endif

      ! Add a copy of all formulas for the other cdump specified by cStr
      ilenstr = len(trim(cStr))
      itof = ntof
      do itof0 = 1, ntof0

        ! Only make a copy of the formula if it does NOT refer to the new state
        isin = 0
        do i = 1, nArg_tof(itof0)
          i2 = len(trim(cArg_tof(itof0,i)))
          i1 = max(1 + i2 - ilenstr, 1)
          if(cArg_tof(itof0,i)(i1:i2) .eq. trim(cStr)) isin = 1
        enddo

cccccc        if(isin .eq. 0) then   ! this lead to problems with missing xxx@^ int tof
        itof = itof + 1
        call add_strs(cFld_tof(itof), cFld_tof(itof0), cStr)
        nDim_tof(itof) = nDim_tof(itof0)
        cPar_tof(itof) = cPar_tof(itof0) 
        cOpr_tof(itof) = cOpr_tof(itof0)
        iOpr_tof(itof) = iOpr_tof(itof0) 
        nBdy_tof(1,itof) = nBdy_tof(1,itof0)
        nBdy_tof(2,itof) = nBdy_tof(2,itof0)
        nBdy_tof(3,itof) = nBdy_tof(3,itof0)
        nArg_tof(itof) = nArg_tof(itof0)
        do i = 1, nArg_tof(itof)
          im = imatch_tof(cArg_tof(itof0,i),ntof0,cFld_tof)
          if(im .lt. 0 .or. isin .eq. 1) then
             cArg_tof(itof,i) = cArg_tof(itof0,i)
          else
            call add_strs(cArg_tof(itof,i), cArg_tof(itof0,i), cStr)
          endif
        enddo
        if(iOpr_tof(itof) .eq. iop_read) then
          cFile_tof(itof) = cDumps(itime1)
        endif
cccccc      endif  ! (isin .eq. 0) if-block
      enddo
      ntof = itof

c        iArg_tof(itof,10)
c        fVal_tof(itof)
c        iRda_tof(itof)
c        itof_Rda(itof)


      if(idebug .gt. 100) then
         write (6,*) "  "
        write (6,998)
998     format("Field       Dim  Operator  nArg  Source")
        do i = 1, ntof
          if(iOpr_tof(i) .eq. iop_read) then
            write (6,999) cFld_tof(i)(1:8), nDim_tof(i),
     1                cOpr_tof(i)(1:9), nArg_tof(i), trim(cFile_tof(i))
999        format(a8, "  ", i5,"  ",a9,i5, "  ", a)
         else
            write (6,997) cFld_tof(i)(1:8), nDim_tof(i),
     1            cOpr_tof(i)(1:9), nArg_tof(i),
     1        (trim(cArg_tof(i,j)),iArg_tof(i,j),j=1,nArg_tof(i))
997        format(a8, "  ", i5,"  ",a9,i5, 10("  ", a))
         endif
         enddo
      endif

      return
      end

C ====================================================================72
      subroutine add_dt_vals(idTime, itime, itime1)
      implicit none
      include 'postproc_info_e3d.h'
      integer idTime, itime, itime1, itp, itm, iv

      nValues = nValues + 1
      fValues(nValues) = fTimes(itime1) - fTimes(itime)
      call setblank(cValues(nValues))
      if(idTime .ge. 0) cValues(nValues)(1:5) = "dTime"
      if(idTime .lt. 0) cValues(nValues)(1:6) = "dTimem"

      itm = -1
      itp = -1
      do iv = 1, nValues
        if(trim(cValues(iv)) .eq. "dTime" ) itp = iv
        if(trim(cValues(iv)) .eq. "dTimem") itm = iv
      enddo
      if(itp .lt. 0  .or.  itm .lt. 0) return

      nValues = nValues + 1
      call setblank(cValues(nValues))
      cValues(nValues)(1:7) = "dTimepm"
      fValues(nValues) = fValues(itp) - fValues(itm)

      return
      end

C ====================================================================72
      subroutine add_strs(a,b,c)
      character*(*) a, b, c
      integer i, ilenbc

      do i = 1, len(a)
        a(i:i) = ' '
      enddo
      
       ilenbc = len(trim(b)) + len(trim(c))
       a(1:ilenbc) = trim(b) // trim(c)

      return
      end

C ====================================================================72
      subroutine count_ops(line, n)
      implicit none
      character*256 line
      integer i, im, n, nlen, i_is_op, im_is_op, ipmtd
      logical boper

      n = 0
      nlen = len(trim(line))
      im_is_op = 0
      do i = 2, nlen-1
        i_is_op = 0
        im = i - 1
c        if(line(i:i) .eq. '+' .or. line(i:i) .eq. '-' .or.
c     1     line(i:i) .eq. '*' .or. line(i:i) .eq. '/'     ) then
        if(ipmtd(line(i:i+1)) .eq. 1) then
          if(line(im:im) .ne. '*') n = n + 1
          i_is_op = 1
        endif
        if(line(i:i) .eq. '(' .and. line(im:im) .ne. ' ' .and.
     1     im_is_op .eq. 0) n =  n + 1
        im_is_op = i_is_op
      enddo

      return
      end

C ====================================================================72
      integer function itd(c)
      character*2 c
      itd = 0
      if(c(1:1) .eq. '*' .or.  c(1:1) .eq. '/') itd = 1
      return
      end

C ====================================================================72
      integer function ipm(c)
      character*2 c
      ipm = 0
      if(c(1:1) .eq. '+' .or.  c(1:1) .eq. '-') ipm = 1
      if(inumb(c(2:2)) .eq. 1) ipm = 0
      
      return
      end

C ====================================================================72
      integer function ipmtd(c)
      character*2 c
      integer ipm, itd
      ipmtd = max(ipm(c), itd(c))
      return
      end

C ====================================================================72
      integer function inumb(c)
      character*1 c
      inumb = 0
      if(ichar('0').le.ichar(c) .and. ichar(c).le.ichar('9')) inumb=1
      return
      end

C ====================================================================72
      subroutine find_epl(cf0, ieq, ip1, ip2, nlen)
      implicit none
      character*256 cf0
      integer i, nlen, ieq, ip1,ip2,npc

      nlen = len(trim(cf0))    ! formula length

      ! Find =
      ieq = -1
      do i = 2, nlen
        if(cf0(i:i) .eq. '=') ieq = i
      enddo
      if(ieq .lt. 0) then
        write (6,*) 'PROBLEM: no = sign in formula'
        write (6,*) '   ', trim(cf0)
        stop
      endif

      ! Find outter most (...)
      npc = 0
      ip1 = -2
      ip2 = -1
      do i = ieq, nlen
        if(cf0(i:i) .eq. '(') then
          if(ip1 .lt. 0) ip1 = i
          npc = npc + 1
        endif
        if(cf0(i:i) .eq. ')') then
          if(npc .eq. 1) ip2 = i
          if(npc .lt. 1) ip2 = -100
          npc = npc - 1
        endif
      enddo
      if(ip2 .le. ip1) then
        write (6,*) 'PROBLEM: parentheses do not balance'
        write (6,*) '   ', trim(cf0)
        stop
      endif

      return
      end

C ====================================================================72
      subroutine trim_paren(cf)
      implicit none
      character*256 cf
      integer i, ieq, ip1, ip2, nlen, is_empty

      call find_epl(cf, ieq, ip1, ip2, nlen)

      if(ip1 .lt. 0) return

      is_empty = 1
      if(ieq+1 .le. ip1-1) then
        do i=ieq+1, ip1-1
          if(cf(i:i) .ne. ' ') is_empty = 0
        enddo
      endif
      if(ip2+1 .le. nlen) then
        do i=ip2+1,nlen
          if(cf(i:i) .ne. ' ') is_empty = 0
        enddo
      endif

      if(is_empty .eq. 1) then
        cf(ip1:ip1) = ' '
        cf(ip2:ip2) = ' '
      endif

      return
      end

C ====================================================================72
      subroutine split_one(cf0, cf1, cf2)
      implicit none
      character*256 cf0, cf1, cf2
      character*10 cHead
      integer i, ieq, ip1, ip2, nlen
      integer nsplit, ipm, itd, ipmtd, iop, ic2, isplita
      common /split_forms/ nsplit

      cHead = "tmp000  = "
      cHead(4:4) = char(ichar('0') +     nsplit / 100    )
      cHead(5:5) = char(ichar('0') + mod(nsplit / 10, 10))
      cHead(6:6) = char(ichar('0') + mod(nsplit     , 10))
c      write (6,*) 'cHead = ', cHead

      do i = 1, 256
        cf1(i:i) = ' '
        cf2(i:i) = ' '
      enddo

      ! Trim off outter parentheses if nothing outside
      call trim_paren(cf0)

      call find_epl(cf0, ieq, ip1, ip2, nlen)

      ! find highest level operator outside of (...)
      iop = -1
      do i = ieq, nlen-1
        if(i .lt. ip1 .or. i .gt. ip2) then
          if(itd(cf0(i:i+1)) .eq. 1) iop=i
        endif
      enddo
      do i = ieq, nlen-1
        if(i .lt. ip1 .or. i .gt. ip2) then
          if(ipm(cf0(i:i+1)) .eq. 1) iop=i
        endif
      enddo
      if(iop .lt. 0) then
        write (6,*) 'PROBLEM: an operator is inside of  (...) in'
        write (6,*) '    ', trim(cf0)
        write (6,*) 'This is not supported yet.'
        stop
      endif

c      write (6,*) 'ip1,ip2,iop=',ip1,ip2,iop
      ! Old formula:    xxx = aaa _iop_ bbb 
      isplita = 0
      if(iop .gt. ip1) isplita = 1
      if(ip1 .lt. 0) then
        ! There's (...) in this formula, but more than one operator
        ! if any operator is in aaa, split it, otherwize bbb
        isplita = 0
        do i = ieq, iop-1
           if(ipmtd(cf0(i:i+1)) .eq. 1) isplita = 1
        enddo
      endif

      if(isplita .eq. 1) then
        ! there are 1 or more ops BEFORE IOP (i.e., in aaa)
        !    tmp = aaa
        !    xxx = tmp _op_ bbb
        cf1(1:10+iop-ieq) = cHead // cf0(ieq+1:iop-1)
        ic2 = ieq+9+nlen-iop
        cf2(1:ic2) = cf0(1:ieq) // " " // cHead(1:7) // cf0(iop:nlen)
      else
        ! there is (...) with ops in it AFTER iop (i.e., bbb has an op)
        !    tmp = bbb
        !    xxx = aaa _op_ tmp
        cf1(1:10+nlen-iop) = cHead // cf0(iop+1:nlen)
        cf2(1:iop+7) = cf0(1:iop) // " " // cHead(1:6)
      endif

      nsplit = nsplit + 1
      return
      end

C ====================================================================72
      subroutine split_formulas()
      implicit none
      include 'postproc_info_e3d.h'
      integer i, j, nfops
      character*32 cStr(100)
      character*256 cf1, cf2
      integer nsplit
      common /split_forms/ nsplit
      nsplit = 0

      j = 1
10    continue
      call count_ops(cFormulas(j), nfops)
      if(idebug .gt. 200) then
        write (6,999) nfops, trim(cFormulas(j))
999     format("nfops=",i4,"       cFormula=",a)
      endif

      if(nfops .gt. 1) then
        call split_one(cFormulas(j), cf1,cf2)

        nFormulas = nFormulas + 1
        if(nFormulas .ge. MAX_FORMULAS) then
          write (6,*) 'PROBLEM: in split_formulas:'
          write (6,*) '   exceeded max number of formulas'
          stop
        endif

        do i = nFormulas, j, -1
           cFormulas(i+1) = cFormulas(i)
        enddo
        cFormulas(j  ) = cf1
        cFormulas(j+1) = cf2
      else
        j = j + 1
      endif
      if(j .le. nFormulas) go to 10

      if(idebug .gt. 100) then
        do j = 1, nFormulas
          call count_ops(cFormulas(j), nfops)
          write (6,999) nfops, trim(cFormulas(j))
        enddo
      endif

      return
      end

C ====================================================================72
      subroutine gen_ops()
      implicit none
c      include 'header_info_e3d.h'
      include 'postproc_info_e3d.h'

      ! Scratch variables
      integer i, j, k, narg, nStr, istr_op, j_last
      character*32 cStr(100)
      character*256 cparam
      integer itof, imatch_tof, imatch_any, iDoNow, ido_tof, istillneed
      integer im, nwrk, iwrk, imdone, iarg, iop, icount, ioff, isindex
      integer iDid_tof(MAX_TOF)  ! (0/1): operation (has not/has) been done to produce field
      integer nbdy(3), jop,jDim,jArg, nbdy_sum, ntof0
      integer in_tof_args, ndvar, icc, is_a_number
      real val

      call get_vars_nams(MAX_TOF, ndvar, cFld_tof, cNam_tof)

      call split_formulas()

      ! Generate Full Table Of Formulas (w/o work array assigments)
      do i=1,ndvar
        nDim_tof(i)  = 1
        iOpr_tof(i)  = iop_read
        cOpr_tof(i)  = "read"
        nArg_tof(i)  = 0
        cFile_tof(i) = cDumpFile0 
      enddo
      ntof = ndvar

      do j = 1, nFormulas
        call parse_formula(cFormulas(j), nStr, cStr, istr_op, cparam)
        i = ndvar + j
        cFld_tof(i)  = cStr(1)
        cPar_tof(i) = cparam
        if(istr_op .gt. 0) then
           cOpr_tof(i) = cStr(istr_op)
          call get_op_prop(cOpr_tof(i),iop,nDim_tof(i),narg,nbdy,
     1                     icomp_tof(i))
c
c dhp: must modify nDim_tof(i) here
c
          iOpr_tof(i) = iop
          if(iOpr_tof(i) .lt. 0) then
            write (6,*) "Problem: do not recognize operation"
            write (6,*) "  ", trim(cStr(istr_op))
            write (6,*) "in formula"
            write (6,*) "   ", trim(cFormulas(j))
            stop
          endif
          nArg_tof(i)  = nstr - 2
          if(narg .ne. nArg_tof(i)) then
            write (6,*) "Problem: wrong number of arguments in formula"
            write (6,*) "   ", trim(cFormulas(j))
          endif
        else
          nArg_tof(i) = nstr - 1
          nDim_tof(i) = nstr - 1
          if(is_a_number(cstr(2),val) .eq. 0) then
            iOpr_tof(i) = iop_noop
            cOpr_tof(i)(1:32) = "no-op                           "
          else
            iOpr_tof(i) = iop_setval
            cOpr_tof(i)(1:32) = "set_val                         "
            fVal_tof(i) = val
          endif
        endif
        narg = 0
        do k = 2, nStr
          if(k .ne. istr_op) then
            narg = narg + 1
            cArg_tof(i,narg) = cStr(k)
          endif
        enddo
      enddo
      ntof = ndvar + nFormulas

       ! Handle sectional arguments as needed
       call do_sectional_args()

       ! add fields to tables for other dumps as needed
       ntof0 = ntof
       if(in_tof_args("@^") .gt. 0) call add_tofs("@^",ntof0)
       ntof0 = ntof
       if(in_tof_args("@v") .gt. 0) call add_tofs("@v",ntof0)
       ntof0 = ntof
       if(in_tof_args("@o") .gt. 0) call add_tofs("@o",ntof0)
c       if(ntof0 .gt. 0) stop

      do i = 1, ntof
        iWrk_tof(i)   = -1
        iArg_tof(i,1) = -1
        nDim_wrk(i) = 0   ! Default is no field starts here
      enddo
c      write (6,*) "Just before 1st call to imatch"
c      write (6,*) trim(cFieldLabel)
      im = imatch_tof(cField, ntof, cFld_tof)
      if(im .lt. 0) then
        write (6,*) "Problem: Target Field does not match any TOF field"
        write (6,*) "  ", trim(cField)
        stop
      endif

c      write (6,*) trim(cField), im
      iWrk_tof(im) = 1 + icomp_tof(im)
c
      call gen_dims()   ! Generate nDim_tof for math operations
c dhp: must modify nDim_tof before here    narg
c          iArg_tof
c
      nwrk = nDim_tof(im)
      if(iOpr_tof(im) .eq. iop_comp) nwrk = 0
10    continue
      imdone = 1
      do i = 1, ntof
        if(iWrk_tof(i) .gt. 0) then
          if(nArg_tof(i) .gt. 0 .and. iArg_tof(i,1) .lt. 0) then
c            write (6,*) '========================'
c            write (6,*) 'iFld=', i, "   cFld=",trim(cFld_tof(i))
            imdone = 0
            iwrk = nwrk + 1
            if(iOpr_tof(i) .eq. iop_noop) iwrk = iWrk_tof(i)
            do iarg = 1, nArg_tof(i)
              im = imatch_any(cArg_tof(i,iarg), val, ioff, isindex)
              if(im .lt. 0) then
                !!!! This if-block is For number/constantns values
                fVal_tof(i) = val
                iArg_tof(i,iarg) = 0
                if(iOpr_tof(i) .eq. iop_times) iOpr_tof(i) = iop_scale
                if(iOpr_tof(i) .eq. iop_plus) iOpr_tof(i) = iop_add
                if(iOpr_tof(i).eq.iop_divide  .and.  iarg.eq.2) then
                  fVal_tof(i) = 1.0 / val
                  iOpr_tof(i) = iop_scale
                endif
                if(iOpr_tof(i) .eq. iop_minus) then
                  if(iarg .eq. 1) then
                    fVal_tof(i) = val
                    iOpr_tof(i) = iop_cminus
                  else
                    fVal_tof(i) = -val
                    iOpr_tof(i) = iop_add
                  endif
                endif
c                write (6,*) 'i, val = ', i, val
              else
                !!!! This if-block is for fields in wrk array
                if(iWrk_tof(im) .lt. 0) then
                  iWrk_tof(im) = iwrk
c                   write (6,899) iarg,trim(cArg_tof(i,iarg)),iwrk
c899                format("iarg=",i3,"    carg=",a,"    iwrk=",i3)
                  iwrk = iwrk + nDim_tof(im)
                  if(iOpr_tof(i) .ne. iop_noop .or.
     1               iOpr_tof(im) .eq. iop_comp) nwrk=nwrk+nDim_tof(im)
                endif
                iArg_tof(i,iarg) = iWrk_tof(im) + ioff
c      write (6,*) "Have set iArg_tof for i,iarg=",i,iarg
c      write (6,*) "iWrk_tof(im)     = ", iWrk_tof(im)
c      write (6,*) "im               = ", im
c      write (6,*) "ioff             = ", ioff
c      write (6,*) "iArg_tof(i,iarg) = ", iArg_tof(i,iarg)
                !!!! This if-block is for fields in wrk array
              endif
            enddo
          endif
        endif
      enddo
      if(imdone .ne. 1) go to 10
      nfv = nwrk
       
      call gen_wrk_dims()   ! Generate nDim of work arrays

      if(idebug .gt. 100) then
         write (6,*) "  "
        write (6,998)
998     format("Field      iWrk  Dim iOpr nArg  DumpFile")
        do i = 1, ntof
          if(iOpr_tof(i) .eq. iop_read) then
            write (6,999) cFld_tof(i)(1:8), iWrk_tof(i), nDim_tof(i),
     1                   iOpr_tof(i), nArg_tof(i), trim(cFile_tof(i))
999        format(a8, "  ", 4i5, "  ", a)
         else
            write (6,997) cFld_tof(i)(1:8), iWrk_tof(i), nDim_tof(i),
     1            iOpr_tof(i), nArg_tof(i),
     1        (trim(cArg_tof(i,j)),iArg_tof(i,j),j=1,nArg_tof(i))
997        format(a8, "  ", 4i5, 10("  ", a, i3))
         endif
         enddo
         write (6,*) "Number of Work arrays: ", nfv
         write (6,*) "  "
      endif


      ! Generate list of operations
      do i = 1, MAX_TOF
        iDid_tof(i) = 0
      enddo
      nops = 0
      icount = 0
20    icount = icount + 1
      ido_tof = -1
      istillneed = -1
      ! Find 1st field in TOF which is needed, and is ready, to be produced
      do itof = ntof, 1, -1
         iDoNow = 1
         if(iWrk_tof(itof) .lt. 0) iDoNow = 0                ! Not needed
         if(iDid_tof(itof) .eq. 1) iDoNow = 0                ! Already done
         if(iDoNow .eq. 1) then
           istillneed = itof                                 ! This field still isneeded
           if(nArg_tof(itof) .gt. 0) then
             do iarg = 1, nArg_tof(itof)
                im = imatch_tof(cArg_tof(itof,iarg),ntof,cFld_tof) ! source field of argument
                if(im .gt. 0) then
                  if(iDid_tof(im) .le. 0) iDoNow = 0               ! Do not have argurment yet
                 endif
             enddo
           endif
         endif

         if(iDoNow .eq. 1) ido_tof = itof
      enddo
      if(ido_tof .gt. 0) then
        iDid_tof(ido_tof) = 1
        nops = nops + 1
        itof_ops(nops) = ido_tof
        if(idebug .gt. 100) then
        write (6,981) nops, ido_tof, trim(cFld_tof(ido_tof))
981     format("iop=",i4,"   ido_tof=",i3,"   Fld=",a)
        endif
      endif
      if(ido_tof .gt. 0  .and.  icount .lt. 10000) go to 20

      if(istillneed .gt. 0) then
        write (6,*) 'PROBLEM: no way to produce: ',
     1                 trim(cFld_tof(istillneed))
c         write (6,*) 'itof=',istillneed
c         do iarg = 1, nArg_tof(istillneed)
c           itof = iArg_tof(istillneed,iarg)
c           write (6,983) itof,trim(cFld_tof(itof)),iDid_tof(itof)
c983        format("itof,Fld,iDid:", i2,"  ",a,"  ",i2)
c         enddo
         stop
      endif

      ! Generate depth of bounaries needed at each operation
      if(idebug .ge. 31) write (6,*) "Generating boundary depth"
      itof = itof_ops(nops)                                 ! Table index of field produced by operation
      do j   = 1, 3
        nBdy_tof(j,itof) = nbdy_output
      enddo

      do iop = nops, 1, -1
        itof = itof_ops(iop)                                 ! Table index of field produced by operation
        call get_op_prop(cOpr_tof(itof),jop,jDim,jArg,nbdy,icc) ! nbdy needed for operation
        if(nArg_tof(itof) .gt. 0) then
          do iarg=1,narg_tof(itof)
c dhp_was            im = imatch_tof(cArg_tof(itof,iarg), ntof, cFld_tof)   ! source field
            im = imatch_any(cArg_tof(itof,iarg),val,ioff,isindex)  !  source field
            if(im .gt. 0) then
              do j=1,3
                nbdy_sum = nBdy_tof(j,itof) + nbdy(j)       ! total boundary depth needed for source field
                nBdy_tof(j,im) = max(nBdy_tof(j,im), nbdy_sum)
              enddo
            endif
          enddo
        endif
      enddo

cccccccccccccccccccccccccccccccccccc
c dhp check only
c      ! Eliminate no-ops from list of operations (itof_ops(*))
c      j = 0
c      j_last = 0
c      do i = 1, nops
c        if(iOpr_tof(itof_ops(i)) .ne. iop_noop) j = j+1
c        if(j .gt. j_last) then
c          if(i .gt. j) itof_ops(j) = itof_ops(i)
c          j_last = j
c        endif
c      enddo
c      nops = j_last
cccccccccccccccccccccccccccccccccccc

      ! Find max boundary depth
      do j=1,3
         nbdy_max(j) = 0
      enddo
      do iop = 1, nops
        i = itof_ops(iop)
        do j=1,3
           nbdy_max(j) = max(nbdy_max(j), nbdy_tof(j,i))
        enddo
      enddo

      ! Assigne read arrays
      nrda = 0
      do iop = 1, nops
        i = itof_ops(iop)
        iRda_tof(i) = 0
        if(iOpr_tof(i) .eq. iop_read) then
          nrda = nrda + 1
          iRda_tof(i) = nrda
          itof_rda(nrda) = i
        endif
      enddo

      ! Support 3-vector strip output
      !   scalar strip output: maxcols should be 2 (which is the defaul)
      !   vector strip output: maxcols should be nDim_tof(itof_ops(nops))
      if(maxcols .eq. 2) then
        if(maxcols .lt. 1+nDim_tof(itof_ops(nops))) then
          maxcols = 1+nDim_tof(itof_ops(nops))
          if(idebug .gt. 1000) write (6,*) 'reset: maxcols = ', maxcols
        endif
      endif

      if(idebug .gt. 10) then
c        write (6,*) 'Number of operations: ', nops
c        write (6,*) 'Number read fields  : ', nrda
c        write (6,*) 'Max boundary depth: : ', (nbdy_max(j),j=1,3)
c        write (6,*) "  "

        write (6,*) "  "
        write (6,993)  nops, nrda,  (nbdy_max(j),j=1,3)
993     format("# Nops=",i3,"   Nreads=",i3,"   Boundary depth=",3i2)

        write (6,990)
990     format("Field               iWrk  Dim Op    nArg  nBdyXYZ",
     1         "  Input(s)          iRda")

        do iop = 1, nops
          i = itof_ops(iop)
          if(iOpr_tof(i) .eq. iop_read) then
            write (6,991) cFld_tof(i)(1:17), iWrk_tof(i), nDim_tof(i),
     1         cOpr_tof(i)(1:9),nArg_tof(i),(nBdy_tof(j,i),j=1,3),
     1         trim(cFile_tof(i)), iRda_tof(i)
991        format(a17, "  ", 2i5," ",a5,i5, 3i3, "  ", a, i8)
         else
            write (6,992) cFld_tof(i)(1:17), iWrk_tof(i), nDim_tof(i),
     1         cOpr_tof(i)(1:9),nArg_tof(i),(nBdy_tof(j,i),j=1,3),
     1         (trim(cArg_tof(i,j)),iArg_tof(i,j),j=1,nArg_tof(i))
992        format(a17, "  ", 2i5," ",a5,i5, 3i3, 10("  ", a, i3))
         endif
         enddo
        
      endif

      call gen_output_label()

      return
      end

C ====================================================================72
      integer function imatch_tof(str, ntof, cFld_tof)
      character*(*) str
      integer ntof, i, imatch
      character*32 cFld_tof(ntof)     ! Name of field

      imatch = -1
      do i = 1, ntof
        if(trim(str) .eq. trim(cFld_tof(i))) imatch = i
      enddo

       imatch_tof = imatch

       return
       end

C ====================================================================72
      integer function is_a_number(str, fValue)
      implicit none
      character*(*) str
      real fvalue
      integer i

      is_a_number = 0
      i = 1
      if(str(i:i) .eq. '-') i = i + 1
      if(str(i:i) .eq. '.') i = i + 1
      if(str(i:i) .eq. '0') is_a_number = 1
      if(str(i:i) .eq. '1') is_a_number = 1
      if(str(i:i) .eq. '2') is_a_number = 1
      if(str(i:i) .eq. '3') is_a_number = 1
      if(str(i:i) .eq. '4') is_a_number = 1
      if(str(i:i) .eq. '5') is_a_number = 1
      if(str(i:i) .eq. '6') is_a_number = 1
      if(str(i:i) .eq. '7') is_a_number = 1
      if(str(i:i) .eq. '8') is_a_number = 1
      if(str(i:i) .eq. '9') is_a_number = 1
      if(is_a_number .eq. 1) read (str,*) fValue

      return
      end

C ====================================================================72
      integer function imatch_any(str, fValue, ioff,isindex)
      implicit none
      include 'postproc_info_e3d.h'
      character*(*) str
      character*32 str0
      integer i, imatch, is_a_number, ioff, imatch_tof, isindex
      real fValue

      isindex = 0
      ioff = 0
      imatch = 0
      do i = 1, ntof
        if(trim(str) .eq. trim(cFld_tof(i))) imatch = i
      enddo
      if(imatch .ne. 0) then
        imatch_any = imatch
        return
      endif

      ! check to see if it's a number or constant
      is_a_number = 0
      if(str(1:1) .eq. '-') is_a_number = 1
      if(str(1:1) .eq. '0') is_a_number = 1
      if(str(1:1) .eq. '1') is_a_number = 1
      if(str(1:1) .eq. '2') is_a_number = 1
      if(str(1:1) .eq. '3') is_a_number = 1
      if(str(1:1) .eq. '4') is_a_number = 1
      if(str(1:1) .eq. '5') is_a_number = 1
      if(str(1:1) .eq. '6') is_a_number = 1
      if(str(1:1) .eq. '7') is_a_number = 1
      if(str(1:1) .eq. '8') is_a_number = 1
      if(str(1:1) .eq. '9') is_a_number = 1
      if(is_a_number .eq. 1) then
         read (str,*) fValue
         imatch_any = -1
         return
      endif

      ! Check to see if it is a defined constant value
      if(nValues .gt. 0) then
        imatch = 0
        do i = 1, nValues
          if(trim(str) .eq. trim(cValues(i))) imatch = i
        enddo
        if(imatch .gt. 0) then
           fValue = fValues(imatch)
           imatch_any = -1
           return
        endif
      endif

      if(nIndexValues .gt. 0) then
        call match_indexes(str, str0, ioff)
        if(ioff .ge. 0) imatch = imatch_tof(str0, ntof, cFld_tof)
        if(imatch .gt. 0) then
          imatch_any = imatch
          isindex = 1
          return
        endif
      endif

      write (6,*) "Problem: string does not match anything"
      write (6,*) "  ", trim(str)
      if(idebug .eq. 29) call write_formulas()

      stop
      end

C ====================================================================72
      subroutine  match_indexes(str, str0, ioff)
      implicit none
      include 'postproc_info_e3d.h'
      character*(*) str, str0
      integer i, i1, n, idel, ioff, iv

      ioff = 0
      n = len(trim(str))
      if(n .lt. 2) return
      idel = 1
      do i = 1, n-1
        do iv = 1, nIndexValues
          if(str(i:i+1) .eq. cIndexValues(iv)) then
            ioff = ioff + idel * (iv-1)
            if(idel .eq. 1) i1 = i
            idel = idel * nIndexValues
          endif
        enddo
      enddo

      if(idel .eq. 1) then
        ioff = -1
       return
      endif
      
      call setblank(str0)
      str0(1:i1-1) = str(1:i1-1)
      return

      end

C ====================================================================72
      integer function  ndim_of_carg(carg)
      implicit none
      include 'postproc_info_e3d.h'
      character*(*) carg
      integer ndim, im, imatch_any, ioff, isindex
      real val

      im = imatch_any(carg, val, ioff, isindex)
cccccccccccccccc
c      if(ioff .ne. 0) then
c        write (6,*) "PROBLEM: ioff != 0 in ndim_of_carg"
c        write (6,*) "Field string = ", trim(carg)
c        stop
c      endif
cccccccccccccccc
      ndim = 1       ! Default: operand is a scalar
      if(im .gt. 0 .and. isindex .eq. 0) ndim = nDim_tof(im)
      ndim_of_carg = ndim

      return
      end

C ====================================================================72
      subroutine gen_dims()
      implicit none
      include 'postproc_info_e3d.h'
      logical notdone
      integer ic, itof, im, nDim1, nDim2, ndim_of_carg

      ic = 0
10    continue
      notdone = .false.
      do 100 itof = 1, ntof
        if(nDim_tof(itof) .eq. -111) then
          ndim1 = ndim_of_carg(cArg_tof(itof,1))
          if(nDim1 .gt. 0) nDim_tof(itof) = nDim1
          notdone = .true.
        endif
        if(nDim_tof(itof) .eq. -999) then
          ndim1 = ndim_of_carg(cArg_tof(itof,1))
          if(nDim1 .gt. 0) nDim_tof(itof) = nDim1 / 9
          notdone = .true.
        endif
        if(nDim_tof(itof) .eq. -333) then
          ndim1 = ndim_of_carg(cArg_tof(itof,1))
          if(nDim1 .gt. 0) nDim_tof(itof) = 3*nDim1
          notdone = .true.
        endif
        if(nDim_tof(itof) .eq. -1000) then
          ndim1 = ndim_of_carg(cArg_tof(itof,1))
          ndim2 = ndim_of_carg(cArg_tof(itof,2))
          if(nDim1 .gt. 0) then
            nDim_tof(itof) = (nDim1/3) * (nDim2/3)
            if(nDim_tof(itof) .lt. 1) then
              write (6,*) "Problem in gen_dims, dot product:"
              write (6,*) "Operatnds must be vectors or of higher rank"
              write (6,*) "But in formula for ", trim(cFld_tof(itof))
              write (6,*) "    dimension of 1st operand is ", nDim1
              write (6,*) "    dimension of 2nd operand is ", nDim2
              call print_tof()
              stop
            endif
          endif
          notdone = .true.
        endif

        if(nArg_tof(itof) .ne. 2) go to 100     ! only for +-*/ ops
        if(nDim_tof(itof) .ge. 0) go to 100     ! already set
        notdone = .true.
        ndim1 = ndim_of_carg(cArg_tof(itof,1))
        ndim2 = ndim_of_carg(cArg_tof(itof,2))
        if(ndim1 .le. 0  .or.  ndim2 .le. 0) go to 100  ! can't set yet

        if(nDim_tof(itof) .eq. -1) then   ! "a * b"
          nDim_tof(itof) = nDim1 * nDim2
        endif
        if(nDim_tof(itof) .eq. -3) then   ! "a + b" or "a - b"
          notdone = .true.
          nDim_tof(itof) = nDim1
          if(nDim2 .ne. nDim1) then
            write (6,*) "Problem in gen_dims:"
            write (6,*) "In + or -: dims of operands must the same."
            write (6,*) "However in formula for ", trim(cFld_tof(itof))
            write (6,*) "dimension of 1st operand is ", nDim1
            write (6,*) "dimension of 2nd operand is ", nDim2
            call print_tof()
            stop
          endif
        endif
        if(nDim_tof(itof) .eq. -2) then   ! "a / b"
          notdone = .true.
          nDim_tof(itof) = nDim1    ! and nDim2 had better be 1
          if(nDim2 .ne. 1) then
            write (6,*) "Problem in gen_dims:"
            write (6,*) "Divisor must be a scalar."
            write (6,*) "However in formula for ", trim(cFld_tof(itof))
            write (6,*) "dimension of 2nd operand is ", nDim2
            call print_tof()
            stop
          endif
        endif
100   continue    ! end itof = 1, ntof  do-loop
      ic = ic + 1
      if(ic .gt. ntof+1) then
        write (6,*) "PROBLEM In gen_dim: can not set dimensions"
        call print_tof()
        stop
      endif
      if(notdone) go to 10

      return
      end

C ====================================================================72
      subroutine gen_wrk_dims()
      implicit none
      include 'postproc_info_e3d.h'
      integer itof, iarg, im, ioff, isindex, imatch_any
      real val

      do itof = 1, ntof
        if(nDIm_tof(itof) .gt. 0 .and. iwrk_tof(itof) .gt. 0) then
          nDim_wrk(iwrk_tof(itof)) = nDim_tof(itof)
        endif
      enddo
      do itof = 1, ntof
c      write (6,*) "cStr_tof(itof): ", trim(cFld_tof(itof))
        do iarg = 1, nArg_tof(itof)
          if(iArg_tof(itof,iarg) .gt. 0) then
c      write (6,999) itof,iarg,iArg_tof(itof,iarg)
c999   format("itof,iarg,iArg_tof(itof,iarg)",3i5)
            ! only do this if iArg_tof has been assigned (= needed)
            iArg_dim(itof,iarg) = nDim_wrk(iArg_tof(itof,iarg))
            im = imatch_any(cArg_tof(itof,iarg), val, ioff, isindex)
            if(isindex .eq. 1) iArg_dim(itof,iarg) = 1
          endif
        enddo
      enddo

      return
      end
C ====================================================================72
      subroutine print_tof()
      implicit none
      include 'postproc_info_e3d.h'
      integer i, j

      write (6,*) "==============================================="
      write (6,*) "  "
      write (6,998)
998   format("Field      iWrk  Dim iOpr nArg  Inputs")
      do i = 1, ntof
        if(iOpr_tof(i) .eq. iop_read) then
          write (6,999) cFld_tof(i)(1:8), iWrk_tof(i), nDim_tof(i),
     1                 iOpr_tof(i), nArg_tof(i), trim(cFile_tof(i))
999      format(a8, "  ", 4i5, "  ", a)
       else
          write (6,997) cFld_tof(i)(1:8), iWrk_tof(i), nDim_tof(i),
     1          iOpr_tof(i), nArg_tof(i),
     1      (trim(cArg_tof(i,j)),iArg_tof(i,j),j=1,nArg_tof(i))
997      format(a8, "  ", 4i5, 10("  ", a, i3))
       endif
       enddo
      write (6,*) "==============================================="

      return
      end
C ====================================================================72
      subroutine split_line_on_char(line0, c, line1, line2)
      implicit none
      character*(*) line0, line1, line2
      character*1 c
      integer n0, n1, n2, ic, n, imatch_str_in_str
      call setblank(line1)
      call setblank(line2)
      n0 = len(trim(line0))    ! length of trimmined input string
      n1 = len(line1)           ! space availablei in line1
      n2 = len(line2)           ! space availablei in line2
      ic = imatch_str_in_str(c, line0)
      if(ic .le. 0) then
        ! no match: all goes in 1st output string
        n = min(n0, n1)
        line1(1:n) = line0(1:n)
        return
      endif
      if(ic.gt. 1) then
        n = min(ic-1, n1)
        line1(1:n) = line0(1:n)
      endif
      if(ic.lt.n0) then
        n = min(n0-ic, n2)
        line2(1:n) = line0(ic+1:ic+n) 

cc        write (6,*) "==============================================="
cc        write (6,*) "ic   = ", ic
cc        write (6,*) "line0: ", trim(line0)
cc        write (6,*) "line1: ", trim(line1)
cc        write (6,*) "line2: ", trim(line2)
cc        write (6,*) "==============================================="
      endif

      return
      end

C ====================================================================72
      subroutine parse_xyz_pars(str, n1, n2, n3)
      implicit none
      include 'postproc_info_e3d.h'
      character*(*) str
      integer n1, n2, n3, i1a,i2a,i3a,i4a, i1b,i2b,i3b
      integer ifind_in_str, is_next
      i1a = ifind_in_str(cxcoord(1:1), str)
      i2a = ifind_in_str(cycoord(1:1), str)
      i3a = ifind_in_str(czcoord(1:1), str)
      i4a = len(trim(str)) + 1
      i1b = is_next(i1a, i2a, i3a, i4a)
      i2b = is_next(i2a, i1a, i3a, i4a)
      i3b = is_next(i3a, i1a, i2a, i4a)

      n1 = 0
      n2 = 0
      n3 = 0
      if(i1a .gt. 0) read(str(i1a+1:i1b-1),*) n1
      if(i2a .gt. 0) read(str(i2a+1:i2b-1),*) n2
      if(i3a .gt. 0) read(str(i3a+1:i3b-1),*) n3

      return
      end

      integer function is_next(i0,ia,ib,ic)
      implicit none
      integer i0,ia,ib,ic, i1
      i1 = max(ia,ib,ic)
      if(i0 .lt. ia  .and.  ia .lt. i1) i1 = ia
      if(i0 .lt. ib  .and.  ib .lt. i1) i1 = ib
      if(i0 .lt. ic  .and.  ic .lt. i1) i1 = ic
      is_next = i1
      return
      end

C ====================================================================72
      subroutine parse_formula(line0, nStr, str, istr_op, cparam)
      implicit none
      character*256 line0, line, cparam
      integer nStr, i, j, nlength, istr, istr_op, is_operator,ipmtd
      character*32 str(100)
      integer is_in_string, was_in_string, jdebug

      jdebug = 0

      call split_line_on_char(line0, ";", line, cparam)
      
      do i = 1, 100
        str(i) = "                                "
      enddo

      nlength = len(trim(line))
c      write (6,999) nlength, trim(line)
c999   format("length=",i3, "     ", a)

      istr_op = -1
      istr = 0
      was_in_string = 0
      do i = 1, nlength
        is_in_string = 1
        is_operator = 0
        if(line(i:i) .eq. ' ' .or.  line(i:i) .eq. ',' .or.
     1     line(i:i) .eq. '=' .or.  line(i:i) .eq. '(' .or.
     1     line(i:i) .eq. ')'        ) is_in_string = 0

c        if(line(i:i) .eq. '+' .or.  line(i:i) .eq. '-' .or.
c     1     line(i:i) .eq. '*' .or.  line(i:i) .eq. '/' ) then
        if(i+1 .le. nlength) then
        if(ipmtd(line(i:i+1)) .eq. 1) then
               is_operator = 1
               is_in_string = 0
         endif
         endif

        if(was_in_string .eq. 0) then
           if(is_in_string .eq. 1) then
             istr = istr + 1
             j = 1
             str(istr)(j:j) = line(i:i)
             was_in_string = 1
           endif
         else
           if(is_in_string .eq. 1) then
             j = j + 1
             str(istr)(j:j) = line(i:i)
           else
              was_in_string = 0
           endif
         endif

         if(is_operator .eq. 1) then
             istr = istr + 1
             str(istr)(1:1) = line(i:i)
             istr_op = istr
         endif

         if(line(i:i) .eq. '(') then
           istr_op = istr
         endif
      enddo
      nstr = istr

      if(jdebug .gt. 10) then
      write (6,*) 'nstr    = ', nstr
      write (6,*) 'istr_op = ', istr_op
      do i = 1, nstr
        write (6,*) i, trim(str(i))
      enddo
      write (6,*)  " "
      endif
c   length= 17     LogRHO = log(RHO)
c   length= 19     V3     = Vx, Vy, Vz
c   length= 17     Vort3  = curl(V3)
c   length= 20     Vort   = norm(Vort3)


      return
      end
C ====================================================================72
      subroutine report_time()
      implicit none
      include 'postproc_info_e3d.h'
      real time
      call get_time(time)
      write (6,999) trim(cDumpFile0), time
999   format(a,"  ", 1pe15.6)
      stop
      end

C ====================================================================72
      subroutine read_formulas_file()
      implicit none
      include 'postproc_info_e3d.h'
      integer i, j, inext

      ! read all lines in formulas file 
      ! Should include definitions for all relevant fields

      open(unit=17, file='formulas.e3d', form='formatted',
     1     status='old', err=80)

      inext = 1
10    i = inext
      do j = 1, 256
        cFormulas(i)(j:j) = ' '
      enddo
      read (17,990,end=90) cFormulas(i)
      inext = i + 1
      if(cFormulas(i)(1:1) .eq. '#') inext = i
      if(cFormulas(i)(1:12) .eq. 'index_values') then
         call set_index_values(cFormulas(i)(13:100))
         inext = i
      endif
      if(cFormulas(i)(1:15) .eq. 'index_variables') then
         call set_index_vars(cFormulas(i)(16:100))
         inext = i
      endif
      if(inext.gt.i .and. nIndexVars.gt.0) call do_indexes(inext)
990   format(a256)
      goto 10

80    continue
      write (6,*) "Problem: open failed on formulas.e3d"
      stop

90    continue
      close(17)

       nFormulas = i - 1

      if(idebug .gt. 100) call write_formulas()

      return
      end

C ====================================================================72
      subroutine write_formulas()
      implicit none
      include 'postproc_info_e3d.h'
      integer i

      write (6,*) ' '
      write (6,*) 'nformulas = ', nformulas
      do i = 1, nformulas
        write (6,999) trim(cFormulas(i))
999     format('Formula: ', a)
      enddo
      write (6,*) ' '

      return
      end
C ====================================================================72
      subroutine split_string(cpat,s12,s1,s2)
      implicit none
      character*(*) cpat,s12,s1,s2
      integer ic, n, i1, ifind_in_str
      ic = ifind_in_str(cpat, s12)
      n = len(trim(s12))
      call setblank(s1)
      call setblank(s2)
      s1(1:ic-1) = s12(1:ic-1)
      i1 = ic + len(trim(cpat))
      s2(1:(n-i1)+1) = s12(i1:n)
      return
      end
      
C ====================================================================72
      integer function ifind_in_str(cPat, cStr)
      implicit none
      character*(*) cPat, cStr
      integer nPat, nStr, i, ifind
      nPat = len(trim(cPat))
      nStr = len(trim(cStr))
      ifind = -1
      do i = 1, nStr+1-nPat
         if(cpat(1:nPat) .eq. cStr(i:i+nPat-1)) ifind = i
      enddo
      ifind_in_str = ifind

      return
      end

C ====================================================================72
      subroutine replace_a_with_b_in_str(a, b, s)
      implicit none
      character*(*) a, b, s
      integer na, ns, i, j
      ! assumes a and b are the same length
      na = len(trim(a))
      ns = len(trim(s))
      do i = 1, ns+1-na
        j = i + na - 1
        if(a(1:na) .eq. s(i:j)) s(i:j) = b(1:na)
      enddo

      return
      end

C ====================================================================72
      subroutine do_indexes(inext)
      implicit none
      include 'postproc_info_e3d.h'
      integer inext, i0, ieq, ifind_in_str, neqn, ntoexpand, ival
      integer iivar, im, j, j1, i1, id, i, ic1, ic2, ic3, ic4, nterm
      character*(1000) ceqn, cterm
      character*2 ca, cb

      ! Return if no index variables are in current formula
      i0 = inext -1
      ic1 = -1
      do j = 1, nIndexVars
        ic1 = max(ic1, ifind_in_str(cIndexVars(j),cFormulas(i0)))
      enddo
      if(ic1 .lt. 0) return
      
      ! Enumerate equations on index vars in RLS of equation
      ceqn = cFormulas(i0)
      ieq = ifind_in_str("=", ceqn)       ! Location of "="
      do iivar = 1, nIndexVars
        im = ifind_in_str(cIndexVars(iivar), cFormulas(i0)(1:ieq))
        if(im .gt. 0) then
c          write (6,*) "Expanding LHS ", trim(cIndexVars(iivar))
          ntoexpand = inext - i0
c          write (6,*) "ntoexpand   = ", ntoexpand
          do j = 1, nIndexValues
            j1 = mod(j, nIndexValues)
            i1 = i0 + j1*ntoexpand
            do id = 0, ntoexpand-1
              i = i1 + id
              if(j1 .gt. 0) cFormulas(i) = cFormulas(i0+id)
              ca = cIndexVars(iivar)
              cb = cIndexValues(j1+1)
              call replace_a_with_b_in_str(ca,cb,cFormulas(i))
c              write (6,*) "i,form ",i,trim(cFormulas(i))
            enddo
          enddo
          inext = inext + ntoexpand * (nindexvalues-1)
        endif
      enddo

       ! If LHS index variables generated an N-tuple of formulas
       ! then add one more formula to point to the N-tuble
       if(inext-i0 .gt. 1) then
         call setblank(cFormulas(inext))
         cFormulas(inext)(1:ieq) = ceqn(1:ieq)
         do iivar = 1, nIndexVars
           ca = cIndexVars(iivar)
           call replace_a_with_b_in_str(ca,"  ",cFormulas(inext))
        enddo
        do i = i0, inext-1
          ic1 = 2 + ieq  + (  i-i0)*(ieq-1)
          ic2 = min(256, 1 + ieq  + (1+i-i0)*(ieq-1))
          cFormulas(inext)(ic1:ic2) = cFormulas(i)(1:ieq-1)
        enddo
         inext = i + 1
       endif

      ! Sum on RHS of equation on inedx vars remaining (only on RHS)
      do i = i0, inext-1
      do iivar = 1, nIndexVars
        ic1 = ieq + 1
        ic2 = len(trim(cFormulas(i)))
        ic3 = ieq + 1
        nterm = 1 + ic2 - ic1
        cterm(1:nterm) = cFormulas(i)(ic1:ic2)   ! template term with varible index

        ca = cIndexVars(iivar)
        if(ifind_in_str(ca,cFormulas(i)) .gt. 0) then
          do ival = 1, nIndexValues
            ic4 = ic3 + nterm - 1
            cFormulas(i)(ic3:ic4) = cterm(1:nterm)
            cb = cIndexValues(ival)
            call replace_a_with_b_in_str(ca,cb,cFormulas(i)(ic3:ic4))
            ic3 = ic3 + nterm
            if(ival .lt. nIndexValues) then
              cFormulas(i)(ic4+1:ic4+2) = " +"
              ic3 = ic3 + 2
            endif
          enddo
        endif
      enddo
      enddo

      if(idebug .ge. 1000) then
        nformulas = inext - 1
        call write_formulas()
         stop
      endif

      return
      end

C ====================================================================72
      subroutine set_index_values(str)
      implicit none
      include 'postproc_info_e3d.h'
      character*(*) str
      character*100 cvals
      integer i, n

      read (str,*) cvals
      n = len(trim(cvals))
      nIndexValues = n - 1

      if(nIndexValues .gt. 8) then
        write (6,*) "PROBLEM: to many index alues:"
        write (6,*) "index value string = ", str(1:30)
        write (6,*) "nIndexValues = ", nIndexValues
        write (6,*) "cvals       = ", cvals(1:n)
        stop
      endif

      do i = 1, nIndexValues
        cIndexValues(i)(1:1) = cvals(1:1)
        cIndexValues(i)(2:2) = cvals(1+i:1+i)
      enddo

      return
      end

C ====================================================================72
      subroutine set_index_vars(str)
      implicit none
      include 'postproc_info_e3d.h'
      character*(*) str
      character*100 cvars
      integer i, n

      read (str,*) cvars
      n = len(trim(cvars))
      nIndexVars   = n - 1

      if(nIndexVars .gt. 20) then
        write (6,*) "PROBLEM: to many index variables:"
        write (6,*) "index variable string = ", str(1:30)
        write (6,*) "nIndexVars   = ", nIndexVars
        write (6,*) "cvars       = ", cvars(1:n)
        stop
      endif

      do i = 1, nIndexVars
        cIndexVars(i)(1:1) = cVars(1:1)
        cIndexVars(i)(2:2) = cVars(1+i:1+i)
      enddo

      return
      end

C ====================================================================72
      subroutine read_input_file()
      implicit none
      include 'postproc_info_e3d.h'
      character*256 cline, ckey
      integer i, iStrLen

      ! Parse input file for dump file, field, and output, output directory, ...
      !
      open(unit=17, file='input.e3d', form='formatted',
     1     status='old', err=80)
10    continue
      do i = 1, 256
        cline(i:i) = ' '
      enddo
      read (17,990,end=90) cline
990   format(a256)
      if(cline(1:10).eq.'dumpfile  ') read (cline,*) ckey,cDumpFile0
      if(cline(1:10).eq.'field     ') read (cline,*) ckey,cField
      if(cline(1:10).eq.'outdir    ') read (cline,*) ckey,cOutDir
      if(cline(1:10).eq.'nfldxyz   ') read (cline,*) ckey,nfx,nfy,nfz
      if(cline(1:10).eq.'minmax    ') read (cline,*) ckey,hvmin,hvmax
      if(cline(1:10).eq.'threshold ') read (cline,*) ckey,threshold
      if(cline(1:10).eq."xrange    ") read (cline,*) ckey,xmin,xmax
      if(cline(1:10).eq."yrange    ") read (cline,*) ckey,ymin,ymax
      if(cline(1:10).eq."zrange    ") read (cline,*) ckey,zmin,zmax
      if(cline(1:10).eq."Rrange    ") read (cline,*) ckey,xmin,xmax
      if(cline(1:10).eq."Zrange    ") read (cline,*) ckey,ymin,ymax
c      if(cline(1:10).eq.'nblend    ') read (cline,*) ckey,nblend
      if(cline(1:10).eq.'xcoord    ') read (cline,*) ckey,cxcoord
      if(cline(1:10).eq.'ycoord    ') read (cline,*) ckey,cycoord
      if(cline(1:10).eq.'zcoord    ') read (cline,*) ckey,czcoord
      if(cline(1:10).eq.'c2label   ') read (cline,*) ckey,c2label
      if(cline(1:10).eq.'c2map     ') read (cline,*) ckey,c2map
      if(cline(1:10).eq.'setconst  ') call get_const(cline)
      if(cline(1:10).eq.'output    ') then
        read (cline,*) ckey,cOutput
c        if(cOutput(1:10).eq.'xstrip    ') read (cline,*) ckey,cOutput
c        if(cOutput(1:10).eq.'ystrip    ') read (cline,*) ckey,cOutput
c        if(cOutput(1:10).eq.'xyslice   ') read (cline,*) ckey,cOutput
        if(cOutput(1:10).eq.'statsz    ') then
          read (cline,*) ckey,cOutput,z00
        endif
      endif
      goto 10

80    continue
      write (6,*) "Problem: open failed on input.e3d"
      stop

90    continue
      close(17)

      cxgrad(6:6) = cxcoord(1:1)
      cygrad(6:6) = cycoord(1:1)
      czgrad(6:6) = czcoord(1:1)

      if(x00 .ne. finvalid) then
         xmin = x00
         xmax = x00
      endif
      if(y00 .ne. finvalid) then
         ymin = y00
         ymax = y00
      endif
      if(z00 .ne. finvalid) then
         zmin = z00
         zmax = z00
      endif

      ! Generate root name of output file based on input file, field
      ! Extension (.*) of output file will depend on output format type
      iStrLen = len(trim(cDumpFile0)) - 3                    ! trim off the last 3 digits
      cOutFile(1:iStrLen) = cDumpFile0(1:iStrLen)
      iStr0 = 1 + iStrLen
      iStrLen = len(trim(cField))
      cOutFile(iStr0:iStr0+iStrLen-1) = cField(1:iStrLen)
      iStr0 = iStr0 + iStrLen

      if(idebug .gt. 100) then
        write (6,*) 'Dump File = ', trim(cDumpFile0)
        write (6,*) 'Field     = ', trim(cField)
        write (6,*) 'Output    = ', trim(cOutput)
        write (6,*) ' '
      endif

      return
      end

C ====================================================================72
      subroutine get_const(cline)
      implicit none
      include 'postproc_info_e3d.h'
      character*256 cline, ckey
      character*40 cstr
      integer i, iStrLen

      cstr = "                                        "
      read (cline,*) ckey,cstr
      if(cstr(1:6) .eq. "debug=") read (cstr(7:40),*) idebug
      if(cStr(2:2) .eq. "=") then
        if(cStr(1:1).eq. cxcoord(1:1)) read (cStr(3:40),*) x00
        if(cStr(1:1).eq. cycoord(1:1)) read (cStr(3:40),*) y00
        if(cStr(1:1).eq. czcoord(1:1)) read (cStr(3:40),*) z00
      endif
      if(cstr(1:5) .eq. "part=") read (cstr(6:40),*) cPart
      if(cstr(1:11) .eq. "cbar_steps=") read (cstr(12:40),*) ncbarsteps
      if(cstr(1:5) .eq. "cbar=") read (cstr(6:40),*) cbarString
      if(cstr(1:5) .eq. "curl=") read (cstr(6:40),*) cCurlCalc
      if(cstr(1:4) .eq. "div=") read (cstr(5:40),*) threshold
      if(cstr(1:8) .eq. "logscale") ilogscale = 1
      if(cstr(1:12) .eq. "logsizescale") isizelogscale = 1
      if(cstr(1:10) .eq. "xboundary=") read (cstr(11:40),*) iboundary(1)
      if(cstr(1:10) .eq. "yboundary=") read (cstr(11:40),*) iboundary(2)
      if(cstr(1:10) .eq. "zboundary=") read (cstr(11:40),*) iboundary(3)
      if(cstr(1:5) .eq. "flip1") iflip1 = 1
      if(cstr(1:10) .eq. "maxrefine=") read (cstr(11:40),*) maxrefine
      if(cstr(1:8) .eq. "maxboxs=") read (cstr(9:40),*) maxboxs_db
      if(cstr(1:5) .eq. "view=") read (cstr(6:40),*) cView
      if(cstr(1:5) .eq. "region=") read (cstr(6:40),*) cRegion
      if(cstr(1:7) .eq. "nblend=") read (cstr(8:40),*) nblend
      if(cstr(1:7) .eq. "maxres=") read (cstr(8:40),*) maxres

      return
      end
C ====================================================================72
      subroutine set_xyz_minmax_from_endpoints()
      include 'postproc_info_e3d.h'

      xmin = min(endpt0(1), endpt1(1))
      xmax = max(endpt0(1), endpt1(1))

      ymin = min(endpt0(2), endpt1(2))
      ymax = max(endpt0(2), endpt1(2))

      zmin = min(endpt0(3), endpt1(3))
      zmax = max(endpt0(3), endpt1(3))

c      write (6,*) " "
c      write (6,*) "xmin,ymin,zmin",xmin,ymin,zmin
c      write (6,*) "xmax,ymax,zmax",xmax,ymax,zmax

      return
      end
C ====================================================================72
      subroutine set_end_points()
      implicit none
      include 'postproc_info_e3d.h'
      real fov, eye(3), center(3), up(3), fnear, ffar, ec(3), pi
      real vecnorm, ecmag, dist
      integer j

      endpt0(1) = xmin
      endpt1(1) = xmax
      endpt0(2) = 0.5*(ymin+ymax)
      endpt1(2) = 0.5*(ymin+ymax)
      endpt0(3) = 0.5*(zmin+zmax)
      endpt1(3) = 0.5*(zmin+zmax)

      if(cView(1:1) .ne. ' ') then
        if(idebug .ge. 21)
     1      write (6,*) "Setting endpoints from view: ", trim(cView)
        call read_view(fov,eye,center,up,fnear,ffar,up)
        call vecsub(ec,eye,center)
        ecmag = sqrt(ec(1)**2 + ec(2)**2 + ec(3)**2)
        call normalize(ec)
        call cross(endpt_vec,up,ec)
        call normalize(endpt_vec)
        pi = 4.0 * atan(1.0)
        dist = ecmag
        if(fnear .gt. 0.333*ecmag) dist = fnear
        extent = dist * tan(0.5*fov*pi/180.0)  ! extent of veiw at dist from eye
        do j = 1, 3
          endpt0(j) = eye(j) - dist*ec(j) - extent*endpt_vec(j)   ! Sample line is on near clipping plane
          endpt1(j) = eye(j) - dist*ec(j) + extent*endpt_vec(j)

          vview(j,1) = endpt_vec(j)
          vview(j,2) = up(j)
          vview(j,3) = ec(j)
        enddo
      endif
      if(idebug .ge. 21) then
        write (6,*) " "
        write (6,*) "extent   :", extent
        write (6,*) "center   :",(      center(j),j=1,3)
        write (6,*) "up       :",(          up(j),j=1,3)
        write (6,*) "ec       :",(          ec(j),j=1,3)
        write (6,*) "endpt_vec:",(   endpt_vec(j),j=1,3)
        write (6,*) "endpt0   :",(      endpt0(j),j=1,3)
        write (6,*) "endpt1   :",(      endpt1(j),j=1,3)
        write (6,*) " "
      endif

      return
      end

C ====================================================================72
      subroutine read_view(fov,eye,center,up,fnear,ffar)
      implicit none
      include 'postproc_info_e3d.h'
      character*16 cstr
      real fov, eye(3), center(3), up(3), fnear, ffar
      integer j

      if(idebug .gt. 20) write (6,*) "IN read_view: ", trim(cView)
      open(unit=17,file=cView,status='old',form='formatted',err=1000)
      read (17,*,err=1001) cstr, fov
      read (17,*,err=1001) cstr, (   eye(j),j=1,3)
      read (17,*,err=1001) cstr, (center(j),j=1,3)
      read (17,*,err=1001) cstr, (    up(j),j=1,3)
      read (17,*,err=1001) cstr, fnear, ffar
      close(17)
      return

1000  write (6,*) "PROBLEM in read_view: could not open ", trim(cView)
      stop
1001  write (6,*) "PROBLEM in read_view: could not read ", trim(cView)
      stop
      end


C ====================================================================72
      subroutine write_view(cView,fov,eye,center,up,fnear,ffar,
     1                      xmin,xmax,ymin,ymax,zmin,zmax)
      implicit none
      character*(*) cView
      real xmin,xmax,ymin,ymax,zmin,zmax
      real fov, eye(3), center(3), up(3), fnear, ffar
      integer k

      open(unit=11,file=cView,form="formatted")
      write (11,901) "fov         ", fov
      write (11,901) "eye         ", (eye(k),     k=1,3)
      write (11,901) "center      ", (center(k), k=1,3)
      write (11,901) "up          ", (up(k), k=1,3)
      write (11,901) "clipping    ", fnear, ffar
901   format(a12, 3f12.6)

      do k=0,5
        write (11,902) "ifauxclip   ", k, 2
902     format(a12, 2i2)
      enddo

      write (11,903) 0, -1.0,  0.0,  0.0, xmin
      write (11,903) 1,  1.0,  0.0,  0.0, xmax
      write (11,903) 2,  0.0, -1.0,  0.0, ymin
      write (11,903) 3,  0.0,  1.0,  0.0, ymax
      write (11,903) 4,  0.0,  0.0, -1.0, zmin
      write (11,903) 5,  0.0,  0.0,  1.0, zmax
903   format("auxclip     ", i2, 3f5.1, f12.6)
      close(11)

      return
      end

C ====================================================================72
      subroutine set_postproc_defaults()
      implicit none
      include 'postproc_info_e3d.h'
      include 'amr_e3d.h'
      integer i, j

      ! Initialize & set defaults
      !
      do i = 1, 256
        cDumpFile0(i:i) = ' '
        cField(i:i)    = ' '
        cFieldLabel(i:i)    = ' '
        cOutput(i:i)   = ' '
        cOutSymb(i:i)  = ' '

        cOutHost(i:i) = ' '
        cOutDir(i:i)  = ' '
        cOutFile(i:i) = ' '
        cCurlCalc(i:i) = ' '

        cInVar(i:i) = ' '
        c_add_tof_problem(i:i) = ' '
        cView(i:i) = ' '
      enddo
      call setblank(c2label)
      call setblank(c2map)
      cxcoord = "x               "
      cycoord = "y               "
      czcoord = "z               "
      cxgrad  = "grad_x"
      cygrad  = "grad_y"
      czgrad  = "grad_z"
      cCurlCalc(1:4) = "cell"

      nfx0      = -1
      nfy0      = -1
      nfz0      = -1
      nfbndx    = 0
      nfbndy    = 0
      nfbndz    = 0
      nblend    = 1
      cOutHost(1:9) = 'localdisk'
      cOutDir(1:1)  = '.'
      i64output  = 0
       hvmin     = 1.0
       hvmax     = 0.0
       small     = 1.0e-30
       xmin      = 1.0
       ymin      = 1.0
       zmin      = 1.0
       xmax      = 0.0
       ymax      = 0.0
       zmax      = 0.0
       threshold = 0.0
       ilogscale = 0
       isizelogscale = 0
       qcmixfac  = 0.5
       qcdivlim0 = -1000000.0
       qcdivlim  = 2.0*qcdivlim0
       maxres    = 2000000000
      nvalues    = 0
      nDumps     = 0
      nadvect    = 6
      c_add_tof_problem(1:4) = 'none'
      finvalid = sqrt(987654321.0)
      x00 = finvalid
      y00 = finvalid
      z00 = finvalid
      cPart = "full                            "
      isRZ = 0
      ncbarsteps = -1
      call setblank(cbarString)
      do j = 1, 3
        iboundary(j) = 0                     ! 0=periodic
      enddo
      nIndexValues = 0
      nIndexVars   = 0
      iflip1 = 0

      nboxs     = 0         ! For AMR: value if NOT AMR
      nlevels   = 1         ! For AMR: value if NOT AMR
      nrep_used = 1          ! For AMR: replication (refinment) factor used
      nboxsused = 0          ! For AMR: number of read bricks (read bricks fit in boxs)
      maxrefine = 1000000000 ! For AMR: max refinment level used
      icurrent_box = -1      ! For AMR: current box number
      maxboxs_db   = -1      ! For AMR: this if for debugging mainly
      ntslab = 32            ! old default value for Z-thickness

      ! Default view cirections are along XYZ
      do j = 1,3
        do i = 1,3
          vview(i,j) = 0.0
          if(i .eq. j) vview(j,j) = 1.0
        enddo
      enddo

      return
      end

C ====================================================================72
      subroutine get_op_prop(cOpr, iOpr, nDim, nArg, nbdy_xyz, icomp)
      implicit none
      include 'postproc_info_e3d.h'
      character*32 copr
      integer nDim, iOpr, nArg, nbdy_xyz(3), nchar, icomp, n1,n2,n3

      iOpr = -1
      nbdy_xyz(1) = 0
      nbdy_xyz(2) = 0
      nbdy_xyz(3) = 0
      icomp = 0
      if(trim(copr) .eq. "read") then
        iOpr = iop_read
        nDim = 1
        nArg = 0
      endif

      if(trim(copr) .eq. "log") then
        iOpr = iop_log
        nDim = 1
        nArg = 1
      endif

      if(trim(copr) .eq. "log10") then
        iOpr = iop_log10
        nDim = 1
        nArg = 1
      endif
      if(trim(copr) .eq. "sqrt") then
        iOpr = iop_sqrt
        nDim = 1
        nArg = 1
      endif
      if(trim(copr) .eq. "step") then
        iOpr = iop_step
        nDim = 1
        nArg = 2
      endif
      if(trim(copr) .eq. "cos") then
        iOpr = iop_cos
        nDim = 1
        nArg = 1
      endif
      if(trim(copr) .eq. "sin") then
        iOpr = iop_sin
        nDim = 1
        nArg = 1
      endif
      if(trim(copr) .eq. "exp") then
        iOpr = iop_exp
        nDim = 1
        nArg = 1
      endif
      if(trim(copr) .eq. "abs") then
        iOpr = iop_abs
        nDim = 1
        nArg = 1
      endif
      if(trim(copr) .eq. "set_constant") then
        iOpr = iop_const
        nDim = 1
        nArg = 1
      endif

      if(trim(copr) .eq. "curl") then
        iOpr = iop_curl
        nDim = 3
        nArg = 1
        nbdy_xyz(1) = 1
        nbdy_xyz(2) = 1
        nbdy_xyz(3) = 1
      endif
      if(trim(copr) .eq. "curld") then
        iOpr = iop_curld
        nDim = 3
        nArg = 1
        nbdy_xyz(1) = 1
        nbdy_xyz(2) = 1
        nbdy_xyz(3) = 1
      endif
      if(trim(copr) .eq. "eigenv_s") then
        iOpr = iop_igenv_s
        nDim = 3
        nArg = 1
        nbdy_xyz(1) = 0
        nbdy_xyz(2) = 0
        nbdy_xyz(3) = 0
      endif

      if(trim(copr) .eq. "grad") then
        iOpr = iop_grad
        nDim = -333  ! means nDim = 3*nDim(1st and only operand)
        nArg = 1
        nbdy_xyz(1) = 1
        nbdy_xyz(2) = 1
        nbdy_xyz(3) = 1
      endif
      if(copr(1:6) .eq. cxgrad) then
        iOpr = iop_grad_1
        nDim = 1
        nArg = 1
        nbdy_xyz(1) = 1
        nbdy_xyz(2) = 0
        nbdy_xyz(3) = 0
      endif
      if(copr(1:6) .eq. cygrad) then
        iOpr = iop_grad_2
        nDim = 1
        nArg = 1
        nbdy_xyz(1) = 0
        nbdy_xyz(2) = 1
        nbdy_xyz(3) = 0
      endif
      if(copr(1:6) .eq. czgrad) then
        iOpr = iop_grad_3
        nDim = 1
        nArg = 1
        nbdy_xyz(1) = 0
        nbdy_xyz(2) = 0
        nbdy_xyz(3) = 1
      endif

      if(trim(copr) .eq. "norm") then
        iOpr = iop_norm
        nDim = 1
        nArg = 1
      endif
      if(copr(1:4) .eq. "getc") then
        iOpr = iop_getc
        nDim = 1
        nArg = 1
      endif
      if(trim(copr) .eq. "radc") then
        iOpr = iop_radc
        nDim = 1
        nArg = 1
      endif


      if(trim(copr) .eq. "subtrace") then
        iOpr = iop_trsub
        nDim = 9
        nArg = 1
      endif
      if(trim(copr) .eq. "trace") then
        iOpr = iop_trace
        nDim = -999    ! means nDim = NDim1 - 9
        nArg = 1
      endif
      if(trim(copr) .eq. "det") then
        iOpr = iop_det
        nDim = 1
        nArg = 1
      endif

      if(copr(1:7) .eq. "sect_av") then
        iOpr = iop_secsum
        nDim = 1
        nArg = 1
      endif
      if(copr(1:8) .eq. "sect_sum") then
        iOpr = iop_secsum
        nDim = 1
        nArg = 1
      endif

      if(copr(1:9) .eq. "component") then
        nchar = len(trim(copr))
        read (copr(10:nchar),*) icomp
        iOpr = iop_comp
        nDim = 1
        nArg = 1
      endif

      if(copr(1:7) .eq. "smooth_") then
        nchar = len(trim(copr))
        call parse_xyz_pars(copr(8:nchar), n1,n2,n3)
        if(max(n1,n2,n3) .le. 0) then
          read (copr(8:nchar),*) nsmooth
          n1 = nsmooth
          n2 = nsmooth
          n3 = nsmooth
        endif
        iOpr = iop_smooth
        nDim = -111  ! means nDim = nDim(1st and only operand)
        nArg = 1
        nbdy_xyz(1) = n1
        nbdy_xyz(2) = n2
        nbdy_xyz(3) = n3
      endif
      if(copr(1:11) .eq. "ballfilter_") then
        nchar = len(trim(copr))
        read (copr(12:nchar),*) nsmooth
        iOpr = iop_ballflt
        nDim = -111  ! means nDim = nDim(1st and only operand)
        nArg = 1
        nbdy_xyz(1) = nsmooth
        nbdy_xyz(2) = nsmooth
        nbdy_xyz(3) = nsmooth
      endif

      if(trim(copr) .eq. "+") then
        iOpr = iop_plus
        nDim = -3    ! means nDim = nDim(1st operand) == ndim(2nd)
        nArg = 2
      endif
      if(trim(copr) .eq. "-") then
        iOpr = iop_minus
        nDim = -3
        nArg = 2
      endif
      if(trim(copr) .eq. "*") then
        iOpr = iop_times
        nDim = -1  ! -1 means nDim = nDim(1st operand) * nDim(2nd operand)
        nArg = 2
      endif
      if(trim(copr) .eq. "/") then
        iOpr = iop_divide
	nDim = -2   ! -2 means nDim = nDim(1st operand) & nDim(2nd) == 1
	nArg = 2
      endif

      if(trim(copr) .eq. "cross") then
        iOpr = iop_cross
        nDim = 3
        nArg = 2
      endif

      if(trim(copr) .eq. "dot") then
        iOpr = iop_dot
        nDim = -1000
        nArg = 2
      endif
      ! why is this repeated???
      if(trim(copr) .eq. "curl") then
        iOpr = iop_curl
        nDim = 3
        nArg = 1
        nbdy_xyz(1) = 1
        nbdy_xyz(2) = 1
        nbdy_xyz(3) = 1
      endif
      if(trim(copr) .eq. "div") then
        iOpr = iop_div
        nDim = 1
        nArg = 1
        nbdy_xyz(1) = 1
        nbdy_xyz(2) = 1
        nbdy_xyz(3) = 1
      endif
      if(trim(copr) .eq. "a.grad_b") then
        iOpr = iop_agradb
        nDim = 3
        nArg = 2
        nbdy_xyz(1) = 1
        nbdy_xyz(2) = 1
        nbdy_xyz(3) = 1
      endif
      if(copr(1:5) .eq. "coord") then
        if(copr(6:6) .eq. "1") iOpr = iop_coord1
        if(copr(6:6) .eq. "2") iOpr = iop_coord2
        if(copr(6:6) .eq. "3") iOpr = iop_coord3
        nDim = 1
        nArg = 1
        nbdy_xyz(1) = 0
        nbdy_xyz(2) = 0
        nbdy_xyz(3) = 0
      endif

      if(trim(copr) .eq. "pullback") then
        iOpr = iop_pullb
        nDim = 1
        nArg = 4
        nbdy_xyz(1) = nadvect
        nbdy_xyz(2) = nadvect
        nbdy_xyz(3) = nadvect
      endif

      if(isRZ .eq. 1) nbdy_xyz(3) = 0

      return
      end

C ====================================================================72
      subroutine check_xyzminmax()
      implicit none
      include 'postproc_info_e3d.h'
      real av

c >>>dhp
      if("strip    " .eq. cOutPut(2:10)) then

        if(cxcoord.eq.cOutPut(1:1) .or. czcoord.eq.cOutPut(1:1)) then
          if(ymin .ne. ymax) then
            av = 0.5 * (ymin + ymax)
            ymin = av
            ymax = av
          endif
        endif

        if(cxcoord.eq.cOutPut(1:1) .or. cycoord.eq.cOutPut(1:1)) then
          if(zmin .ne. zmax) then
            av = 0.5 * (zmin + zmax)
            zmin = av
            zmax = av
          endif
        endif

      endif

      return
      end

C ====================================================================72
      subroutine gen_output_label()
      implicit none
      include 'postproc_info_e3d.h'
      integer i, ilensymb, ilenname
      cFieldLabel = cField   ! ultimate fall back
      do i = 1, ntof
        if(trim(cField) .eq.  trim(cFld_tof(i))) then
           ! (cNam_tof(i)(1:1) .ne. ' ')
           ilensymb = len(trim(cField))
           ilenname = len(trim(cNam_tof(i)))
           if(ilensymb .eq. ilenname) then
              cFieldLabel(1:32) = cNam_tof(i)(1:32)
           endif
        endif
      enddo

      return
      end
C ====================================================================72
      subroutine gen_output_par()
      implicit none
      include 'postproc_info_e3d.h'
      real x0,x1,y0,y1,z0,z1, deltamax
      integer i, k, imax, iStr1
      logical b_file_exists
      real*8 dd8
      character*256 cval

      call get_full_dimensions(nxfull, nyfull, nzfull)
      call get_full_XYZrange(x0,x1,y0,y1,z0,z1)

      ! for output bins and derivivatives
      dx = (x1-x0) / float(nxfull)
      dy = (y1-y0) / float(nyfull)
      dz = (z1-z0) / float(nzfull)
      dx2i = 1.0 / (2.0 * dx)
      dy2i = 1.0 / (2.0 * dy)
      dz2i = 1.0 / (2.0 * dz)
      fLimits(1) = x0
      fLimits(2) = x1
      fLimits(3) = y0
      fLimits(4) = y1
      fLimits(5) = z0
      fLimits(6) = z1

      call calc_nblend_from_maxres()

      if(xmin .gt. xmax) then
        xmin = x0 
        xmax = x1 
      endif
      if(ymin .gt. ymax) then
        ymin = y0 
        ymax = y1 
      endif
      if(zmin .gt. zmax) then
        zmin = z0 
        zmax = z1 
      endif

      izgbllo = 1         ! Default min of Z-range needed for output
      izgblhi = nzfull    ! Default max of Z-range needed for output

      if(cxcoord(1:1) .eq. "R" .and. nzfull .eq. 1) isRZ = 1

      nbdy_output = 0     ! Default number of boundaryies needed for output


c >>>

      if("xyzavs    " .eq. cOutput(1:10) .or.
     1   "correlate " .eq. cOutPut(1:10) .or.
     1   "genview   " .eq. cOutPut(1:10) .or.
     1   "dist      " .eq. cOutPut(1:10) .or.
     1   "vrz       " .eq. cOutPut(1:10) .or.
     1   "vry       " .eq. cOutPut(1:10) .or.
     1   "vrx       " .eq. cOutPut(1:10) .or.
     1    "prof     " .eq. cOutPut(2:10) .or.
     1   "strip     " .eq. cOutPut(1:10) .or.
     1    "strip    " .eq. cOutPut(2:10) .or.
     1    "sect"      .eq. cOutPut(2: 5) .or.
     1    "charge"    .eq. cOutPut(2: 7)      ) then
c        call gen_mesh_range(1,ixout1,ixout2,iyout1,iyout2,izout1,izout2)

      if("strip     " .eq. cOutPut(1:10)) then
        nbdy_output = 1
        call set_end_points()
        call set_xyz_minmax_from_endpoints()
      endif

        call check_xyzminmax()

        call calc_nrep_used(2048)
        call gen_mesh_range(nrep_used,ixout1,ixout2,
     1                                iyout1,iyout2,
     1                                izout1,izout2 )

        izgbllo = izout1
        izgblhi = izout2

        nxout  = (1 + ixout2 - ixout1)
        nyout  = (1 + iyout2 - iyout1)
        nzout  = (1 + izout2 - izout1)

        if(idebug .gt. 30) then
          write (6,*) "izgbllo = ", izgbllo
          write (6,*) "izgblhi = ", izgblhi
          write (6,*) "gen_output_par: iyout1   = ", iyout1
          write (6,*) "gen_output_par: iyout2   = ", iyout2
        endif

c For AMR: support view (=xyzavs): generate:
c DONE      0) put data AMR parameters in thier own common block and header file amr_e3d.h
c DONE         a) amr.h incluided in read_adump_e3d.h (repaces what's in there now)
c DONE         b) amr_e3d.h inlcuded in e3d.f as needed
c DONE     1) dedermine highest replication factor used (nrep_used);
c DONE        a) subroutine calc_nrep_used(iMaxResOut)
c DONE        b) rescale n[xyz]out for nrep_used (nxout *= nrep_used);
c DONE     2) "fitted" box list for read list
c DONE     3)  gen_read_bricks
c DONE     4) work subvolumes from each of the fitted boxs in read list.
c      5) implement  read_fields_from_one_box  in read_adump_e3d.f
c         a) this may need to be an option in read_fld call (so that does not need a stub in other read_*_e3d.f source files)
c      6) Explicitly switche levels as outter loop: reset as needed: d[xyz], mappings for output; ...
c      7) modify fill_avs to replicate nrep_used/nrep(ilevel(ibox))
      endif
      
      iout_rprof   = 0
      if("radprof " .eq.  cOutput(1:8)) then
        iout_rprof = 1
        cOutFile(iStr0:iStr0+7) = ".radprof"
        maxprof = nxfull / 2
        maxcols = 7
        profcoordmin = 0.0
        delbin1 = dx
      endif
      
c      iout_xprof   = 0
c      if("prof " .eq.  cOutput(2:6)  .and.
c     1   cxcoord(1:1) .eq.  cOutput(1:1)) then
c        iout_xprof = 1
c        cOutFile(iStr0:iStr0+5) = ".xprof"
c        cOutFile(iStr0+1:iStr0+1) = cOutput(1:1)
c        maxprof = nxfull
c        maxcols = 7
c        profcoordmin = x0
c        delbin1 = dx
c      endif
      
      iout_xprof   = 0
      if("prof " .eq.  cOutput(2:6)  .and.
     1   cxcoord(1:1) .eq.  cOutput(1:1)) then
        iout_xprof = 1
        cOutFile(iStr0:iStr0+5) = ".xprof"
        cOutFile(iStr0+1:iStr0+1) = cOutput(1:1)
        maxprof = nxout
        maxcols = 7
        profcoordmin = x0 + dx * float(ixout1 - 1)
        delbin1 = dx
      endif

      iout_yprof   = 0
      if("prof " .eq.  cOutput(2:6)  .and.
     1   cycoord(1:1) .eq.  cOutput(1:1)) then
        iout_yprof = 1
        cOutFile(iStr0:iStr0+5) = ".yprof"
        cOutFile(iStr0+1:iStr0+1) = cOutput(1:1)
        maxprof = nyout
        write (6,*) "gen_output_par: nyout   = ", nyout
        write (6,*) "gen_output_par: maxprof = ", maxprof
        maxcols = 7
        profcoordmin = y0 + dy * float(iyout1 - 1)
        delbin1 = dy
      endif

      iout_zprof   = 0
      if("prof " .eq.  cOutput(2:6)  .and.
     1   czcoord(1:1) .eq.  cOutput(1:1)) then
        iout_zprof = 1
        cOutFile(iStr0:iStr0+5) = ".zprof"
        cOutFile(iStr0+1:iStr0+1) = cOutput(1:1)
        maxprof = nzout
        maxcols = 7
        profcoordmin = z0 + dz * float(izout1 - 1)
        delbin1 = dz
      endif

      iout_ysect   = 0
      if(("sect  "  .eq.  cOutput(2:6)  .or.
     1    "charge " .eq.  cOutput(2:8)       )  .and.
     1   (cycoord(1:1).eq.cOutput(1:1) .or. 'H'.eq.cOutput(1:1))) then
        iout_ysect = 1
        cOutFile(iStr0:iStr0+5) = ".ysect"
        if( "c" .eq.  cOutput(2:2)) cOutFile(iStr0:iStr0+7)=".ycharge"
        cOutFile(iStr0+1:iStr0+1) = cOutput(1:1)
        maxprof = nyfull
        maxcols = 1000
        profcoordmin = y0
        delbin1 = dy
      endif

      iout_xsect   = 0
      if("sect " .eq.  cOutput(2:6)  .and.
     1   cxcoord(1:1) .eq.  cOutput(1:1)) then
        iout_xsect = 1
        cOutFile(iStr0:iStr0+5) = ".xsect"
        cOutFile(iStr0+1:iStr0+1) = cOutput(1:1)
        maxprof = nxfull
        maxcols = 1000
        profcoordmin = x0
        delbin1 = dx
      endif

      if(iout_xprof+iout_yprof+iout_zprof+iout_rprof .gt. 0) then
         call write_prof_plt()
         if(b_file_exists(cOutFIle)) stop
      endif

      iout_xyzavs  = 0
      if("xyzavs    " .eq.  cOutput(1:10)) then
        iout_xyzavs = 1
        cOutFile(iStr0:iStr0+5) = ".xysof"
      endif

      iout_vrx  = 0
      if("vrx       " .eq.  cOutput(1:10)) then
        iout_vrx = 1
        cOutFile(iStr0:iStr0+3) = ".vrx"
      endif

      iout_vry  = 0
      if("vry       " .eq.  cOutput(1:10)) then
        iout_vry = 1
        cOutFile(iStr0:iStr0+3) = ".vry"
      endif


      iout_vrz  = 0
      if("vrz       " .eq.  cOutput(1:10)) then
        iout_vrz = 1
        cOutFile(iStr0:iStr0+3) = ".vrz"
      endif

c >>>>
c <<<<<<

      iout_hv = 0
      if("hv        " .eq.  cOutput(1:10)) then
        iout_hv = 1
        ntslab      = 32  ! ntslab should be a multipble of readbrick internal Z-dim
        itslab_zoff = 0
        cOutFile(iStr0:iStr0+3) = ".hv"
        if(hvmin .ge. hvmax) then
          write (6,*) 'Must set valid hv range (hvmin < hvmax) with'
          write (6,*) 'minmax    <hvmin> <hvmax>'
          write (6,*) 'in input.e3d file'
          stop
        endif
        hvscale = 255.99 / (hvmax - hvmin)
      endif

      iout_spc3v = 0
      if("spc3v " .eq.  cOutput(1:6)) then
        iout_spc3v = 1
        cOutFile(iStr0:iStr0+6) = ".spc3v"
        ntslab = 32  ! ntslab should be as large as possible
        nxbof   = nxfull
        nybof   = nyfull
        nzbof   = nzfull
        nvec = 3
      endif

      iout_spcnv = 0
      if("spc1s " .eq.  cOutput(1:6)) then
        iout_spcnv = 1
        cOutFile(iStr0:iStr0+6) = ".spc1s"
        ntslab = 32  ! ntslab should be as large as possible
        nxbof   = nxfull
        nybof   = nyfull
        nzbof   = nzfull
        nvec = 1
      endif


      ndisp_sets = 0
      iout_strfn = 0
      if("strfn3v   " .eq.  cOutput(1:10) .or.
     1   "strfnsc   " .eq.  cOutput(1:10)       ) then
        if(hvmin .ge. hvmax) then
          write (6,*) 'Must set valid min max range (hvmin<hvmax) with'
          write (6,*) 'minmax    <hvmin> <hvmax>'
          write (6,*) 'in f3dinput file'
          stop
        endif
        ! Modify range (specifyed by hvmin & hvmax for differences
        ! Include a factor of 2 extra range
        deltamax = 2.0 * (hvmax - hvmin)
        hvmin = -deltamax
        hvmax =  deltamax
        iout_strfn = 1
        if("strfn3v   " .eq.  cOutput(1:10)) then
          cOutFile(iStr0:iStr0+7) = ".strfn3v"
          nvec = 3
        else
          cOutFile(iStr0:iStr0+7) = ".strfnsc"
          nvec = 1
        endif
        call get_displacements(ndisp_sets, cdisp_names, ndisp_in_set,
     1                       ndiff_sampl,ixdisp,iydisp)
        del_vari  = float(ndiff_sampl) / max(abs(hvmax), abs(hvmin))
        del_vari2 = float(ndiff_sampl) / max(abs(hvmax), abs(hvmin))**2
c >>>
        ntslab = 32  ! ntslab should be as large as possible
        nxbof   = nxfull
        nybof   = nyfull
        nzbof   = nzfull

         ! limit displacements to be no more than half way across mesh
         max_disp_in_set = 0
         do k=1,ndisp_sets
           imax = ndisp_in_set(k)
           do i=ndisp_in_set(k), 2, -1
             if(2*ixdisp(i,k) .gt. nxfull) imax = i-1
             if(2*iydisp(i,k) .gt. nyfull) imax = i-1
           enddo
           ndisp_in_set(k) = imax
           max_disp_in_set = max(imax, max_disp_in_set)
         enddo

         write (6,*) "========== begin debug in main ================"
         write (6,*) 'max_disp_in_set = ', max_disp_in_set
         do k=1,ndisp_sets
           write (6,*) 'set, setname: ', k, cdisp_names(k)(1:20)
           do i=1,ndisp_in_set(k)
             write (6,*) 'i,dispXY: ', i, ixdisp(i,k), iydisp(i,k)
           enddo
         enddo
         write (6,*) "============ end debug in main ================"

      endif


      iout_bof = 0
      if("bof       " .eq.  cOutput(1:10)) then
        iout_bof = 1
        itslab_zoff = 0
c        ntslab = 32  ! ntslab should be as large a power of 2 as possible
        ntslab = nblend * (32 / nblend)
        nxbof   = nxfull / nblend
        nybof   = nyfull / nblend
        nzbof   = nzfull / nblend
        nzbslab =  ntslab /nblend
        write (6,*) "ntslab  = ", ntslab
        write (6,*) "nblend  = ", nblend
        write (6,*) "nzbslab = ", nzbslab
c        if(mod(ntslab,nblend) .ne. 0  .or.
        if(mod(nzfull,nblend) .ne. 0  .or.
     1     mod(nxfull,nblend) .ne. 0  .or.
     1     mod(nyfull,nblend) .ne. 0        ) then
          write (6,*) 'Dims are not a multiple of ntblend'
          write (6,*) 'nxfull,nyfull = ', nxfull,nyfull
          write (6,*) 'nblend = ', nblend
          stop
        endif
        cOutFile(iStr0:iStr0+4) = ".bof"
        if(maxres .lt. 2000000000) then
          cOutFile(iStr0:iStr0+6) = ".mrbof"
        else
cccccccccccccccccccccccccccccccccccccc
c will be
          if(nblend .gt.  1) then
            call gen_trimmed_int_string(nblend, cval)
            iStr1 = iStr0+5+len(trim(cval))
            cOutFile(iStr0:iStr1) = ".b" // trim(cval) // "bof"
            write (6,*) "cOutFile = ", trim(cOutFile)
          endif
cccccccccccccccccccccccccccccccccccccc
c was
        if(nblend .gt.  1) then
          if(nblend .le. 9) write (cOutFile(iStr0:iStr0+10),901) nblend
901       format(".b",i1,"bof")
          if(nblend .gt. 9) write (cOutFile(iStr0:iStr0+10),902) nblend
902       format(".b",i2,"bof")
        endif
cccccccccccccccccccccccccccccccccccccc

        endif
        n64slabsize = 4 * nxbof * nybof * nzbslab
      endif

      iout_strip   = 0
      if("strip " .eq.  cOutput(1:6)) then
        delbin1 = 0.25 * dx / float(nrep_used)
        profcoordmin = -0.25 * dx * float(int(extent/delbin1))
        maxprof = 2*int(extent/delbin1)
        maxcols = 2
        cOutFile(iStr0:iStr0) = "."
        cOutFile(iStr0+1:iStr0+5) = "strip"
        iout_strip   = 1
      endif

      iout_xstrip   = 0
      iout_ystrip   = 0
      if("strip " .eq.  cOutput(2:7)) then

         ! dhp my range
         if(x00 .eq. finvalid) x00 = 0.5*(xmin + xmax)
         if(y00 .eq. finvalid) y00 = 0.5*(ymin + ymax)
         if(z00 .eq. finvalid) z00 = 0.5*(zmin + zmax)
        ix00 = min(nxfull, max(1, int(0.5 + (x00 - x0)/dx)))
        iy00 = min(nyfull, max(1, int(0.5 + (y00 - y0)/dy)))
        iz00 = min(nzfull, max(1, int(0.5 + (z00 - z0)/dz)))

        if(cOutput(1:1) .ne. czcoord(1:1)) then
          if(nboxsused .eq. 0) then
            izgbllo = iz00
            izgblhi = iz00
          endif
          nslab    = 1
          ntslab   = 1
        endif
        cOutFile(iStr0:iStr0) = "."
        cOutFile(iStr0+1:iStr0+6) = cOutput(1:6)
        maxcols = 2

        if(cxcoord(1:1) .eq.  cOutput(1:1)) then
          iout_xstrip = 1
          maxprof = nxout
          profcoordmin = max(xmin, x0)
          delbin1 = dx / float(nrep_used)
        endif
        if(cycoord(1:1) .eq.  cOutput(1:1)) then
          iout_ystrip = 1
          maxprof = nyout
          profcoordmin = max(ymin, y0)
          delbin1 = dy / float(nrep_used)
        endif
        if(idebug .gt. 30) then
          write (6,*) "iout_xstrip = ", iout_xstrip
          write (6,*) "iout_ystrip = ", iout_ystrip
        endif
      endif

      iout_xyslice   = 0
      if("xyslice   " .eq.  cOutput(1:10)) then
        iout_xyslice = 1
        cOutFile(iStr0:iStr0+7) = ".xyslice"
        write (6,*) 'z00 = ', z00
        z00 = (z00-z0)/dz
        iz00 = min(nzfull, max(1, int(0.5 + z00)))
        izgbllo = iz00
        izgblhi = iz00

        if(hvmin .ge. hvmax) then
          write (6,*) 'Must set valid range (hvmin < hvmax) with'
          write (6,*) 'minmax    <hvmin> <hvmax>'
          write (6,*) 'in input.e3d file'
          stop
        endif
        hvscale = 255.99 / (hvmax - hvmin)
      endif

      iout_corr   = 0
      if("correlate " .eq.  cOutput(1:10)) then
        iout_corr   = 1
        cOutFile(iStr0:iStr0+4) = ".corr"
        nbin1 = 256
        nbin2 = 256
        call split_string(":",cField,cFld1,cFld2)
        call get_minmax(cFld1, bin10,bin11)
        call get_minmax(cFld2, bin20,bin21)
      endif

      iout_view   = 0
      if("genview   " .eq.  cOutput(1:10)) iout_view   = 1

      iout_dist   = 0
      if("dist      " .eq.  cOutput(1:10)) then
        if(hvmin .ge. hvmax) then
          write (6,*) 'PDF must set valid hv range (hvmin<hvmax) with'
          write (6,*) 'hvmin = ', hvmin
          write (6,*) 'hvmax = ', hvmax
          write (6,*) 'Set with:'
          write (6,*) 'minmax      <hvmin> <hvmax>'
          write (6,*) 'in f3dinput file'
          stop
        endif
        iout_dist = 1
        cOutFile(iStr0:iStr0+4) = ".dist"
        maxprof  = 256 * 21
        maxcols  = 2
        binscale = 256.0 / (hvmax - hvmin)
        binmin   = hvmin - 10.0*(hvmax - hvmin)
        profcoordmin = binmin
        delbin1      = 1.0 / binscale
        write (6,*) "maxprof = ", maxprof
        write (6,*) "binmin  = ", binmin
      endif

      return
      end

C ====================================================================72
      logical function box_in_vol(i)
      implicit none
      include 'postproc_info_e3d.h'
      include 'amr_e3d.h'
      integer i, irep, k ,ix1,ix2,iy1,iy2,iz1,iz2
      logical bOverlap

      !  TRUE if bounding limits of amr box i overlaps those requested

      if(i .lt. 1  .or.  i .gt. nboxs) then
        box_in_vol = .false.
        return
      endif

      irep = iboxrefine(i)
      call gen_mesh_range(irep,ix1,ix2,iy1,iy2,iz1,iz2)

      bOverlap = .true.
      if(ix1 .gt. iboxoff(1,i)+nboxmesh(1)) bOverlap = .false.
      if(ix2 .lt. iboxoff(1,i)+1          ) bOverlap = .false.

      if(iy1 .gt. iboxoff(2,i)+nboxmesh(2)) bOverlap = .false.
      if(iy2 .lt. iboxoff(2,i)+1          ) bOverlap = .false.

      if(iz1 .gt. iboxoff(3,i)+nboxmesh(3)) bOverlap = .false.
      if(iz2 .lt. iboxoff(3,i)+1          ) bOverlap = .false.

c >>> <<<
c      real xyz(2,3)
c      call gen_box_limits(i,xyz)
c      if(xmin.gt.xyz(2,1) .or. xmax.le.xyz(1,1)) bOverlap = .false.
c      if(ymin.gt.xyz(2,2) .or. ymax.le.xyz(1,2)) bOverlap = .false.
c      if(zmin.gt.xyz(2,3) .or. zmax.le.xyz(1,3)) bOverlap = .false.

      box_in_vol = boverlap

      return
      end

c      write (6,*) "BIV: box,irep,overlap:", i,irep,bOverlap
c      write (6,*) "boxmin:", (xyz(1,k), k=1,3)
c      write (6,*) "boxmax:", (xyz(2,k), k=1,3)
c      write (6,*) "  "
c      write (6,*) "volmin:", xmin, ymin, zmin
c      write (6,*) "volmax:", xmax, ymax, zmax
c      write (6,*) "  "

C ====================================================================72
      subroutine calc_nblend_from_maxres()
      implicit none
      include 'postproc_info_e3d.h'
      real    fmx3, fxyz, fsize
      integer i, ixyzmax, iuse
      fmx3 = 1.1 * float(maxres)**3
      fxyz = float(nxfull) * float(nyfull) * float(nzfull)
      if(fxyz .le. fmx3) return

      if(idebug .gt. 10) then
        write (6,*) "nxfull = ", nxfull
        write (6,*) "nyfull = ", nxfull
        write (6,*) "nzfull = ", nxfull
        write (6,*) "maxres = ", maxres
        write (6,*) "fxyz   = ", fxyz
        write (6,*) "fmx3   = ", fmx3
        write (6,*) " "
        write (6,*) "   i modx mody modz   fsize"
      endif
      ixyzmax = max(nxfull, nyfull, nzfull)
      iuse = -1
      i = 0
10    i = i + 1
      fsize = fxyz / (float(i)**3)
      if(mod(nxfull,i).eq.0 .and. mod(nyfull,i).eq.0 .and.
     1   mod(nzfull,i).eq.0 .and. fsize.le.fmx3            ) iuse = i
      if(idebug .gt. 10) then
        write (6,999) i,mod(nxfull,i),mod(nxfull,i),mod(nxfull,i),fsize
999     format(4i5, f12.1)
      endif
      if(iuse.lt.0  .and.  i.lt.ixyzmax) go to 10
      if(iuse .lt. 0) then
        write (6,*) "================================================"
        write (6,*) "PROBLEM in calc_nblend_from_maxres"
        write (6,*) "Can not find value for nblend that"
        write (6,*) "  1) divides evenly into full dimensions, and"
        write (6,*) "  2) is <= maxres*3 mesh points"
        write (6,*) "nxfull = ", nxfull
        write (6,*) "nyfull = ", nxfull
        write (6,*) "nzfull = ", nxfull
        write (6,*) "maxres = ", maxres
        write (6,*) "================================================"
        stop
      endif

      
      nblend = iuse

      return
      end

C ====================================================================72
      subroutine gen_zcell_wts(irep, iz1, iz2, izoff, wts)
      implicit none
      include 'postproc_info_e3d.h'
      include 'amr_e3d.h'
      integer irep, iz1, iz2, iz, izoff, nzmax
      real wts(iz1:iz2), dzr, z0, z1

      if(zmin .lt. zmax) then
        nzmax = nzfull * irep
        dzr   = (flimits(6)-flimits(5)) / float(nzmax)
        do iz = iz1, iz2
          z0 = max(zmin, dzr*float(izoff+iz-1) + flimits(5))
          z1 = min(zmax, dzr*float(izoff+iz  ) + flimits(5))
          wts(iz) = 0.0
          if(z0 .lt. z1) wts(iz) = (z1 - z0) / (zmax - zmin)
        enddo
      else
        do iz = iz1, iz2
          wts(iz) = 1.0
        enddo
      endif

      if(idebug .gt. 1000) then
        write (6,*) "=============================================="
        write (6,*) "irep = ", irep
        do iz = iz1, iz2
          write (6,999) izoff+iz, wts(iz)
999       format("=====> iz, wt: ", i5, f10.5)
        enddo
        write (6,*) "=============================================="
c        stop
      endif

      return
      end

C ====================================================================72
      subroutine gen_mesh_range(irep,ix1,ix2,iy1,iy2,iz1,iz2)
      implicit none
      include 'postproc_info_e3d.h'
      include 'amr_e3d.h'
      integer irep,ix1,ix2,iy1,iy2,iz1,iz2
      integer nxmax, nymax, nzmax
      real dxr, dyr, dzr

      nxmax = nxfull * irep
      nymax = nyfull * irep
      nzmax = nzfull * irep
      dxr   = (flimits(2)-flimits(1)) / float(nxmax)
      dyr   = (flimits(4)-flimits(3)) / float(nymax)
      dzr   = (flimits(6)-flimits(5)) / float(nzmax)

      ix1 = min(nxmax, max(1, 1+int((xmin-flimits(1))/dxr)))
      ix2 = min(nxmax, max(1, 1+int((xmax-flimits(1))/dxr)))

      iy1 = min(nymax, max(1, 1+int((ymin-flimits(3))/dyr)))
      iy2 = min(nymax, max(1, 1+int((ymax-flimits(3))/dyr)))

      iz1 = min(nzmax, max(1, 1+int((zmin-flimits(5))/dzr)))
      iz2 = min(nzmax, max(1, 1+int((zmax-flimits(5))/dzr)))

      return
      end

C ====================================================================72
      subroutine calc_nrep_used(iMaxResOut)
      implicit none
      include 'postproc_info_e3d.h'
      include 'amr_e3d.h'
      integer iMaxResOut, i, irep, ix1,ix2,iy1,iy2,iz1,iz2
      integer nx1,ny1,nz1
      logical box_in_vol

      nrep_used = 1   ! this variable is in postproc_info_e3d.h
      if(nlevels .le. 1) return

      call gen_mesh_range(1,ix1,ix2,iy1,iy2,iz1,iz2)
      nx1 = 1 + ix2 - ix1
      ny1 = 1 + iy2 - iy1
      nz1 = 1 + iz2 - iz1
      if(idebug .ge. 10) then
        write (6,*) "iMaxResOut :", iMaxResOut
        write (6,*) "nx1,ny1,nz1:", nx1, ny1, nz1
      endif

      ! Find max level of any box intersecting bounding box that is not
      ! at to high a resolution
      do i = 1, nboxs
        irep = iboxrefine(i)
        ! Only consider box if its res will fit.
        if(max(irep*nx1,irep*ny1,irep*nz1).le.iMaxResOut) then
          ! use level if box also intersects the volume requested
          if(box_in_vol(i)) nrep_used = max(nrep_used, irep)
        endif
      enddo
      nrep_used = min(nrep_used, maxrefine)

      if(idebug .ge. 10) write (6,*) "nrep_used = ", nrep_used

c      nrep_used = 1
c      write (6,*) "for development: have reset nrep_used back to 1"

      return
      end

      
      
C ====================================================================72
      subroutine init_frgb(fxy,fxz,fyz)
      implicit none
      include 'postproc_info_e3d.h'
      real fxy(nxout,nyout,3)
      real fxz(nxout,nzout,3)
      real fyz(nyout,nzout,3)
      integer iv,ix,iy,iz, mefirst
      common /comfillavs/ mefirst

      do iv = 1, 3
        do iy = 1, nyout
        do ix = 1, nxout
          fxy(ix,iy,iv) = 0.0
        enddo
        enddo
        do iz = 1, nzout
        do ix = 1, nxout
          fxz(ix,iz,iv) = 0.0
        enddo
        enddo
        do iz = 1, nzout
        do iy = 1, nyout
          fyz(iy,iz,iv) = 0.0
        enddo
        enddo
      enddo
      mefirst = 0

      return
      end

C ====================================================================72
      subroutine init_avs(fxy,fyz,fxz)
      implicit none
      include 'postproc_info_e3d.h'
      real fxy(nxout,nyout,3), fyz(nyout,nzout,3), fxz(nxout,nzout,3)
      integer iv,ix,iy,iz, mefirst
      common /comfillavs/ mefirst
      data mefirst/1/

        do iv = 1, 3
        do iy = 1, nyout
        do ix = 1, nxout
          fxy(ix,iy,iv) = 0.0
        enddo
        enddo
        do iz = 1, nzout
        do iy = 1, nyout
          fyz(iy,iz,iv) = 0.0
        enddo
        enddo
        do ix = 1, nxout
        do iz = 1, nzout
          fxz(ix,iz,iv) = 0.0
        enddo
        enddo
        enddo
        mefirst = 0

      return
      end

C ====================================================================72
      subroutine fxy_to_prof(fxy,profile)
      implicit none
      include 'postproc_info_e3d.h'
      include 'amr_e3d.h'
      real fxy(nxout,nyout,3)
      real*8  profile(maxprof,maxcols)
      integer i

      if(nxout .gt. nyout) then
        do i = 1,nxout
          profile(i,2) = fxy(i,1,1)
        enddo
      else
        do i = 1,nyout
          profile(i,2) = fxy(1,i,1)
        enddo
      endif

      return
      end

C ====================================================================72
      subroutine fill_vr(fxy,fxz,fyz, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      include 'amr_e3d.h'
      real fxy(nxout,nyout,3)
      real fxz(nxout,nzout,3)
      real fyz(nyout,nzout,3)
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max),fninv
      integer irep, i_wb, mefirst
      common /comfillavs/ mefirst

      if(mefirst .eq. 1) call init_frgb(fxy,fxz,fyz)

      if(nboxsused .le. 0) then
        call fill_vr_uniform(fxy,fxz,fyz, i_wb, flds)
      else
c        irep = iboxrefine(icurrent_box)
c        if(irep .lt. nrep_used) call fill_avs_xy_m(fxy, i_wb, flds)
c        if(irep .eq. nrep_used) call fill_avs_xy_0(fxy, i_wb, flds)
      endif

      return
      end
C ====================================================================72
      subroutine fill_avs(fxy,fyz,fxz, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      include 'amr_e3d.h'
      real fxy(nxout,nyout,3), fyz(nyout,nzout,3), fxz(nxout,nzout,3)
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max),fninv
      integer irep, i_wb, mefirst
      common /comfillavs/ mefirst

      if(mefirst .eq. 1) call init_avs(fxy, fyz, fxz)

      if(nboxsused .le. 0) then
        call fill_avs_uniform(fxy,fyz,fxz, i_wb, flds)
      else
        irep = iboxrefine(icurrent_box)
        if(irep .lt. nrep_used) call fill_avs_xy_m(fxy, i_wb, flds)
        if(irep .eq. nrep_used) call fill_avs_xy_0(fxy, i_wb, flds)
      endif

      return
      end

C ====================================================================72
      subroutine fill_avs_xy_m(fxy, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      include 'amr_e3d.h'
      real fxy(nxout,nyout,3)
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max),fninv
      real wt, wtz
      integer itof, iwrk, ix,iy,iz, ixd,iyd,izd, i_wb, nn,n,m
      integer ixd1,ixd2,iyd1,iyd2,iz1,iz2, nrat, ixreloff, iyreloff
      integer irep, ixo1,ixo2,iyo1,iyo2,izo1,izo2, nDim, iv, iw
      allocatable wtz(:)

      itof = itof_ops(nops)
      iwrk = iWrk_tof(itof)
      nDim = nDim_tof(itof)

      wt = 1.0
c     if(nzout .gt. 1) then
c        write (6,*) "fill_avs_xy_levp: nzout>1 not supported yet"
c        stop
c      endif

      irep = iboxrefine(icurrent_box)
      call gen_mesh_range(irep,ixo1,ixo2, iyo1,iyo2, izo1,izo2)

      iz1 = max(1            , izo1 - ioff_wb(3,i_wb))
      iz2 = min(n0_wb(3,i_wb), izo2 - ioff_wb(3,i_wb))
      if(iz1 .gt. iz2) return

      nrat = nrep_used / irep
      ixreloff = nrat*ioff_wb(1,i_wb) - (ixout1-1)
      iyreloff = nrat*ioff_wb(2,i_wb) - (iyout1-1)
      ixd1 = max(1    , ixreloff + 1)
      iyd1 = max(1    , iyreloff + 1)
      ixd2 = min(nxout, ixreloff + nrat*n0_wb(1,i_wb))
      iyd2 = min(nyout, iyreloff + nrat*n0_wb(2,i_wb))
      if(ixd1.gt.ixd2  .or.  iyd1.gt.iyd2) return

      allocate (wtz(iz1:iz2))
      call gen_zcell_wts(irep, iz1, iz2, ioff_wb(3,i_wb), wtz)

      do iz = iz1, iz2
      izd = iz + ioff_wb(3,i_wb) + 1 - izo1

      do iyd = iyd1, iyd2
      iy =  1 + (iyd-1 - iyreloff)/nrat

      do ixd = ixd1, ixd2
      ix =  1 + (ixd-1 - ixreloff)/nrat

        do iv = 1, nDim
        iw = iwrk + iv - 1
          fxy(ixd,iyd,iv) = fxy(ixd,iyd,iv) + flds(ix,iy,iz,iw)*wtz(iz)
        enddo
      enddo
      enddo
      enddo

      deallocate (wtz)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccc
c      write (6,*) " "
c      write (6,999) "iwxlo,iwxhi2", iwxlo,iwxhi,"ixd1,ixd2", ixd1,ixd2
c      write (6,999) "iwylo,iwyhi2", iwylo,iwyhi,"iyd1,iyd2", iyd1,iyd2
c      write (6,*) "iyout2,nrat",iyout2,nrat
c999   format("fill_avs_xy_m: ",a,2i4,"       ",a,2i4)
c      write (6,*) "ioff_wb(2,i_wb),n0_wb(2,i_wb))",
c     1             ioff_wb(2,i_wb),n0_wb(2,i_wb)
c      write (6,*) " "
cccccccccccccccccccccccccccccccccccccccccccccccc

C ====================================================================72
      subroutine fill_avs_xy_0(fxy, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      include 'amr_e3d.h'
      real fxy(nxout,nyout,3)
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max),fninv
      real wt, wtz
      integer itof, iwrk, ix,iy,iz, ixd,iyd,izd, i_wb, nn,n,m
      integer ix1,ix2,iy1,iy2,iz1,iz2, nmin
      integer irep, ixo1,ixo2,iyo1,iyo2,izo1,izo2, nDim, iv, iw
      allocatable wtz(:)

      itof = itof_ops(nops)
      iwrk = iWrk_tof(itof)
      nDim = nDim_tof(itof)

      wt = 1.0
c      if(nzout .gt. 1) then
c        write (6,*) "fill_avs_xy_levp: nzout>1 not supported yet"
c        stop
c      endif

      irep = iboxrefine(icurrent_box)
      call gen_mesh_range(irep,ixo1,ixo2, iyo1,iyo2, izo1,izo2)

      iz1 = max(1            , izo1 - ioff_wb(3,i_wb))
      iz2 = min(n0_wb(3,i_wb), izo2 - ioff_wb(3,i_wb))
      if(iz1 .gt. iz2) return

      iy1 = max(1            , iyo1 - ioff_wb(2,i_wb))
      iy2 = min(n0_wb(2,i_wb), iyo2 - ioff_wb(2,i_wb))
      if(iy1 .gt. iy2) return

      ix1 = max(1            , ixo1 - ioff_wb(1,i_wb))
      ix2 = min(n0_wb(1,i_wb), ixo2 - ioff_wb(1,i_wb))
      if(ix1 .gt. ix2) return

      allocate (wtz(iz1:iz2))
      call gen_zcell_wts(irep, iz1, iz2, ioff_wb(3,i_wb), wtz)

      do iz = iz1, iz2
      izd = iz + ioff_wb(3,i_wb) + 1 - izo1

      do iy = iy1, iy2
      iyd = iy + ioff_wb(2,i_wb) + 1 - iyo1

      do ix = ix1, ix2
      ixd = ix + ioff_wb(1,i_wb) + 1 - ixo1

        do iv = 1, nDim
        iw = iwrk + iv - 1
          fxy(ixd,iyd,iv) = fxy(ixd,iyd,iv) + flds(ix,iy,iz,iw)*wtz(iz)
        enddo
      enddo
      enddo
      enddo

      deallocate (wtz)

      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      write (6,999) iwxlo,iwxhi,ix1,ix2
c999   format("fill_avs_xy_0: iwxlo,iwxhi=",2i5,"      ix1,ix2=", 2i6)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



C ====================================================================72
      subroutine fill_vr_uniform(fxy, fxz,fyz, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real fxy(nxout,nyout,3)
      real fxz(nxout,nzout,3)
      real fyz(nyout,nzout,3)
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max),fninv
      integer itof, iwrk, ix,iy,iz, ixd,iyd,izd, i_wb, nn,n,m
      integer ix1,ix2,iy1,iy2,iz1,iz2
      integer jx1,jx2,jy1,jy2,jz1,jz2, nDim, iv, iw, ic
      character*64 cc, crange
      real alpha, rgb(3,20), scale, a, ea, eam1, c(3)

      itof = itof_ops(nops)
      iwrk = iWrk_tof(itof)
      nDim = nDim_tof(itof)

      iz1 = max(1            , izout1 - ioff_wb(3,i_wb))
      iz2 = min(n0_wb(3,i_wb), izout2 - ioff_wb(3,i_wb))
      if(iz1 .gt. iz2) return

      iy1 = max(1            , iyout1 - ioff_wb(2,i_wb))
      iy2 = min(n0_wb(2,i_wb), iyout2 - ioff_wb(2,i_wb))
      if(iy1 .gt. iy2) return

      ix1 = max(1            , ixout1 - ioff_wb(1,i_wb))
      ix2 = min(n0_wb(1,i_wb), ixout2 - ioff_wb(1,i_wb))
      if(ix1 .gt. ix2) return

      ! Set rgb array
      do iv = 1, 20
      do ic = 1, 3
        rgb(ic,iv) = 0.0
      enddo
      enddo

      ! set default colors & scaling
      rgb(1,1) = 1.0    ! white
      rgb(2,1) = 1.0
      rgb(3,1) = 1.0
      rgb(1,2) = 1.0    ! red
      rgb(3,3) = 1.0    ! blue
      rgb(2,4) = 1.0    ! green
      scale = 0.25      ! ~4 zones shall be ordier unity opaque

      do iz = iz1, iz2
      izd = iz + ioff_wb(3,i_wb) + 1 - izout1

      do iy = iy1, iy2
      iyd = iy + ioff_wb(2,i_wb) + 1 - iyout1

      do ix = ix1, ix2
      ixd = ix + ioff_wb(1,i_wb) + 1 - ixout1

        alpha = 0.0
        do ic = 1, 3
          c(ic) = 0.0
        enddo
        do iv = 1, nDim
           iw = iwrk + iv - 1
           a = max(0.0, scale * flds(ix,iy,iz,iw))
           alpha = alpha + a
           if(a .lt. 88.0) then
             ea = 1.0 - exp(-a)
           else
             ea = 1.0
           endif
           do ic = 1, 3
             c(ic) = c(ic) + ea * rgb(ic,iv)
           enddo
        enddo
        if(alpha .lt. 88.0) then
          eam1 = exp(-alpha)
        else
          eam1 = 0.0
        endif

        do ic = 1, 3
          fxy(ixd,iyd,ic) = eam1 * fxy(ixd,iyd,ic) + c(ic)
          fxz(ixd,izd,ic) = eam1 * fxz(ixd,izd,ic) + c(ic)
          fyz(iyd,izd,ic) = eam1 * fyz(iyd,izd,ic) + c(ic)
        enddo
      enddo
      enddo
      enddo

      return
      end


C ====================================================================72
      subroutine fill_avs_uniform(fxy,fyz,fxz, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real fxy(nxout,nyout,3), fyz(nyout,nzout,3), fxz(nxout,nzout,3)
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max),fninv
      integer itof, iwrk, ix,iy,iz, ixd,iyd,izd, i_wb, nn,n,m
      integer ix1,ix2,iy1,iy2,iz1,iz2
      integer jx1,jx2,jy1,jy2,jz1,jz2, nDim, iv, iw
      character*64 cc, crange

      itof = itof_ops(nops)
      iwrk = iWrk_tof(itof)
      nDim = nDim_tof(itof)

      iz1 = max(1            , izout1 - ioff_wb(3,i_wb))
      iz2 = min(n0_wb(3,i_wb), izout2 - ioff_wb(3,i_wb))
      if(iz1 .gt. iz2) return

      iy1 = max(1            , iyout1 - ioff_wb(2,i_wb))
      iy2 = min(n0_wb(2,i_wb), iyout2 - ioff_wb(2,i_wb))
      if(iy1 .gt. iy2) return

      ix1 = max(1            , ixout1 - ioff_wb(1,i_wb))
      ix2 = min(n0_wb(1,i_wb), ixout2 - ioff_wb(1,i_wb))
      if(ix1 .gt. ix2) return

      do iz = iz1, iz2
      izd = iz + ioff_wb(3,i_wb) + 1 - izout1

      do iy = iy1, iy2
      iyd = iy + ioff_wb(2,i_wb) + 1 - iyout1

      do ix = ix1, ix2
      ixd = ix + ioff_wb(1,i_wb) + 1 - ixout1

        do iv = 1, nDim
        iw = iwrk + iv - 1
          fxy(ixd,iyd,iv) = fxy(ixd,iyd,iv) + flds(ix,iy,iz,iw)
          fyz(iyd,izd,iv) = fyz(iyd,izd,iv) + flds(ix,iy,iz,iw)
          fxz(ixd,izd,iv) = fxz(ixd,izd,iv) + flds(ix,iy,iz,iw)
        enddo
      enddo
      enddo
      enddo

      return
      end



C ====================================================================72
      subroutine rescale_avs(fxy,fyz,fxz)
      implicit none
      include 'postproc_info_e3d.h'
      real fxy(nxout,nyout,3), fyz(nyout,nzout,3), fxz(nxout,nzout,3)
      integer itof, nDim
      integer  ix,iy,iz,iv
      real fninv

      itof = itof_ops(nops)
      nDim = nDim_tof(itof)

      do iv = 1, nDim
          fninv = 1.0 / float(nzout)
          do iy = 1, nyout
          do ix = 1, nxout
            fxy(ix,iy,iv) = fxy(ix,iy,iv) * fninv
          enddo
          enddo
          fninv = 1.0 / float(nxout)
          do iz = 1, nzout
          do iy = 1, nyout
            fyz(iy,iz,iv) = fyz(iy,iz,iv) * fninv
          enddo
          enddo
          fninv = 1.0 / float(nyout)
          do ix = 1, nxout
          do iz = 1, nzout
            fxz(ix,iz,iv) = fxz(ix,iz,iv) * fninv
          enddo
          enddo
      enddo

      return
      end

C ====================================================================72
      subroutine write_vr(frgb,n1,n2)
      implicit none
      include 'postproc_info_e3d.h'
      integer n1, n2
      real frgb(n1,n2,3)
      integer jx1,jx2,jy1,jy2,jz1,jz2, itof, nDim
      integer  ix,iy,ic
      byte crgb
      allocatable crgb(:,:,:)
      character*256 cscFile
      character*32 cxdim,cydim
      integer n

      if(idebug .ge. 0) write (6,*) "write_vr: n1,n2 = ", n1, n2

      allocate(crgb(3,n1,n2))

      do iy = 1, n2
      do ix = 1, n1
      do ic = 1, 3
        crgb(ic,ix,iy) = max(0,min(255,int(255.999 * frgb(ix,iy,ic))))
      enddo
      enddo
      enddo
#ifndef ISGFORTRAN
      open(unit=11, file=cOutFile, form='binary')
#else
      open(unit=11, file=cOutFile, form='unformatted',access='stream')
#endif
      write (11) crgb
      close(11)
      deallocate(crgb)

C convert -size $1x$2 -depth 8 rgb:$3.raw -flip $3.jpg

      call setblank(cscFile)
      call gen_trimmed_int_string(n1, cxdim)
      call gen_trimmed_int_string(n2, cydim)
      n = len(trim(cOutFile))
      cscFile(1:n+3) = trim(cOutFile) // ".sh"
      write (6,*) "Raw to png conversion script: ", trim(cscFile)
      open(unit=11,file=cscFile,form="formatted")
c998   format(a,a,a,a,a,a,a,a)
999   format("convert -size ",a,"x",a," -flip -depth 8 rgb:",a,
     1        "  ", a,".jpg")
      write (11,999) trim(cxdim), trim(cydim), trim(cOutFile),
     1               trim(cOutFile)
      close(11)

c nxout
      return
      end

C ====================================================================72
      subroutine gen_trimmed_int_string(ival,cval)
      implicit none
      character*(*) cval
      integer ival, i, j, n1, n2
      call setblank(cval)
      write (cval,*) ival
      n2 = len(trim(cval))
      n1 = 1
      do i = 1, n2
        if(cval(i:i) .eq. " ") n1 = i + 1
      enddo
      do i = 1, n2
        j = i + n1 - 1
        if(j .le. n2) then
          cval(i:i) = cval(j:j)
        else
          cval(i:i) = " "
        endif
      enddo

      return
      end
      
C ====================================================================72
      subroutine write_avs(fxy,fyz,fxz)
      implicit none
      include 'postproc_info_e3d.h'
      real fxy(nxout,nyout,3), fyz(nyout,nzout,3), fxz(nxout,nzout,3)
      integer jx1,jx2,jy1,jy2,jz1,jz2, itof, nDim
      integer  ix,iy,iz,iv

c      call rescale_avs(fxy, fyz, fxz)

      itof = itof_ops(nops)
      nDim = nDim_tof(itof)

         ! Trim & mirror axes as needed
        call tm_axes(nxfull,fLimits(1),nxout,xmin,xmax,jx1,jx2)
        call tm_axes(nyfull,fLimits(3),nyout,ymin,ymax,jy1,jy2)
        call tm_axes(nzfull,fLimits(5),nzout,zmin,zmax,jz1,jz2)

        call do_sof(nxout,nyout,nzout,nzfull,cxcoord,cycoord,czcoord,
     1        nDim,jx1,jx2,jy1,jy2, xmin,xmax,ymin,ymax,zmin,zmax,fxy)

        call do_sof(nyout,nzout,nxout,nxfull,cycoord,czcoord,cxcoord,
     1        nDim,jy1,jy2,jz1,jz2, ymin,ymax,zmin,zmax,xmin,xmax,fyz)

        call do_sof(nxout,nzout,nyout,nyfull,cxcoord,czcoord,cycoord,
     1        nDim,jx1,jx2,jz1,jz2, xmin,xmax,zmin,zmax,ymin,ymax,fxz)

      return
      end

c        write (6,*) "==================================="
c        write (6,*) "Debug tm_axes"
c        write (6,*) "nxout = ", nxout
c        write (6,*) "nyout = ", nyout
c        write (6,*) "nzout = ", nzout
c        write (6,*) "jx1,jx2 = ", jx1,jx2
c        write (6,*) "jy1,jy2 = ", jy1,jy2
c        write (6,*) "jz1,jz2 = ", jz1,jz2
c        do iy = 1, nyout
c          write (6,*) "iy,fxy(5,iy) = ", iy,fxy(5,iy)
c        enddo
c        write (6,*) "==================================="

C ====================================================================72
      subroutine do_sof(n1,n2,n3,n3f,c1,c2,c3,nDim,i1,i2,j1,j2,
     1                   amin,amax,bmin,bmax,cmin,cmax,ff)
      implicit none
      include 'postproc_info_e3d.h'
      integer n1,n2,n3,nmin, m,i,nmax, nDim, i1,i2,j1,j2,n3f,nlabel
      real ff(n1,n2,3), amin,amax,bmin,bmax,cmin,cmax
      character*6 ctail
      character*16 c1,c2,c3
      character*256 cc,crange

      nmin = min(n1,n2)
      nmax = max(n1,n2,n3)
      if(nmin .le. 1 .or. 20*nmin .le. nmax) return

      call setblank(cc)
      call setblank(crange)

      cOutFile(iStr0:iStr0+5) = "." // c1(1:1) // c2(1:1) // "sof"
      nlabel = len(trim(cFieldLabel))
      if(n3 .eq. n3f) then
c         m = iStr0 + 2 + len(trim(c1))
c         cc(1:m) = trim(c1) // " (" // cOutFile(1:iStr0-1) // ")"
        m = len(trim(c1)) + nlabel + 15
        cc(1:m) = trim(c1) // " (" // cOutFile(1:10) // "  " //
     1            trim(cFieldLabel)  // ")"
      else
        if(n3 .eq. 1) write (crange,999) c3(1:1), cmin
        if(n3 .ne. 1) write (crange,998) cmin, c3(1:1), cmax
999     format(a1,"="f5.2)
998     format(f5.2"<",a1,"<"f5.2)
c        m = iStr0 + 4 + len(trim(c1)) + len(trim(crange))
c        cc(1:m) = trim(c1) // " (" // cOutFile(1:iStr0-1) // ", "
c     1                             // trim(crange)      // ")"
        m = len(trim(c1)) + nlabel + len(trim(crange)) + 17
        cc(1:m) = trim(c1) // " (" // cOutFile(1:10) // ", " //
     1          trim(cFieldLabel) // "  "  // trim(crange) // ")"
      endif

      if(nDim .eq. 1) then
       ! scalar field only
       call sofplot(ff,n1,n2,nDim,amin,amax,bmin,bmax,cOutFile,
     1         trim(cc),trim(c2), hvmin, hvmax,i1,i2,j1,j2,
     1          ncbarsteps, cbarString, -1.0)
      else if(nDim .eq. 2) then
       ! vector field only
       call gen_vecfield_file(ff(1,1,1),n1,n2,amin,amax,bmin,bmax,i1,j1)
c      call vecplot(...)
      else
       ! both scalar & vector fields
       call gen_vecfield_file(ff(1,1,2),n1,n2,amin,amax,bmin,bmax,i1,j1)
       call sofplot(ff,n1,n2,nDim,amin,amax,bmin,bmax,cOutFile,
     1         trim(cc),trim(c2), hvmin, hvmax,i1,i2,j1,j2,
     1          ncbarsteps, cbarString, -1.0)
      endif

      return
      end
C ====================================================================72
      subroutine tm_axes(nfull,fLim,nout,dmin,dmax,i1,i2)
      implicit none
      integer nfull, nout, i1, i2
      real    fLim(2),dmin,dmax,dx, fLimMin

      if(dmax .gt. fLim(2)) dmax = fLim(2)
      fLimMin = min(fLim(1), -fLim(2))
      if(dmin .lt. fLimMin) dmin = fLimMin
      dx = (fLim(2)-fLim(1)) / float(nfull)
      i1 = 1
      i2 = nout
      if(dmin .lt. fLim(1)) then
        i1 = 1 - int(0.5 + (fLim(1)-dmin) / dx)
      endif

c      i2 = nfull
c dhp: this is the problems with range
c nfull,nout=         104          20
c fLim(1),fLim(2)= -0.1585432       10.83284    
c dmin,dmax= -0.1585000       2.000000    
c dx=  0.1056864    
c i1,i2=           1         104

c      if(nfull .eq. 104) then
c         write (6,*) "nfull,nout=",nfull,nout
c         write (6,*) "fLim(1),fLim(2)=",fLim(1),fLim(2)
c         write (6,*) "dmin,dmax=",dmin,dmax
c         write (6,*) "dx=",dx
c         write (6,*) "i1,i2=",i1,i2
c      endif
c      if(nfull .ge. 0) stop
      return
      end
C ====================================================================72
      subroutine fill_hv(bout, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      character*1 bout(nxfull,nyfull,ntslab)
      integer itof, iwrk, ix,iy,iz, ixd,iyd,izd, i_wb, ival, nzmax
      integer iFullDim(3), iCurDim(3), iCurOff(3), iBlock, imap

      itof = itof_ops(nops)
      iwrk = iWrk_tof(itof)

      do iz = 1, n0_wb(3,i_wb)
      izd = iz + ioff_wb(3,i_wb) - itslab_zoff

      do iy = 1, n0_wb(2,i_wb)
      iyd = iy + ioff_wb(2,i_wb)

      do ix = 1, n0_wb(1,i_wb)
      ixd = ix + ioff_wb(1,i_wb)

        ival = int(hvscale*(flds(ix,iy,iz,iwrk) - hvmin))
        bout(ixd,iyd,izd) = char(min(254, max(1, ival)))
      enddo
      enddo
      enddo

      nzmax = min(nzfull, itslab_zoff+ntslab)
      if(ioff_wb(1,i_wb)+n0_wb(1,i_wb) .ge. nxfull .and.
     1   ioff_wb(2,i_wb)+n0_wb(2,i_wb) .ge. nyfull .and.
     1   ioff_wb(3,i_wb)+n0_wb(3,i_wb) .ge. nzmax         ) then

      iFullDim(1) = nxfull
      iFullDim(2) = nyfull
      iFullDim(3) = nzfull
      iCurDim(1)  = nxfull
      iCurDim(2)  = nyfull
      iCurDim(3)  = ntslab
      iCurOff(1)  = 0
      iCurOff(2)  = 0
      iCurOff(3)  = itslab_zoff
      iBlock      = 128 - 2  ! hack
      imap = 0

      write (6,999) nxfull,nyfull,nzmax,
     1       ioff_wb(1,i_wb)+n0_wb(1,i_wb),
     1       ioff_wb(2,i_wb)+n0_wb(2,i_wb),
     1       ioff_wb(3,i_wb)+n0_wb(3,i_wb)
999   format("XYZslab:", 3i5, "     XYZmax_wb:", 3i5)

      call block_tree(cOutFile, iFullDim, iCurDim, iCurOff,
     1                iBlock, fLimits, imap, bout)

        itslab_zoff = itslab_zoff+ntslab
      endif

      return
      end

C ====================================================================72
      subroutine fill_bof(fout, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real fout(nxbof,nybof,nzbslab)
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, iwrk, ix,iy,iz, ixd,iyd,izd, i_wb, ival
      integer iFullDim(3), iCurDim(3), iCurOff(3), iBlock, imap
      real x0,x1,y0,y1,z0,z1, bfac
      integer nnz, nbytes, nwrote, sw_write, izfull

      itof = itof_ops(nops)
      iwrk = iWrk_tof(itof)

      if(idebug .gt. 1000) then
        write (6,*) "nblend          = ", nb lend
        write (6,*) "n0_wb(3,i_wb)   = ", n0_wb(3,i_wb)
        write (6,*) "ioff_wb(3,i_wb) = ", ioff_wb(3,i_wb)
        write (6,*) "itslab_zoff     = ", itslab_zoff
      endif
      do iz = 1, n0_wb(3,i_wb)
      izd = iz + ioff_wb(3,i_wb) - itslab_zoff

      do iy = 1, n0_wb(2,i_wb)
      iyd = iy + ioff_wb(2,i_wb)

      do ix = 1, n0_wb(1,i_wb)
      ixd = ix + ioff_wb(1,i_wb)

      if(nblend .eq. 1) then
        fout(ixd,iyd,izd) = flds(ix,iy,iz,iwrk)
      else
        ixd = 1 + (ix + ioff_wb(1,i_wb)               - 1) / nblend
        iyd = 1 + (iy + ioff_wb(2,i_wb)               - 1) / nblend
        izd = 1 + (iz + ioff_wb(3,i_wb) - itslab_zoff - 1) / nblend
        fout(ixd,iyd,izd) = fout(ixd,iyd,izd) + flds(ix,iy,iz,iwrk)

c        if((ix + ioff_wb(1,i_wb) - 1).eq.3 .and.
c     1     (iy + ioff_wb(2,i_wb) - 1).eq.3      )
c     1     write (6,902) itslab_zoff, ixd,iyd,izd, flds(ix,iy,iz,iwrk)
c902     format("zoff,ixYZd,val:", 4i5, f7.3)
      endif

      enddo
      enddo
      enddo

      izfull = min(nzfull, itslab_zoff+ntslab)
      if(ioff_wb(1,i_wb)+n0_wb(1,i_wb) .ge. nxfull .and.
     1   ioff_wb(2,i_wb)+n0_wb(2,i_wb) .ge. nyfull .and.
     1   ioff_wb(3,i_wb)+n0_wb(3,i_wb) .ge. izfull       ) then

          bfac = 1.0 / float(nblend**3)
          do izd = 1, nzbslab
          do iyd = 1, nybof
          do ixd = 1, nxbof
            fout(ixd,iyd,izd) = bfac*fout(ixd,iyd,izd)
          enddo
          enddo
          enddo

          nnz = min(ntslab, nzfull - itslab_zoff)
c           Full Dimensions:    2052    486   1134
c           blockdim             342     81    189

          nbytes = 4 * nxbof * nybof * (nnz/nblend)
          nwrote = sw_write(cOutHost,cOutDir,cOutFile, i64output,
     1                      nbytes, fout)

          write (6,991) itslab_zoff+1, itslab_zoff+nnz, nbytes
991       format("bof writing: ", 3i8)

          if(nblend .gt. 1) then
            do izd = 1, nzbslab
            do iyd = 1, nybof
            do ixd = 1, nxbof
              fout(ixd,iyd,izd) = 0.0
            enddo
            enddo
            enddo
          endif
          i64output   = i64output   + n64slabsize
          itslab_zoff = itslab_zoff + ntslab

         if(idebug .ge. 10) then
           write (6,*) 'BOF: nz-planes written=', nnz
         endif

         if(itslab_zoff .ge. nzfull) then
           call get_full_XYZrange(x0,x1,y0,y1,z0,z1)
           write (6,*) "nXYZbof = ", nxbof, nybof, nzbof
           open(unit=12,file="fieldformat", form="formatted")
           write (12,992) "xcoord     ", x0, x1, nxbof
           write (12,992) "ycoord     ", y0, y1, nybof
           write (12,992) "zcoord     ", z0, z1, nzbof
992        format(a11, 2f14.6, i6)
           write (12,993) "blockdim   ",nxbof,nybof,nzbof
993        format(a11, 3i6)
           close(12)
           if(idebug .ge. 1) then
             write (6,*) ' '
             write (6,*) 'Output file: ', trim(cOutFile)
           endif
         endif
      endif

      return
      end





C ====================================================================72
      subroutine fill_spc3v(ft3v, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real ft3v(nxbof,nybof,ntslab,3)
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, iwrk, ix,iy,iz, ixd,iyd,izd, i_wb, ival
      integer iFullDim(3), iCurDim(3), iCurOff(3), iBlock, imap
      real x0,x1,y0,y1,z0,z1, bfac
      integer nnz, nbytes, nwrote, sw_write, iv, iwrk0
      real*8 spc(0:4096,3)
      integer icounts(0:4096)

      itof = itof_ops(nops)
      iwrk = iWrk_tof(itof)

      do iv = 1, 3
      iwrk0 = iwrk + iv - 1
      do iz = 1, n0_wb(3,i_wb)
      izd = iz + ioff_wb(3,i_wb) - itslab_zoff

      do iy = 1, n0_wb(2,i_wb)
      iyd = iy + ioff_wb(2,i_wb)

      do ix = 1, n0_wb(1,i_wb)
      ixd = ix + ioff_wb(1,i_wb)

        ft3v(ixd,iyd,izd,iv) = flds(ix,iy,iz,iwrk0)

      enddo
      enddo
      enddo
      enddo

      if(ioff_wb(1,i_wb)+n0_wb(1,i_wb) .ge. nxfull .and.
     1   ioff_wb(2,i_wb)+n0_wb(2,i_wb) .ge. nyfull .and.
     1   ioff_wb(3,i_wb)+n0_wb(3,i_wb) .ge. itslab_zoff+ntslab) then

            call fftxy(nxbof,nybof,ntslab,nzbof,nvec,itslab_zoff,ft3v)

          itslab_zoff = itslab_zoff + ntslab

         if(idebug .gt. 50) then
           write (6,*) 'SPC3V: itslab_zoff = ', itslab_zoff
         endif

         if(itslab_zoff .eq. nzfull) then
           call do_fft_z_spc(nzfull, nvec, icounts, spc)
           call write_spc3v(icounts, spc)
           call get_full_XYZrange(x0,x1,y0,y1,z0,z1)
         endif
      endif

      return
      end

C ====================================================================72
      subroutine write_spc3v(icounts, spc)
      implicit none
      include 'postproc_info_e3d.h'
      real*8 spc(0:4096,3), iT0,ECsum, ETsum
      real*8 pi, fkav, fkmin, fkmax, fac, ET0
      integer icounts(0:4096), imax, i, j
      real time

      imax = nzfull/2 - 1
      write (6,*) 'Openning file ', cOutFile(1:20)
      open(unit=11, file=cOutFile, form='formatted')

      ET0 = spc(0,3)
      ECsum = 0.0d0
      ETsum = 0.0d0
      do i = 1, imax
        ETsum = ETsum + spc(i,3)
        ECsum = ECsum + spc(i,1)
      enddo
      call get_time(time)
      write (11,954) time,ECsum,ETsum-ECsum,ETsum,ET0
954   format("# time = ", f10.5, /,
     1       "# Ecomp, Esoln, Etotal = ", 3e15.6, /,
     1       "# E0 = ", e15.6)

      write (11,955)
955   format("#    kav     EC_norm     ES_norm     ET_norm",
     1                   "     EC_orig     ES_orig     ET_orig")

      pi = 4.0d0 * atan(1.0)
      do i = 1, imax
        fkav  = dble(i)
        fkmin = fkav - 0.5d0
        fkmax = fkav + 0.5d0
        fac = (2.0d0*pi/3.0d0)*(fkmax**3 - fkmin**3)/dble(icounts(i))
        spc(i,2) = spc(i,3) - spc(i,1)
        write (11,958) fkav, (fac*spc(i,j), j=1,3), (spc(i,j), j=1,3)
958     format(f8.2, 1p6e12.4)
      enddo
      close(11)

      write (6,*) ' '
      write (6,*) 'Output:  ', trim(cOutFile)
      write (6,*) ' '

      return
      end

C ====================================================================72
      subroutine fill_spcnv(ftnv, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real ftnv(nxbof,nybof,ntslab,nvec)
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, iwrk, ix,iy,iz, ixd,iyd,izd, i_wb, ival
      integer iFullDim(3), iCurDim(3), iCurOff(3), iBlock, imap
      real x0,x1,y0,y1,z0,z1, bfac
      integer nnz, nbytes, nwrote, sw_write, iv, iwrk0
      real*8 spc(0:4096,3)
      integer icounts(0:4096)

      itof = itof_ops(nops)
      iwrk = iWrk_tof(itof)

      do iv = 1, nvec
      iwrk0 = iwrk + iv - 1
      do iz = 1, n0_wb(3,i_wb)
      izd = iz + ioff_wb(3,i_wb) - itslab_zoff

      do iy = 1, n0_wb(2,i_wb)
      iyd = iy + ioff_wb(2,i_wb)

      do ix = 1, n0_wb(1,i_wb)
      ixd = ix + ioff_wb(1,i_wb)

        ftnv(ixd,iyd,izd,iv) = flds(ix,iy,iz,iwrk0)

      enddo
      enddo
      enddo
      enddo

      if(ioff_wb(1,i_wb)+n0_wb(1,i_wb) .ge. nxfull .and.
     1   ioff_wb(2,i_wb)+n0_wb(2,i_wb) .ge. nyfull .and.
     1   ioff_wb(3,i_wb)+n0_wb(3,i_wb) .ge. itslab_zoff+ntslab) then

            call fftxy(nxbof,nybof,ntslab,nzbof,nvec,itslab_zoff,ftnv)

          itslab_zoff = itslab_zoff + ntslab

         if(idebug .gt. 50) then
           write (6,*) 'spcnv: itslab_zoff, nvec = ', itslab_zoff,nvec
         endif

         if(itslab_zoff .eq. nzfull) then
           call do_fft_z_spc(nzfull, nvec, icounts, spc)
           if(nvec .eq. 1) call write_spc1s(icounts, spc)
           if(nvec .eq. 3) call write_spc3v(icounts, spc)
           call get_full_XYZrange(x0,x1,y0,y1,z0,z1)
         endif
      endif

      return
      end

C ====================================================================72
      subroutine write_spc1s(icounts, spc)
      implicit none
      include 'postproc_info_e3d.h'
      real*8 spc(0:4096,3), iT0,ECsum, ETsum
      real*8 pi, fkav, fkmin, fkmax, fac, ET0
      integer icounts(0:4096), imax, i, j
      real time

      imax = nzfull/2 - 1
      write (6,*) 'Openning file ', cOutFile(1:20)
      open(unit=11, file=cOutFile, form='formatted')

      ET0 = spc(0,3)
      ETsum = 0.0d0
      do i = 1, imax
        ETsum = ETsum + spc(i,3)
      enddo
      call get_time(time)
      write (11,954) time,ETsum,ET0
954   format("# time = ", f10.5, /,
     1       "# Etotal = ", e15.6, /,
     1       "# E0     = ", e15.6)

      write (11,955)
955   format("#    kav     ET_norm",
     1                   "     ET_orig")

      pi = 4.0d0 * atan(1.0)
      do i = 1, imax
        fkav  = dble(i)
        fkmin = fkav - 0.5d0
        fkmax = fkav + 0.5d0
        fac = (2.0d0*pi/3.0d0)*(fkmax**3 - fkmin**3)/dble(icounts(i))
        write (11,958) fkav, fac*spc(i,3), spc(i,3)
958     format(f8.2, 1p2e12.4)
      enddo
      close(11)

      write (6,*) ' '
      write (6,*) 'Output:  ', trim(cOutFile)
      write (6,*) ' '

      return
      end

C ====================================================================72
      subroutine fill_strfn(fnv, i_wb, flds, idispsumS)
      implicit none
      include 'postproc_info_e3d.h'
      real fnv(nxbof,nybof,ntslab,nvec)
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer idispsumS(-ndiff_sampl:ndiff_sampl,
     1                   max_disp_in_set, ndisp_sets)
      integer itof, iwrk, ix,iy,iz, ixd,iyd,izd, i_wb, ival
      integer iFullDim(3), iCurDim(3), iCurOff(3), iBlock, imap
      real x0,x1,y0,y1,z0,z1, bfac
      integer nnz, nbytes, nwrote, sw_write, iv, iwrk0
      real*8 spc(0:4096,3)
      integer icounts(0:4096)

      itof = itof_ops(nops)
      iwrk = iWrk_tof(itof)

      do iv = 1, nvec
      iwrk0 = iwrk + iv - 1
      do iz = 1, n0_wb(3,i_wb)
      izd = iz + ioff_wb(3,i_wb) - itslab_zoff

      do iy = 1, n0_wb(2,i_wb)
      iyd = iy + ioff_wb(2,i_wb)

      do ix = 1, n0_wb(1,i_wb)
      ixd = ix + ioff_wb(1,i_wb)

        fnv(ixd,iyd,izd,iv) = flds(ix,iy,iz,iwrk0)

      enddo
      enddo
      enddo
      enddo

      if(ioff_wb(1,i_wb)+n0_wb(1,i_wb) .ge. nxfull .and.
     1   ioff_wb(2,i_wb)+n0_wb(2,i_wb) .ge. nyfull .and.
     1   ioff_wb(3,i_wb)+n0_wb(3,i_wb) .ge. itslab_zoff+ntslab) then


c >>> do structure function sums here
         if(nvec .eq. 1) call sum_strfnsc(fnv, idispsumS)
c         if(nvec .eq. 3) call sum_strfn3v(icounts, spc)

          itslab_zoff = itslab_zoff + ntslab

         if(idebug .gt. 50) then
           write (6,*) 'fill_strfn: itslab_zoff,nvec=',itslab_zoff,nvec
         endif

         if(itslab_zoff .eq. nzfull) then
c >>> do structure function output here
           if(nvec .eq. 1) call write_strfnsc(idispsumS)
c           if(nvec .eq. 3) call write_strfn3v(icounts, spc)
         endif
      endif

      return
      end

C ====================================================================72
      subroutine sum_strfnsc(fnv, idispsumS)
      implicit none
      include 'postproc_info_e3d.h'
      real fnv(nxbof,nybof,ntslab,nvec)
      integer idispsumS(-ndiff_sampl:ndiff_sampl,
     1                   max_disp_in_set, ndisp_sets)
      integer i,k, ixd,iyd, ix,iy,iz, ixp, iyp, idel
      real dispmagi, delxnorm, delynorm, v1, v2, del, delvl


      do k=1,ndisp_sets
      do i=1, ndisp_in_set(k)
      ixd = ixdisp(i,k)
      iyd = iydisp(i,k)
      dispmagi = 1.0 / sqrt(dble(ixd**2 + iyd**2))
      delxnorm = float(ixd) * dispmagi
      delynorm = float(iyd) * dispmagi
      do iz = 1, ntslab
      do iy = 1, nybof
      do ix = 1, nxbof
        v1 = fnv(ix,iy, iz, 1)

        ixp = ix + ixd
        iyp = iy + iyd
        if(ixp .gt. nxbof) ixp = ixp - nxbof  ! assume periodic in X
        if(iyp .gt. nybof) iyp = iyp - nybof  ! assume periodic in Y
        v2 = fnv(ixp, iyp, iz, 1)

        del = v2 - v1
        idel = int(del * del_vari)
        if(delvL .lt. 0.0) idel = idel-1
        idel = max(-ndiff_sampl, min(ndiff_sampl, idel))
        idispsumS(idel,i,k) = idispsumS(idel,i,k) + 1
      enddo
      enddo
      enddo

      enddo
      enddo

      return
      end

C ====================================================================72
      subroutine write_strfnsc(idispsumS)
      implicit none
      include 'postproc_info_e3d.h'
      integer idispsumS(-ndiff_sampl:ndiff_sampl,
     1                   max_disp_in_set, ndisp_sets)
      integer k, i, ilen, nn, ms, ns, nd
      character*256 cOut
      real time

      call get_time(time)

      do k=1,ndisp_sets
        ! Generate root output file name for this set
        do i=1,256
          cOut(i:i) = ' '
        enddo
        ilen = len(trim(cOutFile)) + len(trim(cdisp_names(k))) + 3
        cOut(1:ilen) = trim(cOutFile)//"-"//trim(cdisp_names(k))//"--"
        write (6,*) len(trim(cOutFile)), "  ", trim(cOutFile)
        write (6,*) len(trim(cdisp_names(k))),"  ",trim(cdisp_names(k))
        write (6,*) "ilen = ", ilen
        write (6,*) "cOut = ", trim(cOut)

        ! idispsum* array dimensions
        nn = ndiff_sampl
        ms= max_disp_in_set
        ns = ndisp_sets
        nd = ndisp_in_set(k)
        call strfnout(-nn,nn,ms,ns,k,nd,idispsumS,del_vari,cOut,"S",
     1                   ixdisp, iydisp, dx, time)
      enddo


      return
      end

C ====================================================================72
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

C ====================================================================72
      subroutine get_displacements(ndisp_sets, cdisp_names,
     1           ndisp_in_set,ndiff_bin,ixdisp,iydisp)
      parameter (MAXDISPSETS=10)
      parameter (MAXNDISPINSET=100)
      character*256 cdisp_names(MAXDISPSETS)
      integer ndisp_in_set(MAXDISPSETS)
      integer ixdisp(MAXNDISPINSET,MAXDISPSETS)
      integer iydisp(MAXNDISPINSET,MAXDISPSETS)
      character*256 cline, ckey

      ndisp_sets = 0
      open(unit=12, file="structfn_XYdisp", form="formatted")
10    continue
      do i = 1, 256
        cline(i:i) = ' '
      enddo
      read (12,990,end=90) cline
990   format(a256)
      if(cline(1:10).eq.'ndiff_bin ') read (cline,*) ckey,ndiff_bin
      if(cline(1:10).eq.'disp_set  ') then
        ndisp_sets = ndisp_sets + 1
        k = ndisp_sets
        ndisp_in_set(k) = 0
        do i = 1, 256
          cdisp_names(k)(i:i) = ' '
        enddo
        read (cline,*) ckey, cdisp_names(k)
      endif
      if(cline(1:10).eq.'xy_disp   ') then
        ndisp_in_set(k) = ndisp_in_set(k) + 1
        n = ndisp_in_set(k)
        read (cline,*) ckey, ixdisp(n,k), iydisp(n,k)
      endif
      go to 10
90    continue
      close(12)

      if(ndisp_sets .eq. 0) then
        write (6,*) "ndisp_sets not set in file structfn_XYdisp"
        stop
      endif

      return
      end

C ====================================================================72
c >>>
C ====================================================================72
      subroutine fill_strip(profile, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real*8  profile(maxprof,maxcols), val
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, iwrk, ixd, i_wb, j, ifld, i, ixyz(3)
      real wb_min(3), wb_max(3), s0,s1, x0, s, dxyz(3), fxyz(3),xyz(3)
      real tri_lin_interp, se

      itof = itof_ops(nops)
      iwrk = iWrk_tof(itof)
      do j = 1, 3
        x0 = flimits(2*j - 1)
        wb_min(j) = dx * float(ioff_wb(j,i_wb)) + x0
        wb_max(j) = dx * float(ioff_wb(j,i_wb) + n0_wb(j,i_wb)) + x0
      enddo

      s0 = -1.0
      s1 =  1.0
      do j = 1, 3
      call overlap_1d(s0,s1, wb_min(j),wb_max(j), endpt0(j),endpt1(j))
      enddo
      if(s0 .ge. s1) return
      s0 = s0 * extent
      s1 = s1 * extent

      if(idebug .gt. 50) then
c      write (6,*) "============================"
c      write (6,999) "wb_max : ", (wb_max(j),j=1,3)
c      write (6,999) "wb_min : ", (wb_min(j),j=1,3)
c      write (6,999) "endpt0 : ", (endpt0(j),j=1,3)
c      write (6,999) "endpt1 : ", (endpt1(j),j=1,3)
c999   format(a, 3f7.3, "   ", 3f7.3, "   s = ", f7.3)
c      write (6,*) "s0,s1 = ", s0, s1
        write (6,*) "  "
        write (6,*) "In fill_strip:"
        write (6,999) (wb_min(j),wb_max(j), j=1,3)
999     format("X:", 2f10.5, "    Y:",2f10.5,"    Z:", 2f10.5)
      endif

      do i = 1, maxprof
        s = profile(i,1)
        if(s .ge. s0 .and. s .le. s1) then
        se = s + extent
        do j = 1, 3
           xyz(j) =  endpt0(j) + se * endpt_vec(j)
          dxyz(j) = (xyz(j) - wb_min(j))/dx
          ixyz(j) = int(dxyz(j) + 0.5)
          fxyz(j) = dxyz(j)+0.5 - float(ixyz(j))
        enddo
         
        if(idebug .gt. 50) then
          write (6,998) (xyz(j),dxyz(j),ixyz(j), j=1,3)
998       format(     "Xvdi:",2f7.3,i4,"     Yvdi:",2f7.3,i4,
     1           "     Zvdi:",2f7.3,i4)
        endif
        profile(i,2) =  tri_lin_interp(flds,iwrk,ixyz,fxyz)

c        write (6,999) "profile: ", (profile(i,j), j=1,2)
        endif
      enddo
c      write (6,*) "============================"

      return
      end

c      write (6,*) "============================"
c      write (6,999) (wb_min(j), wb_max(j), j=1,3), iwxlo,iwxhi
c999   format("wb ranges: ", 3(2f7.3, "   "), "   iwxlohi=", 2i4)

c      write (6,999) "wb_max : ", (wb_max(j),j=1,3)
c      write (6,999) "wb_min : ", (wb_min(j),j=1,3)
c      write (6,999) "endpt0 : ", (endpt0(j),j=1,3)
c      write (6,999) "endpt1 : ", (endpt1(j),j=1,3)
c999   format(a, 3f7.3)
c      write (6,*) "s0,s1 = ", s0, s1
c      write (6,*) "  "
c      write (6,*) "============================"

C ====================================================================72
      real function tri_lin_interp(flds,iwrk,ixyz,fxyz)
      implicit none
      include 'postproc_info_e3d.h'
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      real    fxyz(3), fx0,fy0,fz0, fx1,fy1,fz1
      integer ixyz(3), ix0,iy0,iz0, ix1,iy1,iz1, iwrk

      fx1 = fxyz(1)
      fy1 = fxyz(2)
      fz1 = fxyz(3)
      fx0 = 1.0 - fx1
      fy0 = 1.0 - fy1
      fz0 = 1.0 - fz1
      ix0 = ixyz(1)
      iy0 = ixyz(2)
      iz0 = ixyz(3)
      ix1 = ix0 + 1
      iy1 = iy0 + 1
      iz1 = iz0 + 1

      tri_lin_interp = 
     1                fx0*fy0*fz0 * flds(ix0,iy0,iz0, iwrk)
     1               + fx1*fy0*fz0 * flds(ix1,iy0,iz0, iwrk)
     1               + fx0*fy1*fz0 * flds(ix0,iy1,iz0, iwrk)
     1               + fx1*fy1*fz0 * flds(ix1,iy1,iz0, iwrk)
     1               + fx0*fy0*fz1 * flds(ix0,iy0,iz1, iwrk)
     1               + fx1*fy0*fz1 * flds(ix1,iy0,iz1, iwrk)
     1               + fx0*fy1*fz1 * flds(ix0,iy1,iz1, iwrk)
     1               + fx1*fy1*fz1 * flds(ix1,iy1,iz1, iwrk)

      return
      end


C ====================================================================72
      subroutine overlap_1d(s0,s1, b0,b1, e0, e1)
      implicit none
      real s0,s1, b0,b1, e0, e1, emin,emax, de, dsde, sb0,sb1

      de   = e1 - e0
      if(de .eq. 0.0) then
        if(b0 .le. e0  .and.  e0 .lt. b1) return
        s0 = 2.0
        return
      endif
      dsde = 2.0 / de
      sb0 = dsde * (b0 - e0) - 1.0
      sb1 = dsde * (b1 - e0) - 1.0
      s0 = max(s0, min(sb0, sb1))
      s1 = min(s1, max(sb0, sb1))

      return
      end
C ====================================================================72
      subroutine fill_xstrip(profile, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real*8  profile(maxprof,maxcols), val
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, iwrk, ix, iy, iz, ixd, i_wb, j

      itof = itof_ops(nops)
      iwrk = iWrk_tof(itof)
      iz = iz00 - ioff_wb(3,i_wb)
      iy = iy00 - ioff_wb(2,i_wb)
      if(iy .lt. 1  .or.  iy .gt. n0_wb(2,i_wb)) return
      if(iz .lt. 1  .or.  iz .gt. n0_wb(3,i_wb)) return
      if(idebug .gt. 1000) then
        write (6,999) iy00,iz00,iy,iz,(ioff_wb(j,i_wb),j=1,3)
      endif
999   format("iYZ00:", 2i4, "     iYZ:", 2i4, "     off:",3i4)
      do ix = 1, n0_wb(1,i_wb)
        ixd = ix + ioff_wb(1,i_wb) + 1 - ixout1
        if(ixd .ge. 1  .and.  ixd .le. nxout) then
        do j=2,maxcols
          profile(ixd,j) = flds(ix,iy,iz,iwrk+j-2)
        enddo
        endif
      enddo

      return
      end

C ====================================================================72
      subroutine fill_ystrip(profile, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real*8  profile(maxprof,maxcols), val
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, iwrk, ix, iy, iz, iyd, i_wb, j

      itof = itof_ops(nops)
      iwrk = iWrk_tof(itof)
      iz = iz00 - ioff_wb(3,i_wb)
      ix = ix00 - ioff_wb(1,i_wb)
      if(ix .lt. 1  .or.  ix .gt. n0_wb(1,i_wb)) return
      if(iz .lt. 1  .or.  iz .gt. n0_wb(3,i_wb)) return
      if(idebug .gt. 1000) then
        write (6,999) ix00,iz00,ix,iz,(ioff_wb(j,i_wb),j=1,3)
      endif
999   format("iXZ00:", 2i4, "     iXZ:", 2i4, "     off:",3i4)
      do iy = 1, n0_wb(2,i_wb)
        iyd = iy + ioff_wb(2,i_wb) + 1 - iyout1
        if(iyd .ge. 1  .and.  iyd .le. nyout) then
        do j=2,maxcols
          profile(iyd,j) = flds(ix,iy,iz,iwrk+j-2)
        enddo
        endif
      enddo

      return
      end

C ====================================================================72
      subroutine fill_xyslice(bout, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      character bout(nxfull,nyfull)
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, iwrk, ix, iy, iz, ixd, iyd, i_wb, ival

      itof = itof_ops(nops)
      iwrk = iWrk_tof(itof)
      iz = iz00 - ioff_wb(3,i_wb)
      if(iz .lt. 1  .or.  iz .gt. n0_wb(3,i_wb)) return

      do iy = 1, n0_wb(2,i_wb)
      iyd = iy + ioff_wb(2,i_wb)

      do ix = 1, n0_wb(1,i_wb)
      ixd = ix + ioff_wb(1,i_wb)

        ival = int(hvscale*(flds(ix,iy,iz,iwrk) - hvmin))
        bout(ixd,iyd) = char(min(254, max(1, ival)))

      enddo
      enddo

      return
      end

C ====================================================================72
      subroutine fill_rprof(profile, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real*8  profile(maxprof,maxcols), val
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, iwrk, ix, iy, iz, ird, i_wb
      real xx, yy, zz, xyzmid0

      if(idebug .gt. 10) write (6,*) "top of fill_rprof"

      itof = itof_ops(nops)
      iwrk = iWrk_tof(itof)
      xyzmid0 = 0.5 * float(nxfull + 1)
      do iz = 1, n0_wb(3,i_wb)
      zz = float(iz + ioff_wb(3,i_wb)) - xyzmid0
      do iy = 1, n0_wb(2,i_wb)
      yy = float(iy + ioff_wb(2,i_wb)) - xyzmid0
      do ix = 1, n0_wb(1,i_wb)
        xx = float(ix + ioff_wb(1,i_wb)) - xyzmid0
        ird = 1 + int(sqrt(xx**2 + yy**2 + zz**2))
        if(ird .le. nxfull/2) then
          val = flds(ix,iy,iz,iwrk)
          profile(ird,2) = profile(ird,2) + val
          profile(ird,3) = profile(ird,3) + val**2
          profile(ird,7) = profile(ird,7) + 1.0d+00
          if(profile(ird,5) .gt. val) profile(ird,5) = val
          if(profile(ird,6) .lt. val) profile(ird,6) = val
        endif
      enddo
      enddo
      enddo

      if(idebug .gt. 10) write (6,*) "bot of fill_rprof"

      return
      end


C ====================================================================72
      subroutine fill_xprof(profile, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real*8  profile(maxprof,maxcols), val
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, iwrk, ix, iy, iz, ixd, i_wb, ix1,ix2,iy1,iy2,iz1,iz2

      itof = itof_ops(nops)
      iwrk = iWrk_tof(itof)

      ix1 = max(            1, ixout1-ioff_wb(1,i_wb))
      ix2 = min(n0_wb(1,i_wb), ixout2-ioff_wb(1,i_wb))
      iy1 = max(            1, iyout1-ioff_wb(2,i_wb))
      iy2 = min(n0_wb(2,i_wb), iyout2-ioff_wb(2,i_wb))
      iz1 = max(            1, izout1-ioff_wb(3,i_wb))
      iz2 = min(n0_wb(3,i_wb), izout2-ioff_wb(3,i_wb))

      do iz = iz1, iz2
      do iy = iy1, iy2
      do ix = ix1, ix2
        ixd = ix + ioff_wb(1,i_wb) - (ixout1-1)   ! << subtrace x offset & boundary depth here
        val = flds(ix,iy,iz,iwrk)
        profile(ixd,2) = profile(ixd,2) + val
        profile(ixd,3) = profile(ixd,3) + val**2
        profile(ixd,7) = profile(ixd,7) + 1.0d+00
        if(profile(ixd,5) .gt. val) profile(ixd,5) = val
        if(profile(ixd,6) .lt. val) profile(ixd,6) = val
      enddo
      enddo
      enddo

      if(idebug .gt. 1000) then
        write (6,*) 'iwrk = ', iwrk
        write (6,*) 'iz,ixd,flds(1,1,iz,iwrk)'
        do ix = 1, n0_wb(1,i_wb)
        ixd = ix + ioff_wb(1,i_wb)
          write (6,*) ix,ixd,flds(ix,1,1,iwrk)
        enddo
        stop
      endif

      return
      end

C ====================================================================72
      subroutine fill_yprof(profile,  i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real*8  profile(maxprof,maxcols), val
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, iwrk, ix, iy, iz, ixd, iyd, i_wb
      integer  ix1,ix2,iy1,iy2,iz1,iz2
      real*8 area(nxfull), afac

      itof = itof_ops(nops)
      iwrk = iWrk_tof(itof)

      ix1 = max(            1, ixout1-ioff_wb(1,i_wb))
      ix2 = min(n0_wb(1,i_wb), ixout2-ioff_wb(1,i_wb))
      iy1 = max(            1, iyout1-ioff_wb(2,i_wb))
      iy2 = min(n0_wb(2,i_wb), iyout2-ioff_wb(2,i_wb))
      iz1 = max(            1, izout1-ioff_wb(3,i_wb))
      iz2 = min(n0_wb(3,i_wb), izout2-ioff_wb(3,i_wb))

      if(isRZ .eq. 1) then
        afac = 4.0d0 * atan(1.0d0) * dble(dx**2)
        do ixd = 1, nxfull
          area(ixd) = afac * ((dble(ixd))**2 - (dble(ixd-1))**2)
        enddo
      else
        do ixd = 1, nxfull
          area(ixd) = dx * dz
        enddo
      endif

      do iz = iz1, iz2
      do iy = iy1, iy2
      iyd = iy + ioff_wb(2,i_wb) - (iyout1-1)   ! << subtrace x offset & boundary depth here
      do ix = ix1, ix2
        ixd = ix + ioff_wb(1,i_wb) - (ixout1-1) ! << subtrace x offset & boundary depth here
        val = flds(ix,iy,iz,iwrk)
        profile(iyd,2) = profile(iyd,2) + area(ixd) * val
        profile(iyd,3) = profile(iyd,3) + area(ixd) * val**2
        profile(iyd,7) = profile(iyd,7) + area(ixd)
        if(profile(iyd,5) .gt. val) profile(iyd,5) = val
        if(profile(iyd,6) .lt. val) profile(iyd,6) = val
      enddo
      enddo
      enddo

c      write (6,*) "mid of fill_yprof"

      if(idebug .gt. 2000) then
        write (6,*) 'iwrk = ', iwrk
        write (6,*) 'iy,iyd,flds(1,1,iz,iwrk)'
        do iy = 1, n0_wb(2,i_wb)
        iyd = iy + ioff_wb(2,i_wb)
          write (6,*) iz,iyd,flds(1,iy,1,iwrk)
        enddo
        stop
      endif

c      write (6,*) "bot of fill_yprof"
c      write (6,*) "==================================="

      return
      end

C ====================================================================72
      subroutine fill_corr(icorrelate, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      integer icorrelate(nbin1,nbin2), i_wb
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, iwrk1, iwrk2, nDim, ix,iy,iz, ibin1,ibin2
      integer ix1,ix2,iy1,iy2,iz1,iz2
      real s1, s2

      itof  = itof_ops(nops)
      iwrk1 = iWrk_tof(itof)
      iwrk2 = iwrk1 + 1
      nDim  = nDim_tof(itof)
      s1 = float(nbin1) / (bin11 - bin10)
      s2 = float(nbin2) / (bin21 - bin20)

      ix1 = max(            1, ixout1-ioff_wb(1,i_wb))
      ix2 = min(n0_wb(1,i_wb), ixout2-ioff_wb(1,i_wb))
      iy1 = max(            1, iyout1-ioff_wb(2,i_wb))
      iy2 = min(n0_wb(2,i_wb), iyout2-ioff_wb(2,i_wb))
      iz1 = max(            1, izout1-ioff_wb(3,i_wb))
      iz2 = min(n0_wb(3,i_wb), izout2-ioff_wb(3,i_wb))

      do iz = iz1, iz2
      do iy = iy1, iy2
      do ix = ix1, ix2
        ibin1 = int(s1 * (flds(ix,iy,iz,iwrk1) - bin10))
        ibin1 = max(1,min(nbin1, ibin1))
        ibin2 = int(s2 * (flds(ix,iy,iz,iwrk2) - bin20))
        ibin2 = max(1,min(nbin2, ibin2))
        icorrelate(ibin1,ibin2) = icorrelate(ibin1,ibin2) + 1
      enddo
      enddo
      enddo

      return
      end

C ====================================================================72
      subroutine fill_view(sym33, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real*8 sym33(3,3)
      integer i_wb
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      real*8 v1,v2,v3
      integer ic1, ic2, ic3, nDim, ndim3, ix,iy,iz, iv
      integer itof, iwrk0, ix1,ix2,iy1,iy2,iz1,iz2

      itof  = itof_ops(nops)
      iwrk0 = iWrk_tof(itof)
      nDim  = nDim_tof(itof)
      nDim3 = nDim / 3

      ix1 = max(            1, ixout1-ioff_wb(1,i_wb))
      ix2 = min(n0_wb(1,i_wb), ixout2-ioff_wb(1,i_wb))
      iy1 = max(            1, iyout1-ioff_wb(2,i_wb))
      iy2 = min(n0_wb(2,i_wb), iyout2-ioff_wb(2,i_wb))
      iz1 = max(            1, izout1-ioff_wb(3,i_wb))
      iz2 = min(n0_wb(3,i_wb), izout2-ioff_wb(3,i_wb))

      do iv = 0, nDim3-1
        ic1 = iWrk0 + iv
        ic2 = ic1 +   nDim3
        ic3 = ic1 + 2*nDim3
        do iz = iz1, iz2
        do iy = iy1, iy2
        do ix = ix1, ix2
          v1 = flds(ix,iy,iz,ic1)
          v2 = flds(ix,iy,iz,ic2)
          v3 = flds(ix,iy,iz,ic3)
          sym33(1,1) = sym33(1,1) + v1*v1
          sym33(1,2) = sym33(1,2) + v1*v2
          sym33(1,3) = sym33(1,3) + v1*v3
          sym33(2,2) = sym33(2,2) + v2*v2
          sym33(2,3) = sym33(2,3) + v2*v3
          sym33(3,3) = sym33(3,3) + v3*v3
        enddo
        enddo
        enddo
      enddo

      return
      end

C ====================================================================72
      subroutine fill_dist(profile, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real*8  profile(maxprof,maxcols)
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, iwrk, ix, iy, iz,  i_wb, ibin
      integer ix1,ix2,iy1,iy2,iz1,iz2

      itof = itof_ops(nops)
      iwrk = iWrk_tof(itof)

      ix1 = max(            1, ixout1-ioff_wb(1,i_wb))
      ix2 = min(n0_wb(1,i_wb), ixout2-ioff_wb(1,i_wb))
      iy1 = max(            1, iyout1-ioff_wb(2,i_wb))
      iy2 = min(n0_wb(2,i_wb), iyout2-ioff_wb(2,i_wb))
      iz1 = max(            1, izout1-ioff_wb(3,i_wb))
      iz2 = min(n0_wb(3,i_wb), izout2-ioff_wb(3,i_wb))

      do iz = iz1, iz2
      do iy = iy1, iy2
      do ix = ix1, ix2
        ibin = int(binscale * (flds(ix,iy,iz,iwrk) - binmin))
        ibin = max(1,min(maxprof, ibin))
        profile(ibin,2) = profile(ibin,2) + 1.0d0
      enddo
      enddo
      enddo

      return
      end

C ====================================================================72
      subroutine fill_zprof(profile, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real*8  profile(maxprof,maxcols), val
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      integer itof, iwrk, ix, iy, iz, izd,i_wb, ix1,ix2,iy1,iy2,iz1,iz2

      itof = itof_ops(nops)
      iwrk = iWrk_tof(itof)

      ix1 = max(            1, ixout1-ioff_wb(1,i_wb))
      ix2 = min(n0_wb(1,i_wb), ixout2-ioff_wb(1,i_wb))
      iy1 = max(            1, iyout1-ioff_wb(2,i_wb))
      iy2 = min(n0_wb(2,i_wb), iyout2-ioff_wb(2,i_wb))
      iz1 = max(            1, izout1-ioff_wb(3,i_wb))
      iz2 = min(n0_wb(3,i_wb), izout2-ioff_wb(3,i_wb))

      do iz = iz1, iz2
      izd = iz + ioff_wb(3,i_wb) - (izout1-1)   ! << subtrace z boundary depth here
      !  iy = 1, n0_wb(2,i_wb)
      !  ix = 1, n0_wb(1,i_wb)

      do iy = iy1, iy2
      do ix = ix1, ix2
        val = flds(ix,iy,iz,iwrk)
        profile(izd,2) = profile(izd,2) + val
        profile(izd,3) = profile(izd,3) + val**2
        profile(izd,7) = profile(izd,7) + 1.0d+00
        if(profile(izd,5) .gt. val) profile(izd,5) = val
        if(profile(izd,6) .lt. val) profile(izd,6) = val
      enddo
      enddo
      enddo

      if(idebug .gt. 1000) then
        write (6,*) 'iwrk = ', iwrk
        write (6,*) 'iz,izd,flds(1,1,iz,iwrk)'
        do iz = 1, n0_wb(3,i_wb)
        izd = iz + ioff_wb(3,i_wb)
          write (6,*) iz,izd,flds(1,1,iz,iwrk)
        enddo
        stop
      endif

      return
      end

C ====================================================================72
      subroutine gen_sect_bins(jdim)
      implicit none
      include 'postproc_info_e3d.h'
      integer itof, iarg, jdim, ii, isec_value

      itof = itof_ops(nops)
      do iarg=1,nArg_tof(itof)
       ii = isec_value(cArg_tof(itof,iarg),cSecFormat(jdim),cmsigns)
       if(iarg .eq. 1) then
         isect_min = ii
         isect_max = ii
       else
         isect_min = min(isect_min, ii)
         isect_max = max(isect_max, ii)
       endif
      enddo
      if(jdim .eq. 1) isect_min = 1

      do iarg=1,nArg_tof(itof)
        ii = isec_value(cArg_tof(itof,iarg),cSecFormat(jdim),cmsigns)
        isect_bin(iarg) = 2 + ii - isect_min   ! bin 1 is reserved for coofdinate
      enddo

      return
      end

C ====================================================================72
      subroutine fill_xsect(profile, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real*8  profile(maxprof,maxcols), val
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      real f0, f1, sval
      integer itof,iwrk0, ix,iy,iz, ixd,iyd, i_wb, nDim,iarg,jdim,ii
      integer isd, ismin, iy1, iy2
      real*8 area, afac, tot

!          ixout1 ixout1 ixout1 ixout1
      iy1 = max(1            , iyout1 - ioff_wb(2,i_wb))
      iy2 = min(n0_wb(2,i_wb), iyout2 - ioff_wb(2,i_wb))
      if(iy1 .gt. iy2) return


      jdim = 1  ! default is 1st sectional dimension (usually dimaeter)
      itof = itof_ops(nops)
      iwrk0 = iWrk_tof(itof) - 1
      nDim = nDim_tof(itof)
      call gen_sect_bins(1)   ! argument (jdim) is 1 here: y-axis is sectional size with charges summed

      ! Normalize area to produce average in [iyout1,iyout2] rather than total
      tot = 0.0d0
      area = dx * dz
      do iyd = 1, nyfull
        if(iyout1.le.iyd .and.iyd.le.iyout2) tot = tot + area
      enddo
      area = area / tot
c coord1_min

      do iarg=1,nArg_tof(itof)
        isd = isect_bin(iarg)
        do iz = 1, n0_wb(3,i_wb)
        !  iy = 1, n0_wb(2,i_wb)
        do iy = iy1, iy2
        iyd = iy + ioff_wb(2,i_wb)
        do ix = 1, n0_wb(1,i_wb)
          ixd = ix + ioff_wb(1,i_wb)
          val = flds(ix,iy,iz,iwrk0+iarg)
          profile(ixd,isd) = profile(ixd,isd) + area * val
        enddo
        enddo
        enddo
      enddo

      return
      end

C ====================================================================72
      subroutine fill_ysect(profile, i_wb, flds)
      implicit none
      include 'postproc_info_e3d.h'
      real*8  profile(maxprof,maxcols), val
      real flds(iwxlo:iwxhi, iwylo:iwyhi, iwzlo:iwzhi, nfv0_max)
      real f0, f1, sval
      integer itof,iwrk0, ix,iy,iz, ixd,iyd, i_wb, nDim,iarg,jdim,ii
      integer isd, ismin, ix1, ix2
      real*8 area(maxprof), afac, tot

      ix1 = max(1            , ixout1 - ioff_wb(1,i_wb))
      ix2 = min(n0_wb(1,i_wb), ixout2 - ioff_wb(1,i_wb))
      if(idebug .gt. 50) then
         write (6,902) i_wb,ixout1,ixout2,ioff_wb(1,i_wb)
902      format("i_wb,ixout[12],ioff_wb(1,i_wb):",5i5)
      endif
      if(ix1 .gt. ix2) return

      itof = itof_ops(nops)
      iwrk0 = iWrk_tof(itof) - 1
      nDim = nDim_tof(itof)
      if(cOutput(2:2) .eq. 's') jdim = 1    ! y-axis is sectional size with charges summed
      if(cOutput(2:2) .eq. 'c') jdim = 2    ! y-axis is charge with  sectional sizes summed
      call gen_sect_bins(jdim)
      if(idebug .ge. 21) write (6,901) jdim, isect_min, isect_max
901   format("in fill_ysect: jdim = ",i2,"    isect_[min,max] = ", 2i3)

      ! ixout1  ixout2

      tot = 0.0d0
      if(isRZ .eq. 1) then
        afac = 4.0d0 * atan(1.0d0) * dble(dx**2)
        do ixd = 1, nxfull
          area(ixd) = afac * ((dble(ixd))**2 - (dble(ixd-1))**2)
          if(ixout1.le.ixd .and.ixd.le.ixout2) tot = tot + area(ixd)
        enddo
      else
        do ixd = 1, nxfull
          area(ixd) = dx * dz
          if(ixout1.le.ixd .and.ixd.le.ixout2) tot = tot + area(ixd)
        enddo
      endif
      ! Normalize area to produce average in [ixout1,ixout2] rather than total
      do ixd = 1, nxfull
        area(ixd) = area(ixd) / tot
      enddo
c coord1_min

      do iarg=1,nArg_tof(itof)
        isd = isect_bin(iarg)
        do iz = 1, n0_wb(3,i_wb)
        do iy = 1, n0_wb(2,i_wb)
        iyd = iy + ioff_wb(2,i_wb)
        !  ix = 1, n0_wb(1,i_wb)
        do ix = ix1, ix2
          ixd = ix + ioff_wb(1,i_wb)
          val = flds(ix,iy,iz,iwrk0+iarg)
          profile(iyd,isd) = profile(iyd,isd) + area(ixd) * val
        enddo
        enddo
        enddo
      enddo

      return
      end



ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      ! debug
c      ii = 2 + isect_max - isect_min
c      write (6,*) "iyd,(profile(iyd,isd),isd=1,ii)"
c      do iy = 1, n0_wb(2,i_wb)
c        iyd = iy + ioff_wb(2,i_wb)
c        write (6,888) iyd,(profile(iyd,isd),isd=1,ii)
c888     format(i4,1p20e12.4)
c      enddo
c      write (6,*) "fill_ysect is being implemented"
c      if(nDim .ge. -1000) stop
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


C ====================================================================72
      subroutine write_strip(prof)
      implicit none
      include 'postproc_info_e3d.h'
      real*8  prof(maxprof,maxcols)
      real x0,x1,y0,y1,z0,z1
      integer i, j, n
      character*256 cPltFile, ss, s1, c1, c2, c3
      character*1 cc
      character*10 cmax

      if(idebug .gt. 10) then
        write (6,*) ' '
        write (6,*) 'Output File:    ', trim(cOutFile)
      endif
      open(unit=11,file=cOutFile, form='formatted')
      do i = 1, maxprof
        write (11,999) (prof(i,j), j=1,maxcols)
999     format(1p100e15.6)
      enddo
      close(11)
      call setblank(cPltFile)
      n = len(trim(cOutFile))
      cPltFile(1:n+4) = trim(cOutFile) // ".plt"
      write (6,*) "strip plot file: ", trim(cPltFile)
      open(unit=11,file=cPltFile,form="formatted")
1     format(a,a,a,a,a,a,a,a)

      ! xlabel
      call setblank(c1)
      if(iout_xstrip.eq.1) call addstr(c1,cxcoord)
      if(iout_ystrip.eq.1) call addstr(c1,cycoord)
      if(iout_strip.eq.1) call addstr(c1,"S")
      if(iflip1 .eq. 1) then
        cc = c1(1:1)
        call setblank(c1)
        c1(1:6) = cc // "max-" // cc
      endif
22    format(a,"=",f5.2)
      call setblank(c2)
      if(iout_xstrip.eq.1) write (c2,22) cycoord(1:1),y00
      if(iout_ystrip.eq.1) write (c2,22) cxcoord(1:1),x00
      if(iout_strip .eq.1) write (c2,1) "view = ", trim(cView)
      if(nzfull .gt. 1  .and.  iout_strip .ne. 1) then
        call setblank(c3)
        write (c3,22) czcoord(1:1),z00
        call addstr(c2,",")
        call addstr(c2,c3)
      endif
      write (11,1) "set xlabel '", trim(c1),"   (",trim(c2), ")'"
      write (11,1) "set ylabel '", trim(cFieldLabel), "'"
      if(ilogscale .eq. 1) write (11,1) "set logscale y"
      if(iflip1 .eq. 0) then
        write (11,1) "plot '", trim(cOutFile),
     1               "' using 1:2 notitle with linespoints"
      else
        call get_full_XYZrange(x0,x1,y0,y1,z0,z1)
        if(iout_xstrip.eq.1) write (cmax,23) x1
        if(iout_ystrip.eq.1) write (cmax,23) y1
23      format(f10.5)
        write (11,1) "plot '", trim(cOutFile), "' using ",
     1               "(",cmax,"-$1):2 notitle with linespoints"
      endif
      close(11)
      return
      end

C ====================================================================72
      subroutine write_xyslice(bout)
      implicit none
      include 'postproc_info_e3d.h'
      character*1 bout(nxfull,nyfull)
      integer nbytes, nwrote, sw_write

      nbytes = nxfull * nyfull
      nwrote = sw_write(cOutHost,cOutDir,cOutFile, i64output,
     1                  nbytes, bout)

      call w_ffdatadims(nxfull,nyfull,nzfull)

      write (6,*) 'Output file: ', trim(cOutFile)

      return
      end

C ====================================================================72
      subroutine w_ffdatadims(nfx,nfy,nfz)
      open(unit=15,file='ffdatadims',form='formatted')
      call wrt_line_to_15("export FFXDIM=", nfx)
      call wrt_line_to_15("export FFYDIM=", nfy)
      call wrt_line_to_15("export FFZDIM=", nfz)
      close(15)
      return
      end

C=====================================================================72
      subroutine wrt_line_to_15(ckey, nn)
      character*14 ckey
      character*8 c8num
      character*32 cline

      write (c8num,901) nn
      do i=1,8
        if(c8num(i:i) .eq. ' ') i1 = i+1
      enddo
      cline( 1:14   ) = ckey
      cline(15:23-i1) = c8num(i1:8)
      do i=24-i1, 32
        cline(i:i) = ' '
      enddo
      write (15,902) cline

901   format(i8)
902   format(a32)

      return
      end



C ====================================================================72
      subroutine gen_coord123(cxyz, i0, i1, i2)
      implicit none
      include 'postproc_info_e3d.h'
      integer i0, i1, i2
      character*1 cxyz(0:2)
      cxyz(0) = cxcoord(1:1)
      cxyz(1) = cycoord(1:1)
      cxyz(2) = czcoord(1:1)
      if(iout_xprof.eq.1) i0 = 0
      if(iout_yprof.eq.1) i0 = 1
      if(iout_zprof.eq.1) i0 = 2
      i1 = mod(i0+1,3)
      i2 = mod(i0+2,3)
      return
      end

c---------------------------------------------------------------------72
      subroutine write_prof(prof)
      implicit none
      include 'postproc_info_e3d.h'
      real*8  prof(maxprof,maxcols)
      integer i, j, i0, i1, i2
      real*8  av, ms, d2, total
      character*1 cxyz(0:2)
      character*256 cPltFile
      
      write (6,*) 'Output File:    ', trim(cOutFile)
      open(unit=11,file=cOutFile, form='formatted')
      if(iout_rprof .eq. 0) then
        call gen_coord123(cxyz, i0, i1, i2)
        write (11,998) "#             ",cxyz(i0),
     1 "        Average"," RMS Dispersion","         R.M.S.",
     1 "        Minimum","        Maximum","          Total"
      else
        write (11,998) "#        Radius",
     1 "        Average"," RMS Dispersion","         R.M.S.",
     1 "        Minimum","        Maximum","          Total"

      endif
998   format(a,a,a,a,a,a,a,a,a,a)
      ! if(idebug .gt. 10) 
      write (6,*) "write_prof: maxprof =", maxprof
      do i = 1, maxprof
      if(prof(i,7) .gt. 1.0d-30) then
        total = prof(i,2)
        av = prof(i,2) / prof(i,7)
        ms = prof(i,3) / prof(i,7)
        d2 = ms - av**2
        prof(i,2) = av
        prof(i,4) = sqrt(ms)
        prof(i,3) = 0.0
        if(d2 .gt. 0.0) prof(i,3) = sqrt(d2)
        write (11,999) (prof(i,j), j=1,6), total
999     format(1p7e15.6)
      endif
      enddo
      close(11)

      return
      end

c---------------------------------------------------------------------72
      real function value_sect(fsect0, jdim)
      implicit none
      include 'postproc_info_e3d.h'
      integer isect, jdim
      real fsect0, fsect, f0, f1, fval

c cSecMap(n), fSecF0(n), fSecF1(n)
      fsect = fsect0

      f0 = fSecF0(jdim)
      f1 = fSecF1(jdim)
      fval = f0 + f1*fsect                             ! default map is linear
      if(cSecMap(jdim)(1:3) .eq. "exp") fval = f0 * f1**fsect
      value_sect = fval

      return
      end
      
c---------------------------------------------------------------------72
      subroutine set_sect_title(str)
      implicit none
      include 'postproc_info_e3d.h'
      character*(*) str
      integer nlen
      nlen = len(str)
      call setblank(cSecTitle)
      cSecTitle(1:nlen) = str
      return
      end

c---------------------------------------------------------------------72
      subroutine write_Hsect(prof)
      implicit none
      include 'postproc_info_e3d.h'
      real*8  prof(maxprof,maxcols)
      integer i, j, i0, i1, i2,nsectmax, jdim, isect, ncbar, npersect
      real*8  av, ms, d2, total, v1min, v1max, v2min, v2max, dsect
      real value_sect, fsect, fs1, fs0
      character*256 cPltFile, crange, cmin, cmax
      real*8 hist
      allocatable hist(:,:)

      npersect = 1
      allocate (hist(isect_min:isect_max,2))

      if(coutput(2:2) .eq. 's') jdim = 1     ! section dimension
      if(coutput(2:2) .eq. 'c') jdim = 2     ! charge dimension

      if(idebug .ge. 21) write (6,901) jdim, isect_min, isect_max
901   format("jdim = ",i2,"    isect_[min,max] = ", 2i3)

      do j = isect_min, isect_max

        if(jdim .eq. 1) then
          fsect = float(j) - 0.5
          hist(j,1) = value_sect(fsect,jdim)
          fs0 = value_sect(float(j-1),jdim)
          fs1 = value_sect(float(j  ),jdim)
          dsect = fs1 - fs0
        else
          hist(j,1) = value_sect(float(j),jdim)
          dsect = 1.0
          if(isizelogscale.eq.1) dsect = log(fs1) - log(fs0)
        endif

        hist(j,2) = 0.0d0
        do i = iyout1, iyout2
            hist(j,2) = hist(j,2) +  prof(i,j+2-isect_min)
        enddo
        hist(j,2) = hist(j,2) / (dsect * dble(1+iyout2-iyout1))
      enddo

      write (6,*) 'Sectional Histogram Output File: ', trim(cOutFile)
      open(unit=11,file=cOutFile,form="formatted")
      do j = isect_min, isect_max
        write (11,999) (hist(j,i),i=1,2)
999     format(1p2e15.6)
      enddo

996   format(f5.2)
997   format(a," [",e15.6,":"e15.6,"]")
998   format(a,a,a,a,a,a,a,a,a,a)
      call setblank(cPltFile)
      write (cPltFile,998) trim(cOutFile), ".plt"
      open(unit=11,file=cPltFile, form='formatted')
      write (11,998) "set style fill solid 1.0 border -1"
      write (11,998) "set xlabel '",trim(cSecName(jdim)),"  [",
     1               trim(cSecUnit(jdim)),"]'"

      if(coutput(2:2) .eq. 's') then
      if(isizelogscale.eq.0) then
        write (11,998) "set ylabel 'dN/d(Dp)   [1/(nm-cc)]'"
      else
        write (11,998) "set ylabel 'dN/d(log Dp)   [1/cc]'"
        write (11,998) "set logscale x"
        fs0 = 0.95 * value_sect(float(isect_min),jdim)
        fs1 = 1.05 * value_sect(float(isect_max),jdim)
        write (11,997) "set xrange" , fs0, fs1
        if(fs1 .lt. 11.0) then
          write (11,998) "set xtics (0.125,0,25,0.5,1,2,4,8)"
        else if(fs1 .lt. 100.0) then
          write (11,998) "set xtics (0,25,0.5,1,2,4,8,16,32,64)"
        endif
      endif
      else
        write (11,998) "set ylabel 'Particle Density [1/cc]'"
        fs0 = value_sect(float(isect_min),jdim) - 0.7
        fs1 = value_sect(float(isect_max),jdim) + 0.7
        write (11,997) "set xrange" , fs0, fs1
      endif

      if(ilogscale .eq. 1) write (11,998) "set logscale y"
      write (11,998) "plot '", trim(cOutFIle), 
     1               "' using 1:2 notitle with boxes"

      close(11)

c  sof hist

      return
      end

c---------------------------------------------------------------------72
      subroutine write_sect(prof)
      implicit none
      include 'postproc_info_e3d.h'
      real*8  prof(maxprof,maxcols)
      integer i, j, i0, i1, i2,nsectmax, jdim, isect, ncbar, npersect
      real*8  av, ms, d2, total, v1min, v1max, v2min, v2max
      character*256 cPltFile, crange, cmin, cmax
      character*256 cFile, cxlabel, cylabel
      real sof, cbar(4,100)
      real value_sect, fsect, f0,f1, fs1,fs0, dsect
      integer nnnx, nnny
      allocatable sof(:,:)

      if(cOutput(1:1) .eq. 'H') then
        call write_Hsect(prof)
        return
      endif


      ! Generate the sof
      npersect = max(1, 512/isect_max)
      nsectmax = isect_max * npersect
      allocate (sof(0:maxprof,0:nsectmax))
      sof(0,0) = float(maxprof)
      v1max = profcoordmin + delbin1 * float(maxprof)
      do i = 1, maxprof
        if(trim(c2map) .eq. "max_c-c" .and. iout_ysect .eq. 1) then
          sof(i,0) = v1max - prof(i,1)
        else
          sof(i,0) = prof(i,1)
        endif
      enddo
      jdim = 1  ! section dimension
      do j = 1, nsectmax
        fsect = (float(j)-0.5) / float(npersect)
        sof(0,j) = value_sect(fsect,jdim)
      enddo

      ! map prof (filled as sum of concentrations) to distribution
      do isect = 1, isect_max
        f0 = (float(isect)-1.0)
        f1 = (float(isect)-0.0)
        fs0 = value_sect(f0,jdim)
        fs1 = value_sect(f1,jdim)
        dsect = fs1 - fs0
        if(isizelogscale.eq.1) dsect = log(fs1) - log(fs0)
        do i = 1, maxprof
          prof(i,isect+1) = prof(i,isect+1) / dsect
        enddo
      enddo

      ! fill sof with distribution values in prof 
      do j = 1, nsectmax
        isect = 1 + (j-1) / npersect
        do i = 1, maxprof
          sof(i,j) = prof(i,isect+1)
        enddo
      enddo

      ! write the sof
      write (6,*) 'Output File:    ', trim(cOutFile)
#ifndef ISGFORTRAN
      open(unit=11,file=cOutFile, form='binary')
#else
      open(unit=11,file=cOutFile, form='unformatted', access='stream')
#endif
      write (11) sof
      close(11)

      if(ncbarsteps.le.1) then
        call string_to_cbar(cSectCbar, cbar, ncbar)
      else
        call string_to_cbar(cSectCbar, cbar, ncbar)
        call cbar_substeps(ncbar,ncbarsteps,cbar)
      endif

       ! Generate plot limits
      v1min = profcoordmin
      v1max = profcoordmin + delbin1 * float(maxprof)
        if(trim(c2map) .eq. "max_c-c" .and. iout_ysect .eq. 1) then
          v1min = 0.0
          v1max = delbin1 * float(maxprof)
        endif
       v2min = value_sect(0.0,jdim)
       v2max = value_sect(float(isect_max),jdim)
      if(v2min*30.0 .gt. v2max .and. isizelogscale.eq.0) v2min = 0.0

      ! generate cb range if not already set
      if(hvmin .gt. hvmax) then
        hvmin = 1.0
        hvmax = prof(1,2)
        do isect = 1, isect_max
        do i = 1, maxprof
            hvmax = max(hvmax, prof(i,isect+1))
        enddo
        enddo
        hvmax = 10.0**(float(int(1001.0+alog(hvmax)/alog(10.0))-1000))
      endif

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ! Generate plot file

996   format(f5.2)
997   format(a," [",e15.6,":"e15.6,"]")
998   format(a,a,a,a,a,a,a,a,a,a)

      call setblank(cRange)
      call setblank(cmin)
      call setblank(cmax)
      if(iout_ysect .eq. 1) then
        write (cmin,996) xmin
        write (cmax,996) xmax
        if(ixout1 .eq. ixout2) then
          write (cRange,998) "at ", cxcoord(1:1), " = ", trim(cmin)
        else
          write (cRange,998) " averaged over  ", cxcoord(1:1), " in [",
     1                         trim(cmin), ", ", trim(cmax), "]"
        endif
      else
        write (cmin,996) ymin
        write (cmax,996) ymax
        if(iyout1 .eq. iyout2) then
          write (cRange,998) "at ", cycoord(1:1), " = ", trim(cmin)
        else
          write (cRange,998) " averaged over  ", cycoord(1:1), " in [",
     1                         trim(cmin), ", ", trim(cmax), "]"
        endif
      endif

      call setblank(cSecTitle)
      if(isizelogscale .eq. 0) then
        write (cSecTitle,998) "dN/d(Dp): Particle Size Distribution",
     1                        " [1/(cc nm)]  ", trim(cRange)
      else
        write (cSecTitle,998) "dN/d(Log Dp) [1/cc]  ", trim(cRange)
      endif
      call setblank(cPltFile)
      write (cPltFile,998) trim(cOutFile), ".plt"
      open(unit=11,file=cPltFile, form='formatted')
      write (11,998) "set title '", trim(cSecTitle), "'"
      write (11,998) "set pm3d map"
      write (11,997) "set xrange" , v1min, v1max
      write (11,997) "set yrange" , v2min, v2max
      write (11,997) "set cbrange", hvmin, hvmax
      write (11,998) "set pm3d explicit at b"
      call write_cbar(11, ncbar, ncbarsteps, cbar)
      write (11,998) "set clip two"
      call setblank(cxlabel)
      call setblank(cylabel)
      if(iout_ysect .eq. 1) then
        if(c2label(1:1) .eq. " ") then
          write (11,998) "set xlabel '",cycoord(1:1),"'"
          write (cxlabel,998) cycoord(1:1)
        else
          write (11,998) "set xlabel '",trim(c2label),"'"
          write (cxlabel,998) trim(c2label)
        endif
      else
        write (11,998) "set xlabel '",cxcoord(1:1),"'"
        write (cxlabel,998) cxcoord(1:1)
      endif

      write (11,998) "set ylabel '",trim(cSecName(jdim)),"  [",
     1               trim(cSecUnit(jdim)),"]'"
      write (cylabel,998) trim(cSecName(jdim)),"  [",
     1                    trim(cSecUnit(jdim)),"]"
      if(cSecMap(jdim)(1:3) .eq. "exp" .and. v2min .ne. 0.0) then
        write (11,998) "set logscale y"
        if(v2max .lt. 10.0) write (11,998)
     1 "set ytics (0.1,0,2,0.3,0.4,0.5,0.6,0.7,0.8,0.9",
     1             ",1,2,3,4,5,6,7,8,9,10)"
      endif
      if(hvmin .gt. 0.0) write (11,998) "set logscale zcb"
      write (11,998) "splot '", trim(cOutFile),
     1               "' binary notitle with pm3d"
      close(11)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      write (cFile,998) trim(cOutFile), ".pdt"
      call write_tec(cfile,cSecTitle,cxlabel,cylabel, "PSD",
     1               maxprof, nsectmax, sof)
      deallocate(sof)


      return
      end

c---------------------------------------------------------------------72
      subroutine write_corr(icorrelate)
      implicit none
      include 'postproc_info_e3d.h'
      integer icorrelate(nbin1,nbin2), i10max, i,j, i1,i2,n1,n2, nrep
      real ff, fmin, fmax, fval
      character*256 cc
      integer*8 ic
      integer m
      allocatable ff(:,:)

      nrep = 4
      n1 = nrep * nbin1
      n2 = nrep * nbin2
      allocate (ff(n1,n2))

      ic = 0
      fmin = 0.1
      fmax = float(icorrelate(1,1))
      do j = 1, nbin2
      do i = 1, nbin1
        ic = ic + icorrelate(i,j)
        fval = float(icorrelate(i,j))
        fmax = max(fmax, fval)
        do i1 = 1 + (i-1)*nrep, i*nrep
        do i2 = 1 + (j-1)*nrep, j*nrep
          ff(i1,i2) = fval
        enddo
        enddo
      enddo
      enddo
      i10max = int(log(fmax) / log(10.0) + 0.999999)
      fmax  = 10.0 ** i10max
      cbarString(1:20) = "w0grybap            "
      do i = i10max+3, 20
        cbarString(i:i) = ' '
      enddo

      call setblank(cc)
      m   = len(trim(cFld1) // "  (" // cOutFile(1:10) // ")")
      cc(1:m) = trim(cFld1) // "  (" // cOutFile(1:10) // ")"
      call sofplot(ff,n1,n2,1,bin10,bin11,bin20,bin21,cOutFile,
     1         trim(cc),trim(cFld2), fmin, fmax,1,n1,1,n2,
     1         ncbarsteps, cbarString, 1.0)

      deallocate (ff)

      write (6,*) "Output file : ", trim(cOutFile)
      write (6,*) "Cell count  = ", ic

      return
      end

c---------------------------------------------------------------------72
      subroutine write_view_from_sym33(sym33)
      implicit none
      include 'postproc_info_e3d.h'
      real*8 sym33(3,3), ye(3), ee(3,3), dfac
      integer ierr, ie, ic, i,j
      character*256 cViewOut
      real fov, cent(3),eye(3),up(3), fnear, ffar

      sym33(2,1) = sym33(1,2)
      sym33(3,1) = sym33(1,3)
      sym33(3,2) = sym33(2,3)

      dfac = 1.0d0 / dble(nxout*nyout*nzout)
      do j = 1, 3
      do i = 1, 3
        sym33(i,j) = sym33(i,j) * dfac
      enddo
c      write (6,998) (sym33(i,j), i=1,3)
998   format("sym33: ", 3f10.5)
      enddo

      call princip(sym33,ye,ee, ierr)
      if(ierr .lt. 0) then
        write (6,*) "Error in principle called from do_eigenv_s"
        stop
      endif

      do ie = 1, 3
        write (6,999) ye(ie), (ee(ic,ie),ic=1,3)
999     format("Variation: ", f10.5, "    dir: ",3f8.3)
      enddo

      fov = 35.308399
      call gen_view_par(ee,fov,eye,cent,up,fnear,ffar,
     1                   xmin,xmax,ymin,ymax,zmin,zmax)

      call setblank(cViewOut)
      cViewOut(1:6) = "a.view"
      if(idebug .gt. 20) write (6,*) "View file: ", trim(cViewOut)
      call write_view(cViewOut,fov,eye,cent,up,fnear,ffar,
     1                xmin,xmax,ymin,ymax,zmin,zmax)

      return
      end

c---------------------------------------------------------------------72
      subroutine gen_view_par(ee,fov,eye,cent,up,fnear,ffar,
     1                        xmin,xmax,ymin,ymax,zmin,zmax)
      implicit none
      real*8 ee(3,3)
      real   fov, eye(3), cent(3), up(3), fnear, ffar
      real   xmin,xmax,ymin,ymax,zmin,zmax
      real   distance, pi
      integer i

      pi = 4.0 * atan(1.0)
      distance = 0.5*(xmax-xmin) / tan(0.5 * fov*pi/180.0)
      cent(1) = 0.5 * (xmin + xmax)
      cent(2) = 0.5 * (ymin + ymax)
      cent(3) = 0.5 * (zmin + zmax)
      do i = 1, 3
        up(i)    = ee(i,2)
        eye(i) = cent(i) + distance * ee(i,3)
      enddo
      fnear = 0.1 * distance
      ffar  = 3.0 * distance

      return
      end
c---------------------------------------------------------------------72
      subroutine write_dist(prof)
      implicit none
      include 'postproc_info_e3d.h'
      real*8  prof(maxprof,maxcols)
      integer i, j, i0, i1, i2
      real*8  av, ms, d2, total
      character*1 cxyz(0:2)
      character*256 cPltFile
      real*8 totcells, facpdf
      real  value, pdf, cumul, time
      integer imax, imin, icounts
      integer(kind=8) i64sum

      ! Write out bin, counsts, PDF, & Cumulative Probability
c      totcells = dble(nxfull)*dble(nyfull)*dble(nzfull)
      totcells = 0.0d0
      do i = 1, maxprof
      totcells = totcells + dble(prof(i,2))
        if(prof(          i,2) .gt. 0) imax = i
        if(prof(1+maxprof-i,2) .gt. 0) imin = 1+maxprof-i
      enddo
      facpdf = dble(binscale) / totcells

      call get_time(time)

      open(unit=11, file=cOutFile, form='formatted')
#ifndef ISGFORTRAN
      write (11,999) "#time       ", time
999   format(a, f)
#else
      write (11) "#time       ", time
#endif
      write (11,998) "#col_1      ",trim(cFieldLabel)
      write (11,998) "#col_2      ","Counts"
      write (11,998) "#col_3      ","Probability Distribution Function"
      write (11,998) "#col_4      ","Cumulative Probability"
998   format(a, a)
      i64sum = 0
      do i = imin, imax
        value = (0.5 + float(i))/binscale + binmin
        icounts = int(prof(i,2))
        pdf = facpdf * prof(i,2)
        i64sum = i64sum + icounts
        cumul = dble(i64sum) / totcells
        write (11,997) value, icounts, pdf, cumul
997     format(f12.6, i15, 2f20.15)
      enddo
      close(11)

      return
      end

c---------------------------------------------------------------------72
      subroutine write_prof_plt()
      implicit none
      include 'postproc_info_e3d.h'
      integer i0, i1, i2
      real    v1min, v1max
      character*1 cxyz(0:2)
      character*256 cPltFile, cu

      v1min = profcoordmin
      v1max = profcoordmin + delbin1 * float(maxprof)
      
      call setblank(cPltFile)
998   format(a,a,a,a,a,a,a,a,a,a,a,a)
      write (cPltFile,998) trim(cOutFile), ".plt"
      open(unit=11,file=cPltFile,form="formatted")
997   format("set xrange [",f10.5,":"f10.5,"]")

      call setblank(cu)
      if(iflip1 .eq. 0) then
        cu(1:9) = "' using 1"
        write (11,997) v1min, v1max
      else
        write (cu,990) v1max
990     format("' using (", f10.5, "-$1)")
        write (11,997) 0, v1max-v1min
      endif
      if(iout_rprof .eq. 0) then
        call gen_coord123(cxyz, i0, i1, i2)
        if(iflip1 .eq. 0) then
          write (11,998) "set xlabel '",cxyz(i0),"'"
        else
          write (11,998) "set xlabel '",cxyz(i0),"max-",cxyz(i0),"'"
        endif
      else
        write (11,998) "set xlabel 'Radius'"
      endif
      write (11,998) "set ylabel '",trim(cFieldLabel),"'"

      if(cPart(1:4) .eq. "full") write (11,998) "plot '",trim(cOutFile),
     1 trim(cu), ":2:3 title 'RMS Dispersion' with yerrorlines, ",
     1 "'",trim(cOutFile),trim(cu), ":5 title 'Minimum' with lines, ",
     1 "'",trim(cOutFile),trim(cu), ":6 title 'Maximum' with lines"

      if(cPart(1:2) .eq. "av") write (11,998) "plot '",trim(cOutFile),
     1 trim(cu), ":2 title 'Average' with linespoints"

      if(cPart(1:3) .eq. "sum") write (11,998) "plot '",trim(cOutFile),
     1 trim(cu), ":7 title 'Sum' with linespoints"

      close(11)
      return
      end

C ====================================================================72
      subroutine setblank(str)
      character*(*) str
      integer n, i
      n = len(str)
      do i = 1, n
        str(i:i) = ' '
      enddo
      return
       end

C ====================================================================72
      subroutine addstr(str1, str2)
      character*(*) str1, str2
      integer n1, n2
      n1 = len(trim(str1))
      n2 = len(trim(str2))
      str1(n1+1:n1+n2) = trim(str2)
      return
       end

c=====================================================================72
      subroutine get_minmax(cFld, fmin, fmax)
      implicit none
      character*(*) cFld
      real fmin, fmax
      character*256 cFile, ckey
      integer n
      call setblank(cFile)
      n = len(trim(cFld))
      cFile(1:n+7) = trim(cFld) // ".minmax"
      open(unit=91,file=cFile,status="old",form="formatted",err=100)
      read (91,*,err=101) ckey, fmin, fmax
      close(91)
      return

100   write (6,*) "Can not open: ", trim(cFile)
      stop

101   write (6,*) "trouble reading minmax in: ", trim(cFile)
      stop

      end

c=====================================================================72
      logical function b_file_exists(cFile)
      implicit none
      character*(*) cFile
      open(unit=91,file=cFile,status="old",err=100)
      close(91)
      b_file_exists = .true.
      return

100   b_file_exists = .false.
      return

      end
c=====================================================================72
