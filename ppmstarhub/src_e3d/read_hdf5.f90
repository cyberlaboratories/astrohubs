

      subroutine read_hdf5(nx,ny,nz,cfile,cfield,fld)
      use hdf5
      implicit none
      INTEGER nx,ny,nz
      character(255) cfile,cfield
      REAL fld(nx,ny,nz)
      INTEGER hdferr
      INTEGER(HID_T) file, group,dset
      INTEGER(HSIZE_T) dims(3)

      dims(1) = nx
      dims(2) = ny
      dims(3) = nz
      CALL h5open_f(hdferr)
      CALL h5fopen_f(trim(adjustl(cfile)),H5F_ACC_RDONLY_F,file,hdferr)
      CALL h5dopen_f(file,trim(adjustl(cfield)),group,hdferr)
      CALL h5dread_f(group,H5T_NATIVE_REAL,fld,dims,hdferr)
      CALL h5fclose_f(file,hdferr)
      CALL h5close_f(hdferr)

      return
      end

!=======================================================================

      subroutine read_hdf5_sub(nx,ny,nz,iox,ioy,ioz,cfile,cfld,fld)
      use hdf5
      implicit none
      INTEGER nx,ny,nz,iox,ioy,ioz, j
      character(255) cfile,cfld
      REAL fld(nx,ny,nz)
      INTEGER(HSIZE_T) idims(3), ioff(3), ioff_out(3)
      INTEGER(HID_T) file_id, dset_id, memspace,dataspace, error   ! HDF5 handles & pointers
      INTEGER hdferr
      idims(1) = nx
      idims(2) = ny
      idims(3) = nz
      ioff(1) = iox
      ioff(2) = ioy
      ioff(3) = ioz
      ioff_out(1) = 0
      ioff_out(2) = 0
      ioff_out(3) = 0

     ! Open file &  dataset.
      CALL h5open_f(hdferr)
     CALL h5fopen_f (trim(cfile), H5F_ACC_RDONLY_F, file_id, error)
     CALL h5dopen_f(file_id, trim(cfld), dset_id, error)
     CALL h5dget_space_f(dset_id, dataspace, error)             ! Get dataset id
     CALL h5sselect_hyperslab_f(dataspace, H5S_SELECT_SET_F, &
                                ioff, idims, error)            ! Select hyperslab
     CALL h5screate_simple_f(3, idims, memspace, error)        ! Create memory dataspace.
     CALL h5sselect_hyperslab_f(memspace, H5S_SELECT_SET_F, &
                                ioff_out, idims, error)        ! Select hyperslab in memory.
     CALL H5dread_f(dset_id, H5T_NATIVE_REAL, fld, idims, &
                     error, memspace, dataspace)               ! Read hyperslab

     CALL h5sclose_f(dataspace, error)
     CALL h5sclose_f(memspace, error)
     CALL h5dclose_f(dset_id, error)
     CALL h5fclose_f(file_id, error)
      CALL h5close_f(hdferr)


      return
      end

