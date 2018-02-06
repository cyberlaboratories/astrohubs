c
c
      subroutine avminmax(nx0,ny0,nz0,nx1,ny1,nz1,vvv0,vvv1,vmin,vmax)
      character(1) vvv0(nx0*ny0*nz0), vvv1(nx1*ny1*nz1)
      character(1) vmin(nx1*ny1*nz1), vmax(nx1*ny1*nz1)
      integer vv0(4096),vv1(4096),vv2(4096),vv3(4096)
      integer vvmin(4096), vvmax(4096)
      integer ivav(4096), ivmin(4096), ivmax(4096)

      nxy0 = nx0 * ny0;
      do k1=0, nz1-1
      k0 = k1 * 2
      kadd = nxy0;
      if(k0+1 .ge. nz0) kadd = 0;

      do j1=0, ny1-1
      j0 = j1 * 2
      jadd = nx0
      if(j0+1 .ge. ny0) jadd = 0

      jk1off = nx1 * (j1 + ny1 * k1)
      jk0off = nx0 * (j0 + ny0 * k0)

      jk0off0 = jk0off
      jk0off1 = jk0off + jadd
      jk0off2 = jk0off + kadd
      jk0off3 = jk0off + jadd + kadd
      
!DEC$ VECTOR ALWAYS
c!DEC$ VECTOR ALIGNED
ccdir$ swp
      do i=1,nx0
        vv0(i) = ichar(vvv0(jk0off0+i))
        vv1(i) = ichar(vvv0(jk0off1+i))
        vv2(i) = ichar(vvv0(jk0off2+i))
        vv3(i) = ichar(vvv0(jk0off3+i))
        vvmin(i) = min(vv0(i),vv1(i),vv2(i),vv3(i))
        vvmax(i) = max(vv0(i),vv1(i),vv2(i),vv3(i))
      enddo

!DEC$ VECTOR ALWAYS
c!DEC$ VECTOR ALIGNED
ccdir$ swp
      do i1=1,nx1
        i0   = 1 + (i1-1) * 2
        i0p1 = i0 + 1

        ivav(i1) = (vv0(i0) + vv0(i0p1) + vv1(i0) + vv1(i0p1) +
     1              vv2(i0) + vv2(i0p1) + vv3(i0) + vv3(i0p1)   )/8
        ivmax(i1) = max(vvmax(i0), vvmax(i0p1))
        ivmin(i1) = min(vvmin(i0), vvmin(i0p1))
      enddo

!DEC$ VECTOR ALWAYS
c!DEC$ VECTOR ALIGNED
ccdir$ swp
      do i1=1,nx1
        vvv1(jk1off+i1) = char(ivav(i1))
        vmax(jk1off+i1) = char(ivmax(i1))
        vmin(jk1off+i1) = char(ivmin(i1))
      enddo


      enddo
      enddo

c         v0 = vv0[i0  ];
c         v1 = vv0[i0p1];
c         v2 = vv1[i0  ];
c         v3 = vv1[i0p1];
c         v4 = vv2[i0  ];
c         v5 = vv2[i0p1];
c         v6 = vv3[i0  ];
c         v7 = vv3[i0p1];

      return
      end

c---------------------------------------------------------------------72
ccccccccccccccc
c   ixoff = itx * nx;
c   iyoff = ity * ny;
c   for(iz=0; iz<curdim[2]; iz++) {
c   for(iy=0; iy<ny       ; iy++) {
c       idstoff = ixoff + iNx * (iy+iyoff + iNy * iz);
c       isrcoff =  nx * (iy + ny * iz);
c       //memcpy(&ucBlockTreeData[idstoff], &ucSubSlice[isrcoff], nx*4);
c       for(ix=0; ix<nx       ; ix++) {
c         ucBlockTreeData[ix + idstoff] = ucSubSlice[ix + isrcoff];
c       }
c   }}
ccccccccccccccc
      subroutine copyinto(nx,ny,nz,itx,ity,itz,nfx,nfy,nfz,src,dst)
      byte src(nx*ny*nz), dst(nfx*nfy*nfz)

      ixoff = itx * nx
      iyoff = ity * ny

      do iz = 0, nfz-1
      do iy = 0, ny-1
        idstoff = ixoff + nfx * (iy+iyoff + nfy * iz)
        isrcoff = nx * (iy + ny * iz)

!DEC$ VECTOR ALWAYS
c!DEC$ VECTOR ALIGNED
ccdir$ swp
        do ix = 1, nx
          dst(ix+idstoff) = src(ix+isrcoff)
        enddo

      enddo
      enddo

      return
      end
