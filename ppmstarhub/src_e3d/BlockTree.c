//#include "stdafx.h"
#include "StdAfx.h"

#ifdef WIN32
#include <windows.h>
#include <sys/timeb.h>
#endif
//#include <io.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <sys/time.h>
unsigned long GetTickCount()
{
  struct timeval times;
  double timeNow;
  static double time0=0.0;
  gettimeofday( &times, NULL );
  if(time0 == 0.0) { time0 = times.tv_sec; }
  timeNow = (times.tv_sec - time0) + times.tv_usec*1.0e-06;
  //timeNow = (times.tv_sec -1000000000) + times.tv_usec*1.0e-06;
  return (unsigned long)(1000.0 * timeNow);
}


#define MAX_NODES 8192
#define MAX_CHILDREN 64
#define NBLACK 4*1024*1024
#define FILE_BLOCK 64*1024

static unsigned char *ucSlab;
static unsigned char *fieldav[65];
static unsigned char *fieldmin[65];
static unsigned char *fieldmax[65];

void Block_Tree(char *pfile, int fulldim[3], int curdim[3], int curoff[3], int *pblock,
                float fLimits[6], int imap, unsigned char *field)
{
    int nx, ny, nz, nx1, ny1, nz1, iBlock;
    static int midpoint[MAX_NODES][3];
    static __int64 byteoff, tile_byteoff[MAX_NODES];
    static int tile_dim[MAX_NODES][3], tile_offset[MAX_NODES][3], tile_rep[MAX_NODES];
    static int tile_parent[MAX_NODES], tile_error[MAX_NODES];
    static int tile_children[MAX_NODES][MAX_CHILDREN], tile_N_children[MAX_NODES];
    static int irep, ireplast, nlines, iNnodes=0;
    int i, j, k, ih, jh, kh, node, node_root, ichild, more_nodes, maxerrtmp;
    int i0, j0, k0, i1, j1, k1, i0p1;
    __int64 nblocks, leftover;
    static int fd;
    static int iprint = 0;
    int irl, jrl, krl, irh, jrh, krh;
    int itl, jtl, ktl, ith, jth, kth;
    int kil, kih;
    int irn, jrn, krn, irl0, jrl0, krl0, irh0, jrh0, krh0, iOff, nChar;
    int ks, js, ijks, nslab;
    int itn, jtn, ktn, iErr;
    int vmin[1024], vmax[1024], iadd, jadd, kadd, ijk0, ijk1, ijk;
        register int v0, v1, v2, v3, v4, v5, v6, v7;
    register int vmax0,vmax1,vmax2,vmax3,vmax4,vmax5;
        int vv0[2048], vv1[2048], vv2[2048], vv3[2048];
    int nx0, ny0, nz0, nxy0, irep1, irep0, jkoff;
    int ijkd, itlu, ithu, jtlu, jthu, ktlu, kthu, is, kilu, kihu;
    static unsigned char fSlab[NBLACK],  fSlabMin[NBLACK],  fSlabMax[NBLACK];
    static int iDoRep[65];
    int iIntersection, nxyzmax, iBndry;
    __int64 i64Offset;
    char cLine[256];
    long dwSize;
    float dlog10, dlog10i, dsum, dmap[256];
    void TreeScanOffsets(int node, __int64 *pbyteoff, __int64 tile_byteoff[MAX_NODES],
       int irep, int tile_rep[MAX_NODES], int tile_dim[MAX_NODES][3],
       int tile_N_children[MAX_NODES], int tile_children[MAX_NODES][MAX_CHILDREN]);
    void avminmax_(int *nx0, int *ny0, int *nz0, int *nx1,int *ny1, int *nz1,
                   unsigned char *fieldav0, unsigned char *fieldav1,
                   unsigned char *fieldmin, unsigned char *fieldmax);
    unsigned long ulAvStart, ulRefStart, ulStart, ulwrite0Start;
    static float fAvTime, fRefTime, fTotalTime, fwrite0Time;
        int jk1off, jk0off, jk0off0, jk0off1, jk0off2, jk0off3, ijk00, ijk01, ijk02, ijk03;
        unsigned char *fldin, *fldav, *fldmn, *fldmx;
        int itlu0, ithu0;
        float gmin, gmax, gminav, gmaxav;

    ulStart = GetTickCount();

    // Dimensions of full domain, and max there of
    nx = fulldim[0];
    ny = fulldim[1];
    nz = fulldim[2];
    nxyzmax = nx;
    if(ny > nxyzmax) { nxyzmax = ny; }
    if(nz > nxyzmax) { nxyzmax = nz; }
    iBlock = *pblock;

    // Generate tile decomposition layout on 1st pass
    if(iNnodes == 0) {

        fRefTime    = 0.0;
        fAvTime     = 0.0;
        fwrite0Time = 0.0;
        fTotalTime  = 0.0;

        nlines     = 2 * (1+nx/iBlock)*(1+ny/iBlock)*(1+nz/iBlock);         // Max # of lines for hv info
        byteoff = (__int64)(FILE_BLOCK * (1 + (80 * nlines)/(FILE_BLOCK))); // offset to start of byte info

        for(irep=0; irep<65; irep++) { iDoRep[irep] = 0; }
        node = 0;
        irep = 1;
        more_nodes = 2;
        while(more_nodes > 1) {
//            printf("irep = %d\n", irep);
            iDoRep[irep] = 1;

            iBlock = *pblock;
//            if(nxyzmax > irep*(iBlock+2)  ||  nxyzmax <= irep*iBlock)
            if(nxyzmax > irep*(iBlock+2)) {
                iBndry = 1;        // Here, tiles need overlapping boundaries or it fits in one tile even w/ bndrys
            } else {
                iBndry = 0;        // Here, tile spans entire domain, no overlapping boundaries needed
                iBlock += 2;
            }

//            printf("iBlock=%d  nx=%d  ny=%d  nz=%d\n", iBlock, nx, ny, nz);

            more_nodes = 0;
            for(k=0; k<nz; k += iBlock*irep) { kh = k+iBlock*irep;  if(kh > nz) { kh = nz; }
            for(j=0; j<ny; j += iBlock*irep) { jh = j+iBlock*irep;  if(jh > ny) { jh = ny; }
            for(i=0; i<nx; i += iBlock*irep) { ih = i+iBlock*irep;  if(ih > nx) { ih = nx; }
                // midpoint of node will be used in determining parent node
                midpoint[node][0] = (i+ih) / 2;
                midpoint[node][1] = (j+jh) / 2;
                midpoint[node][2] = (k+kh) / 2;

                // xyz size of tile including 2 boundary zones CHANGE
                tile_dim[node][0] = 2*iBndry + (ih - i + irep - 1) / irep;
                tile_dim[node][1] = 2*iBndry + (jh - j + irep - 1) / irep;
                tile_dim[node][2] = 2*iBndry + (kh - k + irep - 1) / irep;  // if any of these are 2, will skip node!

                // offset (in bytes) from beginning of file to start of tile data
                tile_byteoff[node] = byteoff;
                byteoff += (__int64)(FILE_BLOCK*(1+
                    (tile_dim[node][0]*tile_dim[node][1]*tile_dim[node][2] - 1)
                    /(FILE_BLOCK)));

                // replication factor of tile
                tile_rep[node] = irep;

                // xyz offset of tile (including boundary indentation) in units of its own zones
                tile_offset[node][0] = i/irep - iBndry;
                tile_offset[node][1] = j/irep - iBndry;
                tile_offset[node][2] = k/irep - iBndry;

                // Default parent node is -1 (no parent).
                // Set parents of all child nodes based on midpoint of child being inside current node
                // As parents are found for children, tell the parents about it as well.
                tile_parent[node] = -1;
                tile_N_children[node] = 0;
                for(ichild=0; ichild<node; ichild++) {
                    if(2*tile_rep[ichild] == irep) {
                        if(i < midpoint[ichild][0] && midpoint[ichild][0] < ih &&
                           j < midpoint[ichild][1] && midpoint[ichild][1] < jh &&
                           k < midpoint[ichild][2] && midpoint[ichild][2] < kh   )
                        {
                            tile_parent[ichild] = node;
                            tile_children[node][tile_N_children[node]] = ichild;
                            tile_N_children[node]++;
                        }
                    }
                }

                // Initialize max error in using this tile
                tile_error[node]  =  0;

                // By construction, root node is always the last one in the list
                node_root = node;

                node++;
                more_nodes++;
            }}}
//            printf("node=%d   more_nodes=%d\n", node, more_nodes);
            irep *= 2;
        }
        iNnodes = node;
        node_root = node-1;            // By construction, root node is always the last one in the list

        // Now reset byte offsets into file to be in the order of a tree scan for each level
        // Will leave the root of the tree as the last brick, so that if you see that, then you have it all.
        nlines     = iNnodes + 100;                                           // # of lines for hv info + other header info
//        printf("nlines = %d\n", nlines);
        byteoff = (__int64)(FILE_BLOCK * (1 + (80 * nlines) / (FILE_BLOCK)));   // offset to start of byte info

        // For each replication factor used, do a  tree scan to
        //   1) find nodes at this replication factor, and
        //   2) assign byte offsets into the hv file in tree scan order
        for(irep=1; irep<65; irep++) { if(iDoRep[irep] > 0) {
            TreeScanOffsets(node_root, &byteoff, tile_byteoff,
                            irep, tile_rep, tile_dim,
                            tile_N_children, tile_children    );
        }}

        // Allocate space for reduced fields
        fieldav[1] = field;
        irep = 2;
//        printf("nx = %d   ny = %d   nz = %d\n", nx, ny, nz);
        while(irep < 65) {
//            printf("irep=%d   iDoRep[irep]=%d\n", irep, iDoRep[irep]);
            if(iDoRep[irep] > 0) {
                nx1 = (curdim[0] + irep - 1) / irep;
                ny1 = (curdim[1] + irep - 1) / irep;
                nz1 = (curdim[2] + irep - 1) / irep;
                dwSize = (long)nx1 * (long)ny1 * (long)nz1;
#ifndef LINUX
                fieldav[irep]  = (unsigned char *)VirtualAlloc(NULL, dwSize, MEM_COMMIT, PAGE_READWRITE);
                fieldmax[irep] = (unsigned char *)VirtualAlloc(NULL, dwSize, MEM_COMMIT, PAGE_READWRITE);
                fieldmin[irep] = (unsigned char *)VirtualAlloc(NULL, dwSize, MEM_COMMIT, PAGE_READWRITE);
#else
                fieldav[irep]  = (unsigned char *)malloc(dwSize);
                fieldmax[irep] = (unsigned char *)malloc(dwSize);
                fieldmin[irep] = (unsigned char *)malloc(dwSize);
#endif

// DEBUG
//printf("field[av,max,min]: irep = %2d, allocated %d bytes each\n", irep, dwSize);
if(fieldav[irep]  == NULL) { printf("fieldav  NULL: irep=%d  dwSize=%d\n", irep,dwSize); }
if(fieldmin[irep] == NULL) { printf("fieldmin NULL: irep=%d  dwSize=%d\n", irep,dwSize); }
if(fieldmax[irep] == NULL) { printf("fieldmax NULL: irep=%d  dwSize=%d\n", irep,dwSize); }
            }
            irep *= 2;
        }


        // Write out file filled with 0s
        ulwrite0Start = GetTickCount();  // TEST MOD
//printf("Start write of 0s\n");
#ifndef LINUX
        ucSlab = (unsigned char *)VirtualAlloc(NULL, NBLACK, MEM_COMMIT, PAGE_READWRITE);
#else
        ucSlab = (unsigned char *)malloc(NBLACK);
#endif
        // TEST MOD
        for(i=0; i<NBLACK; i++) { ucSlab[i] = (unsigned char)0; }
        // TEST MOD

//printf("byteoff=%lld\n", byteoff);

        leftover = byteoff;
        nblocks = 0;
        while(leftover >= (__int64)NBLACK) { nblocks++; leftover -= (__int64)NBLACK; }
//printf("NBLACK=%d   nblocks=%lld   byteoff=%lld   leftover=%lld\n", NBLACK, nblocks, byteoff, leftover);

        if(nblocks > 100000) { exit(0); }
#ifndef LINUX
        fd = _open(pfile, _O_WRONLY | _O_CREAT | _O_BINARY , _S_IREAD | _S_IWRITE );
#else
        fd = open(pfile,O_RDWR|O_CREAT, 0644 );
//        fd = open(pfile,O_CREAT|O_WRONLY);
#endif
        // TEST MOD
        for(i=0; i<nblocks; i++) { write(fd, ucSlab, NBLACK); }
        if(leftover > 0)         { write(fd, ucSlab, (unsigned int)leftover); }
        // TEST MOD
//        printf("Finish write of 0s\n");
        fwrite0Time += (float)(GetTickCount() - ulwrite0Start) / (float)1000.0;
    }

    if(imap == 1) {
        // i = 50.0 * (1.0 +  log(v)/log(10.0))
        // v  = exp(log(10.0) * (0.02 * (float)i - 1.0));
        dlog10  = (float)log(10.0);
        dlog10i = (float)1.0 / dlog10;
        for(i=0; i<256; i++) {
            dmap[i] = (float)exp(dlog10 * ((float)0.02 * (float)i - (float)1.0));
        }
    }

    ulAvStart = GetTickCount();  // TEST MOD

    // For every subregion, start by generating all of the reduced resolution fields
    irep1 = 2;
    while(irep1 < 65) {
        if(iDoRep[irep1] > 0) {
//            printf("start field average for irep = %2d\n", irep1);
            irep0 = irep1 / 2;
            nx1 = (curdim[0] + irep1 - 1) / irep1;
            ny1 = (curdim[1] + irep1 - 1) / irep1;
            nz1 = (curdim[2] + irep1 - 1) / irep1;
            nx0 = (curdim[0] + irep0 - 1) / irep0;
            ny0 = (curdim[1] + irep0 - 1) / irep0;
            nz0 = (curdim[2] + irep0 - 1) / irep0;

            fldin = fieldav[irep0];
            avminmax_(&nx0,&ny0,&nz0,&nx1,&ny1,&nz1,
              fieldav[irep0],fieldav[irep1],fieldmin[irep1],fieldmax[irep1]);

/*
for(k=0; k<nz1; k++) { for(j=0; j<ny1; j++) { for(i=0; i<nx1; i++) {
  ijk = i + nx1*(j + ny1*k);
  if(fieldav[irep1][ijk] != fieldmin[irep1][ijk]  ||
     fieldav[irep1][ijk] != fieldmax[irep1][ijk]    ) {
  printf("ixyz=%d %d %d   fieldav  = %d\n", i,j,k, (int)fieldav[irep1][ijk]);
  printf("ixyz=%d %d %d   fieldmin = %d\n", i,j,k, (int)fieldmin[irep1][ijk]);
  printf("ixyz=%d %d %d   fieldmax = %d\n", i,j,k, (int)fieldmax[irep1][ijk]);
  exit(0);
  }
}}}
*/

//            printf("field averaged for irep = %2d\n", irep1);
        }
        irep1 *= 2;
    }

    fAvTime += (float)(GetTickCount() - ulAvStart) / (float)1000.0;
    ulRefStart = GetTickCount();

    // For every subregion currently in memory (pointed to by field)
    // for every tile in the Block tree (0 <= node < iNnodes)
    //    find the intersection
    //    if not empty: 1) copy intersecion into subregion of bob-tile and write to disk file
    //                  2) update max error tile
    //
    // Assumption: 1) Calling routine generate a sequence of regions which cover the entire domain.
    //             2) Region limits in each direction are multiples of the max.tile_rep.
    //             3) Region limits in XY span each tile that is updated

    // Used for generating offset in field array from ijk coords
    irn = curdim[0];
    jrn = curdim[1];
    krn = curdim[2];

    for(node=0; node<iNnodes; node++) {

        irep = tile_rep[node];

        fldav = fieldav[irep];
        fldmn = fieldmin[irep];
        fldmx = fieldmax[irep];

        // generate: Low and Hight inclusive limits on xyz (=ijk) ranges for "field" Region
        irl = curoff[0];   irh = curoff[0] + curdim[0] - 1;
        jrl = curoff[1];   jrh = curoff[1] + curdim[1] - 1;
        krl = curoff[2];   krh = curoff[2] + curdim[2] - 1;

        // save actual IJK limits of array field
        irl0 = irl;                        irh0 = irh;
        jrl0 = jrl;                        jrh0 = jrh;
        krl0 = krl;                        krh0 = krh;

        // If Region is at domain edge, then -conceptually- it includes a boundary of depth irep.
        // Values in this boundary region are gotten from the closest interior points of "field"
        if(irl == 0) { irl -= irep; }      if(irh == nx-1) { irh += 2 * irep; }
        if(jrl == 0) { jrl -= irep; }      if(jrh == ny-1) { jrh += 2 * irep; }
        if(krl == 0) { krl -= irep; }      if(krh == nz-1) { krh += 2 * irep; }

        // dimensions of tile
        itn = tile_dim[node][0];
        jtn = tile_dim[node][1];
        ktn = tile_dim[node][2];

        // generate: Low and Hight inclusive limits on xyz (=IJK) ranges for block tree Tile
        itl = irep*tile_offset[node][0];   ith = irep*(tile_offset[node][0]+tile_dim[node][0])-1;
        jtl = irep*tile_offset[node][1];   jth = irep*(tile_offset[node][1]+tile_dim[node][1])-1;
        ktl = irep*tile_offset[node][2];   kth = irep*(tile_offset[node][2]+tile_dim[node][2])-1;

        // generate: Low and Hight inclusive limits on z (=K) range for Intersection of Tile and Region
        if(krl > ktl) { kil = krl; } else { kil = ktl; }
        if(krh < kth) { kih = krh; } else { kih = kth; }

        // if the Low limit of the intersection > the High limit then there is no intersection
        iIntersection = 1;
        if(kil > kih) { iIntersection = 0; }

        // if any part of the tile in XY (=IJ) directions is not contained in the array field Region
        // then dont do tile.  Assume this block-tree tile (=node) is covered in some other call to this subroutine.
        if(itl < irl  ||  ith > irh) { iIntersection = 0; }
        if(jtl < jrl  ||  jth > jrh) { iIntersection = 0; }

        if(iIntersection == 1) {

            // Generate IJK Low and High coords in tile Units
            itlu = tile_offset[node][0];   ithu = tile_offset[node][0] + tile_dim[node][0];
            jtlu = tile_offset[node][1];   jthu = tile_offset[node][1] + tile_dim[node][1];
            ktlu = tile_offset[node][2];   kthu = tile_offset[node][2] + tile_dim[node][2];

            // This if-block is where all the copy to tiles and output is done.

            // 1) Array fField is copied, or averaged if irep>1, as floats into array fSlab
            //    max and min fField values (for each zone) are gathered into fSlabMax, fSlabMin arrays
            //    The three fSlab arrays are always in the tile space.
            kilu =  kil    / irep;
            kihu = (kih+1) / irep;
            nslab = itn * jtn * (kihu - kilu);        // size of slab

            nx0 = (curdim[0] + irep - 1) / irep;
            ny0 = (curdim[1] + irep - 1) / irep;
            nz0 = (curdim[2] + irep - 1) / irep;

            ijkd = 0;
            if(irep == 1) {
                for(k=kilu; k<kihu; k++) {
                    ks = k - curoff[2] / irep;
                    if(ks <     0) { ks =     0; }
                    if(ks > nz0-1) { ks = nz0-1; }
//                    ks = ks - curoff[2] / irep;
                    for(j=jtlu; j<jthu; j++) {
                        js = j;
                        if(js <     0) { js =     0; }
                        if(js > ny0-1) { js = ny0-1; }
                        jkoff = nx0*(js + ny0*ks);

                        itlu0 = itlu;   if(itlu <   0) { itlu0 =   0; }
                        ithu0 = ithu;   if(ithu > nx0) { ithu0 = nx0; }

                        if(itlu < 0) {
                            ucSlab[ijkd] = fldav[itlu0+jkoff];
                            ijkd++;
                        }

                        memcpy(&ucSlab[ijkd], &fldav[itlu0+jkoff], ithu0-itlu0);
                        ijkd += ithu0-itlu0;

                        if(nx0 < ithu) {
                            ucSlab[ijkd] = fldav[ithu0+jkoff];
                            ijkd++;
                        }

/**********************************************************************************
TEST MOD
                        for(i=itlu; i<ithu; i++) {
                            is = i;
                            if(is <     0) { is =     0; }
                            if(is > nx0-1) { is = nx0-1; }
                            ijks = is + jkoff;
                            ucSlab[ijkd] = fieldav[irep][ijks];
                            ijkd++;
                            ijks++;
                        }
TEST MOD
**********************************************************************************/

                    }
                }
            } else {

//if(irep > 1) printf("After 2a: kilu=%d  kihu=%d  jtlu=%d  jthu=%d  curoff[2]=%d\n", kilu,kihu,jtlu,jthu,curoff[2]);

                for(k=kilu; k<kihu; k++) {
                    ks = k - curoff[2] / irep;
                    if(ks <     0) { ks =     0; }
                    if(ks > nz0-1) { ks = nz0-1; }
//                    ks = ks - curoff[2] / irep;
                    for(j=jtlu; j<jthu; j++) {
                        js = j;
                        if(js <     0) { js =     0; }
                        if(js > ny0-1) { js = ny0-1; }
                        jkoff = nx0*(js + ny0*ks);

                        itlu0 = itlu;   if(itlu <   0) { itlu0 =   0; }
                        ithu0 = ithu;   if(ithu > nx0) { ithu0 = nx0; }

                        if(itlu < 0) {
                            ucSlab[ijkd]   = fldav[itlu0+jkoff];
                            fSlabMin[ijkd] = fldmn[itlu0+jkoff];
                            fSlabMax[ijkd] = fldmx[itlu0+jkoff];

//                            and she was (before 2008 Sept. 30th)
//                            fSlabMin[ijkd] = fldmn[ijks];
//                            fSlabMax[ijkd] = fldmx[ijks];
                            ijkd++;
                        }
                        memcpy(&ucSlab[ijkd],   &fldav[itlu0+jkoff], ithu0-itlu0);
                        memcpy(&fSlabMin[ijkd], &fldmn[itlu0+jkoff], ithu0-itlu0);
                        memcpy(&fSlabMax[ijkd], &fldmx[itlu0+jkoff], ithu0-itlu0);
                        ijkd += ithu0-itlu0;

                        if(nx0 < ithu) {
                            ucSlab[ijkd]   = fldav[ithu0+jkoff];
                            fSlabMin[ijkd] = fldmn[ithu0+jkoff];
                            fSlabMax[ijkd] = fldmx[ithu0+jkoff];

//                            and she was (before 2008 Sept. 30th)
//                            fSlabMin[ijkd] = fldmn[ijks];
//                            fSlabMax[ijkd] = fldmx[ijks];

                            ijkd++;
                        }
                    }
                }
            }

            // 3) Max error for node is updated based on min and max values under each voxel
            //    Only do this for tiles of reduced resolution, otherwise error is always 0
            if(irep > 1) {
// if(irep > 1) printf("After 3a\n");
                maxerrtmp = tile_error[node];
// if(irep > 1) printf("After 3b\n");
                for(i=0; i<nslab; i++) {
                    iErr = (int)ucSlab[i] - (int)fSlabMin[i];
                    if(maxerrtmp < iErr) { maxerrtmp = iErr; }

                    iErr = (int)fSlabMax[i] - (int)ucSlab[i];
                    if(maxerrtmp < iErr) { maxerrtmp = iErr; }
                }
// if(irep > 1) printf("After 3c\n");
                tile_error[node] = maxerrtmp;
// if(irep > 1) printf("After 3d\n");
            }

            // 4) ucSlab data is written to disk file
            //    Assume slab spans X- and Y-dims of tile
            i64Offset = tile_byteoff[node] + (__int64)(itn * jtn * (kil-ktl)/irep);

//            printf("XYS-Dim = %3d %3d %2d  nbytes=%2d  val=%3d  seek=%15lld\n",
//                itn, jtn, (kih-kil+1)/irep, nslab, (int)ucSlab[nslab/3], i64Offset);
//            printf("XYS-Dim = %3d %3d %2d  irep=%2d  ki[lh]=%3d %3d  kt[lh]=%3d %3d\n",
//                itn, jtn, (kih-kil+1)/irep, irep, kil, kih, ktl, kth);
#ifndef LINUX
            // TEST MOD
            _lseeki64(fd, i64Offset, SEEK_SET);
            // TEST MOD
#else
            // TEST MOD
            lseek(fd,i64Offset,SEEK_SET);
            // TEST MOD
#endif
            // TEST MOD
            write(fd, ucSlab, nslab);
            // TEST MOD
        }                            // end of iIntersection if-block

    }                                // end of node loop


    fRefTime   += (float)(GetTickCount() - ulRefStart) / (float)1000.0;
    fTotalTime += (float)(GetTickCount() - ulStart) / (float)1000.0;


    // Output Tree
//    printf("z:offset=%d   z:size=%d   z:offset+size=%d\n", curoff[2],curdim[2],curoff[2]+curdim[2]);
/* TEST MOD */
    if(curoff[2]+curdim[2] >= nz) {

        // TEST MOD
        printf("0s,Av,Ref,Total = %f  %f  %f  %f\n", fwrite0Time, fAvTime, fRefTime, fTotalTime);

        iOff = 0;
        for(node=0; node<iNnodes; node++) {
#ifndef LINUX
            sprintf(cLine, "%s %10I64d %2d   %4d %4d %4d   %5d %5d %5d %8d %3d\n",
                pfile,                    tile_byteoff[node],        tile_rep[node],
                tile_dim[node][0],        tile_dim[node][1],        tile_dim[node][2],
                tile_offset[node][0],    tile_offset[node][1],    tile_offset[node][2],
                tile_parent[node],        tile_error[node]                            );
#else
            sprintf(cLine, "%s %10lld %2d   %4d %4d %4d   %5d %5d %5d %8d %3d\n",
                pfile,                    tile_byteoff[node],        tile_rep[node],
                tile_dim[node][0],        tile_dim[node][1],        tile_dim[node][2],
                tile_offset[node][0],    tile_offset[node][1],    tile_offset[node][2],
                tile_parent[node],        tile_error[node]                            );
#endif

            nChar = 0;
            while ( (int)cLine[nChar] != (int)'\n'  &&  nChar < 256 ) { nChar++; }
            nChar++;            
//            if(nChar > 120) { printf("nChar = %d\n", nChar); nChar = 120; }
            for(i=0; i<nChar; i++) { ucSlab[iOff+i] = (unsigned char)cLine[i]; }
            iOff += nChar;
        }
        
        if(fLimits[0] < fLimits[1]) {
            sprintf(cLine, "XYZ_Limits %f %f %f %f %f %f\n",
                fLimits[0], fLimits[1], fLimits[2], fLimits[3], fLimits[4], fLimits[5]);

            nChar = 0;
            while ( (int)cLine[nChar] != (int)'\n'  &&  nChar < 256 ) { nChar++; }
            nChar++;
//            if(nChar > 200) { printf("nChar = %d\n", nChar); nChar = 120; }
            for(i=0; i<nChar; i++) { ucSlab[iOff+i] = (unsigned char)cLine[i]; }
            iOff += nChar;
        }

        // Pad header with a whole bunch of charage returns
        sprintf(cLine, "\n");
        for(i=0; i<256; i++) { ucSlab[iOff+i] = (unsigned char)cLine[0]; }
        nChar += iOff + 256;
#ifndef LINUX
        _lseeki64(fd, 0, SEEK_SET);
#else
        lseek(fd,0,SEEK_SET);
#endif
        write(fd, ucSlab, nChar);
#ifndef LINUX
        _close(fd);
#else
        close(fd);
#endif
    }

/* TEST MOD */  
}

void TreeScanOffsets(int node, __int64 *pbyteoff, __int64 tile_byteoff[MAX_NODES],
                     int irep, int tile_rep[MAX_NODES], int tile_dim[MAX_NODES][3],
                     int tile_N_children[MAX_NODES], int tile_children[MAX_NODES][MAX_CHILDREN])
{
    if(irep == tile_rep[node]) {
        tile_byteoff[node] = *pbyteoff;
        *pbyteoff += (__int64)(FILE_BLOCK*(1 + (tile_dim[node][0]*tile_dim[node][1]*tile_dim[node][2]-1)/
                                               (FILE_BLOCK)));
    } else {
        int ichild;
        for(ichild=0; ichild<tile_N_children[node]; ichild++) {
            TreeScanOffsets(tile_children[node][ichild], pbyteoff, tile_byteoff,
                            irep, tile_rep, tile_dim,
                            tile_N_children, tile_children                       );
        }
    }

    return;
}


/*******************************************************************
 *               WRAPPER FOR CALLS FROM FORTRAN                    *
 *******************************************************************/

#ifdef _WIN32
#define block_tree_ BLOCK_TREE
#endif

#define CARG(t,s,n) i = 0; while(i<n  &&  s[i] != ' ') { t[i] = s[i]; i++; } t[i] = 0

void block_tree_(char *pfile, int fulldim[3], int curdim[3], int curoff[3], int *pblock,
                 float fLimits[6], int *pimap, unsigned char *field, int n1)
{
  int i, imap;
  char str1[1024];
  CARG(str1,pfile,n1);
  imap = *pimap;
  Block_Tree(str1, fulldim, curdim, curoff, pblock, fLimits, imap, field);
/*
void Block_Tree(char *pfile, int fulldim[3], int curdim[3], int curoff[3], int *pblock,
                float fLimits[6], int imap, unsigned char *field)
*/
}

