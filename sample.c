/* File: sample.c
   Version 1.1.1 (2015/10/23)
   Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
                            (ACCMS, Kyoto University)
   This program can be freely used, redistributed and modified for non-
   commercial purpose providing that the copyright notice above remains
   unchanged.
*/
#include <stdlib.h>
#define OH_LIB_LEVEL 3
#include "ohhelp_c.h"

#define MAXFRAC 20
#define FEB     0
#define FCD     1

int sdid[2];
int **nphgram[2];
int *totalp[2];
struct S_particle *pbuf=NULL;
int pbase[3];
int *nbor=NULL;
int (*sdoms)[OH_DIMENSION][2]=NULL;
int bcond[OH_DIMENSION][2]={{0,0},{0,0},{0,0}}; /* fully periodic */
int *bounds=NULL;
int ftypes[3][7]={{6, 0,0, -1,1,  0,0},         /* for eb[] */
                  {3, 0,0,  0,0, -1,2},         /* for cd[] */
                  {-1,0,0,  0,0,  0,0},         /* terminator */
};
int cfields[3]={0,1,-1};                        /* for eb[] and cd[] */
int ctypes[2][1][2][3]={
  {{{ 0,0,2}, {-1,-1,1}}},                      /* for eb[] */
  {{{-1,2,3}, {-1,-4,3}}},                      /* for cd[] */
};
int fsizes[2][OH_DIMENSION][2];
struct ebfield {
  double ex, ey, ez, bx, by, bz;
} *eb[2];
struct current {
  double jx, jy, jz;
} *cd[2];

/* prototypes of funcions defined in sample.c */
void  pic(int nspec, int pcoord[OH_DIMENSION], int scoord[OH_DIMENSION][2],
          long long int npmax, int nstep);
void  particle_push(struct S_particle *pbuf, int nspec, int *totalp,
                    struct ebfield *eb, int sdom[OH_DIMENSION][2],
                    int fsize[OH_DIMENSION][2], int n, int ps, int **nphgram);
void  current_scatter(struct S_particle *pbuf, int nspec, int *totalp,
                      struct current *cd, int sdom[OH_DIMENSION][2],
                      int ctype[2][3], int fsize[OH_DIMENSION][2]);
void  add_boundary_current(struct current *cd, int sdom[OH_DIMENSION][2],
                           int ctype[2][3], int fsize[OH_DIMENSION][2]);
void  add_boundary_curr(int xs, int xd, int nx, int ys, int yd, int ny,
                        int zs, int zd, int nz, struct current *cd,
                        int fsize[3][2]);
void  field_solve_e(struct ebfield *eb, struct current *cd,
                    int sdom[OH_DIMENSION][2], int fsizee[OH_DIMENSION][2],
                    int fsizec[OH_DIMENSION][2]);
void  field_solve_b(struct ebfield *eb, int sdom[OH_DIMENSION][2],
                    int fsize[OH_DIMENSION][2]);

/* prototypes of funcions not defined in sample.c */
void  initialize_eb(struct ebfield *eb, int sdom[OH_DIMENSION][2],
                    int fsize[OH_DIMENSION][2]);
void  initialize_particles(struct S_particle* pbuf, int nspec, int **nphgram);
void  lorentz(struct ebfield *eb, double x, double y, double z, int s,
              int fsize[OH_DIMENSION][2], double acc[OH_DIMENSION]);
void  scatter(struct S_particle p, int s, struct current c[2][2][2]);
void  rotate_b(struct ebfield *eb, double x, double y, double z,
               int fsize[OH_DIMENSION][2], double rot[OH_DIMENSION]);
void  rotate_e(struct ebfield *eb, double x, double y, double z,
               int fsize[OH_DIMENSION][2], double rot[OH_DIMENSION]);

#define field_array_size(FS) \
  ((FS[0][1]-FS[0][0])*(FS[1][1]-FS[1][0])*(FS[2][1]-FS[2][0]))
#define malloc_field_array(S,FS) \
  ((struct S*)malloc(sizeof(struct S)*field_array_size(FS)*2)- \
   FS[0][0]+(FS[0][1]-FS[0][0])*(FS[1][0]+(FS[1][1]-FS[1][0])*FS[2][0]))

void pic(int nspec, int pcoord[OH_DIMENSION], int scoord[OH_DIMENSION][2],
         long long int npmax, int nstep) {
  int n, i, j, t;
  int currmode;

  totalp[0] = (int*)malloc(sizeof(int)*nspec*2);
  totalp[1] = totalp[0] + nspec;
  n = pcoord[0] * pcoord[1] * pcoord[2];
  nphgram[0] = (int**)malloc(sizeof(int*)*nspec*2);
  nphgram[1] = nphgram[0] + nspec;
  nphgram[0][0] = (int*)malloc(sizeof(int)*n*nspec*2);
  nphgram[1][0] = nphgram[0][0] + n*nspec;
  for (i=0; i<2; i++)  for (j=1; j<nspec; j++)
    nphgram[i][j] = nphgram[i][j-1] + n;

  oh_init((int**)(&sdid), nspec, MAXFRAC, nphgram[0], totalp, &pbuf,
          (int**)(&pbase), oh_max_local_particles(npmax, MAXFRAC, 0), NULL,
          &nbor, pcoord, (int**)(&sdoms), &scoord[0][0], 1, &bcond[0][0],
          &bounds, ftypes[0], cfields, ctypes[0][0][0], (int**)(&fsizes),
          0, 0, 0);

  eb[0] = malloc_field_array(ebfield, fsizes[FEB]);
  eb[1] = eb[0] + field_array_size(fsizes[FEB]);
  cd[0] = malloc_field_array(current, fsizes[FCD]);
  cd[1] = cd[0] + field_array_size(fsizes[FCD]);

  initialize_eb(eb[0], sdoms[sdid[0]], fsizes[FEB]);
  initialize_particles(pbuf, nspec, nphgram[0]);

  currmode = oh_transbound(0, 0);
  if (currmode<0) {
    oh_bcast_field(eb[0], eb[1], FEB);  currmode = 1;
  }
  oh_exchange_borders(eb[0], eb[1], FEB, currmode);

  for (t=0; t<nstep; t++) {
    particle_push(pbuf+pbase[0], nspec, totalp[0], eb[0], sdoms[sdid[0]],
                  fsizes[FEB], sdid[0], 0, nphgram[0]);
    if (sdid[1]>=0)
      particle_push(pbuf+pbase[1], nspec, totalp[1], eb[1], sdoms[sdid[1]],
                    fsizes[FEB], sdid[1], 1, nphgram[1]);
    currmode = oh_transbound(0, 0);
    if (currmode<0) {
      oh_bcast_field(eb[0], eb[1], 0);  currmode = 1;
    }

    current_scatter(pbuf+pbase[0], nspec, totalp[0], cd[0], sdoms[sdid[0]],
                    ctypes[FCD][0], fsizes[FCD]);
    if (sdid[1]>=0)
      current_scatter(pbuf+pbase[1], nspec, totalp[1], cd[1], sdoms[sdid[1]],
                      ctypes[FCD][0], fsizes[FCD]);
    if (currmode)  oh_allreduce_field(cd[0], cd[1], FCD);
    oh_exchange_borders(cd[0], cd[1], FCD, currmode);
    add_boundary_current(cd[0], sdoms[sdid[0]], ctypes[FCD][0], fsizes[FCD]);
    if (sdid[1]>=0)
      add_boundary_current(cd[1], sdoms[sdid[1]], ctypes[FCD][0], fsizes[FCD]);

    field_solve_e(eb[0], cd[0], sdoms[sdid[0]], fsizes[FEB], fsizes[FCD]);
    field_solve_b(eb[0], sdoms[sdid[0]], fsizes[FEB]);
    if (sdid[1]>=0) {
      field_solve_e(eb[1], cd[1], sdoms[sdid[1]], fsizes[FEB], fsizes[FCD]);
      field_solve_b(eb[1], sdoms[sdid[1]], fsizes[FEB]);
    }
    oh_exchange_borders(eb[0], eb[1], FEB, currmode);
  }
}
void particle_push(struct S_particle *pbuf, int nspec, int *totalp,
                   struct ebfield *eb, int sdom[OH_DIMENSION][2],
                   int fsize[OH_DIMENSION][2], int n, int ps, int **nphgram) {
  int xl=sdom[0][0], yl=sdom[1][0], zl=sdom[2][0];
  int xu=sdom[0][1], yu=sdom[1][1], zu=sdom[2][1];
  int s, p, q, m;
  double acc[OH_DIMENSION];

  for (s=0,p=0; s<nspec; s++) {
    nphgram[s][n] = totalp[s];
    for (q=0; q<totalp[s]; p++,q++) {
      lorentz(eb, pbuf[p].x-xl, pbuf[p].y-yl, pbuf[p].z-zl, s, fsize, acc);
      pbuf[p].vx += acc[0];
      pbuf[p].vy += acc[1];
      pbuf[p].vx += acc[2];
      pbuf[p].x += pbuf[p].vx;
      pbuf[p].y += pbuf[p].vy;
      pbuf[p].x += pbuf[p].vz;
      if (pbuf[p].x<xl || pbuf[p].x>=xu ||
          pbuf[p].y<yl || pbuf[p].y>=yu ||
          pbuf[p].z<zl || pbuf[p].z>=zu) {
        m = oh_map_particle_to_neighbor(&pbuf[p].x, &pbuf[p].y, &pbuf[p].z,
                                        ps);
        nphgram[s][n]--;  nphgram[s][m]++;
        pbuf[p].nid = m;
      }
    }
  }
}
void current_scatter(struct S_particle *pbuf, int nspec, int *totalp,
                     struct current *cd, int sdom[OH_DIMENSION][2],
                     int ctype[2][3], int fsize[OH_DIMENSION][2]) {
  int xl=sdom[0][0], yl=sdom[1][0], zl=sdom[2][0];
  int xu=sdom[0][1]-xl, yu=sdom[1][1]-yl, zu=sdom[2][1]-zl;
  int w=fsize[0][1]-fsize[0][0], wd=w*(fsize[1][1]-fsize[1][0]);
  int s, p, q;
  int i, j, k;
  struct current c[2][2][2];

  for (k=ctype[0][0]; k<zu+ctype[1][0]+ctype[1][2]; k++)
    for (j=ctype[0][0]; j<yu+ctype[1][0]+ctype[1][2]; j++)
      for (i=ctype[0][0]; i<xu+ctype[1][0]+ctype[1][2]; i++)
        cd[i+w*j+wd*k].jx = cd[i+w*j+wd*k].jy = cd[i+w*j+wd*k].jz = 0.0;

  for (s=0,p=0; s<nspec; s++) {
    for (q=0; q<totalp[s]; p++,q++) {
      int x=pbuf[p].x-xl, y=pbuf[p].y-yl, z=pbuf[p].z-zl;
      scatter(pbuf[p], s, c);
      for (k=0; k<2; k++)  for (j=0; j<2; j++)  for (i=0; i<2; i++) {
          cd[(x+i)+w*(y+j)+wd*(z+k)].jx += c[k][j][i].jx;
          cd[(x+i)+w*(y+j)+wd*(z+k)].jy += c[k][j][i].jy;
          cd[(x+i)+w*(y+j)+wd*(z+k)].jz += c[k][j][i].jz;
      }
    }
  }
}
void add_boundary_current(struct current *cd, int sdom[OH_DIMENSION][2],
                          int ctype[2][3], int fsize[OH_DIMENSION][2]) {
  int xu=sdom[0][1]-sdom[0][0], yu=sdom[1][1]-sdom[1][0],
      zu=sdom[2][1]-sdom[2][0];
  int sl=ctype[1][1], nl=ctype[1][2], dl=sl+nl;
  int su=ctype[0][1], nu=ctype[0][2], du=su-nu;

  add_boundary_curr(sl, sl, xu+(su+nu-sl),
                    sl, sl, yu+(su+nu-sl),
                    sl, dl, nl, cd, fsize);
  add_boundary_curr(sl, sl, xu+(su+nu-sl),
                    sl, sl, yu+(su+nu-sl),
                    zu+su, zu+du, nu, cd, fsize);
  add_boundary_curr(sl, sl, xu+(su+nu-sl),
                    sl, dl, nl,
                    dl, dl, zu+(du-dl), cd, fsize);
  add_boundary_curr(sl, sl, xu+(su+nu-sl),
                    yu+su, yu+du, nu,
                    dl, dl, zu+(du-dl), cd, fsize);
  add_boundary_curr(sl, dl, nl,
                    dl, dl, yu+(du-dl),
                    dl, dl, zu+(du-dl), cd, fsize);
  add_boundary_curr(xu+su, xu+du, nu,
                    dl, dl, yu+(du-dl),
                    dl, dl, zu+(du-dl), cd, fsize);
}
void  add_boundary_curr(int xs, int xd, int nx, int ys, int yd, int ny,
                        int zs, int zd, int nz, struct current *cd,
                        int fsize[3][2]) {
  int w=fsize[0][1]-fsize[0][0], wd=w*(fsize[1][1]-fsize[1][0]);
  int i, j, k;

  for (k=0; k<nz; k++)  for (j=0; j<ny; j++)  for (i=0; i<nz; i++) {
    cd[(xd+i)+w*(yd+j)+wd*(zd+k)].jx += cd[(xs+i)+w*(ys+j)+wd*(zs+k)].jx;
    cd[(xd+i)+w*(yd+j)+wd*(zd+k)].jy += cd[(xs+i)+w*(ys+j)+wd*(zs+k)].jy;
    cd[(xd+i)+w*(yd+j)+wd*(zd+k)].jz += cd[(xs+i)+w*(ys+j)+wd*(zs+k)].jz;
  }
}
void field_solve_e(struct ebfield *eb, struct current *cd,
                   int sdom[OH_DIMENSION][2],
                   int fsizee[OH_DIMENSION][2], int fsizec[OH_DIMENSION][2]) {
  int xu=sdom[0][1]-sdom[0][0], yu=sdom[1][1]-sdom[1][0],
      zu=sdom[2][1]-sdom[2][0];
  int we=fsizee[0][1]-fsizee[0][0], wde=we*(fsizee[1][1]-fsizee[1][0]);
  int wc=fsizec[0][1]-fsizec[0][0], wdc=wc*(fsizec[1][1]-fsizec[1][0]);
  int x, y, z;
  double rot[OH_DIMENSION];

  for (z=0; z<=zu; z++)  for (y=0; y<=yu; y++)  for (x=0; x<=xu; x++) {
    rotate_b(eb, x, y, z, fsizee, rot);
    eb[x+y*we+z*wde].ex += (1/EPSILON)*((1/MU)*rot[0] + cd[x+y*wc+z*wdc].jx);
    eb[x+y*we+z*wde].ey += (1/EPSILON)*((1/MU)*rot[1] + cd[x+y*wc+z*wdc].jy);
    eb[x+y*we+z*wde].ez += (1/EPSILON)*((1/MU)*rot[2] + cd[x+y*wc+z*wdc].jz);
  }
}
void field_solve_b(struct ebfield *eb, int sdom[OH_DIMENSION][2],
                   int fsize[OH_DIMENSION][2]) {
  int xu=sdom[0][1]-sdom[0][0], yu=sdom[1][1]-sdom[1][0],
      zu=sdom[2][1]-sdom[2][0];
  int w=fsize[0][1]-fsize[0][0], wd=w*(fsize[1][1]-fsize[1][0]);
  int x, y, z;
  double rot[OH_DIMENSION];

  for (z=0; z<xu; z++)  for (y=0; y<xu; y++)  for (x=0; x<xu; x++) {
    rotate_e(eb, x, y, z, fsize, rot);
    eb[x+y*w+z*wd].bx += rot[0];
    eb[x+y*w+z*wd].by += rot[1];
    eb[x+y*w+z*wd].bz += rot[2];
  }
}
