/* File: ohhelp3.c
   Version 1.1.1 (2015/10/23)
   Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
                            (ACCMS, Kyoto University)
   This program can be freely used, redistributed and modified for non-
   commercial purpose providing that the copyright notice above remains
   unchanged.
*/
#define EXTERN extern
#include "ohhelp1.h"
#include "ohhelp2.h"
#undef  EXTERN
#define EXTERN
#include "ohhelp3.h"

static void init_subdomain_actively(int (*sd)[OH_DIMENSION][2],
                                    int sc[OH_DIMENSION][2],
                                    int *pcoord, int bc[OH_DIMENSION][2],
                                    int (*bd)[OH_DIMENSION][2], int nb,
                                    int bbase);
static void init_subdomain_passively(int (*sd)[OH_DIMENSION][2],
                                     int (*bd)[OH_DIMENSION][2], int nb,
                                     int bbase);
static int  comp_xyz(const void* aa, const void* bb);
static void init_fields(int (*ft)[OH_FTYPE_N], int *cf, int cfid,
                        int (*ct)[2][OH_CTYPE_N], int nb,
                        int sd[OH_DIMENSION][2], int **fsizes);
static void set_border_exchange(int e, int ps, MPI_Datatype type);
static void set_border_comm(int esize, int f, int *xyz, int *wdh,
                            int (*exti)[2], int (*exto)[2],
                            int (*off)[2], int (*size)[2],
                            int lu, int sr, MPI_Datatype basetype,
                            struct S_borderexc bx[OH_DIMENSION][2]);
static int  transbound3(int currmode, int stats, int level);
static int  map_irregular(double p0, double p1, double p2, int dim, int from,
                          int n);
static int  map_irregular_range(double p, int dim, int from, int to);

void
oh3_init_(int *sdid, int *nspec, int *maxfrac, int *nphgram,
          int *totalp, struct S_particle *pbuf, int *pbase, int *maxlocalp,
          struct S_mycommf *mycomm, int *nbor, int *pcoord,
          int *sdoms, int *scoord, int *nbound, int *bcond, int *bounds,
          int *ftypes, int *cfields, int *ctypes, int *fsizes,
          int *stats, int *repiter, int *verbose) {
  specBase = 1;
  init3(&sdid, *nspec, *maxfrac, &nphgram, &totalp, NULL, NULL, &pbuf, &pbase,
        *maxlocalp, NULL, mycomm, &nbor, pcoord, &sdoms, scoord, *nbound,
        bcond, &bounds, ftypes, cfields, -1, ctypes, &fsizes,
        *stats, *repiter, *verbose, 0);
}
void
oh3_init(int **sdid, int nspec, int maxfrac, int **nphgram,
         int **totalp, struct S_particle **pbuf, int **pbase, int maxlocalp,
         void *mycomm, int **nbor, int *pcoord,
         int **sdoms, int *scoord, int nbound, int *bcond, int **bounds,
         int *ftypes, int *cfields, int *ctypes, int **fsizes,
         int stats, int repiter, int verbose) {
  specBase = 0;
  init3(sdid, nspec, maxfrac, nphgram, totalp, NULL, NULL, pbuf, pbase,
        maxlocalp, (struct S_mycommc*)mycomm, NULL, nbor, pcoord, sdoms,
        scoord, nbound, bcond, bounds, ftypes, cfields, 0, ctypes, fsizes,
        stats, repiter, verbose, 0);
}
void
oh13_init_(int *sdid, int *nspec, int *maxfrac, int *nphgram,
           int *totalp, int *rcounts, int *scounts,
           struct S_mycommf *mycomm, int *nbor, int *pcoord,
           int *sdoms, int *scoord, int *nbound, int *bcond, int *bounds,
           int *ftypes, int *cfields, int *ctypes, int *fsizes,
           int *stats, int *repiter, int *verbose) {
  init3(&sdid, *nspec, *maxfrac, &nphgram, &totalp, &rcounts, &scounts,
        NULL, NULL, 0, NULL, mycomm, &nbor, pcoord, &sdoms, scoord, *nbound,
        bcond, &bounds, ftypes, cfields, -1, ctypes, &fsizes,
        *stats, *repiter, *verbose, 1);
}
void
oh13_init(int **sdid, int nspec, int maxfrac, int **nphgram,
          int **totalp, int **rcounts, int **scounts,
          void *mycomm, int **nbor, int *pcoord,
          int **sdoms, int *scoord, int nbound, int *bcond, int **bounds,
          int *ftypes, int *cfields, int *ctypes, int **fsizes,
          int stats, int repiter, int verbose) {
  init3(sdid, nspec, maxfrac, nphgram, totalp, rcounts, scounts, NULL, NULL,
        0, (struct S_mycommc*)mycomm, NULL, nbor, pcoord, sdoms, scoord,
        nbound, bcond, bounds, ftypes, cfields, 0, ctypes, fsizes,
        stats, repiter, verbose, 1);
}
void
init3(int **sdid, int nspec, int maxfrac, int **nphgram,
      int **totalp, int **rcounts, int **scounts,
      struct S_particle **pbuf, int **pbase, int maxlocalp,
      struct S_mycommc *mycommc, struct S_mycommf *mycommf,
      int **nbor, int *pcoord, int **sdoms, int *scoord,
      int nbound, int *bcond, int **bounds, int *ftypes,
      int *cfields, int cfid, int *ctypes, int **fsizes,
      int stats, int repiter, int verbose, int skip2) {
  int nn;
  int (*sd)[OH_DIMENSION][2]=(int(*)[OH_DIMENSION][2])*sdoms;
  double (*sdf)[OH_DIMENSION][2];
  int (*sc)[2]=(int(*)[2])scoord;
  int (*bc)[2]=(int(*)[2])bcond;
  int (*bd)[OH_DIMENSION][2]=(int(*)[OH_DIMENSION][2])*bounds;
  int (*ft)[OH_FTYPE_N]=(int(*)[OH_FTYPE_N])ftypes;
  int (*ct)[2][OH_CTYPE_N]=(int(*)[2][OH_CTYPE_N])ctypes;
  int d, n, m;

  if (skip2)
    init1(sdid, nspec, maxfrac, nphgram, totalp, rcounts, scounts,
          mycommc, mycommf, nbor, pcoord, stats, repiter, verbose);
  else
    init2(sdid, nspec, maxfrac, nphgram, totalp, pbuf, pbase, maxlocalp,
          mycommc, mycommf, nbor, pcoord, stats, repiter, verbose);
  excludeLevel2 = skip2;
  nn = nOfNodes;

  if (!sd) {
    sd = (int(*)[OH_DIMENSION][2])
         (*sdoms = (int*)mem_alloc(sizeof(int), nn*OH_DIMENSION*2,
                                   "SubDomains"));
    sd[0][OH_DIM_X][OH_LOWER] = 0;  sd[0][OH_DIM_X][OH_UPPER] = -1;
  }
  if (!bd)
    bd = (int(*)[OH_DIMENSION][2])
         (*bounds = (int*)mem_alloc(sizeof(int), nn*OH_DIMENSION*2,
                                    "Boundaries"));

  for (d=0,n=1,m=OH_NEIGHBORS>>1; d<OH_DIMENSION; d++,n*=3) {
    int nl=DstNeighbors[m-n], nu=DstNeighbors[m+n];
    Adjacent[d][OH_LOWER] = nl<0 ? -(nl+1) : nl;
    Adjacent[d][OH_UPPER] = nu<0 ? -(nu+1) : nu;
  }
  if (sd[0][OH_DIM_X][OH_LOWER]>sd[0][OH_DIM_X][OH_UPPER])
    init_subdomain_actively(sd, sc, pcoord, bc, bd, nbound, -cfid);
  else
    init_subdomain_passively(sd, bd, nbound, -cfid);

  SubDomains = (int(*)[OH_DIMENSION][2])
               mem_alloc(sizeof(int), nn*OH_DIMENSION*2, "SubDomains");
  sdf = SubDomainsFloat =
        (double(*)[OH_DIMENSION][2])
        mem_alloc(sizeof(double), nn*OH_DIMENSION*2, "SubDomainsFloat");
  Boundaries = (int(*)[OH_DIMENSION][2])
               mem_alloc(sizeof(int), nn*OH_DIMENSION*2, "Boundaries");
  memcpy(SubDomains, sd, sizeof(int)*nn*OH_DIMENSION*2);
  for (n=0; n<nn; n++)  for (d=0; d<OH_DIMENSION; d++) {
    sdf[n][d][OH_LOWER] = sd[n][d][OH_LOWER];
    sdf[n][d][OH_UPPER] = sd[n][d][OH_UPPER];
  }
  memcpy(Boundaries, bd, sizeof(int)*nn*OH_DIMENSION*2);
  bd = Boundaries;
  if (cfid) {
    for (n=0; n<nn; n++)  for (d=0; d<OH_DIMENSION; d++) {
      bd[n][d][OH_LOWER]--;  bd[n][d][OH_UPPER]--;
    }
  }
  init_fields(ft, cfields, cfid, ct, nbound, sd[myRank], fsizes);
}
static void
init_subdomain_actively(int (*sd)[OH_DIMENSION][2], int sc[OH_DIMENSION][2],
                        int *pcoord, int bc[OH_DIMENSION][2],
                        int (*bd)[OH_DIMENSION][2], int nb, int bbase) {
  int nn=nOfNodes, pqr=1;
  int d, lu, i, j, k, x, y, z, n;

  SubDomainDesc = NULL;
  for (d=0; d<OH_DIMENSION; d++) {
    int lo = Grid[d].coord[OH_LOWER] = sc[d][OH_LOWER];
    int up = Grid[d].coord[OH_UPPER] = sc[d][OH_UPPER];
    int size = up - lo;
    int ave, nl;
    Grid[d].fcoord[OH_LOWER] = lo;  Grid[d].fcoord[OH_UPPER] = up;
    n = Grid[d].n = pcoord[d];
    if (n<=0)
      errstop("# of %c-nodes (%d) should be positive", Message.xyz[d], n);
    if (size<=0)
      errstop("upper edge of %c-coordinate (%d) should be greater than "
              "lower edge (%d)", Message.xyz[d], up, lo);
    ave = Grid[d].light.size = size/n;
    Grid[d].light.rfsize = 1.0/(double)ave;
    Grid[d].light.rfsizeplus = 1.0/(double)(ave+1);
    nl =  Grid[d].light.n = n - size%n;
    Grid[d].light.fthresh = (Grid[d].light.thresh = lo + nl * ave);
    Grid[d].fsize = (Grid[d].size = n==nl ? ave : ave+1);
    Grid[d].gsize = Grid[d].rgsize = 1.0;
    pqr *= n;
  }
  for (; d<3; d++) {
    Grid[d].n = Grid[d].light.n = 1;
    Grid[d].coord[OH_LOWER] = Grid[d].coord[OH_UPPER] = 0;
    Grid[d].fcoord[OH_LOWER] = Grid[d].fcoord[OH_UPPER] = 0.0;
    Grid[d].size = Grid[d].light.size = Grid[d].light.thresh = 0;
    Grid[d].fsize = Grid[d].light.rfsize
                  = Grid[d].light.rfsizeplus = Grid[d].light.fthresh = 0.0;
    Grid[d].gsize = Grid[d].rgsize = 1.0;
  }
  if (pqr!=nn) {
    if (OH_DIMENSION==1)
      errstop("<# of x-nodes>(%d) should be eqal to <# of nodes>(%d)",
              pcoord[0], nn);
    else if (OH_DIMENSION==2)
      errstop("<# of x-nodes>(%d) * <# of y-nodes>(%d) "
              "should be eqal to <# of nodes>(%d)",
              pcoord[0], pcoord[1], nn);
    else
      errstop("<# of x-nodes>(%d) * <# of y-nodes>(%d) * <# of z-nodes>(%d) "
              "should be eqal to <# of nodes>(%d)",
              pcoord[0], pcoord[1], pcoord[2], nn);
  }
  for (d=0; d<OH_DIMENSION; d++) {
    for (lu=OH_LOWER; lu<=OH_UPPER; lu++) {
      if (bc[d][lu]<bbase || bc[d][lu]>=nb+bbase)
        errstop("system's %s boundary condition for %c-coordinate %d is "
                "invalid",
                Message.loup[lu], Message.xyz[d], bc[d][lu]);
    }
  }
  for (i=0,z=Grid[OH_DIM_Z].coord[OH_LOWER],n=0; i<Grid[OH_DIM_Z].n; i++) {
    int bot=z, top=z+Grid[OH_DIM_Z].light.size;
    int bzlo = i==0 && OH_DIMENSION>OH_DIM_Z ?
                 bc[OH_DIM_Z][OH_LOWER] : bbase;
    int bzup = i==Grid[OH_DIM_Z].n-1  && OH_DIMENSION>OH_DIM_Z ?
                 bc[OH_DIM_Z][OH_UPPER] : bbase;
    if (i>=Grid[OH_DIM_Z].light.n)  top++;
    z = top;
    for (j=0,y=Grid[OH_DIM_Y].coord[OH_LOWER]; j<Grid[OH_DIM_Y].n; j++) {
      int south=y, north=y+Grid[OH_DIM_Y].light.size;
      int bylo = j==0 && OH_DIMENSION>OH_DIM_Y ?
                   bc[OH_DIM_Y][OH_LOWER] : bbase;
      int byup = j==Grid[OH_DIM_Y].n-1  && OH_DIMENSION>OH_DIM_Y ?
                   bc[OH_DIM_Y][OH_UPPER] : bbase;
      if (j>=Grid[OH_DIM_Y].light.n)  north++;
      y = north;
      for (k=0,x=Grid[OH_DIM_X].coord[OH_LOWER]; k<Grid[OH_DIM_X].n; k++,n++) {
        int west=x, east=x+Grid[OH_DIM_X].light.size;
        if (k>=Grid[OH_DIM_X].light.n)  east++;
        x = east;
        sd[n][OH_DIM_X][OH_LOWER] = west;
        sd[n][OH_DIM_X][OH_UPPER] = east;
        bd[n][OH_DIM_X][OH_LOWER] = bd[n][OH_DIM_X][OH_UPPER] = bbase;
        if (OH_DIMENSION>OH_DIM_Y) {
          sd[n][OH_DIM_Y][OH_LOWER] = south;
          sd[n][OH_DIM_Y][OH_UPPER] = north;
          bd[n][OH_DIM_Y][OH_LOWER] = bylo;
          bd[n][OH_DIM_Y][OH_UPPER] = byup;
        }
        if (OH_DIMENSION>OH_DIM_Z) {
          sd[n][OH_DIM_Z][OH_LOWER] = bot;
          sd[n][OH_DIM_Z][OH_UPPER] = top;
          bd[n][OH_DIM_Z][OH_LOWER] = bzlo;
          bd[n][OH_DIM_Z][OH_UPPER] = bzup;
        }
      }
      bd[n-Grid[OH_DIM_X].n][OH_DIM_X][OH_LOWER] = bc[OH_DIM_X][OH_LOWER];
      bd[n-1][OH_DIM_X][OH_UPPER] = bc[OH_DIM_X][OH_UPPER];
    }
  }
}
static void
init_subdomain_passively(int (*sd)[OH_DIMENSION][2],
                         int (*bd)[OH_DIMENSION][2], int nb, int bbase) {
  int nn=nOfNodes;
  struct S_subdomdesc *sdd = SubDomainDesc =
    (struct S_subdomdesc*)mem_alloc(sizeof(struct S_subdomdesc), nn,
                                    "SubDomainDesc");
  int min[OH_DIMENSION], max[OH_DIMENSION];
  int smin[OH_DIMENSION], smax[OH_DIMENSION];
  int me=myRank;
  int i, d, dd, lu, l;
  int lo[OH_DIMENSION-1], up[OH_DIMENSION-1], h[OH_DIMENSION-1];

  for (d=0; d<OH_DIMENSION; d++) {
    min[d] = sd[0][d][OH_LOWER];  max[d] = sd[0][d][OH_UPPER];
    smin[d] = smax[d] = max[d] - min[d];
  }
  for (i=0; i<nn; i++) {
    for (d=0; d<OH_DIMENSION; d++) {
      int lo=sd[i][d][OH_LOWER], up=sd[i][d][OH_UPPER], n=up-lo;
      sdd[i].coord[d].fc[OH_LOWER] = (sdd[i].coord[d].c[OH_LOWER] = lo);
      sdd[i].coord[d].fc[OH_UPPER] = (sdd[i].coord[d].c[OH_UPPER] = up);
      sdd[i].coord[d].n = 0;
      if (n<smin[d])  smin[d] = n;
      if (n>smax[d])  smax[d] = n;
      if (lo<min[d])  min[d] = lo;
      if (up>max[d])  max[d] = up;
      if (n<=0)
        errstop("subdomain %d has %c-coordinate lower boundary %d "
                "not less than upper boundary %d", i, Message.xyz[d], lo, up);
      for (lu=OH_LOWER; lu<=OH_UPPER; lu++) {
        if (bd[i][d][lu]<bbase || bd[i][d][lu]>=nb+bbase)
          errstop("rank-%d's %s boundary condition for %c-coordinate %d is "
                  "invalid",
                  i, Message.loup[lu], Message.xyz[d], bd[i][d][lu]);
      }
    }
    sdd[i].id = i;
  }
  for (d=0; d<OH_DIMENSION; d++) {
    Grid[d].fsize = (Grid[d].size = smax[i]);
    Grid[d].light.size = smin[d];
    Grid[d].light.rfsize = Grid[d].light.rfsizeplus = 0.0;
    Grid[d].fcoord[OH_LOWER] = (Grid[d].coord[OH_LOWER] = min[d]);
    Grid[d].fcoord[OH_UPPER] = (Grid[d].coord[OH_UPPER] = max[d]);
    Grid[d].n = Grid[d].light.n = 0;    /* never referred but ... */
    Grid[d].light.thresh = 0;  Grid[d].light.fthresh = 0.0;
    Grid[d].gsize = Grid[d].rgsize = 1.0;
    for (lu=OH_LOWER; lu<=OH_UPPER; lu++) {
      int n=Adjacent[d][lu];
      if (n==nn || bd[me][d][lu]!=bbase) continue;
      for (dd=0; dd<OH_DIMENSION; dd++) {
        if (d==dd) {
          int diff = sd[n][dd][OH_UPPER-lu] - sd[me][dd][lu];
          int dsize = max[dd] - min[dd];
          if (diff!=0 && diff!=dsize && diff!=-dsize)
            local_errstop("rank-%d and its %c-%s neighbor rank-%d have "
                          "incompatible %s/%s boundaries of %c-coordinate "
                          "%d and %d",
                          me, Message.xyz[d], Message.loup[lu], n,
                          Message.loup[lu], Message.loup[OH_UPPER-lu],
                          Message.xyz[dd],
                          sd[me][dd][lu], sd[n][dd][OH_UPPER-lu]);
        } else {
          for (l=OH_LOWER; l<=OH_UPPER; l++) {
            if (sd[n][dd][l]!=sd[me][dd][l])
              local_errstop("rank-%d and its %c-%s neighbor rank-%d have "
                            "incompatible %s boundary of %c-coordinate "
                            "%d and %d",
                            me, Message.xyz[d], Message.loup[lu], n,
                            Message.loup[l], Message.xyz[dd],
                            sd[me][dd][l], sd[n][dd][l]);
          }
        }
      }
    }
  }
  for (; d<3; d++) {
    Grid[d].n = Grid[d].light.n = 1;
    Grid[d].coord[OH_LOWER] = Grid[d].coord[OH_UPPER] = 0;
    Grid[d].fcoord[OH_LOWER] = Grid[d].fcoord[OH_UPPER] = 0.0;
    Grid[d].size = Grid[d].light.size = Grid[d].light.thresh = 0;
    Grid[d].fsize = Grid[d].light.rfsize
                  = Grid[d].light.rfsizeplus = Grid[d].light.fthresh = 0.0;
    Grid[d].gsize = Grid[d].rgsize = 1.0;
  }
  qsort(sdd, nn, sizeof(struct S_subdomdesc), comp_xyz);
  for (d=0; d<OH_DIMENSION-1; d++) {
    sdd[0].coord[d].h = h[d] = 0;
    lo[d] = sdd[0].coord[d].c[OH_LOWER];
    up[d] = sdd[0].coord[d].c[OH_UPPER];
  }
  for (i=1; i<nn; i++) {
    for (d=0; d<OH_DIMENSION-1; d++) {
      if (lo[d]!=sdd[i].coord[d].c[OH_LOWER] ||
          up[d]!=sdd[i].coord[d].c[OH_UPPER]) {
        for (dd=d; dd<OH_DIMENSION-1; dd++) {
          sdd[h[dd]].coord[dd].n = i - h[dd];
          sdd[i].coord[dd].h = h[dd] = i;
          lo[dd] = sdd[i].coord[dd].c[OH_LOWER];
          up[dd] = sdd[i].coord[dd].c[OH_UPPER];
        }
        break;
      } else {
        sdd[i].coord[d].h = h[d];
      }
    }
    sdd[i].coord[OH_DIMENSION-1].n = 1;  sdd[i].coord[OH_DIMENSION-1].h = i;
  }
  for (d=0; d<OH_DIMENSION-1; d++)  sdd[h[d]].coord[d].n = nn - h[d];
}
static int
comp_xyz(const void* aa, const void* bb) {
  struct S_subdomdesc *a=(struct S_subdomdesc*)aa, *b=(struct S_subdomdesc*)bb;
  int d;

  for (d=0; d<OH_DIMENSION; d++) {
    if (a->coord[d].c[OH_LOWER]+a->coord[d].c[OH_UPPER]<
        b->coord[d].c[OH_LOWER]+b->coord[d].c[OH_UPPER])  return(-1);
    if (a->coord[d].c[OH_LOWER]+a->coord[d].c[OH_UPPER]>
        b->coord[d].c[OH_LOWER]+b->coord[d].c[OH_UPPER])  return(1);
    if (a->coord[d].c[OH_LOWER]<b->coord[d].c[OH_LOWER])  return(-1);
    if (a->coord[d].c[OH_LOWER]>b->coord[d].c[OH_LOWER])  return(1);
    if (a->coord[d].c[OH_UPPER]<b->coord[d].c[OH_UPPER])  return(-1);
    if (a->coord[d].c[OH_UPPER]>b->coord[d].c[OH_UPPER])  return(1);
  }
  return(a->id<b->id ? -1 : 1);
}
#if OH_DIMENSION==1
#define Field_Disp(F,X,Y,Z) (FieldDesc[F].esize * (X))
#elif OH_DIMENSION==2
#define Field_Disp(F,X,Y,Z)\
  (FieldDesc[F].esize *\
   ((X) + FieldDesc[F].size[OH_DIM_X] * (Y)))
#else
#define Field_Disp(F,X,Y,Z)\
  (FieldDesc[F].esize *\
   ((X) + FieldDesc[F].size[OH_DIM_X] *\
    ((Y) + FieldDesc[F].size[OH_DIM_Y] * (Z))))
#endif
static void
init_fields(int (*ft)[OH_FTYPE_N], int *cf, int cfid, int (*ct)[2][OH_CTYPE_N],
            int nb, int sd[OH_DIMENSION][2], int **fsizes) {
  struct S_flddesc *fd;
  struct S_borderexc (*bx)[2][OH_DIMENSION][2];
  int (*fs)[OH_DIMENSION][2]=(int(*)[OH_DIMENSION][2])*fsizes;
  int nf, ne;
  int f, e, b, d, lu, i, *tmp;

  nOfBoundaries = nb;
  for (nf=0; ft[nf][OH_FTYPE_ES]>0; nf++);
  nOfFields = nf;
#ifdef OH_POS_AWARE
  for (ne=0; cf[ne]>=0; ne++);
#else
  for (ne=0; cf[ne]+cfid>=0; ne++);
#endif
  nOfExc = ne;

  FieldDesc = fd = (struct S_flddesc*)mem_alloc(sizeof(struct S_flddesc), nf,
                                                "FieldDesc");
#ifndef OH_POS_AWARE
  FieldTypes = (int(*)[OH_FTYPE_N])
               mem_alloc(sizeof(int), nf*OH_FTYPE_N, "FieldTypes");
  BoundaryCommTypes  = (int(*)[2][OH_CTYPE_N])
                       mem_alloc(sizeof(int), ne*nb*2*OH_CTYPE_N,
                                 "BoundaryCommTypes");
  memcpy(FieldTypes, ft, sizeof(int)*nf*OH_FTYPE_N);
  memcpy(BoundaryCommTypes, ct, sizeof(int)*ne*nb*2*OH_CTYPE_N);
  ft = FieldTypes;  ct = BoundaryCommTypes;

  tmp = (int*)mem_alloc(sizeof(int), ne, "BoundaryCommFields");
  for (e=0; e<nOfExc; e++)  tmp[e] = cf[e] + cfid;
  BoundaryCommFields = cf = tmp;
#endif

  if (!fs)
    fs = (int(*)[OH_DIMENSION][2])
         (*fsizes = (int*)mem_alloc(sizeof(int), nf*OH_DIMENSION*2,
                                    "FieldSizes"));
  for (f=0; f<nf; f++) {
    int lo=ft[f][OH_FTYPE_LO], up=ft[f][OH_FTYPE_UP];
    fd[f].esize = ft[f][OH_FTYPE_ES];
    for (lu=OH_FTYPE_BL; lu<OH_FTYPE_RU; lu+=2) {
      int lot=ft[f][lu], upt=ft[f][lu+1];
      if (lot<lo)  lo = lot;
      if (upt>up)  up = upt;
    }
    fd[f].ext[OH_LOWER] = lo;  fd[f].ext[OH_UPPER] = up;
  }
  for (e=0,i=0; e<ne; e++) {
    int f=cf[e];
    int lo, up;
    if (f>=nf)
      errstop("boundary communication #%d cannot be defined for "
              "undefined field #%d", e-cfid, f-cfid);
    lo = fd[f].ext[OH_LOWER];  up = fd[f].ext[OH_UPPER];
    for (b=0; b<nb; b++,i++) {
      int sl=ct[i][OH_LOWER][OH_CTYPE_SIZE];
      int su=ct[i][OH_UPPER][OH_CTYPE_SIZE];
      int lo1=ct[i][OH_LOWER][OH_CTYPE_FROM];
      int lo2=ct[i][OH_UPPER][OH_CTYPE_TO];
      int up1=ct[i][OH_LOWER][OH_CTYPE_TO]   + sl;
      int up2=ct[i][OH_UPPER][OH_CTYPE_FROM] + su;
      if (sl && lo1<lo)  lo = lo1;
      if (su && lo2<lo)  lo = lo2;
      if (sl && up1>up)  up = up1;
      if (su && up2>up)  up = up2;
    }
    fd[f].ext[OH_LOWER] = lo;  fd[f].ext[OH_UPPER] = up;
  }
  for (f=0; f<nf; f++) {
    int lo=fd[f].ext[OH_LOWER], up=fd[f].ext[OH_UPPER];
    for (d=0; d<OH_DIMENSION; d++) {
      fs[f][d][OH_LOWER] = lo;
      fs[f][d][OH_UPPER] = (Grid[d].size+cfid) + up;
      fd[f].size[d] = Grid[d].size + (up - lo);
    }
  }
  for (f=0; f<nf; f++) {
    int bl = ft[f][OH_FTYPE_BL];
    int rl = ft[f][OH_FTYPE_RL];
    fd[f].bc.base  = Field_Disp(f, bl, bl, bl);
    fd[f].red.base = Field_Disp(f, rl, rl, rl);
  }
  set_field_descriptors(ft, sd, 0);

  BorderExc = bx =
    (struct S_borderexc(*)[2][OH_DIMENSION][2])
    mem_alloc(sizeof(struct S_borderexc), ne*2*OH_DIMENSION*2, "BorderExc");

  for (e=0; e<ne; e++) {
    for (d=0; d<OH_DIMENSION; d++) {
      for (lu=0; lu<2; lu++)
        bx[e][1][d][lu].send.deriv = bx[e][1][d][lu].recv.deriv = 0;
    }
#ifdef OH_POS_AWARE
    set_border_exchange(e, 0, e<ne-1 ? MPI_DOUBLE : MPI_LONG_LONG_INT);
#else
    set_border_exchange(e, 0, MPI_DOUBLE);
#endif
  }
  clear_border_exchange();
}
void
set_field_descriptors(int (*ft)[OH_FTYPE_N], int sd[OH_DIMENSION][2], int ps) {

  int nf=nOfFields;
  struct S_flddesc *fd = FieldDesc;
  int size[3] = {0,0,0};
  int d, f;

  for (d=0; d<OH_DIMENSION; d++)  size[d] = sd[d][OH_UPPER] - sd[d][OH_LOWER];
  for (f=0; f<nf; f++) {
    int bu = ft[f][OH_FTYPE_BU] - 1;
    int ru = ft[f][OH_FTYPE_RU] - 1;
    int es = ft[f][OH_FTYPE_ES];
    fd[f].bc.size[ps] =
      Field_Disp(f, size[OH_DIM_X]+bu, size[OH_DIM_Y]+bu, size[OH_DIM_Z]+bu) -
      fd[f].bc.base + es;
    fd[f].red.size[ps] =
      Field_Disp(f, size[OH_DIM_X]+ru, size[OH_DIM_Y]+ru, size[OH_DIM_Z]+ru) -
      fd[f].red.base + es;
  }
}
static void
set_border_exchange(int e, int ps, MPI_Datatype type) {
  struct S_borderexc (*bx)[2] = BorderExc[e][ps];
  int f = BoundaryCommFields[e];
  int nb = nOfBoundaries;
  int (*bt)[2][OH_CTYPE_N] = &BoundaryCommTypes[e*nb];
  int (*bd)[2] = Boundaries[RegionId[ps]];
  int (*sd)[2] = SubDomains[RegionId[ps]];
  struct S_flddesc *fd = &FieldDesc[f];
  int esize = fd->esize;
  int fext = fd->ext[OH_UPPER] - fd->ext[OH_LOWER];
  int xyz[3] = {
    sd[OH_DIM_X][OH_UPPER]-sd[OH_DIM_X][OH_LOWER],
    OH_DIMENSION>OH_DIM_Y ? sd[OH_DIM_Y][OH_UPPER]-sd[OH_DIM_Y][OH_LOWER] : 0,
    OH_DIMENSION>OH_DIM_Z ? sd[OH_DIM_Z][OH_UPPER]-sd[OH_DIM_Z][OH_LOWER] : 0
  };
  int *wdh = fd->size;
  int exti[OH_DIMENSION][2], exto[OH_DIMENSION][2];
  int soff[OH_DIMENSION][2], roff[OH_DIMENSION][2];
  int ssize[OH_DIMENSION][2], rsize[OH_DIMENSION][2];
  int d, lu;

  for (d=0; d<OH_DIMENSION; d++) {
    int blo=bd[d][OH_LOWER], bup=bd[d][OH_UPPER];
    exti[d][OH_LOWER] = bt[blo][OH_LOWER][OH_CTYPE_FROM];
    exti[d][OH_UPPER] =
      bt[bup][OH_UPPER][OH_CTYPE_FROM] + bt[bup][OH_UPPER][OH_CTYPE_SIZE];
    exto[d][OH_LOWER] = bt[blo][OH_UPPER][OH_CTYPE_TO];
    exto[d][OH_UPPER] =
      bt[bup][OH_LOWER][OH_CTYPE_TO] + bt[bup][OH_LOWER][OH_CTYPE_SIZE];
    for (lu=OH_LOWER; lu<=OH_UPPER; lu++) {
      int sb=bd[d][lu], rb=bd[d][1-lu];
      soff[d][lu] = bt[sb][lu][OH_CTYPE_FROM];
      roff[d][lu] = bt[rb][lu][OH_CTYPE_TO];
      ssize[d][lu] = bt[sb][lu][OH_CTYPE_SIZE];
      rsize[d][lu] = bt[rb][lu][OH_CTYPE_SIZE];
    }
  }
  for (lu=OH_LOWER; lu<=OH_UPPER; lu++) {
    set_border_comm(esize, f, xyz, wdh, exti, exto, soff, ssize, lu, 0, type,
                    bx);
    set_border_comm(esize, f, xyz, wdh, exti, exto, roff, rsize, lu, 1, type,
                    bx);
  }
}
static void
set_border_comm(int esize, int f, int *xyz, int *wdh,
                int exti[OH_DIMENSION][2], int exto[OH_DIMENSION][2],
                int off[OH_DIMENSION][2], int size[OH_DIMENSION][2],
                int lu, int sr, MPI_Datatype basetype,
                struct S_borderexc bx[OH_DIMENSION][2]) {
  int bl[2]={1,1};
  MPI_Datatype tmptype[2]={MPI_DATATYPE_NULL,MPI_UB};
  int w=wdh[OH_DIM_X], wd=w*esize;
  int dp=OH_DIMENSION==1 ? 1 : wdh[OH_DIM_Y];
  MPI_Aint dispz[2]={0, wd*dp*sizeof(double)};
  struct S_bcomm *bcx, *bcy, *bcz;
  int xexto, yexti, yexto, zexti;
  int s;
  int lower = sr ? lu==OH_UPPER : lu==OH_LOWER;

  bcx = (sr==0) ? &bx[OH_DIM_X][lu].send : &bx[OH_DIM_X][lu].recv;
  bcx->deriv = 0;
  xexto = xyz[OH_DIM_X] + exto[OH_DIM_X][OH_UPPER] - exto[OH_DIM_X][OH_LOWER];
  if (OH_DIMENSION>OH_DIM_Y) {
    bcy = (sr==0) ? &bx[OH_DIM_Y][lu].send : &bx[OH_DIM_Y][lu].recv;
    bcy->deriv = 0;
    yexti = xyz[OH_DIM_Y] +
            exti[OH_DIM_Y][OH_UPPER] - exti[OH_DIM_Y][OH_LOWER];
    yexto = xyz[OH_DIM_Y] +
            exto[OH_DIM_Y][OH_UPPER] - exto[OH_DIM_Y][OH_LOWER];
  }
  if (OH_DIMENSION>OH_DIM_Z) {
    bcz = (sr==0) ? &bx[OH_DIM_Z][lu].send : &bx[OH_DIM_Z][lu].recv;
    bcz->deriv = 0;
    zexti = xyz[OH_DIM_Z] +
            exti[OH_DIM_Z][OH_UPPER] - exti[OH_DIM_Z][OH_LOWER];
  }
  if (OH_DIMENSION==1) {
    if ((s=size[OH_DIM_X][lu])==0) {
      bcx->buf = bcx->count = 0;  bcx->type = MPI_DATATYPE_NULL;
    } else {
      bcx->type = basetype;
      bcx->count = s * esize;
      bcx->buf =
        Field_Disp(f,
                   lower ? off[OH_DIM_X][lu] : xyz[OH_DIM_X]+off[OH_DIM_X][lu],
                   0, 0);
    }
  } else if (OH_DIMENSION==2) {
    if ((s=size[OH_DIM_X][lu])==0) {
      bcx->buf = bcx->count = 0;  bcx->type = MPI_DATATYPE_NULL;
    } else {
      MPI_Type_vector(yexti, s*esize, wd, basetype, &(bcx->type));
      MPI_Type_commit(&(bcx->type));  bcx->deriv = 1;
      bcx->count = 1;
      bcx->buf =
        Field_Disp(f,
                   lower ? off[OH_DIM_X][lu] : xyz[OH_DIM_X]+off[OH_DIM_X][lu],
                   exti[OH_DIM_Y][OH_LOWER], 0);
    }
    if ((s=size[OH_DIM_Y][lu])==0) {
      bcy->buf = bcy->count = 0;  bcy->type = MPI_DATATYPE_NULL;
    } else {
      if (xexto==w) {
        bcy->type = basetype;
        bcy->count = s * wd;
      } else {
        MPI_Type_vector(s, xexto*esize, wd, basetype, &(bcy->type));
        MPI_Type_commit(&(bcy->type));  bcy->deriv = 1;
        bcy->count = 1;
      }
      bcy->buf =
        Field_Disp(f, exto[OH_DIM_X][OH_LOWER],
                   lower ? off[OH_DIM_Y][lu] : xyz[OH_DIM_Y]+off[OH_DIM_Y][lu],
                   0);
    }
  } else {
    if ((s=size[OH_DIM_X][lu])==0) {
      bcx->buf = bcx->count = 0;  bcx->type = MPI_DATATYPE_NULL;
    } else {
      MPI_Type_vector(yexti, s*esize, wd, basetype, tmptype);
      MPI_Type_struct(2, bl, dispz, tmptype, &(bcx->type));
      MPI_Type_commit(&(bcx->type));  bcx->deriv = 1;
      bcx->count = zexti;
      bcx->buf =
        Field_Disp(f,
                   lower ? off[OH_DIM_X][lu] : xyz[OH_DIM_X]+off[OH_DIM_X][lu],
                   exti[OH_DIM_Y][OH_LOWER], exti[OH_DIM_Z][OH_LOWER]);
    }
    if ((s=size[OH_DIM_Y][lu])==0) {
      bcy->buf = bcy->count = 0;  bcy->type = MPI_DATATYPE_NULL;
    } else {
      if (xexto==w) {
        MPI_Type_vector(zexti, s*wd, wd*dp, basetype, &(bcy->type));
        bcy->count = 1;
      } else {
        MPI_Type_vector(s, xexto*esize, wd, basetype, tmptype);
        MPI_Type_struct(2, bl, dispz, tmptype, &(bcy->type));
        bcy->count = zexti;
      }
      MPI_Type_commit(&(bcy->type));  bcy->deriv = 1;
      bcy->buf =
        Field_Disp(f, exto[OH_DIM_X][OH_LOWER],
                   lower ? off[OH_DIM_Y][lu] : xyz[OH_DIM_Y]+off[OH_DIM_Y][lu],
                   exti[OH_DIM_Z][OH_LOWER]);
    }
    if ((s=size[OH_DIM_Z][lu])==0) {
      bcz->buf = bcz->count = 0;  bcz->type = MPI_DATATYPE_NULL;
    } else {
      if (xexto==w && yexto==dp) {
        bcz->type = basetype;
        bcz->count = s * wd * dp;
      } else {
        if (xexto==w) {
          MPI_Type_vector(s, wd*yexto, wd*dp, basetype, &(bcz->type));
          bcz->count = 1;
        } else if (yexto==dp) {
          MPI_Type_vector(s*yexto, xexto*esize, wd, basetype, &(bcz->type));
          bcz->count = 1;
        } else {
          MPI_Type_vector(yexto, xexto*esize, wd, basetype, tmptype);
          MPI_Type_struct(2, bl, dispz, tmptype, &(bcz->type));
          bcz->count = s;
        }
        MPI_Type_commit(&(bcz->type));  bcz->deriv = 1;
      }
      bcz->buf =
        Field_Disp(f, exto[OH_DIM_X][OH_LOWER], exto[OH_DIM_Y][OH_LOWER],
                   lower ? off[OH_DIM_Z][lu] :
                           xyz[OH_DIM_Z]+off[OH_DIM_Z][lu]);
    }
  }
}
void
clear_border_exchange() {
  int ne=nOfExc, e, d, lu;
  struct S_borderexc (*bx)[2][OH_DIMENSION][2] = BorderExc;

  for (e=0; e<ne; e++) {
    for (d=0; d<OH_DIMENSION; d++) {
      for (lu=OH_LOWER; lu<=OH_UPPER; lu++) {
        if (bx[e][1][d][lu].send.deriv)
          MPI_Type_free(&bx[e][1][d][lu].send.type);
        if (bx[e][1][d][lu].recv.deriv)
          MPI_Type_free(&bx[e][1][d][lu].recv.type);
        bx[e][1][d][lu].send.buf =   bx[e][1][d][lu].recv.buf = 0;
        bx[e][1][d][lu].send.count = bx[e][1][d][lu].recv.count = -1;
        bx[e][1][d][lu].send.deriv = bx[e][1][d][lu].recv.deriv = 0;
        bx[e][1][d][lu].send.type =  bx[e][1][d][lu].recv.type =
          MPI_DATATYPE_NULL;
      }
    }
  }
}
void
oh3_grid_size_(double size[OH_DIMENSION]) {
  oh3_grid_size(size);
}
void
oh3_grid_size(double size[OH_DIMENSION]) {
  int d, n, nn=nOfNodes;
  for (d=0; d<OH_DIMENSION; d++) {
    double s = (Grid[d].gsize = size[d]);
    Grid[d].rgsize = 1 / s;
    for (n=0; n<nn; n++) {
      SubDomainsFloat[n][d][OH_LOWER] *= s;
      SubDomainsFloat[n][d][OH_UPPER] *= s;
    }
    Grid[d].fcoord[OH_LOWER] *= s;
    Grid[d].fcoord[OH_UPPER] *= s;
    Grid[d].fsize *= s;
    Grid[d].light.rfsize /= s;
    Grid[d].light.rfsizeplus /= s;
    Grid[d].light.fthresh *= s;
    if (SubDomainDesc) {
      for (n=0; n<nn; n++) {
        SubDomainDesc[n].coord[d].fc[OH_LOWER] *= s;
        SubDomainDesc[n].coord[d].fc[OH_UPPER] *= s;
      }
    }
  }
}
int
oh3_transbound_(int *currmode, int *stats) {
  return(transbound3(*currmode, *stats, 3));
}
int
oh3_transbound(int currmode, int stats) {
  return(transbound3(currmode, stats, 3));
}
static int
transbound3(int currmode, int stats, int level) {
  int oldp=RegionId[1], newp;

  currmode = excludeLevel2 ? transbound1(currmode, stats, 1) :
                             transbound2(currmode, stats, level);
  newp = RegionId[1];
  if (oldp!=newp) {
    if (oldp>=0)  clear_border_exchange();
    if (newp>=0)  set_field_descriptors(FieldTypes, SubDomains[newp], 1);
  }
  return(currmode);
}
#define Map_Particle_To_Neighbor(XYZ,RID,DIM,N,INC) {\
  double xyz=*XYZ;\
  if (xyz<SubDomainsFloat[RID][DIM][OH_LOWER]) {\
    N -= INC;\
    if (xyz<Grid[DIM].fcoord[OH_LOWER]) {\
      if (Boundaries[RID][DIM][OH_LOWER])  return(-1);\
      *XYZ += Grid[DIM].fcoord[OH_UPPER] - Grid[DIM].fcoord[OH_LOWER];\
    }\
  } else if (xyz>=SubDomainsFloat[RID][DIM][OH_UPPER]) {\
    N += INC;\
    if (xyz>=Grid[DIM].fcoord[OH_UPPER]) {\
      if (Boundaries[RID][DIM][OH_UPPER])  return(-1);\
      *XYZ -= Grid[DIM].fcoord[OH_UPPER] - Grid[DIM].fcoord[OH_LOWER];\
    }\
  }\
}
#define Neighbor_Id(N) ((n=(N))<0 ? ((n=-n-1)<nOfNodes ? n : -1) : n)
#if OH_DIMENSION==1
int
oh3_map_region_to_adjacent_node_(double *x, int *ps) {
  return(oh3_map_particle_to_neighbor(x, *ps));
}
int
oh3_map_particle_to_neighbor_(double *x, int *ps) {
  return(oh3_map_particle_to_neighbor(x, *ps));
}
int
oh3_map_particle_to_neighbor(double *x, int ps) {
  int rid=RegionId[ps], n=OH_NEIGHBORS>>1;

  Map_Particle_To_Neighbor(x, rid, OH_DIM_X, n, 1);
  return(Neighbor_Id(Neighbors[ps][n]));
}
#elif OH_DIMENSION==2
int
oh3_map_region_to_adjacent_node_(double *x, double *y, int *ps) {
  return(oh3_map_particle_to_neighbor(x, y, *ps));
}
int
oh3_map_particle_to_neighbor_(double *x, double *y, int *ps) {
  return(oh3_map_particle_to_neighbor(x, y, *ps));
}
int
oh3_map_particle_to_neighbor(double *x, double *y, int ps) {
  int rid=RegionId[ps], n=OH_NEIGHBORS>>1;

  Map_Particle_To_Neighbor(x, rid, OH_DIM_X, n, 1);
  Map_Particle_To_Neighbor(y, rid, OH_DIM_Y, n, 3);
  return(Neighbor_Id(Neighbors[ps][n]));
}
#else
int
oh3_map_region_to_adjacent_node_(double *x, double *y, double *z, int *ps) {
  return(oh3_map_particle_to_neighbor(x, y, z, *ps));
}
int
oh3_map_particle_to_neighbor_(double *x, double *y, double *z, int *ps) {
  return(oh3_map_particle_to_neighbor(x, y, z, *ps));
}
int
oh3_map_particle_to_neighbor(double *x, double *y, double *z, int ps) {
  int rid=RegionId[ps], n=OH_NEIGHBORS>>1;

  Map_Particle_To_Neighbor(x, rid, OH_DIM_X, n, 1);
  Map_Particle_To_Neighbor(y, rid, OH_DIM_Y, n, 3);
  Map_Particle_To_Neighbor(z, rid, OH_DIM_Z, n, 9);
  return(Neighbor_Id(Neighbors[ps][n]));
}
#endif
#define Map_Particle_To_Subdomain(XYZ,DIM,SDOM) {\
  double thresh = Grid[DIM].light.fthresh;\
  if (XYZ<Grid[DIM].fcoord[OH_LOWER] || XYZ>=Grid[DIM].fcoord[OH_UPPER])\
    return(-1);\
  if (XYZ<thresh)\
    SDOM = (XYZ - Grid[DIM].fcoord[OH_LOWER]) * Grid[DIM].light.rfsize;\
  else  SDOM = (int)((XYZ - thresh) * Grid[DIM].light.rfsizeplus) + \
               Grid[DIM].light.n;\
}
#define Adjust_Subdomain(XYZ,DIM,SDOM,INC) {\
  if (XYZ<SubDomainsFloat[SDOM][DIM][OH_LOWER])  SDOM-=INC;\
  else if (XYZ>=SubDomainsFloat[SDOM][DIM][OH_UPPER])  SDOM+=INC;\
}
#if OH_DIMENSION==1
int
oh3_map_region_to_node_(double *x) {
  return(oh3_map_particle_to_subdomain(*x));
}
int
oh3_map_particle_to_subdomain_(double *x) {
  return(oh3_map_particle_to_subdomain(*x));
}
int
oh3_map_particle_to_subdomain(double x) {
  int sdx;

  if (SubDomainDesc)  return(map_irregular_subdomain(x, 0.0, 0.0));
  Map_Particle_To_Subdomain(x, OH_DIM_X, sdx);
  Adjust_Subdomain(x, OH_DIM_X, sdx, 1);
  return(sdx);
}
#elif OH_DIMENSION==2
int
oh3_map_region_to_node_(double *x, double *y) {
  return(oh3_map_particle_to_subdomain(*x, *y));
}
int
oh3_map_particle_to_subdomain_(double *x, double *y) {
  return(oh3_map_particle_to_subdomain(*x, *y));
}
int
oh3_map_particle_to_subdomain(double x, double y) {
  int sdx, sdy, sd, nx=Grid[OH_DIM_X].n;

  if (SubDomainDesc)  return(map_irregular_subdomain(x, y, 0.0));
  Map_Particle_To_Subdomain(x, OH_DIM_X, sdx);
  Map_Particle_To_Subdomain(y, OH_DIM_Y, sdy);
  sd = sdx + nx * sdy;
  Adjust_Subdomain(x, OH_DIM_X, sd, 1);
  Adjust_Subdomain(y, OH_DIM_Y, sd, nx);
  return(sd);
}
#else
int
oh3_map_region_to_node_(double *x, double *y, double *z) {
  return(oh3_map_particle_to_subdomain(*x, *y, *z));
}
int
oh3_map_particle_to_subdomain_(double *x, double *y, double *z) {
  return(oh3_map_particle_to_subdomain(*x, *y, *z));
}
int
oh3_map_particle_to_subdomain(double x, double y, double z) {
  int sdx, sdy, sdz, sd, nx=Grid[OH_DIM_X].n, nxy=nx*Grid[OH_DIM_Y].n;

  if (SubDomainDesc)  return(map_irregular_subdomain(x, y, z));
  Map_Particle_To_Subdomain(x, OH_DIM_X, sdx);
  Map_Particle_To_Subdomain(y, OH_DIM_Y, sdy);
  Map_Particle_To_Subdomain(z, OH_DIM_Z, sdz);
  sd = sdx + nx * sdy + nxy * sdz;
  Adjust_Subdomain(x, OH_DIM_X, sd, 1);
  Adjust_Subdomain(y, OH_DIM_Y, sd, nx);
  Adjust_Subdomain(z, OH_DIM_Z, sd, nxy);
  return(sd);
}
#endif
int
map_irregular_subdomain(double x, double y, double z) {
  return(map_irregular(x, y, z, OH_DIM_X, 0, nOfNodes));
}
static int
map_irregular(double p0, double p1, double p2, int dim, int from, int n) {
  double size=Grid[dim].fsize;
  int to=from+n, lo, up, i;
  struct S_subdomdesc *sd = SubDomainDesc;

  lo = map_irregular_range(p0*2.0-size, dim, from, to);
  up = map_irregular_range(p0*2.0+size, dim, lo, to);
  for (i=lo; i<up; ) {
    int n = sd[i].coord[dim].n;
    if (p0>=sd[i].coord[dim].fc[OH_LOWER] &&
        p0< sd[i].coord[dim].fc[OH_UPPER]) {
      if (dim<OH_DIMENSION-1) {
        int ret = map_irregular(p1, p2, 0.0, dim+1, i, n);
        if (ret>=0)  return(ret);
      }
      else
        return(sd[i].id);
    }
    i += n;
  }
  return(-1);
}
static int
map_irregular_range(double p, int dim, int from, int to) {
  struct S_subdomdesc *sd = SubDomainDesc;
  int i;

  if (from==to) return(to);
  if (p<sd[from].coord[dim].fc[OH_LOWER]+sd[from].coord[dim].fc[OH_UPPER])
    return(from);
  if (p>=sd[to-1].coord[dim].fc[OH_LOWER]+sd[to-1].coord[dim].fc[OH_UPPER])
    return(to);
  for (i=(from+to)>>1; from<i; i=(from+to)>>1) {
    if (p<sd[i].coord[dim].fc[OH_LOWER]+sd[i].coord[dim].fc[OH_UPPER])
      to = i;
    else
      from = i;
  }
  return(to);
}
void
oh3_bcast_field_(void *pfld, void *sfld, int *ftype) {
  int base=FieldDesc[*ftype-1].bc.base;
  int *size=FieldDesc[*ftype-1].bc.size;

  oh1_broadcast((double*)pfld+base, (double*)sfld+base, size[0], size[1],
                MPI_DOUBLE, MPI_DOUBLE);
}
void
oh3_bcast_field(void *pfld, void *sfld, int ftype) {
  int base=FieldDesc[ftype].bc.base;
  int *size=FieldDesc[ftype].bc.size;

  oh1_broadcast((double*)pfld+base, (double*)sfld+base, size[0], size[1],
                MPI_DOUBLE, MPI_DOUBLE);
}
void
oh3_reduce_field_(void *pfld, void *sfld, int *ftype) {
  int base=FieldDesc[*ftype-1].red.base;
  int *size=FieldDesc[*ftype-1].red.size;

  oh1_reduce((double*)pfld+base, (double*)sfld+base, size[0], size[1],
             MPI_DOUBLE, MPI_DOUBLE, MPI_SUM, MPI_SUM);
}
void
oh3_reduce_field(void *pfld, void *sfld, int ftype) {
  int base=FieldDesc[ftype].red.base;
  int *size=FieldDesc[ftype].red.size;

  oh1_reduce((double*)pfld+base, (double*)sfld+base, size[0], size[1],
             MPI_DOUBLE, MPI_DOUBLE, MPI_SUM, MPI_SUM);
}
void
oh3_allreduce_field_(void *pfld, void *sfld, int *ftype) {
  int base=FieldDesc[*ftype-1].red.base;
  int *size=FieldDesc[*ftype-1].red.size;

  oh1_all_reduce((double*)pfld+base, (double*)sfld+base, size[0], size[1],
                 MPI_DOUBLE, MPI_DOUBLE, MPI_SUM, MPI_SUM);
}
void
oh3_allreduce_field(void *pfld, void *sfld, int ftype) {
  int base=FieldDesc[ftype].red.base;
  int *size=FieldDesc[ftype].red.size;

  oh1_all_reduce((double*)pfld+base, (double*)sfld+base, size[0], size[1],
                 MPI_DOUBLE, MPI_DOUBLE, MPI_SUM, MPI_SUM);
}
void
oh3_exchange_borders_(void *pfld, void *sfld, int *ctype, int *bcast) {
  oh3_exchange_borders(pfld, sfld, *ctype-1, *bcast);
}
void
oh3_exchange_borders(void *pfld, void *sfld, int ctype, int bcast) {
  MPI_Status st;
  int d, lu;
  double *pf=(double*)pfld, *sf=(double*)sfld;

  for (d=0; d<OH_DIMENSION; d++) {
    for (lu=OH_LOWER; lu<=OH_UPPER; lu++) {
      int dst=Adjacent[d][lu], src=Adjacent[d][1-lu];
      struct S_borderexc *bx=&BorderExc[ctype][0][d][lu];
      int scount=bx->send.count;
      int rcount=bx->recv.count;
      if (scount && rcount)
        MPI_Sendrecv(pf+bx->send.buf, scount, bx->send.type, dst, 0,
                     pf+bx->recv.buf, rcount, bx->recv.type, src, 0,
                     MCW, &st);
      else if (scount)
        MPI_Send(pf+bx->send.buf, scount, bx->send.type, dst, 0, MCW);
      else if (rcount)
        MPI_Recv(pf+bx->recv.buf, rcount, bx->recv.type, src, 0, MCW, &st);
    }
  }
  if (Mode_PS(currMode) && bcast) {
    if (RegionId[1]>=0 &&
        BorderExc[ctype][1][OH_DIM_X][OH_LOWER].send.count<0)
      set_border_exchange(ctype, 1, MPI_DOUBLE);
    for (d=0; d<OH_DIMENSION; d++) {
      for (lu=OH_LOWER; lu<=OH_UPPER; lu++) {
        struct S_borderexc *bxp=&BorderExc[ctype][0][d][lu];
        struct S_borderexc *bxs=&BorderExc[ctype][1][d][lu];
        oh1_broadcast(pf+bxp->recv.buf, sf+bxs->recv.buf,
                      bxp->recv.count, bxs->recv.count,
                      bxp->recv.type, bxs->recv.type);
      }
    }
  }
}
