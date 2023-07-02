/* File: ohhelp4s.c
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
#include "ohhelp3.h"
#undef  EXTERN
#define EXTERN
#include "ohhelp4s.h"

static void init4s(int **sdid, const int nspec, const int maxfrac,
                   const dint npmax, const int minmargin, const int maxdensity,
                   int **totalp, int **pbase, int *maxlocalp, int *cbufsize,
                   struct S_mycommc *mycommc, struct S_mycommf *mycommf,
                   int **nbor, int *pcoord, int **sdoms, int *scoord,
                   const int nbound, int *bcond, int **bounds, int *ftypes,
                   int *cfields, const int cfid, int *ctypes, int **fsizes,
                   int **zbound,
                   const int stats, const int repiter, const int verbose);
static int  transbound4s(int currmode, int stats, const int level);
static int  try_primary4s(const int currmode, const int level,
                          const int stats);
static int  try_stable4s(const int currmode, const int level, const int stats);
static void rebalance4s(const int currmode, const int level, const int stats);
static void exchange_particles4s(int currmode, const int nextmode,
                                 const int level, int reb, int oldp, int newp,
                                 const int stats);
static void count_population(const int nextmode, const int psnew,
                             const int stats);
static void exchange_population(const int currmode);
static void reduce_population();
static void add_population(dint *npd, const int xl, const int xu,
                           const int yl, const int yu, const int zl,
                           const int zu, const int src);
static void make_recv_list(const int currmode, const int level, const int reb,
                           const int oldp, const int newp, const int stats);
static void sched_recv(const int reb, const int get, const int stay,
                       const int nid, const int tag,
                       struct S_recvsched_context *context);
static void make_send_sched(const int reb, const int pcode, const int oldp,
                            const int newp, struct S_commlist *rlist[2],
                            int *rlidx[2], int *nacc, int *nsendptr);
static int  make_send_sched_body(const int ps, const int n, const int sdid,
                                 struct S_commlist *rlist);
static void make_send_sched_self(const int psor2, struct S_commlist *rlist,
                                 int *naccptr);
static void make_send_sched_hplane(const int psor2, const int z, int *naccptr,
                                   int *np, int *buf);
static void update_descriptors(const int oldp, const int newp);
static void update_neighbors(const int ps);
static void set_grid_descriptor(const int idx, const int nid);
static void adjust_field_descriptor(const int ps);
static void update_real_neighbors(const int mode, const int dosec,
                                  const int oldp, const int newp);
static void upd_real_nbr(const int root, const int psp, const int pss,
                         const int nbr, const int dosec, struct S_node *node,
                         struct S_realneighbor rnbrptr[2], int *occur[2]);
static void exchange_xfer_amount(const int trans, const int psnew,
                                 const int nextmode);
static void make_bxfer_sched(const int trans, const int psnew,
                             struct S_commlist *rlist[2], int *rlidx[2]);
static void make_bsend_sched(const int psor2, const int n, const int nx,
                             const int ny, struct S_commlist *rlist,
                             int *nsendptr, int *vpptr);
static void make_brecv_sched(const int psor2, const int n, const int nx,
                             const int ny, struct S_commlist *rlist,
                             int *nrecvptr, int vpidx);
static void move_to_sendbuf_4s(const int nextmode, const int psold,
                               const int psnew, const int trans,
                               const int oldp, const int *nacc,
                               const int nsend, const int stats);
static void move_to_sendbuf_uw4s(const int ps, const int mysd, const int cbase,
                                 const int nbase);
static void move_to_sendbuf_dw4s(const int ps, const int mysd, const int ctail,
                                 const int ntail);
static void sort_particles(const int nextmode, const int psnew,
                           const int stats);
static void move_and_sort(const int nextmode, const int psold, const int psnew,
                          const int oldp, const int *nacc, const int stats);
static void sort_received_particles(const int nextmode, const int psnew,
                                    const int stats);
static void set_sendbuf_disps4s(const int nextmode, const int trans);
static void xfer_particles(const int trans, const int psnew,
                           const int nextmode, struct S_particle *sbuf);
static void xfer_boundary_particles_v(const int psnew, const int pcode,
                                      const int d);
static void xfer_boundary_particles_h(const int psnew);
static void exchange_border_data_v(void *buf, void *sbuf, void *rbuf,
                                   MPI_Datatype type, const MPI_Aint esize,
                                   const int d);
static void exchange_border_data_h(void *buf, MPI_Datatype type,
                                   const MPI_Aint esize);

#define If_Dim(D, ET, EF)  (ET)
#define For_Y(LINIT, LCONT, LNEXT) for(LINIT; LCONT; LNEXT)
#define For_Z(LINIT, LCONT, LNEXT) for(LINIT; LCONT; LNEXT)
#define Do_Y(ACT) ACT
#define Do_Z(ACT) ACT
#define Coord_To_Index(GX, GY, GZ, W, DW)  ((GX) + (GY)*(W) + (GZ)*(DW))
#define Decl_For_All_Grid()\
  int fag_x1, fag_y1, fag_z1;\
  int fag_xidx, fag_yidx, fag_zidx, fag_gx, fag_gy, fag_gz;\
  int fag_w, fag_dw;
#define For_All_Grid(PS, X0, Y0, Z0, X1, Y1, Z1)\
  For_Z((fag_zidx=(Z0), fag_x1=GridDesc[PS].x+(X1),\
         fag_y1=GridDesc[PS].y+(Y1), fag_z1=GridDesc[PS].z+(Z1),\
         fag_w=GridDesc[PS].w, fag_dw=GridDesc[PS].dw,\
         fag_gz=Coord_To_Index(X0,Y0,Z0,fag_w,fag_dw)),\
        (fag_zidx<fag_z1), (fag_zidx++,fag_gz+=fag_dw))\
    For_Y((fag_yidx=(Y0), fag_gy=fag_gz),\
          (fag_yidx<fag_y1), (fag_yidx++,fag_gy+=fag_w))\
      for (fag_xidx=(X0),fag_gx=fag_gy; fag_xidx<fag_x1; fag_xidx++,fag_gx++)
#define For_All_Grid_Abs(PS, X0, Y0, Z0, X1, Y1, Z1)\
  For_Z((fag_zidx=(Z0), fag_x1=(X1), fag_y1=(Y1), fag_z1=(Z1),\
         fag_w=GridDesc[PS].w, fag_dw=GridDesc[PS].dw,\
         fag_gz=Coord_To_Index(X0,Y0,Z0,fag_w,fag_dw)),\
        (fag_zidx<fag_z1), (fag_zidx++,fag_gz+=fag_dw))\
    For_Y((fag_yidx=(Y0), fag_gy=fag_gz),\
          (fag_yidx<fag_y1), (fag_yidx++,fag_gy+=fag_w))\
      for (fag_xidx=(X0),fag_gx=fag_gy; fag_xidx<fag_x1; fag_xidx++,fag_gx++)
#define The_Grid()  (fag_gx)
#define Grid_X()  (fag_xidx)
#define Grid_Y()  (fag_yidx)
#define Grid_Z()  (fag_zidx)

#define URN_PRI 0
#define URN_SEC 1
#define URN_TRN 2

void
oh4s_init_(int *sdid, const int *nspec, const int *maxfrac, const dint *npmax,
           const int *minmargin, const int *maxdensity, int *totalp,
           int *pbase, int *maxlocalp, int *cbufsize, struct S_mycommf *mycomm,
           int *nbor, int *pcoord, int *sdoms, int *scoord, const int *nbound,
           int *bcond, int *bounds, int *ftypes, int *cfields, int *ctypes,
           int *fsizes, int *zbound,
           const int *stats, const int *repiter, const int *verbose) {
  specBase = 1;
  init4s(&sdid, *nspec, *maxfrac, *npmax, *minmargin, *maxdensity, &totalp,
         &pbase, maxlocalp, cbufsize, NULL, mycomm, &nbor, pcoord, &sdoms,
         scoord, *nbound, bcond, &bounds, ftypes, cfields, -1, ctypes, &fsizes,
         &zbound,
         *stats, *repiter, *verbose);
}
void
oh4s_init(int **sdid, const int nspec, const int maxfrac, const dint npmax,
          const int minmargin, const int maxdensity, int **totalp,
          int **pbase, int *maxlocalp, int *cbufsize, void *mycomm,
          int **nbor, int *pcoord, int **sdoms, int *scoord, const int nbound,
          int *bcond, int **bounds, int *ftypes, int *cfields, int *ctypes,
          int **fsizes, int **zbound,
          const int stats, const int repiter, const int verbose) {
  specBase = 0;
  init4s(sdid, nspec, maxfrac, npmax, minmargin, maxdensity, totalp,
         pbase, maxlocalp, cbufsize, (struct S_mycommc*)mycomm, NULL, nbor,
         pcoord, sdoms, scoord, nbound, bcond, bounds, ftypes, cfields, 0,
         ctypes, fsizes, zbound,
         stats, repiter, verbose);
}
#define Allocate_NOfPGrid(BODY, NPG, TYPE, SIZE, MSG) {\
  const int ns2 = nOfSpecies<<1;\
  const int gridsize = SIZE;\
  TYPE *npg = BODY;\
  TYPE **npgp = (TYPE**)mem_alloc(sizeof(TYPE*), ns2, MSG);\
  int s, g, exto=OH_PGRID_EXT*3;\
  const int base = Coord_To_Index(exto, exto, exto,\
                                  GridDesc[0].w, GridDesc[0].dw);\
  if (!npg)\
    BODY = npg = (TYPE*)mem_alloc(sizeof(TYPE), ns2*gridsize, MSG) + base;\
  for (s=0; s<ns2; s++,npg+=gridsize) {\
    npgp[s] = npg;\
    for (g=0; g<gridsize; g++)  npg[g-base] = 0;\
  }\
  NPG[0] = npgp;  NPG[1] = npgp + nOfSpecies;\
}
static int nOfLocalPLimitShadow = -1;
static void
init4s(int **sdid, const int nspec, const int maxfrac, const dint npmax,
       const int minmargin, const int maxdensity, int **totalp, int **pbase,
       int *maxlocalp, int *cbufsize, struct S_mycommc *mycommc,
       struct S_mycommf *mycommf, int **nbor, int *pcoord, int **sdoms,
       int *scoord, const int nbound, int *bcond, int **bounds, int *ftypes,
       int *cfields, const int cfid, int *ctypes, int **fsizes, int **zbound,
       const int stats, const int repiter, const int verbose) {
  int nn, me, nnns, nnns2, n;
  int (*ft)[OH_FTYPE_N] = (int(*)[OH_FTYPE_N])ftypes;
  int *cf = cfields;
  int (*ct)[2][OH_CTYPE_N] = (int(*)[2][OH_CTYPE_N])ctypes;
  int nf, ne, c, b, size, ps, s, tr, i, x, y, z;
  int *nphgram = NULL;
  int *rnbr, *iptr;
  dint *npgdummy = NULL,  *npgtdummy = NULL;
  int loggrid;
  dint idmax;
  const int ext = OH_PGRID_EXT, ext2 = ext<<1, ext3 = ext*3;
  struct S_particle pbufdummy, *pbufdummyptr = &pbufdummy;
  dint npl;

  if (OH_DIMENSION!=3)
    errstop("dimension size %d is not 3 which level-4s extension requires.",
            OH_DIMENSION);
  if (OH_PGRID_EXT!=1)
    errsotp("boundary plane thickness %d is not 1 which level-4s extension "
            "requires.", OH_PGRID_EXT);

  MPI_Comm_size(MCW, &nn);  nnns = nn * nspec;  nnns2 = nnns << 1;
  TempArray = (int*)mem_alloc(sizeof(int), nn<<2, "TempArray");

  for (nf=0; ft[nf][OH_FTYPE_ES]>0; nf++);
  for (ne=0; cf[ne]+cfid>=0; ne++);
  FieldTypes = (int(*)[OH_FTYPE_N])
               mem_alloc(sizeof(int), (nf+2)*OH_FTYPE_N, "FieldTypes");
  BoundaryCommFields = cf =
    (int(*))mem_alloc(sizeof(int), ne+2, "BoundaryCommFields");
  BoundaryCommTypes = (int(*)[2][OH_CTYPE_N])
                      mem_alloc(sizeof(int), (ne+1)*nbound*2*OH_CTYPE_N,
                                "BoundaryCommTypes");
  memcpy(FieldTypes, ft, sizeof(int)*nf*OH_FTYPE_N);
  for (c=0; c<ne; c++)  cf[c] = cfields[c] + cfid;
  memcpy(BoundaryCommTypes, ct, sizeof(int)*ne*nbound*2*OH_CTYPE_N);
  ft = FieldTypes;  ct = BoundaryCommTypes + ne * nbound;
  ft[nf][OH_FTYPE_ES] = 1;
  ft[nf][OH_FTYPE_LO] = 0;  ft[nf][OH_FTYPE_UP] = 0;
  ft[nf][OH_FTYPE_BL] = -ext;  ft[nf][OH_FTYPE_BU] = ext;
  ft[nf][OH_FTYPE_RL] = -ext;  ft[nf][OH_FTYPE_RU] = ext;
  ft[nf+1][OH_FTYPE_ES] = -1;
  cf[ne] = nf;  cf[ne+1] = -1;
  ct[0][OH_LOWER][OH_CTYPE_FROM] = -ext;
  ct[0][OH_LOWER][OH_CTYPE_TO]   = ext;
  ct[0][OH_LOWER][OH_CTYPE_SIZE] = ext2;
  ct[0][OH_UPPER][OH_CTYPE_FROM] = -ext;
  ct[0][OH_UPPER][OH_CTYPE_TO]   = -ext3;
  ct[0][OH_UPPER][OH_CTYPE_SIZE] = ext2;
  for (b=1; b<nbound; b++)
    ct[b][OH_LOWER][OH_CTYPE_FROM] =
      ct[b][OH_LOWER][OH_CTYPE_TO]   =
      ct[b][OH_LOWER][OH_CTYPE_SIZE] =
      ct[b][OH_UPPER][OH_CTYPE_FROM] =
      ct[b][OH_UPPER][OH_CTYPE_TO]   =
      ct[b][OH_UPPER][OH_CTYPE_SIZE] = 0;

  init3(sdid, nspec, maxfrac, &nphgram, totalp, NULL, NULL, &pbufdummyptr,
        pbase, 0, mycommc, mycommf, nbor, pcoord, sdoms, scoord, nbound,
        bcond, bounds, (int*)ft, cf, cfid, (int*)BoundaryCommTypes, fsizes,
        stats, repiter, verbose, 0);

  size =
    ((Grid[OH_DIM_X].size+ext2)*(Grid[OH_DIM_Y].size+ext2)*
     (Grid[OH_DIM_Z].size+ext2)-
     Grid[OH_DIM_X].size*Grid[OH_DIM_Y].size*Grid[OH_DIM_Z].size) +
    Grid[OH_DIM_X].size*Grid[OH_DIM_Y].size;
  npl = (dint)oh2_max_local_particles(npmax, maxfrac, minmargin) +
    2 * maxdensity * size;
  if (npl>INT_MAX) mem_alloc_error("Particles", 0);
  nOfLocalPLimitShadow = *maxlocalp = npl;
  size =
    2 * maxdensity * Grid[OH_DIM_Z].size *
    ((Grid[OH_DIM_X].size+ext2)*(Grid[OH_DIM_Y].size+ext2)-
     Grid[OH_DIM_X].size*Grid[OH_DIM_Y].size);
  BoundarySendBuf =
    (struct S_particle*)mem_alloc(sizeof(struct S_particle), size,
                                  "BoundarySendBuf");
  size = Grid[OH_DIM_X].size + ext2;
  if (size<Grid[OH_DIM_Y].size)  size = Grid[OH_DIM_Y].size;
  *cbufsize = 2 * maxdensity * Grid[OH_DIM_Z].size * size;

  me = myRank;
  PbufIndex = NULL;
  set_grid_descriptor(0, me);
  size = GridDesc[0].dw * GridDesc[0].h;
  Allocate_NOfPGrid(npgdummy, NOfPGrid, dint, size, "NOfPGrid");
  Allocate_NOfPGrid(npgtdummy, NOfPGridTotal, dint, size, "NOfPGridTotal");
  NOfPGridZ = (dint*)mem_alloc(sizeof(dint), Grid[OH_DIM_Z].size,
                               "NOfPGridZ");

  size = Coord_To_Index(Grid[OH_DIM_X].size-1,
                        If_Dim(OH_DIM_Y, Grid[OH_DIM_Y].size-1, 0),
                        If_Dim(OH_DIM_Z, Grid[OH_DIM_Z].size-1, 0),
                        GridDesc[0].w, GridDesc[0].dw);
  for (loggrid=0; size; loggrid++,size>>=1);
  idmax = (dint)(((nn+OH_NEIGHBORS)<<1)-1)<<loggrid;
  if (idmax>INT_MAX && sizeof(OH_nid_t)==sizeof(int)) {
    const int ext6 = ext3<<1;
    errstop("local grid size (%d+%d)*(%d+%d)*(%d+%d) times number of nodes %d "
            "is too large for OH_nid_t=int and thus OH_BIG_SPACE should be "
            "defined.",
            GridDesc[0].w-ext6, ext6,
            GridDesc[0].d-ext6, ext6,
            GridDesc[0].h-ext6, ext6, nn);
  }
  logGrid = loggrid;  gridMask = (1 << loggrid) - 1;
  adjust_field_descriptor(0);

  iptr = (int*)mem_alloc(sizeof(int), 2*2*4*nspec, "HPlane");
  for (ps=0; ps<2; ps++)  for (i=OH_LOWER; i<=OH_UPPER; i++) {
    HPlane[ps][i].nsend = iptr;  iptr += nspec;
    HPlane[ps][i].nrecv = iptr;  iptr += nspec;
    HPlane[ps][i].sbuf  = iptr;  iptr += nspec;
    HPlane[ps][i].rbuf  = iptr;  iptr += nspec;
    HPlane[ps][i].nbor  = MPI_PROC_NULL;
  }
  size = 2*nn + 2*2 + 2;
  VPlane = (struct S_vplane*)mem_alloc(sizeof(struct S_vplane), size,
                                       "VPlane");
  VPlaneHead[0] = VPlaneHead[1] = VPlaneHead[2] = VPlaneHead[3] =
    VPlaneHead[4] = VPlaneHead[5] = VPlaneHead[6] = VPlaneHead[7] =
    VPlaneHead[8] = 0;

  iptr = *zbound;
  if (!iptr)  iptr = *zbound = mem_alloc(sizeof(int), 4, "ZBound");
  ZBoundShadow = (int(*)[2])iptr;
  ZBoundShadow[0][OH_LOWER] = 0;  ZBoundShadow[0][OH_UPPER] = GridDesc[0].z;
  ZBoundShadow[1][OH_UPPER] = ZBoundShadow[1][OH_UPPER] = 0;

  InteriorParts = mem_alloc(sizeof(struct S_interiorp), nspec*2,
                            "InteriorParts");

  MPI_Type_vector(nspec, 1, nn, MPI_INT, &T_Hgramhalf);
  MPI_Type_commit(&T_Hgramhalf);
  for (n=0; n<nnns2; n++)  NOfSend[n] = 0;

  for (z=0,n=0; z<3; z++) {
    int (*bd)[OH_DIMENSION][2] = Boundaries;
    const int nonpz = z!=1 && bd[me][OH_DIM_Z][z>>1];
    for (y=0; y<3; y++) {
      const int nonpy = nonpz || (y!=1 && bd[me][OH_DIM_Y][y>>1]);
      for (x=0; x<3; x++,n++) {
        int dnbr = DstNeighbors[n];
        const int nrev = OH_NEIGHBORS - 1 - n;
        if (nonpy || (x!=1 && bd[me][OH_DIM_X][x>>1]))
          DstNeighbors[n] = SrcNeighbors[nrev] = -(nn+1);
        else if (dnbr<0 && dnbr>=-nn)
          DstNeighbors[n] = SrcNeighbors[nrev] = -(dnbr+1);
        else
          SrcNeighbors[nrev] = dnbr;
      }
    }
  }
  for (i=0; i<nn; i++)  TempArray[i] = 0;
  for (n=0; n<OH_NEIGHBORS; n++) {
    const int dnbr = DstNeighbors[n],  snbr = SrcNeighbors[n];
    int *sfirst = TempArray + nn;
    if (dnbr>=0) {
      if (TempArray[dnbr]&1)  DstNeighbors[n] = -(dnbr+1);
      else                    TempArray[dnbr] |= 1;
    }
    if (snbr>=0) {
      if (TempArray[snbr]&2) {
        SrcNeighbors[n] = -(snbr+1);
        FirstNeighbor[n] = sfirst[snbr];
      } else {
        FirstNeighbor[n] = sfirst[snbr] = n;
        TempArray[snbr] |= 2;
      }
    } else
      FirstNeighbor[n] = n;
    PrimaryRLIndex[n] = n;
  }
  update_neighbors(0);
  rnbr = (int*)mem_alloc(sizeof(int), nn*2*2*2, "RealNeighbors");
  for (tr=0; tr<2; tr++)  for (ps=0; ps<2; ps++,rnbr+=nn) {
    RealDstNeighbors[tr][ps].n = RealSrcNeighbors[tr][ps].n = 0;
    RealDstNeighbors[tr][ps].nbor = rnbr;
    RealSrcNeighbors[tr][ps].nbor = rnbr + nn*2*2;
  }
  update_real_neighbors(URN_PRI, 0, -1, -1);

  if (!SubDomainDesc)
    memcpy(BoundaryCondition, bcond, sizeof(int)*OH_DIMENSION*2);
}
void
oh4s_particle_buffer_(const int *maxlocalp, struct S_particle *pbuf) {
  oh4s_particle_buffer(*maxlocalp, &pbuf);
}
void
oh4s_particle_buffer(const int maxlocalp, struct S_particle **pbuf) {

  if (nOfLocalPLimitShadow<0)
    errstop("oh4s_particle_buffer() has to be called after oh4s_init()");
  else if (maxlocalp<nOfLocalPLimitShadow)
    errstop("argument maxlocalp %d given to oh4s_particle_buffer() is less "
            "than that calculated by oh4s_init() %d",
            maxlocalp, nOfLocalPLimitShadow);
  if (*pbuf)
    Particles = *pbuf;
  else
    Particles = *pbuf =
      (struct S_particle*)mem_alloc(sizeof(struct S_particle),
                                    maxlocalp<<1, "Particles");
  SendBuf = Particles + maxlocalp;
  nOfLocalPLimit = totalParts = maxlocalp;
}
void
oh4s_per_grid_histogram_(int *pghgram, int *pgindex) {
  oh4s_per_grid_histogram(&pghgram, &pgindex);
}
void
oh4s_per_grid_histogram(int **pghgram, int **pgindex) {
  int *npgo=NULL, *npgi=NULL;
  const int size = GridDesc[0].dw*GridDesc[0].h;
  Allocate_NOfPGrid(npgo, NOfPGridOut, int, size, "NOfPGridOut");
  Allocate_NOfPGrid(*pghgram, NOfPGridOutShadow, int, size,
                    "NOfPGridOutShadow");
  Allocate_NOfPGrid(npgi, NOfPGridIndex, int, size, "NOfPGridIndex");
  Allocate_NOfPGrid(*pgindex, NOfPGridIndexShadow, int, size,
                    "NOfPGridIndexShadow");
}
int
oh4s_transbound_(int *currmode, int *stats) {
  return(transbound4s(*currmode, *stats, 4));
}
int
oh4s_transbound(int currmode, int stats) {
  return(transbound4s(currmode, stats, 4));
}
static int
transbound4s(int currmode, int stats, const int level) {
  int ret=MODE_NORM_SEC;
  const int nn=nOfNodes, ns=nOfSpecies, ns2=ns<<1, nnns2=nn*ns2;
  struct S_particle *tmp;
  int i, ps, s, tp;
  Decl_For_All_Grid();

  stats = stats && statsMode;
  currmode = transbound1(currmode, stats, level);

  ZBound[0][OH_LOWER] = ZBound[0][OH_UPPER] = 0;
  ZBound[1][OH_UPPER] = ZBound[1][OH_UPPER] = 0;
  if (try_primary4s(currmode, level, stats))  ret = MODE_NORM_PRI;
  else if (!Mode_PS(currmode) || !try_stable4s(currmode, level, stats)) {
    rebalance4s(currmode, level, stats);  ret = MODE_REB_SEC;
  }
  if (!PbufIndex)
    PbufIndex = (int*)mem_alloc(sizeof(int), ns2+1, "PbufIndex");
  for (i=0; i<nnns2; i++) NOfPLocal[i] = 0;
  for (s=0,tp=0; s<ns2; s++) {
    TotalP[s] = TotalPNext[s];  PbufIndex[s] = tp;  tp += TotalPNext[s];
  }
  PbufIndex[s] = totalParts = *totalLocalParticles = tp;  nOfInjections = 0;
  for (s=0; s<ns2; s++)  InjectedParticles[s] = 0;

  for (ps=0; ps<=Mode_PS(ret); ps++) {
    const int extio = (ps==1 && ret<0) ? OH_PGRID_EXT*3 : OH_PGRID_EXT;
    for (s=0; s<ns; s++) {
      dint *npg = NOfPGrid[ps][s];
      For_All_Grid(ps, -extio, -extio, -extio, extio, extio, extio)
        npg[The_Grid()] = 0;
    }
  }
  ZBoundShadow[0][0] = ZBound[0][0];    ZBoundShadow[0][1] = ZBound[0][1];
  ZBoundShadow[1][0] = ZBound[1][0];    ZBoundShadow[1][1] = ZBound[1][1];
  tmp = Particles;  Particles = SendBuf;  SendBuf = tmp;
  currMode = ret<0 ? -ret : ret;
  return(ret);
}
static int
try_primary4s(const int currmode, const int level, const int stats) {
  const int oldp = RegionId[1];

  if (!try_primary1(currmode, level, stats)) return(FALSE);
  exchange_particles4s(currmode, 0, level, 0, oldp, -1, stats);
  if (Mode_PS(currmode))  update_real_neighbors(URN_PRI, 0, -1, -1);
  return(TRUE);
}
static int
try_stable4s(const int currmode, const int level, const int stats) {
  if (!try_stable1(currmode, (Mode_Acc(currmode) ? level : -level), stats))
    return(FALSE);
  exchange_particles4s(currmode, 1, level, 0, RegionId[1], RegionId[1], stats);
  return(TRUE);
}
static void
rebalance4s(const int currmode, const int level, const int stats) {
  const int me=myRank, ns=nOfSpecies;
  const int oldp = RegionId[1],  amode = Mode_Acc(currmode);
  const int ninj = nOfInjections;
  int s, n, newp;

  rebalance1(currmode, (amode ? level : -level), stats);
  newp = amode ? Nodes[me].parentid : NodesNext[me].parentid;
  if (ninj && amode && oldp!=newp) {
    int *sinj = InjectedParticles + ns;
    const int sbase=specBase;
    int i;
    struct S_particle *p;
    Decl_Grid_Info();
    for (s=0; s<ns; s++)  sinj[s] = 0;
    if (newp>=0) {
      for (i=0,p=Particles+totalParts; i<ninj; i++,p++) {
        const OH_nid_t nid = p->nid;
        int sdid;
        if (Secondary_Injected(nid)) {
          Primarize_Id(p, sdid);  Secondarize_Id(p);
          if (sdid==newp)  sinj[Particle_Spec(p->spec-sbase)]++;
        }
      }
    }
  }
  exchange_particles4s(currmode, 1, level, 1, oldp, newp, stats);
  if (!amode) {
    set_grid_descriptor(1, newp);
    for (n=0; n<OH_NEIGHBORS; n++)  Neighbors[1][n] = Neighbors[2][n];
    update_neighbors(1);
  }
}
#define Parent_Old(PCODE)       ((PCODE) & 4)
#define Parent_New(PCODE)       ((PCODE) & 2)
#define Parent_New_Same(PCODE)  (((PCODE) & 3) == 3)
#define Parent_New_Diff(PCODE)  (((PCODE) & 3) == 2)
static void
exchange_particles4s(int currmode, const int nextmode, const int level,
                     int reb, int oldp, int newp, const int stats) {
  const int ns=nOfSpecies, exti=OH_PGRID_EXT;
  const int trans = !Mode_Acc(currmode) && reb ? 1 : 0;
  int pcode =
    (oldp>=0 ? 4 : 0) + (newp>=0 ? 2 : 0) + (oldp==newp ? 1 : 0);
  int ps, psold, psnew, s;
  int nacc[2], nsend, tp;
  struct S_commlist *rlist[2];
  int *rlidx[2];
  Decl_For_All_Grid();

  if (Mode_Acc(currmode)) {
    if (nextmode) {
      int i;
      const int nnns2 = nOfNodes * nOfSpecies * 2;
      if (reb) {
        exchange_particles(SecRList, SecRLSize, oldp, 0, currmode, stats);
        update_descriptors(oldp, newp);
        set_grid_descriptor(1, newp);
        update_neighbors(1);
        update_real_neighbors(URN_SEC, 0, -1, newp);
      }
      else
        exchange_particles(CommList+SLHeadTail[1], SecSLHeadTail[0], oldp, 0,
                           currmode, stats);
      for (i=0; i<nnns2; i++)  NOfSend[i] = 0;
    } else {
      move_to_sendbuf_primary(Mode_PS(currmode), stats);
      exchange_primary_particles(currmode, stats);
    }
    count_population(nextmode, (Parent_New(pcode) ? 1 : 0), 0);
    currmode = Mode_Set_Any(nextmode);
    reb = 0;  oldp = newp;  pcode = newp>=0 ? 7 : 0;
  }
  exchange_population(currmode);
  psold = Parent_Old(pcode) ? 1 : 0;
  psnew = Parent_New(pcode) ? 1 : 0;
  if (nextmode) {
    make_recv_list(currmode, level, reb, oldp, newp, stats);
    rlist[0] = CommList;  rlist[1] = SecRList;
    rlidx[0] = RLIndex;   rlidx[1] = SecRLIndex;
  } else {
    rlist[0] = PrimaryCommList[0];  rlist[1] = PrimaryCommList[1];
    rlidx[0] = rlidx[1] = PrimaryRLIndex;
  }
  make_send_sched(reb, pcode, oldp, newp, rlist, rlidx, nacc, &nsend);
  exchange_xfer_amount(trans, psnew, nextmode);

  for (ps=0,tp=0; ps<=psnew; ps++) {
    const int psor2 = ps ? trans + 1 : 0;
    const int sb = specBase;
    for (s=0; s<ns; s++) {
      dint *npgt=NOfPGridTotal[ps][s];
      int *npgo=NOfPGridOut[ps][s], *npgos=NOfPGridOutShadow[ps][s];
      int *npgi=NOfPGridIndex[ps][s], *npgis=NOfPGridIndexShadow[ps][s];
      For_All_Grid(psor2, -exti, -exti, -exti, exti, exti, exti) {
        const int g = The_Grid(),  np = npgo[g];
        npgos[g] = np;  npgt[g] = npgi[g] = tp;  npgis[g] = tp + sb;  tp += np;
      }
    }
  }
  if (trans || (dint)nacc[1]+(dint)nsend>(dint)nOfLocalPLimit) {
    move_to_sendbuf_4s(nextmode, psold, psnew, trans, oldp, nacc, nsend,
                       stats);
    xfer_particles(trans, psnew, nextmode, SendBuf);
    make_bxfer_sched(trans, psnew, rlist, rlidx);
    sort_particles(nextmode, psnew, stats);
  } else {
    make_bxfer_sched(0, psnew, rlist, rlidx);
    move_and_sort(nextmode, psold, psnew, oldp, nacc, stats);
    xfer_particles(trans, psnew, nextmode, SendBuf+nacc[1]);
    sort_received_particles(nextmode, psnew, stats);
  }
  xfer_boundary_particles_v(psnew, trans, 0);
  xfer_boundary_particles_v(psnew, trans, 1);
  xfer_boundary_particles_h(psnew);
}
static void
count_population(const int nextmode, const int psnew, const int stats) {
  int ps, s, t, i, j, tp;
  const int ns=nOfSpecies, exti=OH_PGRID_EXT;
  Decl_For_All_Grid();
  Decl_Grid_Info();

  if (stats) oh1_stats_time(STATS_TB_SORT, nextmode);
  for (ps=0,t=0,j=0,tp=0; ps<=psnew; ps++) {
    for (s=0; s<ns; s++,t++) {
      dint *npgs = NOfPGrid[ps][s];
      const int tpn = TotalP[t] = TotalPNext[t];
      tp += tpn;
      For_All_Grid(ps, -exti, -exti, -exti, exti, exti, exti)
        npgs[The_Grid()]=0;
      for (i=0; i<tpn; i++,j++) {
        const int g = Grid_Position(Particles[j].nid);
        npgs[g]++;
        Particles[j].nid = Combine_Subdom_Pos(OH_NBR_SELF, g);
      }
    }
    if (ps==0)  primaryParts = tp;
  }
  totalParts = tp;  nOfInjections = 0;
}
static void
exchange_population(const int currmode) {
  const int ns=nOfSpecies;
  int s, zz;
  dint **npg = NOfPGridTotal[0];
  const int ct=nOfExc-1;
  const int ext = OH_PGRID_EXT,  ext2 = ext<<1,  ext3 = ext*3;
  const int x = GridDesc[0].x,  y = GridDesc[0].y,  z = GridDesc[0].z;
  const int w = GridDesc[0].w,  dw = GridDesc[0].dw;
  Decl_For_All_Grid();

  if (Mode_PS(currmode))  reduce_population();
  else {
    for (s=0; s<ns; s++) {
      dint *npgs = NOfPGrid[0][s],  *npgt = npg[s];
      For_All_Grid(0, -ext, -ext, -ext, ext, ext, ext)
        npgt[The_Grid()] = npgs[The_Grid()];
    }
  }
  for (zz=0; zz<z; zz++)  NOfPGridZ[zz] = 0;
  for (s=0; s<ns; s++) {
    dint *npgt = npg[s];
    oh3_exchange_borders(npgt, NULL, ct, 0);
    add_population(npgt, -ext3, x+ext3, -ext3, y+ext3, -ext, ext,  -dw*ext2);
    add_population(npgt, -ext3, x+ext3, -ext3, y+ext3, z-ext, z+ext, dw*ext2);
    add_population(npgt, -ext3, x+ext3, -ext, ext, -ext, z+ext, -w*ext2);
    add_population(npgt, -ext3, x+ext3, y-ext, y+ext, -ext, z+ext,  w*ext2);
    add_population(npgt, -ext,  ext,   -ext, y+ext, -ext, z+ext, -ext2);
    add_population(npgt, x-ext, x+ext, -ext, y+ext, -ext, z+ext,  ext2);

    For_All_Grid(0, 0, 0, 0, 0, 0, 0)
      NOfPGridZ[Grid_Z()] += npgt[The_Grid()];
  }
}
static void
reduce_population() {
  const int ft=nOfFields-1;
  const int base = FieldDesc[ft].red.base;
  const int *size = FieldDesc[ft].red.size;

  if (MyComm->black) {
    if (MyComm->prime!=MPI_COMM_NULL)
      MPI_Reduce(NOfPGrid[0][0]+base, NOfPGridTotal[0][0]+base, size[0],
                 MPI_LONG_LONG_INT, MPI_SUM, MyComm->rank, MyComm->prime);
    if (MyComm->sec!=MPI_COMM_NULL)
      MPI_Reduce(NOfPGrid[1][0]+base, NOfPGridTotal[1][0]+base, size[1],
                 MPI_LONG_LONG_INT, MPI_SUM, MyComm->root, MyComm->sec);
  } else {
    if (MyComm->sec!=MPI_COMM_NULL)
      MPI_Reduce(NOfPGrid[1][0]+base, NOfPGridTotal[1][0]+base, size[1],
                 MPI_LONG_LONG_INT, MPI_SUM, MyComm->root, MyComm->sec);
    if (MyComm->prime!=MPI_COMM_NULL)
      MPI_Reduce(NOfPGrid[0][0]+base, NOfPGridTotal[0][0]+base, size[0],
                 MPI_LONG_LONG_INT, MPI_SUM, MyComm->rank, MyComm->prime);
  }
  if (MyComm->prime==MPI_COMM_NULL)
    memcpy(NOfPGridTotal[0][0]+base, NOfPGrid[0][0]+base,
           size[0]*sizeof(dint));
}
static void
add_population(dint *npd, const int xl, const int xu, const int yl,
               const int yu, const int zl, const int zu,
               const int src) {
  dint *nps=npd+src;
  Decl_For_All_Grid();

  For_All_Grid_Abs(0, xl, yl, zl, xu, yu, zu)
    npd[The_Grid()] += nps[The_Grid()];
}
static void
make_recv_list(const int currmode, const int level, const int reb,
               const int oldp, const int newp, const int stats) {
  const int me = myRank, ns=nOfSpecies, nn=nOfNodes, nnns=nn*ns;
  const int nn2 = nn<<1;
  struct S_node *nodes = reb ? NodesNext : Nodes;
  struct S_node *mynode = nodes + me;
  struct S_node *ch;
  struct S_recvsched_context
    context = {0, 0, 0, CommList};
  int rlsize, rlidx;
  const int ft=nOfFields-1;
  const int npgbase = FieldDesc[ft].bc.base;
  const int *npgsize = FieldDesc[ft].bc.size;
  const int zmax = GridDesc[0].z-1;
  struct S_commlist *lastrl;
  int i;

  for (ch=mynode->child; ch; ch=ch->sibling)
    sched_recv(reb, ch->get.sec, ch->stay.sec, ch->id, nnns, &context);
  sched_recv(0, mynode->get.prime, mynode->stay.prime, me, 0, &context);

  rlidx = rlsize = context.cptr - CommList;  lastrl = context.cptr - 1;
  if (rlsize==0) {
    struct S_commlist *rl = CommList;
    rl->rid = me;  rl->tag = 0;  rl->sid = 0;  rl->count = 0;
    rl->region = zmax;
    rlidx = rlsize = 1;
  } else
    lastrl->region = zmax;

  for (i=0; i<OH_NEIGHBORS; i++) {
    const int dst=DstNeighbors[i], src=SrcNeighbors[i];
    int rc;
    MPI_Status st;
    if (dst==me) {
      RLIndex[i] = 0;  continue;
    }
    if (src>=0) {
      RLIndex[i] = rlidx;
      if (dst>=0)
        MPI_Sendrecv(CommList, rlsize, T_Commlist, dst, 0,
                     CommList+rlidx, nn2, T_Commlist, src, 0, MCW, &st);
      else
        MPI_Recv(CommList+rlidx, nn2, T_Commlist, src, 0, MCW, &st);
      MPI_Get_count(&st, T_Commlist, &rc);  rlidx += rc;
    } else {
      if (dst>=0)
        MPI_Send(CommList, rlsize, T_Commlist, dst, 0, MCW);
      RLIndex[i] = (src<-nn) ? rlidx : RLIndex[FirstNeighbor[i]];
    }
  }
  RLIndex[OH_NEIGHBORS] = rlidx;  SecRLIndex[OH_NEIGHBORS] = 0;
  AltSecRList = SecRList = CommList + rlidx;  AltSecRLIndex[OH_NEIGHBORS] = 0;
  if (Mode_PS(currmode)) {
    oh1_broadcast(RLIndex, SecRLIndex, OH_NEIGHBORS+1, OH_NEIGHBORS+1,
                  MPI_INT, MPI_INT);
    oh1_broadcast(CommList, SecRList, rlidx,
                  SecRLIndex[OH_NEIGHBORS], T_Commlist, T_Commlist);
    AltSecRList += SecRLIndex[OH_NEIGHBORS];
  }
  if (reb) {
    build_new_comm(currmode, -level, 2, stats);
    update_descriptors(oldp, newp);
    set_grid_descriptor(2, newp);
    update_real_neighbors(URN_TRN, Mode_PS(currmode), oldp, newp);
    oh1_broadcast(RLIndex, AltSecRLIndex, OH_NEIGHBORS+1, OH_NEIGHBORS+1,
                  MPI_INT, MPI_INT);
    oh1_broadcast(CommList, AltSecRList, RLIndex[OH_NEIGHBORS],
                  AltSecRLIndex[OH_NEIGHBORS], T_Commlist, T_Commlist);
  }
  oh1_broadcast(NOfPGridTotal[0][0]+npgbase, NOfPGridTotal[1][0]+npgbase,
                npgsize[0], npgsize[1], MPI_LONG_LONG_INT, MPI_LONG_LONG_INT);
}
static void
sched_recv(const int reb, const int get, const int stay, const int nid,
           const int tag, struct S_recvsched_context *context) {
  const int z0=context->z;
  dint nptotal=context->nptotal;
  dint nplimit=context->nplimit;
  struct S_commlist *cptr=context->cptr;
  const int ns=nOfSpecies;
  int z;
  const int zz = GridDesc[0].z;

  if (reb)
    nplimit += get;
  else
    nplimit += get + stay;

  context->nplimit = nplimit;
  if (nptotal>=nplimit)  return;
  cptr->rid = nid;  cptr->tag = tag;  cptr->sid = 0;  cptr->count = 0;
  for (z=z0; z<zz; z++) {
    nptotal += NOfPGridZ[z];
    if (nptotal>=nplimit) {
      cptr->region = z;  context->z = z + 1;
      context->nptotal = nptotal;  context->cptr = cptr + 1;
      return;
    }
  }
  local_errstop("per-plane histogram total %d is less than the total particle "
                "population %d up to node %d",
                nptotal, nplimit, nid);
}
static void
make_send_sched(const int reb, const int pcode, const int oldp,
                const int newp, struct S_commlist *rlist[2],
                int *rlidx[2], int *nacc, int *nsendptr) {
  const int psold = Parent_Old(pcode) ? 1 : 0;
  const int psnew = Parent_New(pcode) ? 1 : 0;
  const int ns = nOfSpecies, ns2 = ns<<1,  nn = nOfNodes;
  const int tagt = OH_NBR_TCC * ns,  tagb = OH_NBR_BCC * ns;
  const int tag1 = OH_NEIGHBORS * ns;
  int s, ps, n;
  int nsend=0;

  for (s=0; s<ns2; s++)  TotalPNext[s] = 0;
  for (ps=0; ps<=psold; ps++) {
    const int root = ps ? oldp : myRank;
    for (n=0; n<OH_NEIGHBORS; n++) {
      const int nrev = OH_NEIGHBORS - 1 - n;
      int sdid = Neighbors[ps][n];
      if (sdid<0)  sdid = -(sdid+1);
      if (sdid<nn && (n==OH_NBR_SELF || sdid!=root))
        nsend += make_send_sched_body(ps, n, sdid, rlist[ps]+rlidx[ps][nrev]);
    }
  }
  nacc[0] = nacc[1] = 0;
  for (ps=0; ps<=psnew; ps++) {
    int psor2;
    int *nbors, *ri;
    struct S_commlist *rl;
    struct S_hplane *hp = HPlane[ps];
    if (ps && Parent_New_Diff(pcode)) {
      psor2 = 2;  nbors = Neighbors[2];  rl = AltSecRList;  ri = AltSecRLIndex;
    } else {
      psor2 = ps;  nbors = Neighbors[ps];  rl = rlist[ps];  ri = rlidx[ps];
    }
    make_send_sched_self(psor2, rl+ri[OH_NBR_SELF], nacc+ps);
    if (ZBound[ps][OH_UPPER]==0)  continue;
    if (hp[OH_LOWER].nbor==nn) {
      int sdid = nbors[OH_NBR_BCC];
      if (sdid<0)  sdid = -(sdid+1);
      if (sdid<nn) {
        const int zmax = (SubDomains[sdid][OH_DIM_Z][OH_UPPER] -
                          SubDomains[sdid][OH_DIM_Z][OH_LOWER]) - 1;
        struct S_commlist *rlb = rl + ri[OH_NEIGHBORS-1-OH_NBR_BCC];
        int rlz = rlb->region;
        while (rlz<zmax)  rlz = (++rlb)->region;
        hp[OH_LOWER].nbor = rlb->rid;
        hp[OH_LOWER].stag = (rlb->tag) ? tagb + tag1 : tagb;
      } else {
        hp[OH_LOWER].nbor = MPI_PROC_NULL;
        hp[OH_LOWER].stag = tagb;
      }
    }
    if (hp[OH_UPPER].nbor==nn) {
      int sdid = nbors[OH_NBR_TCC];
      struct S_commlist *rlt = rl + ri[OH_NEIGHBORS-1-OH_NBR_TCC];
      if (sdid<0)  sdid = -(sdid+1);
      if (sdid<nn) {
        hp[OH_UPPER].nbor = rlt->rid;
        hp[OH_UPPER].stag = (rlt->tag) ?  tagt + tag1 : tagt;
      } else {
        hp[OH_UPPER].nbor = MPI_PROC_NULL;
        hp[OH_UPPER].stag = tagt;
      }
    }
    if (!ps)  nacc[1] = nacc[0];
  }
  *nsendptr = nsend;
}
#define For_All_Grid_Z(PS, X0, Y0, Z0, X1, Y1, Z1)\
  For_Z((fag_zidx=(Z0), fag_x1=GridDesc[PS].x+(X1),\
         fag_y1=GridDesc[PS].y+(Y1), fag_z1=GridDesc[PS].z+(Z1),\
         fag_w=GridDesc[PS].w, fag_dw=GridDesc[PS].dw,\
         fag_gz=Coord_To_Index(X0,Y0,Z0,fag_w,fag_dw)),\
        (fag_zidx<fag_z1), (fag_zidx++,fag_gz+=fag_dw))
#define For_All_Grid_XY(PS, X0, Y0, X1, Y1)\
  For_Y((fag_yidx=(Y0), fag_gy=fag_gz),\
        (fag_yidx<fag_y1), (fag_yidx++,fag_gy+=fag_w))\
    for (fag_xidx=(X0),fag_gx=fag_gy; fag_xidx<fag_x1; fag_xidx++,fag_gx++)

#define Grid_Exterior_Boundary(N, GS, PL, PU) {\
  const int e = OH_PGRID_EXT;\
  if (N==0)      { PL = -e;    PU = -(GS); }\
  else if (N==1) { PL = 0;     PU = 0; }\
  else           { PL = (GS);  PU = e; }\
}
#define Grid_Interior_Boundary(N, GS, PL, PU) {\
  const int e = OH_PGRID_EXT;\
  if (N==0)      { PL = 0;       PU = -(GS)+e; }\
  else if (N==1) { PL = 0;       PU = 0; }\
  else           { PL = (GS)-e;  PU = 0; }\
}

#define Make_Send_Sched_Body(MYSELF) {\
  int s, nofsidx=nofsbase;\
  for (s=0; s<ns; s++,nofsidx+=nn) {\
    dint *npg = NOfPGrid[ps][s];\
    int nsendofs=0;\
    For_All_Grid_XY(ps, xl, yl, xu, yu) {\
      const int g = The_Grid();\
      if (MYSELF)  npg[g] = 0;\
      else {\
        nsendofs += npg[g];  npg[g] = nofsidx + 1;\
      }\
    }\
    nsend += nsendofs;  NOfSend[nofsidx] += nsendofs;\
  }\
}
static int
make_send_sched_body(const int ps, const int n, const int sdid,
                     struct S_commlist *rlist) {
  const int me=myRank, ns=nOfSpecies, nn=nOfNodes;
  const int nx = n % 3, ny = n/3 % 3, nz = n/9;
  int xl, xu, yl, yu, zl, zu, zn;
  int rlz = rlist->region, rid, ridp=-1, ridn=-1, nofsbase;
  int nsend = 0;
  const int zmax = (SubDomains[sdid][OH_DIM_Z][OH_UPPER] -
                    SubDomains[sdid][OH_DIM_Z][OH_LOWER]) - 1;
  Decl_For_All_Grid();

  Grid_Exterior_Boundary(nx, GridDesc[ps].x, xl, xu);
  Grid_Exterior_Boundary(ny, GridDesc[ps].y, yl, yu);
  Grid_Exterior_Boundary(nz, GridDesc[ps].z, zl, zu);
  zn = (nz==0) ? zmax + 1 - OH_PGRID_EXT : 0;
  while (rlz<zn)  rlz = (++rlist)->region;
  rid = rlist->rid;  nofsbase = rlist->tag + rid;

  For_All_Grid_Z(ps, xl, yl, zl, xu, yu, zu) {
    if (n==OH_NBR_SELF && rid==me) {
      Make_Send_Sched_Body(1);
    }
    else {
      Make_Send_Sched_Body(0);
    }
    if (++zn>rlz && zn<=zmax) {
      rlz = (++rlist)->region;  rid = rlist->rid;  nofsbase = rlist->tag + rid;
    }
  }
  return(nsend);
}
static void
make_send_sched_self(const int psor2, struct S_commlist *rlist, int *naccptr) {
  const int me=myRank, nn=nOfNodes, ns=nOfSpecies;
  const int tag1 = OH_NEIGHBORS * ns;
  const int tagt = OH_NBR_TCC * ns,  tagb = OH_NBR_BCC * ns;
  const int ps = psor2==0 ? 0 : 1,  rtag = ps ? tag1 : 0;
  const int exti = OH_PGRID_EXT;
  const int zmax = GridDesc[psor2].z - 1;
  int rlz = -1, rid = nn, ridp = -1, ridn, stag = 0;
  struct S_hplane *hp = HPlane[ps];
  int *zb = ZBound[ps];
  int np = *naccptr,  *tpn = TotalPNext + (ps ? ns : 0);
  int s;
  Decl_For_All_Grid();

  hp[OH_LOWER].nbor = hp[OH_UPPER].nbor = MPI_PROC_NULL;
  For_All_Grid_Z(psor2, -exti, -exti, -exti, exti, exti, exti) {
    const int z = Grid_Z();
    ridn = (z==rlz) ? (z<zmax ? rlist->rid : nn) : -1;
    if (ridp==me) {
      zb[OH_UPPER] = z;  hp[OH_UPPER].nbor = rid;
      hp[OH_UPPER].stag = stag + tagt;
      hp[OH_UPPER].rtag = rtag + tagb;
    } else if (ridn==me) {
      zb[OH_LOWER] = z + 1;  hp[OH_LOWER].nbor = rid;
      hp[OH_LOWER].stag = stag + tagb;
      hp[OH_LOWER].rtag = rtag + tagt;
    }
    if (rid==me) {
      if (ridp>=0) {
        make_send_sched_hplane(psor2, z, naccptr,
                               hp[OH_LOWER].nsend, hp[OH_LOWER].sbuf);
        if (ridn>=0) {
          for (s=0; s<ns; s++) {
            hp[OH_UPPER].nsend[s] = hp[OH_LOWER].nsend[s];
            hp[OH_UPPER].sbuf[s] = hp[OH_LOWER].sbuf[s];
          }
        }
      }
      else if (ridn>=0)
        make_send_sched_hplane(psor2, z, naccptr,
                               hp[OH_UPPER].nsend, hp[OH_UPPER].sbuf);
      else
        make_send_sched_hplane(psor2, z, naccptr, NULL, NULL);
    } else {
      if (ridp==me)
        make_send_sched_hplane(psor2, z, naccptr,
                               hp[OH_UPPER].nrecv, hp[OH_UPPER].rbuf);
      else if (ridn==me)
        make_send_sched_hplane(psor2, z, naccptr,
                               hp[OH_LOWER].nrecv, hp[OH_LOWER].rbuf);
      else {
        for (s=0; s<ns; s++) {
          int *npgo=NOfPGridOut[ps][s];
          For_All_Grid_XY(psor2, -exti, -exti, exti, exti)
            npgo[The_Grid()] = 0;
        }
      }
    }
    ridp = -1;
    if (z==rlz) {
      ridp = rid;
      if (z<zmax) {
        rlz = rlist->region;  rid = rlist->rid;
        stag = rlist->tag ? tag1 : 0;  rlist++;
      } else {
        rlz++;  rid = nn;  stag = 0;
      }
    }
  }
  for (s=0; s<ns; s++) {
    hp[OH_LOWER].sbuf[s] += np;  hp[OH_LOWER].rbuf[s] += np;
    hp[OH_UPPER].sbuf[s] += np;  hp[OH_UPPER].rbuf[s] += np;
    np += tpn[s];
  }
}
#define For_All_Grid_XY_At_Z(PS, X0, Y0, X1, Y1, Z0)\
  For_Z((fag_zidx=(Z0), fag_x1=GridDesc[PS].x+(X1),\
         fag_y1=GridDesc[PS].y+(Y1), fag_z1=(Z0)+1,\
         fag_w=GridDesc[PS].w, fag_dw=GridDesc[PS].dw,\
         fag_gz=Coord_To_Index(X0,Y0,Z0,fag_w,fag_dw)),\
        (fag_zidx<fag_z1), (fag_zidx++,fag_gz+=fag_dw))\
    For_Y((fag_yidx=(Y0), fag_gy=fag_gz),\
          (fag_yidx<fag_y1), (fag_yidx++,fag_gy+=fag_w))\
      for (fag_xidx=(X0),fag_gx=fag_gy; fag_xidx<fag_x1; fag_xidx++,fag_gx++)
static void
make_send_sched_hplane(const int psor2, const int z, int *naccptr,
                        int *np, int *buf) {
  const int ns=nOfSpecies, exti=OH_PGRID_EXT;
  const int ps = psor2==0 ? 0 : 1,  nsor0 = ps ? ns : 0;
  int nacc = *naccptr, s;
  Decl_For_All_Grid();

  for (s=0; s<ns; s++) {
    dint *npgt = NOfPGridTotal[ps][s];
    int *npgo = NOfPGridOut[ps][s];
    int npofs = 0;
    if (buf)  buf[s] = TotalPNext[nsor0+s];
    For_All_Grid_XY_At_Z(psor2, -exti, -exti, exti, exti, z) {
      const int g = The_Grid();
      npofs += (npgo[g] = npgt[g]);
    }
    nacc += npofs;  TotalPNext[nsor0+s] += npofs;
    if (np)  np[s] = npofs;
  }
  *naccptr = nacc;
}
static void
update_descriptors(const int oldp, const int newp) {
  int n;

  if (oldp!=newp) {
    if (oldp>=0)  clear_border_exchange();
    if (newp>=0) {
      set_field_descriptors(FieldTypes, SubDomains[newp], 1);
      adjust_field_descriptor(1);
    }
  }
}
#define Neighbor_Grid_Offset(PS, N, SD, D, XYZ)\
  (N==0 ? 0 : (N<0 ? SubDomains[SD][D][OH_LOWER]-SubDomains[SD][D][OH_UPPER] :\
                     GridDesc[ps].XYZ))
static void
update_neighbors(const int ps) {
  int n, nx, ny, nz;
  const int nn = nOfNodes;
  struct S_commlist *cl = PrimaryCommList[ps];

  for (nz=-1,n=0; nz<2; nz++) {
    for (ny=-1; ny<2; ny++) {
      for (nx=-1; nx<2; nx++,n++) {
        int nbr = Neighbors[ps][n];
        const int nrev = OH_NEIGHBORS - 1 - n;
        nbr = AbsNeighbors[ps][n] = nbr<0 ? -(nbr+1) : nbr;
        cl[nrev].rid = nbr;  cl[nrev].tag = cl[nrev].sid = cl[nrev].count = 0;
        if (nbr>=nn) {
          GridOffset[ps][n] = 0;  cl[nrev].region = 0;
        } else {
          GridOffset[ps][n] =
            Coord_To_Index(Neighbor_Grid_Offset(ps, nx, nbr, OH_DIM_X, x),
                           Neighbor_Grid_Offset(ps, ny, nbr, OH_DIM_Y, y),
                           Neighbor_Grid_Offset(ps, nz, nbr, OH_DIM_Z, z),
                           GridDesc[0].w, GridDesc[0].dw);
          cl[nrev].region = SubDomains[nbr][OH_DIM_Z][OH_UPPER] -
                            SubDomains[nbr][OH_DIM_Z][OH_LOWER] - 1;
        }
      }
    }
  }
}
static void
set_grid_descriptor(const int idx, const int nid) {
  const int exti6 = OH_PGRID_EXT*6;
  const int w = GridDesc[idx].w = Grid[OH_DIM_X].size+(exti6);
  const int d = GridDesc[idx].d =
                If_Dim(OH_DIM_Y, Grid[OH_DIM_Y].size+(exti6), 1);

  GridDesc[idx].h = If_Dim(OH_DIM_Z, Grid[OH_DIM_Z].size+(exti6), 1);
  GridDesc[idx].dw = d * w;
  if (nid>=0) {
    GridDesc[idx].x = SubDomains[nid][OH_DIM_X][OH_UPPER] -
                      SubDomains[nid][OH_DIM_X][OH_LOWER];
    GridDesc[idx].y = If_Dim(OH_DIM_Y,
                             SubDomains[nid][OH_DIM_Y][OH_UPPER] -
                             SubDomains[nid][OH_DIM_Y][OH_LOWER], 0);
    GridDesc[idx].z = If_Dim(OH_DIM_Z,
                             SubDomains[nid][OH_DIM_Z][OH_UPPER] -
                             SubDomains[nid][OH_DIM_Z][OH_LOWER], 0);
  } else {
    GridDesc[idx].x = GridDesc[idx].y = GridDesc[idx].z = -exti6;
    /* to ensure, e.g., x+3*(OH_PGRID_EXT)<=-3*(OH_PGRID_EXT) */
  }
}
static void
adjust_field_descriptor(const int ps) {
  const int f = nOfFields - 1,  ns = nOfSpecies;
  int d, fs;

  for (d=0,fs=1; d<OH_DIMENSION; d++)  fs *= FieldDesc[f].size[d];
  fs *= ns-1;
  FieldDesc[f].bc.size[ps] += fs;    FieldDesc[f].red.size[ps] += fs;
}
static void
update_real_neighbors(const int mode, const int dosec, const int oldp,
                      const int newp) {
  const int me=myRank, nn=nOfNodes, nn4=nn<<2;
  const int dosec0 = mode != URN_PRI;
  int i, nbridx, ps, *doccur[2], *soccur[2];

  for (i=0; i<nn4; i++)  TempArray[i] = 0;
  doccur[0] = TempArray;       doccur[1] = doccur[0] + nn;
  soccur[0] = doccur[1] + nn;  soccur[1] = soccur[0] + nn;

  if (mode==URN_TRN) {
    int *tmp = RealSrcNeighbors[1][0].nbor;
    RealSrcNeighbors[1][0].n = RealSrcNeighbors[0][0].n;
    RealSrcNeighbors[1][0].nbor = RealSrcNeighbors[0][0].nbor;
    RealSrcNeighbors[0][0].nbor = tmp;
  }
  RealDstNeighbors[0][0].n = RealDstNeighbors[0][1].n = 0;
  RealSrcNeighbors[0][0].n = RealSrcNeighbors[0][1].n = 0;
  upd_real_nbr(me, 0, 1, 0, dosec0, Nodes, RealDstNeighbors[0], doccur);
  upd_real_nbr(me, 0, 0, 0, dosec0, Nodes, RealSrcNeighbors[0], soccur);
  if (mode==URN_PRI)  return;

  nbridx = mode==URN_TRN ? 2 : 1;
  upd_real_nbr(newp, 0, 1, nbridx, 1, Nodes, RealDstNeighbors[0], doccur);
  upd_real_nbr(newp, 1, 1, nbridx, 1, Nodes, RealSrcNeighbors[0], soccur);
  if (mode!=URN_TRN)  return;

  for (ps=0; ps<2; ps++) {
    const int nd = RealDstNeighbors[0][ps].n;
    const int ns = RealSrcNeighbors[0][ps].n;
    for (i=0; i<nd; i++)  doccur[ps][RealDstNeighbors[0][ps].nbor[i]] = 0;
    for (i=0; i<ns; i++)  soccur[ps][RealSrcNeighbors[0][ps].nbor[i]] = 0;
  }
  RealDstNeighbors[1][0].n = RealDstNeighbors[1][1].n = 0;
  RealSrcNeighbors[1][1].n = 0;
  upd_real_nbr(me,   0, 1, 0, 1,     Nodes,     RealDstNeighbors[1], doccur);
  upd_real_nbr(oldp, 0, 1, 1, 1,     Nodes,     RealDstNeighbors[1], doccur);
  upd_real_nbr(newp, 1, 1, 2, dosec, NodesNext, RealSrcNeighbors[1], soccur);
}
static void
upd_real_nbr(const int root, const int psp, const int pss,
             const int nbr, const int dosec, struct S_node *nodes,
             struct S_realneighbor rnbrptr[2], int *occur[2]) {
  const int me=myRank;
  struct S_realneighbor *pnbr = rnbrptr+psp,  *snbr = rnbrptr+pss;
  int *poccur = occur[psp],  *soccur = occur[pss];
  int i;

  if (root<0)  return;
  if (root!=me && !poccur[root]) {
    pnbr->nbor[pnbr->n++] = root;  poccur[root] = 1;
  }
  if (dosec) {
    struct S_node *ch;
    for (ch=nodes[root].child; ch; ch=ch->sibling) {
      const int nid = ch->id;
      if (nid!=me && !soccur[nid]) {
        snbr->nbor[snbr->n++] = nid;  soccur[nid] = 1;
      }
    }
  }
  for (i=0; i<OH_NEIGHBORS; i++) {
    const int nid = Neighbors[nbr][i];
    struct S_node *ch;
    if (nid<0 || nid==root)  continue;
    if (!poccur[nid]) {
      pnbr->nbor[pnbr->n++] = nid;  poccur[nid] = 1;
    }
    if (dosec) {
      for (ch=nodes[nid].child; ch; ch=ch->sibling) {
        const int cid = ch->id;
        if (!soccur[cid]) {
          snbr->nbor[snbr->n++] = cid;  soccur[cid] = 1;
        }
      }
    }
  }
}
static void
exchange_xfer_amount(const int trans, const int psnew, const int nextmode) {
  const struct S_realneighbor *snbr = RealSrcNeighbors[trans];
  const struct S_realneighbor *dnbr = RealDstNeighbors[trans];
  const int nnns = nOfNodes * nOfSpecies;
  int ps, tag, req;

  for (ps=0,tag=0,req=0; ps<=psnew; ps++,tag+=nnns) {
    const int n = snbr[ps].n;
    const int *nbor = snbr[ps].nbor;
    int i,  *nrbase = NOfRecv + tag;
    for (i=0; i<n; i++,req++) {
      const int nid = nbor[i];
      MPI_Irecv(nrbase+nid, 1, T_Hgramhalf, nid, tag, MCW, Requests+req);
    }
  }
  for (ps=0,tag=0; ps<=nextmode; ps++,tag+=nnns) {
    const int n = dnbr[ps].n;
    const int *nbor = dnbr[ps].nbor;
    int i,  *nsbase = NOfSend + tag;
    for (i=0; i<n; i++,req++) {
      const int nid = nbor[i];
      MPI_Isend(nsbase+nid, 1, T_Hgramhalf, nid, tag, MCW, Requests+req);
    }
  }
  MPI_Waitall(req, Requests, Statuses);
}
static void
make_bxfer_sched(const int trans, const int psnew, struct S_commlist *rlist[2],
                 int *rlidx[2]) {
  const int nn = nOfNodes;
  int d, ps, du;
  int (*vph)[2][2] = (int(*)[2][2])VPlaneHead;
  int nsendib=0, nrecveb=0, vpidx=0;

  for (d=OH_DIM_X; d<=OH_DIM_Y; d++) {
    for (ps=0; ps<2; ps++) {
      const int psor2 = ps ? trans+1 : 0;
      struct S_commlist *rl;
      int *ri;
      if (ps>psnew || ZBound[ps][OH_UPPER]==0) {
        vph[d][ps][0] = vph[d][ps][1] = vpidx;
        continue;
      }
      if (psor2==2) {
        rl = AltSecRList;  ri = AltSecRLIndex;
      } else {
        rl = rlist[ps];  ri = rlidx[ps];
      }
      for (du=OH_LOWER; du<=OH_UPPER; du++) {
        const int vpisave = vpidx;
        const int nx = (d==OH_DIM_X) ? du<<1 : 1;
        const int ny = (d==OH_DIM_X) ? 1 : du<<1;
        const int n = 3*3 + 3*ny + nx;
        const int nrev = OH_NEIGHBORS - 1 - n;
        int nbor = Neighbors[psor2][n];
        vph[d][ps][du] = vpidx;
        if (nbor<0)  nbor = -(nbor+1);
        if (nbor<nn) {
          make_bsend_sched(psor2, n, nx, ny, rl+ri[nrev], &nsendib, &vpidx);
          make_brecv_sched(psor2, n, nx, ny, rl+ri[nrev], &nrecveb, vpisave);
        }
      }
    }
  }
  vph[2][0][0] = vpidx;
}
#define Add_Pillar_Voxel(I)     (((dint)(I))<<32)
#define Is_Pillar_Voxel(V)      (V>=((dint)1)<<32)
#define Pillar_Lower(V)         (V&INT_MAX)
#define Pillar_Upper(V)         ((V)>>32)

static void
make_bsend_sched(const int psor2, const int n, const int nx, const int ny,
                 struct S_commlist *rlist, int *nsendptr, int *vpptr) {
  const int ns=nOfSpecies;
  const int ps = psor2==0 ? 0 : 1;
  const int tag1 = OH_NEIGHBORS;
  const int stag = n;
  const int rtag = ((ps ? OH_NEIGHBORS : 0) + (OH_NEIGHBORS - 1 - n));
  struct S_commlist *rl = rlist;
  int rlz = rl->region;
  int nsend = *nsendptr, nsendsave = nsend,  vpidx = *vpptr;
  int xl, xu, yl, yu;
  const int zbl = ZBound[ps][OH_LOWER];
  const int zbu = ZBound[ps][OH_UPPER] - GridDesc[psor2].z;
  const int zmax = ZBound[ps][OH_UPPER] - 1;
  const int xtop = GridDesc[psor2].x;
  int s;
  Decl_For_All_Grid();

  if (ny==1) {
    Grid_Interior_Boundary(nx, GridDesc[psor2].x, xl, xu);
  } else {
    xl = -OH_PGRID_EXT;  xu = OH_PGRID_EXT;
  }
  Grid_Interior_Boundary(ny, GridDesc[psor2].y, yl, yu);

  while (rlz<zbl)  rlz = (++rl)->region;
  VPlane[vpidx].nbor = rl->rid;
  VPlane[vpidx].stag = (rl->tag ? tag1 : 0) + stag;
  VPlane[vpidx].rtag = rtag;
  VPlane[vpidx].sbuf = nsend;
  For_All_Grid_Z(psor2, xl, yl, zbl, xu, yu, zbu) {
    const int z = Grid_Z();
    for (s=0; s<ns; s++) {
      dint *npg = NOfPGrid[ps][s];
      int *npgo = NOfPGridOut[ps][s];
      For_All_Grid_XY(psor2, xl, yl, xu, yu) {
        const int g = The_Grid();
        const dint dst = npg[g];
        if (Grid_X()<0 || Grid_X()>=xtop)
          npg[g] += Add_Pillar_Voxel(nsend+1);
        else if (dst>=0)
          npg[g] = -(nsend+1);
        else
          npg[g] -= Add_Pillar_Voxel(nsend+1);
        nsend += npgo[g];
      }
    }
    if (z==rlz && z<zmax) {
      VPlane[vpidx++].nsend = nsend - nsendsave;
      rlz = (++rl)->region;
      VPlane[vpidx].nbor = rl->rid;
      VPlane[vpidx].stag = (rl->tag ? tag1 : 0) + stag;
      VPlane[vpidx].rtag = rtag;
      VPlane[vpidx].sbuf = nsendsave = nsend;
    }
  }
  VPlane[vpidx].nsend = nsend - nsendsave;
  *nsendptr = nsend;  *vpptr = vpidx + 1;
}
static void
make_brecv_sched(const int psor2, const int n, const int nx, const int ny,
                 struct S_commlist *rlist, int *nrecvptr, int vpidx) {
  const int ns=nOfSpecies;
  const int ps = psor2==0 ? 0 : 1;
  int nrecv = *nrecvptr,  nrecvsave = nrecv;
  struct S_commlist *rl = rlist;
  int rlz = rl->region;
  int xl, xu, yl, yu;
  const int zbl = ZBound[ps][OH_LOWER];
  const int zbu = ZBound[ps][OH_UPPER] - GridDesc[psor2].z;
  const int zmax = ZBound[ps][OH_UPPER] - 1;
  int s;
  Decl_For_All_Grid();

  if (ny==1) {
    Grid_Exterior_Boundary(nx, GridDesc[psor2].x, xl, xu);
  } else {
    xl = -OH_PGRID_EXT;  xu = OH_PGRID_EXT;
  }
  Grid_Exterior_Boundary(ny, GridDesc[psor2].y, yl, yu);

  while (rlz<zbl)  rlz = (++rl)->region;
  VPlane[vpidx].rbuf = nrecv;
  For_All_Grid_Z(psor2, xl, yl, zbl, xu, yu, zbu) {
    const int z = Grid_Z();
    for (s=0; s<ns; s++) {
      int *npgo = NOfPGridOut[ps][s];
      For_All_Grid_XY(psor2, xl, yl, zu, yu)
        nrecv += npgo[The_Grid()];
    }
    if (z==rlz && z<zmax) {
      VPlane[vpidx++].nrecv = nrecv - nrecvsave;
      rlz = (++rl)->region;
      VPlane[vpidx].rbuf = nrecvsave = nrecv;
    }
  }
  VPlane[vpidx].nrecv = nrecv - nrecvsave;  *nrecvptr = nrecv;
}
#define Local_Grid_Position(G, NID, PS)  ((G) + GridOffset[PS][NID>>loggrid])

#define Move_Or_Do(P, PS, MYSD, TOSB, ACT, PIL) {\
  const OH_nid_t nid = P->nid;\
  int g = Grid_Position(nid);\
  int sdid;\
  dint dst;\
  if (nid<0)  continue;\
  sdid = Neighbor_Subdomain_Id(nid, PS);\
  if (sdid!=(MYSD)) g = Local_Grid_Position(g, nid, PS);\
  dst = npg[g];\
  if (dst==0)  { ACT; }\
  else if (!PIL) {\
    if (TOSB)  sb[NOfSend[dst-1]++] = *P;\
  }\
  else if (dst>0)\
    sb[NOfSend[Pillar_Lower(dst)-1]++] = *P;\
  else {\
    ACT;\
    const dint bsbidx = -dst;\
    if (!Is_Pillar_Voxel(bsbidx)) {\
      BoundarySendBuf[bsbidx-1] = *P;  npg[g] = dst - 1;\
    }\
    else {\
      BoundarySendBuf[Pillar_Lower(bsbidx)-1] =\
        BoundarySendBuf[Pillar_Upper(bsbidx)-1] = *P;\
      npg[g] = dst - (Add_Pillar_Voxel(1) + 1);\
    }\
  }\
}
static void
move_to_sendbuf_4s(const int nextmode, const int psold, const int psnew,
                   const int trans, const int oldp, const int *nacc,
                   const int nsend, const int stats) {
  const int me=myRank, ns=nOfSpecies, nn=nOfNodes, sbase=specBase;
  const int ninj=nOfInjections, nplim=nOfLocalPLimit;
  int ps, s, t, i;
  int *nofr;
  int ninjp=0, ninjs=nplim;
  struct S_particle *sb = SendBuf,  *p;
  Decl_Grid_Info();

  if (stats) oh1_stats_time(STATS_TB_MOVE, nextmode);
  set_sendbuf_disps4s(nextmode, trans);

  for (ps=0,t=0,nofr=NOfRecv; ps<2; ps++) {
    const int nnbr = RealSrcNeighbors[trans][ps].n;
    const int *rnbr = RealSrcNeighbors[trans][ps].nbor;
    if (ps<=psnew) {
      for (s=0; s<ns; s++,t++,nofr+=nn) {
        int n, nrec;
        for (n=0,nrec=0; n<nnbr; n++)  nrec += nofr[rnbr[n]];
        InteriorParts[t].size = nrec;
      }
    } else
      for (s=0; s<ns; s++,t++)  InteriorParts[t].size = 0;
  }
  for (i=0,p=Particles+totalParts; i<ninj; i++,p++) {
    const int s = Particle_Spec(p->spec-sbase);
    const OH_nid_t nid = p->nid;
    const int ps = Secondary_Injected(nid) ? 1 : 0;
    dint *npg = NOfPGrid[ps][s];
    if (nid<0) continue;
    if (ps) {
      Primarize_Id_Only(p);
      Move_Or_Do(p, ps, oldp, 1,
                 (sb[--ninjs]=*p, InteriorParts[ns+s].size++), 0);
    } else
      Move_Or_Do(p, ps, me, 1,
                 (sb[nsend+ninjp++]=*p, InteriorParts[s].size++), 0);
  }
  move_to_sendbuf_uw4s(0, me, 0, 0);
  if (psold) {
    move_to_sendbuf_uw4s(1, oldp, primaryParts, nacc[0]);
    move_to_sendbuf_dw4s(1, oldp, totalParts, nacc[1]);
  } else {
    struct S_particle *rbb=Particles+nacc[0];
    int s;
    for (s=0; s<ns; s++) {
      RecvBufBases[ns+s] = rbb;  InteriorParts[ns+s].head = rbb - Particles;
      rbb += TotalPNext[ns+s];
    }
  }
  move_to_sendbuf_dw4s(0, me, primaryParts, nacc[0]);

  for (i=0,p=SendBuf+nsend; i<ninjp; i++,p++)
    *(RecvBufBases[Particle_Spec(p->spec-sbase)]++) = *p;
  for (i=ninjs,p=SendBuf+ninjs; i<nplim; i++,p++)
    *(RecvBufBases[Particle_Spec(p->spec-sbase)+ns]++) = *p;

  primaryParts = *secondaryBase = nacc[0];
}
static void
move_to_sendbuf_uw4s(const int ps, const int mysd, const int cbase,
                     const int nbase) {
  const int ns=nOfSpecies;
  const int nsor0 = ps ? ns : 0;
  const int *ctp = TotalP + nsor0,  *ntp = TotalPNext + nsor0;
  struct S_interiorp *ip = InteriorParts + nsor0;
  struct S_particle *p,  **rbb = RecvBufBases + nsor0,  *sb = SendBuf;
  int s, c, d, cn, dn;
  Decl_Grid_Info();

  for (s=0,c=cbase,d=nbase; s<ns; s++,c=cn,d=dn) {
    dint *npg = NOfPGrid[ps][s];
    cn = c + ctp[s];  dn  = d + ntp[s];
    ip[s].head = d;
    if (d<=c) {
      for (p=Particles+c; c<cn; c++,p++)
        Move_Or_Do(p, ps, mysd, 1, (Particles[d++]=*p), 0);
      rbb[s] = Particles + d;  ip[s].size += d - ip[s].head;
    } else if (dn<=cn) {
      const int cb = c;
      int cm, dm;
      for (p=Particles+c; c<d; c++,p++)
        Move_Or_Do(p, ps, mysd, 0, (d++), 0);
      cm = c - 1;  dm = d - 1;
      for (p=Particles+c; c<cn; c++,p++)
        Move_Or_Do(p, ps, mysd, 1, (Particles[d++]=*p), 0);
      rbb[s] = Particles + d;  ip[s].size += d - ip[s].head;
      for (c=dm,d=dm,p=Particles+c; c>=cb; c--,p--)
        Move_Or_Do(p, ps, mysd, 1, (Particles[d--]=*p), 0);
    }
  }
}
static void
move_to_sendbuf_dw4s(const int ps, const int mysd, const int ctail,
                     const int ntail) {
  const int ns=nOfSpecies;
  const int nsor0 = ps ? ns : 0;
  const int *ctp = TotalP + nsor0,  *ntp = TotalPNext + nsor0;
  struct S_interiorp *ip = InteriorParts + nsor0;
  struct S_particle *sb = SendBuf,  *p,  **rbb = RecvBufBases + nsor0;
  int s, c, d, cn, dn;
  Decl_Grid_Info();

  cn = ctail;  dn = ntail;
  for (s=ns-1,c=cn-1,d=dn-1; s>=0; s--,c=cn-1,d=dn-1) {
    dint *npg = NOfPGrid[ps][s];
    const int dd = d;
    cn -= ctp[s];  dn -= ntp[s];
    if (c>=d || cn>=dn)  continue;
    for (p=Particles+c; c>=cn; c--,p--)
      Move_Or_Do(p, ps, mysd, 1, (Particles[d--]=*p), 0);
    ip[s].head = d - ip[s].size + 1;  ip[s].size += dd - d;
    rbb[s] = Particles + ip[s].head;
  }
}
#define Sort_Particle(P) {\
  const int g = Grid_Position(P->nid);\
  const dint dst = npg[g];\
  SendBuf[npgt[g]++] = *P;\
  if (dst<0) {\
    const dint bsbidx = -dst;\
    if (!Is_Pillar_Voxel(bsbidx)) {\
      BoundarySendBuf[bsbidx-1] = *P;  npg[g] = dst - 1;\
    }\
    else {\
      BoundarySendBuf[Pillar_Lower(bsbidx)-1] =\
        BoundarySendBuf[Pillar_Upper(bsbidx)-1] = *P;\
      npg[g] = dst - (Add_Pillar_Voxel(1) + 1);\
    }\
  }\
}
static void
sort_particles(const int nextmode, const int psnew, const int stats) {
  const int ns=nOfSpecies;
  struct S_particle *p;
  int ps, s, t, i;
  Decl_Grid_Info();

  if (stats) oh1_stats_time(STATS_TB_SORT, nextmode);
  for (ps=0,t=0; ps<=psnew; ps++) {
    for (s=0; s<ns; s++,t++) {
      dint *npg = NOfPGrid[ps][s], *npgt = NOfPGridTotal[ps][s];
      const int ips = InteriorParts[t].size;
      for (i=0,p=Particles+InteriorParts[t].head; i<ips; i++,p++)
        Sort_Particle(p);
    }
  }
}
static void
move_and_sort(const int nextmode, const int psold, const int psnew,
              const int oldp, const int *nacc, const int stats) {
  const int me=myRank, ns=nOfSpecies, nn=nOfNodes, sbase=specBase;
  const int mysubdom[2] = {me, oldp},  ninj = nOfInjections;
  struct S_particle *p,  *rbb,  *sb = SendBuf + nacc[1];
  int *nofr;
  int ps, s, t, i;
  Decl_For_All_Grid();
  Decl_Grid_Info();

  if (stats) oh1_stats_time(STATS_TB_MOVE, nextmode);
  set_sendbuf_disps4s(nextmode, 0);
  for (ps=0,t=0,nofr=NOfRecv,rbb=Particles; ps<=psnew; ps++) {
    const int nnbr = RealSrcNeighbors[0][ps].n;
    const int *rnbr = RealSrcNeighbors[0][ps].nbor;
    for (s=0; s<ns; s++,t++,nofr+=nn) {
      int n, nrec;
      for (n=0,nrec=0; n<nnbr; n++)  nrec += nofr[rnbr[n]];
      RecvBufBases[t] = rbb;  rbb += nrec;
    }
  }
  RecvBufBases[t] = rbb;

  for (ps=0,p=Particles,t=0; ps<=psold; ps++) {
    const int mysd = mysubdom[ps];
    for (s=0; s<ns; s++,t++) {
      dint *npg = NOfPGrid[ps][s],  *npgt = NOfPGridTotal[ps][s];
      const int itail = TotalP[t];
      for (i=0; i<itail; i++,p++)
        Move_Or_Do(p, ps, mysd, 1, (SendBuf[npgt[g]++]=*p), 1);
    }
  }
  for (i=0; i<ninj; i++,p++) {
    const int s = Particle_Spec(p->spec-sbase);
    const OH_nid_t nid = p->nid;
    const int ps = Secondary_Injected(nid) ? 1 : 0;
    const int mysd = mysubdom[ps];
    dint *npg = NOfPGrid[ps][s],  *npgt = NOfPGridTotal[ps][s];
    if (nid<0) continue;
    if (ps)  Primarize_Id_Only(p);
    Move_Or_Do(p, ps, mysd, 1, (SendBuf[npgt[g]++]=*p), 1);
  }
  primaryParts = *secondaryBase = nacc[0];
}
static void
sort_received_particles(const int nextmode, const int psnew, const int stats) {
  const int ns=nOfSpecies;
  int ps, s;
  struct S_particle *p = Particles, **rbb = RecvBufBases+1;
  Decl_Grid_Info();

  if (stats) oh1_stats_time(STATS_TB_SORT, nextmode);
  for (ps=0; ps<=psnew; ps++) {
    for (s=0; s<ns; s++,rbb++) {
      dint *npg = NOfPGrid[ps][s], *npgt = NOfPGridTotal[ps][s];
      const struct S_particle *rbtail = *rbb;
      for (; p<rbtail; p++)  Sort_Particle(p);
    }
  }
}
static void
set_sendbuf_disps4s(const int nextmode, const int trans) {
  const int nn=nOfNodes, ns=nOfSpecies;
  int ps, s, i, np, *sbd;

  for (ps=0,sbd=NOfSend,np=0; ps<=nextmode; ps++) {
    const int n = RealDstNeighbors[trans][ps].n;
    const int *nbor = RealDstNeighbors[trans][ps].nbor;
    for (s=0; s<ns; s++,sbd+=nn) {
      for (i=0; i<n; i++) {
        const int nid = nbor[i];
        const int nsend = sbd[nid];
        sbd[nid] = np;  np += nsend;
      }
    }
  }
}
static void
xfer_particles(const int trans, const int psnew, const int nextmode,
               struct S_particle *sbuf) {
  const int nn=nOfNodes, ns=nOfSpecies;
  int ps, s, t, i, req, sdisp, *nofr, *nofs;

  for (ps=0,t=0,nofr=NOfRecv,req=0; ps<=psnew; ps++) {
    const int n = RealSrcNeighbors[trans][ps].n;
    const int *nbor = RealSrcNeighbors[trans][ps].nbor;
    for (s=0; s<ns; s++,t++,nofr+=nn) {
      struct S_particle *rbuf = RecvBufBases[t];
      for (i=0; i<n; i++) {
        const int nid = nbor[i];
        const int nrecv = nofr[nid];
        if (nrecv) {
          MPI_Irecv(rbuf, nrecv, T_Particle, nid, t, MCW, Requests+req++);
          rbuf += nrecv;
        }
      }
    }
  }
  for (ps=0,t=0,sdisp=0,nofs=NOfSend; ps<=nextmode; ps++) {
    const int n = RealDstNeighbors[trans][ps].n;
    const int *nbor = RealDstNeighbors[trans][ps].nbor;
    for (s=0; s<ns; s++,t++,nofs+=nn) {
      for (i=0; i<n; i++) {
        const int nid = nbor[i];
        const int sdnxt = nofs[nid];
        const int nsend = sdnxt - sdisp;
        nofs[nid] = 0;
        if (nsend) {
          MPI_Isend(sbuf+sdisp, nsend, T_Particle, nid, t, MCW,
                    Requests+req++);
        }
        sdisp = sdnxt;
      }
    }
  }
  MPI_Waitall(req, Requests, Statuses);
}
static void
xfer_boundary_particles_v(const int psnew, const int trans, const int d) {
  const int ns=nOfSpecies;
  int vphi=d*2*2;
  const int vphead=VPlaneHead[vphi], vptail=VPlaneHead[vphi+2*2];
  int i, s, req=0, ps;
  struct S_vplane *vp;
  struct S_particle *p;
  Decl_For_All_Grid();

  if (vphead==vptail)  return;

  for (i=vphead,vp=VPlane+vphead; i<vptail; i++,vp++) {
    const int nrecv = vp->nrecv;
    if (nrecv)
      MPI_Irecv(Particles+vp->rbuf, nrecv, T_Particle, vp->nbor, vp->rtag,
                MCW, Requests+req++);
  }
  for (i=vphead,vp=VPlane+vphead; i<vptail; i++,vp++) {
    const int nsend = vp->nsend;
    if (nsend)
      MPI_Isend(BoundarySendBuf+vp->sbuf, nsend, T_Particle, vp->nbor,
                vp->stag, MCW, Requests+req++);
  }
  if (req==0)  return;
  MPI_Waitall(req, Requests, Statuses);

  p = Particles + VPlane[vphead].rbuf;
  for (ps=0; ps<=psnew; ps++) {
    const int psor2 = ps ? trans + 1 : 0;
    const int zl = ZBound[ps][OH_LOWER];
    const int zu = ZBound[ps][OH_UPPER] - GridDesc[psor2].z;
    int du;
    for (du=OH_LOWER; du<=OH_UPPER; du++,vphi++) {
      int ny;
      int xl, yl, xu, yu;
      if (VPlaneHead[vphi]==VPlaneHead[vphi+1])  continue;
      if (d==OH_DIM_X) {
        ny = 1;
        Grid_Exterior_Boundary(du<<1, GridDesc[psor2].x, xl, xu);
      } else {
        ny = du<<1;
        xl = -OH_PGRID_EXT;  xu = OH_PGRID_EXT;
      }
      Grid_Exterior_Boundary(ny, GridDesc[psor2].y, yl, yu);
      For_All_Grid_Z(psor2, xl, yl, zl, xu, yu, zu) {
        for (s=0; s<ns; s++) {
          dint *npg=NOfPGrid[ps][s];
          int *npgo=NOfPGridOut[ps][s], *npgi=NOfPGridIndex[ps][s];
          For_All_Grid_XY(psor2, xl, yl, xu, yu) {
            const int g = The_Grid(),  tail = npgi[g] + npgo[g];
            const dint dst = npg[g];
            int i;
            if (Is_Pillar_Voxel(dst)) {
              struct S_particle *q = p;
              int j = Pillar_Upper(dst) - 1;
              for (i=npgi[g]; i<tail; i++)  BoundarySendBuf[j++] = *q++;
            }
            for (i=npgi[g]; i<tail; i++) {
              SendBuf[i] = *p++;  SendBuf[i].nid = -2;
            }
          }
        }
      }
    }
  }
}
static void
xfer_boundary_particles_h(const int psnew) {
  const int ns=nOfSpecies;
  int ps, ud, s, req=0;

  for (ps=0; ps<=psnew; ps++) {
    for (ud=OH_LOWER; ud<=OH_UPPER; ud++) {
      struct S_hplane *hp = HPlane[ps] + ud;
      int *nrecv = hp->nrecv,  *rbuf = hp->rbuf;
      const int nbor = hp->nbor,  tag = hp->rtag;
      if (nbor!=MPI_PROC_NULL) {
        for (s=0; s<ns; s++) {
          if (nrecv[s])
            MPI_Irecv(SendBuf+rbuf[s], nrecv[s], T_Particle, nbor, tag+s, MCW,
                      Requests+req++);
        }
      }
    }
  }
  for (ps=0; ps<=psnew; ps++) {
    for (ud=OH_LOWER; ud<=OH_UPPER; ud++) {
      struct S_hplane *hp = HPlane[ps] + ud;
      int *nsend = hp->nsend,  *sbuf = hp->sbuf;
      const int nbor = hp->nbor,  tag = hp->stag;
      if (nbor!=MPI_PROC_NULL) {
        for (s=0; s<ns; s++) {
          if (nsend[s])
            MPI_Isend(SendBuf+sbuf[s], nsend[s], T_Particle, nbor, tag+s, MCW,
                      Requests+req++);
        }
      }
    }
  }
  if (req==0)  return;
  MPI_Waitall(req, Requests, Statuses);

  for (ps=0; ps<=psnew; ps++) {
    for (ud=OH_LOWER; ud<=OH_UPPER; ud++) {
      struct S_hplane *hp = HPlane[ps] + ud;
      int *nrecv = hp->nrecv,  *rbuf = hp->rbuf;
      if (hp->nbor!=MPI_PROC_NULL) {
        for (s=0; s<ns; s++) {
          const int tail = rbuf[s] + nrecv[s];
          int i;
          for (i=rbuf[s]; i<tail; i++)  SendBuf[i].nid = -2;
        }
      }
    }
  }
}
void
oh4s_exchange_border_data_(void *buf, void *sbuf, void *rbuf, int *type) {
  oh4s_exchange_border_data(buf, sbuf, rbuf, MPI_Type_f2c(*type));
}
void
oh4s_exchange_border_data(void *buf, void *sbuf, void *rbuf,
                          MPI_Datatype type) {
  MPI_Aint esize, lb;

  MPI_Type_get_extent(type, &lb, &esize);
  exchange_border_data_v(buf, sbuf, rbuf, type, esize, 0);
  exchange_border_data_v(buf, sbuf, rbuf, type, esize, 1);
  exchange_border_data_h(buf, type, esize);
}
static void
exchange_border_data_v(void *buf, void *sbuf, void *rbuf, MPI_Datatype type,
                       const MPI_Aint esize, const int d) {
  char *b = (char*)buf, *sb = (char*)sbuf,  *rb = (char*)rbuf;
  const int ns=nOfSpecies, pscurr=RegionId[1]<0 ? 0 : 1, vphi=d*2*2;
  const int vphead=VPlaneHead[vphi], vptail=VPlaneHead[vphi+2*2];
  struct S_vplane *vp;
  int ps, s, i, req=0;
  Decl_For_All_Grid();

  if (vphead == vptail)  return;

  for (ps=0,i=vphi; ps<=pscurr; ps++) {
    int du;
    const int zl = ZBound[ps][OH_LOWER];
    const int zu = ZBound[ps][OH_UPPER] - GridDesc[ps].z;
    for (du=OH_LOWER; du<=OH_UPPER; du++,i++) {
      int ny;
      int xl, yl, xu, yu;
      if (VPlaneHead[i]==VPlaneHead[i+1])  continue;
      if (d==OH_DIM_X) {
        ny = 1;
        Grid_Interior_Boundary(du<<1, GridDesc[ps].x, xl, xu);
      } else {
        ny = du<<1;
        xl = -OH_PGRID_EXT;  xu = OH_PGRID_EXT;
      }
      Grid_Interior_Boundary(ny, GridDesc[ps].y, yl, yu);
      For_All_Grid_Z(ps, xl, yl, zl, xu, yu, zu) {
        for (s=0; s<ns; s++) {
          int *npgo=NOfPGridOut[ps][s], *npgi=NOfPGridIndex[ps][s];
          For_All_Grid_XY(ps, xl, yl, xu, yu) {
            const int g = The_Grid(),  nbyte = npgo[g]*esize;
            memcpy(sb, b+npgi[g]*esize, nbyte);
            sb += nbyte;
          }
        }
      }
    }
  }
  rb -= VPlane[vphead].rbuf * esize;
  for (i=vphead,vp=VPlane+vphead; i<vptail; i++,vp++) {
    const int nrecv = vp->nrecv;
    if (nrecv)
      MPI_Irecv(rb+vp->rbuf*esize, nrecv, type, vp->nbor, vp->rtag, MCW,
                Requests+req++);
  }
  sb = (char*)sbuf - VPlane[vphead].sbuf * esize;
  for (i=vphead,vp=VPlane+vphead; i<vptail; i++,vp++) {
    const int nsend = vp->nsend;
    if (nsend)
      MPI_Isend(sb+vp->sbuf*esize, nsend, type, vp->nbor, vp->stag, MCW,
                Requests+req++);
  }
  if (req==0)  return;
  MPI_Waitall(req, Requests, Statuses);

  rb = (char*)rbuf;
  for (ps=0,i=vphi; ps<=pscurr; ps++) {
    const int zl = ZBound[ps][OH_LOWER];
    const int zu = ZBound[ps][OH_UPPER] - GridDesc[ps].z;
    int du;
    for (du=OH_LOWER; du<=OH_UPPER; du++,i++) {
      int ny;
      int xl, yl, xu, yu;
      if (VPlaneHead[i]==VPlaneHead[i+1])  continue;
      if (d==OH_DIM_X) {
        ny = 1;
        Grid_Exterior_Boundary(du<<1, GridDesc[ps].x, xl, xu);
      } else {
        ny = du<<1;
        xl = -OH_PGRID_EXT;  xu = OH_PGRID_EXT;
      }
      Grid_Exterior_Boundary(ny, GridDesc[ps].y, yl, yu);
      For_All_Grid_Z(ps, xl, yl, zl, xu, yu, zu) {
        for (s=0; s<ns; s++) {
          int *npgo=NOfPGridOut[ps][s], *npgi=NOfPGridIndex[ps][s];
          For_All_Grid_XY(ps, xl, yl, xu, yu) {
            const int g = The_Grid(),  nbyte = npgo[g]*esize;
            memcpy(b+npgi[g]*esize, rb, nbyte);
            rb += nbyte;
          }
        }
      }
    }
  }
}
static void
exchange_border_data_h(void *buf, MPI_Datatype type, const MPI_Aint esize) {
  char *b=(char*)buf;
  const int ns=nOfSpecies, pscurr=RegionId[1]<0 ? 0 : 1;
  int ps, ud, s, req=0;
  Decl_For_All_Grid();

  for (ps=0; ps<=pscurr; ps++) {
    for (ud=OH_LOWER; ud<=OH_UPPER; ud++) {
      struct S_hplane *hp = HPlane[ps] + ud;
      int *nrecv = hp->nrecv,  *rbuf = hp->rbuf;
      const int nbor = hp->nbor,  tag = hp->rtag;
      if (nbor!=MPI_PROC_NULL) {
        for (s=0; s<ns; s++) {
          if (nrecv[s])
            MPI_Irecv(b+rbuf[s]*esize, nrecv[s], type, nbor, tag+s, MCW,
                      Requests+req++);
        }
      }
    }
  }
  for (ps=0; ps<=pscurr; ps++) {
    for (ud=OH_LOWER; ud<=OH_UPPER; ud++) {
      struct S_hplane *hp = HPlane[ps] + ud;
      int *nsend = hp->nsend,  *sbuf = hp->sbuf;
      const int nbor = hp->nbor,  tag = hp->stag;
      if (nbor!=MPI_PROC_NULL) {
        for (s=0; s<ns; s++) {
          if (nsend[s])
            MPI_Isend(b+sbuf[s]*esize, nsend[s], type, nbor, tag+s, MCW,
                      Requests+req++);
        }
      }
    }
  }
  if (req)  MPI_Waitall(req, Requests, Statuses);
}
#ifndef OH_NO_CHECK
#define Check_Particle_Location(P, PS, S, NS, INJ) {\
  const int t = (PS) ? (S)+(NS) : (S);\
  const int pidx = (P) - Particles;\
  if ((PS)<0 || (PS)>1 || (S)<0 || (S)>=(NS) ||\
      (PbufIndex && ((INJ) ?\
                     (((PS)&&RegionId[1]<0) ||\
                         pidx>=totalParts+nOfInjections) :\
                     (pidx<PbufIndex[t] || pidx>=PbufIndex[t+1]))))\
    local_errstop("'part' argument pointing %c%d%c of the particle buffer is "\
                  "inconsistent with 'ps'=%d and 's'=%d",\
                  specBase?'(':'[', pidx+specBase,\
                  specBase?')':']', PS, (S)+specBase);\
}
#else
#define Check_Particle_Location(P, PS, S, NS, INJ)
#endif
#define Map_Particle_To_Neighbor(P, XYZ, DIM, MYSD, K, INC, UB, G, IDX) {\
  const double xyz = XYZ;\
  const double gsize = Grid[DIM].gsize;\
  const double lb = Grid[DIM].fcoord[OH_LOWER];\
  const double gf =\
    (G = (xyz-lb)*Grid[DIM].rgsize + Grid[DIM].coord[OH_LOWER]) * gsize;\
  if (gf>xyz) G--;\
  else if (gf+gsize<=xyz) G++;\
  G -= SubDomains[MYSD][DIM][OH_LOWER];  IDX += G;\
  if (G<0) {\
    K -= INC;\
    if (xyz<lb) {\
      if (Boundaries[MYSD][DIM][OH_LOWER]) { P->nid = -1;  return(-1); }\
      XYZ += Grid[DIM].fcoord[OH_UPPER] - lb;\
    }\
    if (G<-OH_PGRID_EXT)  K = -OH_NEIGHBORS;\
  } else if (G>=UB) {\
    double ub = Grid[DIM].fcoord[OH_UPPER];\
    K += INC;\
    if (xyz>=ub) {\
      if (Boundaries[MYSD][DIM][OH_UPPER]) { P->nid = -1;  return(-1); }\
      XYZ -= ub - lb;\
    }\
    G-=UB;\
    if (G>=OH_PGRID_EXT)  K = -OH_NEIGHBORS;\
  }\
}
#define Adjust_Neighbor_Grid(G, N, DIM)\
  if (G<0) G += SubDomains[N][DIM][OH_UPPER]-SubDomains[N][DIM][OH_LOWER];
int
oh4s_map_particle_to_neighbor_(struct S_particle *part, const int *ps,
                               const int *s) {
  return(oh4s_map_particle_to_neighbor(part, *ps, *s-1));
}
int
oh4s_map_particle_to_neighbor(struct S_particle *part, const int ps,
                              const int s) {
  const int ns = nOfSpecies,  inj = part>=Particles+totalParts;
  int x, y, z, w, d, dw, mysd;
  const int psnn = ps ? (s+nOfSpecies)*nOfNodes : s*nOfNodes;
  int k = OH_NBR_SELF,  idx = 0;
  int gz, gy, gx;
  int sd;
  Decl_Grid_Info();

  Check_Particle_Location(part, ps, s, ns, inj);
  x = GridDesc[ps].x;  y = GridDesc[ps].y;  z = GridDesc[ps].z;
  w = GridDesc[ps].w;  d = GridDesc[ps].d;  dw = GridDesc[ps].dw;
  mysd = RegionId[ps];
  Do_Z(Map_Particle_To_Neighbor(part, part->z, OH_DIM_Z, mysd, k, 9, z, gz,
                                idx));
  Do_Z(idx *= d);
  Do_Y(Map_Particle_To_Neighbor(part, part->y, OH_DIM_Y, mysd, k, 3, y, gy,
                                idx));
  Do_Y(idx *= w);
  Map_Particle_To_Neighbor(part, part->x, OH_DIM_X, mysd, k, 1, x, gx, idx);

  if (k==OH_NBR_SELF) {
    NOfPGrid[ps][s][idx]++;
    NOfPLocal[psnn+mysd]++;
    part->nid = Combine_Subdom_Pos(k, idx);
    if (inj) {
      if (ps) {
        InjectedParticles[ns+s]++;  Secondarize_Id(part);
      } else {
        InjectedParticles[s]++;
      }
    }
    return(mysd);
  } else if (k<0)
    return(oh4s_map_particle_to_subdomain(part, ps, s));
  sd = AbsNeighbors[ps][k];
  if (sd>=nOfNodes) {
    part->nid = -1;  return(-1);
  }
  Adjust_Neighbor_Grid(gx, sd, OH_DIM_X);
  Do_Y(Adjust_Neighbor_Grid(gy, sd, OH_DIM_Y));
  Do_Z(Adjust_Neighbor_Grid(gz, sd, OH_DIM_Z));
  NOfPLocal[psnn+sd]++;

  if (sd==mysd) {
    idx = Coord_To_Index(gx, gy, gz, w, dw);
    NOfPGrid[ps][s][idx]++;
    part->nid = Combine_Subdom_Pos(OH_NBR_SELF, idx);
    if (inj)  InjectedParticles[ps ? ns+s : s]++;
  } else {
    NOfPGrid[ps][s][idx]++;
    part->nid = Combine_Subdom_Pos(k, Coord_To_Index(gx, gy, gz, w, dw));
  }
  if (inj && ps)  Secondarize_Id(part);
  return(sd);
}
#define Map_To_Grid(P, PXYZ, XYZ, DIM, GG, LG) {\
  const double gsize = Grid[DIM].gsize;\
  const double lb = Grid[DIM].fcoord[OH_LOWER];\
  const double ub = Grid[DIM].fcoord[OH_UPPER];\
  double gf;\
  XYZ = PXYZ;\
  LG = 0;\
  if (XYZ<lb) {\
    if (BoundaryCondition[DIM][OH_LOWER]) { P->nid = -1;  return(-1); }\
    XYZ += (ub - lb);  PXYZ = XYZ;\
    LG = Grid[DIM].coord[OH_LOWER] - Grid[DIM].coord[OH_UPPER];\
  }\
  else if (XYZ>=ub) {\
    if (BoundaryCondition[DIM][OH_UPPER]) { P->nid = -1;  return(-1); }\
    XYZ -= (ub - lb);  PXYZ = XYZ;\
    LG = Grid[DIM].coord[OH_UPPER] - Grid[DIM].coord[OH_LOWER];\
  }\
  GG = (XYZ-lb)*Grid[DIM].rgsize + Grid[DIM].coord[OH_LOWER];\
  gf = GG * gsize;\
  if (gf>XYZ) GG--;\
  else if (gf+gsize<=XYZ) GG++;\
  LG += GG;\
}
#define Map_Particle_To_Subdomain(XYZ, DIM, SDOM) {\
  double thresh = Grid[DIM].light.thresh;\
  if (XYZ<thresh)\
    SDOM = (XYZ - Grid[DIM].coord[OH_LOWER]) / Grid[DIM].light.size;\
  else\
    SDOM = (XYZ - thresh)/ (Grid[DIM].light.size + 1) + Grid[DIM].light.n;\
}
#define Local_Coordinate(N, MYSD, GG, LG, DIM, K, INC, AA) {\
  GG -= SubDomains[N][DIM][OH_LOWER];\
  if (N==MYSD)  LG = GG;\
  else {\
    const int ub = SubDomains[MYSD][DIM][OH_UPPER];\
    if (LG>=ub+OH_PGRID_EXT)  AA = 1;\
    else {\
      const int inc = LG<ub ? 0 : INC;\
      LG -= SubDomains[MYSD][DIM][OH_LOWER];\
      if (LG<-OH_PGRID_EXT)  AA = 1;\
      k += LG<0 ? -INC : inc;\
    }\
  }\
}
int
oh4s_map_particle_to_subdomain_(struct S_particle *part, const int *ps,
                                const int *s) {
  return(oh4s_map_particle_to_subdomain(part, *ps, *s-1));
}
int
oh4s_map_particle_to_subdomain(struct S_particle *part, const int ps,
                               const int s) {
  const int ns = nOfSpecies,  inj = part>=Particles+totalParts;
  const int nx  = Grid[OH_DIM_X].n;
  const int nxy = If_Dim(OH_DIM_Y, nx*Grid[OH_DIM_Y].n, 0);
  const int t = ps ? ns + s : s;
  int w, dw, mysd;
  int sd;
  double x, y, z;
  int px, py, pz;
  int gx, gy, gz;
  int lx, ly, lz;
  int k = OH_NBR_SELF,  aacc = 0;
  Decl_Grid_Info();

  Check_Particle_Location(part, ps, s, ns, inj);
  w = GridDesc[ps].w;  dw = GridDesc[ps].dw;  mysd = RegionId[ps];
  Map_To_Grid(part, part->x, x, OH_DIM_X, gx, lx);
  Do_Y(Map_To_Grid(part, part->y, y, OH_DIM_Y, gy, ly));
  Do_Z(Map_To_Grid(part, part->z, z, OH_DIM_Z, gz, lz));
  if (SubDomainDesc) {
    sd = map_irregular_subdomain(x, If_Dim(OH_DIM_Y ,y, 0),
                                    If_Dim(OH_DIM_Z, z, 0));
    if (sd<0) { part->nid = -1;  return(-1); }
  } else {
    Map_Particle_To_Subdomain(gx, OH_DIM_X, px);
    Do_Y(Map_Particle_To_Subdomain(gy, OH_DIM_Y, py));
    Do_Z(Map_Particle_To_Subdomain(gz, OH_DIM_Z, pz));
    sd = Coord_To_Index(px, py, pz, nx, nxy);
  }
  Local_Coordinate(sd, mysd, gx, lx, OH_DIM_X, k, 1, aacc);
  Do_Y(Local_Coordinate(sd, mysd, gy, ly, OH_DIM_Y, k, 3, aacc));
  Do_Z(Local_Coordinate(sd, mysd, gz, lz, OH_DIM_Z, k, 9, aacc));
  NOfPLocal[t*nOfNodes+sd]++;
  if (aacc) {
    currMode = Mode_Set_Any(currMode);
    part->nid = Combine_Subdom_Pos(sd+OH_NEIGHBORS,
                                   Coord_To_Index(gx, gy, gz, w, dw));
  } else {
    NOfPGrid[ps][s][Coord_To_Index(lx, ly, lz, w, dw)]++;
    part->nid = Combine_Subdom_Pos(k, Coord_To_Index(gx, gy, gz, w, dw));
  }
  if (inj) {
    if (sd==mysd)  InjectedParticles[t]++;
    if (ps)  Secondarize_Id(part);
  }
  return(sd);
}
int
oh4s_inject_particle_(const struct S_particle *part, const int *ps) {
  return(oh4s_inject_particle(part, *ps));
}
int
oh4s_inject_particle(const struct S_particle *part, const int ps) {
  const int ns = nOfSpecies;
  int inj = totalParts + nOfInjections++;
  struct S_particle *p = Particles + inj;
  int s = Particle_Spec(part->spec - specBase);
  int sd;

#ifndef OH_HAS_SPEC
  if (ns!=1)
    local_errstop("particles cannot be injected when S_particle does not "
                  "have 'spec' element and you have two or more species");
#endif
  if (inj>=nOfLocalPLimit)
    local_errstop("injection causes local particle buffer overflow");
  *p = *part;
  sd = oh4s_map_particle_to_neighbor(p, ps, s);
  if (sd<0)  nOfInjections--;
  return(sd);
}
void
oh4s_remove_mapped_particle_(struct S_particle *part, const int *ps,
                             const int *s) {
  oh4s_remove_mapped_particle(part, *ps, *s-1);
}
void
oh4s_remove_mapped_particle(struct S_particle *part, const int ps,
                            const int s) {
  const int nn = nOfNodes, ns = nOfSpecies, inj = part>=Particles+totalParts;
  OH_nid_t nid = part->nid;
  int sd, g, psreal=ps, mysd, t;
  Decl_Grid_Info();

  Check_Particle_Location(part, psreal, s, ns, inj);
  if (nid<0)  return;
  sd = Subdomain_Id(nid, psreal);
  g = Grid_Position(nid);
  if (sd>=nn) {
    psreal = 1;  Primarize_Id(part, sd);  nid = part->nid;
  }
  mysd = RegionId[psreal];
  part->nid = -1;
  t = psreal ? ns+s : s;
  NOfPLocal[t*nn+sd]--;
  if (inj && sd==mysd)  InjectedParticles[t]--;
  if (Mode_Acc(currMode))  return;
  if (sd!=mysd)  g = Local_Grid_Position(g, nid, psreal);
  NOfPGrid[psreal][s][g]--;
}
int
oh4s_remap_particle_to_neighbor_(struct S_particle *part, const int *ps,
                                 const int *s) {
  return(oh4s_remap_particle_to_neighbor(part, *ps, *s-1));
}
int
oh4s_remap_particle_to_neighbor(struct S_particle *part, const int ps,
                                const int s) {
  oh4s_remove_mapped_particle(part, ps, s);
  return(oh4s_map_particle_to_neighbor(part, ps, s));
}
int
oh4s_remap_particle_to_subdomain_(struct S_particle *part, const int *ps,
                                  const int *s) {
  return(oh4s_remap_particle_to_subdomain(part, *ps, *s-1));
}
int
oh4s_remap_particle_to_subdomain(struct S_particle *part, const int ps,
                                 const int s) {
  oh4s_remove_mapped_particle(part, ps, s);
  return(oh4s_map_particle_to_subdomain(part, ps, s));
}
