/* File: ohhelp4s.h
   Version 1.1.1 (2015/10/23)
   Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
                            (ACCMS, Kyoto University)
   This program can be freely used, redistributed and modified for non-
   commercial purpose providing that the copyright notice above remains
   unchanged.
*/
#define OH_PGRID_EXT    1
#define OH_NBR_SELF     (OH_NEIGHBORS>>1)
#define OH_NBR_BCC      (1+1*3+0*3*3)
#define OH_NBR_TCC      (1+1*3+2*3*3)

#define Grid_Position(ID)         ((ID)&gridmask)
#define Combine_Subdom_Pos(ID, G) (((OH_nid_t)(ID)<<loggrid) + (G))
#define Primarize_Id_Only(P) \
  (P)->nid -= (OH_nid_t)(nOfNodes+OH_NEIGHBORS)<<loggrid
#define Secondarize_Id(P) \
  (P)->nid += (OH_nid_t)(nOfNodes+OH_NEIGHBORS)<<loggrid
#define Secondary_Injected(ID) \
  ((ID>>loggrid)>=nOfNodes+OH_NEIGHBORS)
#define Neighbor_Subdomain_Id(ID, PS) \
  AbsNeighbors[PS][(ID)>>loggrid]

EXTERN int *PbufIndex;                                  /* [2*ns+1] */
EXTERN dint **NOfPGrid[2], **NOfPGridTotal[2];          /* [2][ns][z][y][x] */
EXTERN int **NOfPGridOut[2], **NOfPGridOutShadow[2];    /* [2][ns][z][y][x] */
EXTERN int **NOfPGridIndex[2], **NOfPGridIndexShadow[2];/* [2][ns][z][y][x] */
EXTERN dint *NOfPGridZ;                                 /* [z] */
EXTERN int ZBound[2][2], (*ZBoundShadow)[2];
struct S_hplane {
  int nbor, stag, rtag;
  int *nsend, *nrecv, *sbuf, *rbuf;                     /* [ns] */
};
EXTERN struct S_hplane HPlane[2][2];                    /* [2][2] */
struct S_vplane {
  int nbor, stag, rtag;
  int nsend, nrecv, sbuf, rbuf;
};
EXTERN struct S_vplane *VPlane;                         /* [2*nn+6] */
EXTERN int VPlaneHead[2*2*2+1];
EXTERN struct S_particle *BoundarySendBuf;

struct S_griddesc {
  int x, y, z, w, d, h, dw;
};
EXTERN struct S_griddesc GridDesc[3];

struct S_interiorp {
  int head, size;
};
EXTERN struct S_interiorp *InteriorParts;

EXTERN struct S_commlist *AltSecRList, PrimaryCommList[2][OH_NEIGHBORS];
EXTERN int SecRLIndex[OH_NEIGHBORS+1], AltSecRLIndex[OH_NEIGHBORS+1];
EXTERN int PrimaryRLIndex[OH_NEIGHBORS];

struct S_recvsched_context {
  int z;
  dint nptotal, nplimit;
  struct S_commlist *cptr;
};
EXTERN MPI_Datatype T_Hgramhalf;

EXTERN int FirstNeighbor[OH_NEIGHBORS], GridOffset[2][OH_NEIGHBORS];
struct S_realneighbor {
  int n, *nbor;
};
EXTERN struct S_realneighbor RealDstNeighbors[2][2], RealSrcNeighbors[2][2];

EXTERN int BoundaryCondition[OH_DIMENSION][2];

/* Prototypes for the functions called from simulator code */
void oh4s_init(int **sdid, const int nspec, const int maxfrac,
               const dint npmax, const int minmargin,
               const int maxdensity, int **totalp, int **pbase,
               int *maxlocalp, int *cbufsize, void *mycomm, int **nbor,
               int *pcoord, int **sdoms, int *scoord, const int nbound,
               int *bcond, int **bounds, int *ftypes, int *cfields,
               int *ctypes, int **fsizes, int **zbound,
               const int stats, const int repiter, const int verbose);
void oh4s_particle_buffer(const int maxlocalp, struct S_particle **pbuf);
void oh4s_per_grid_histogram(int **pghgram, int **pgindex);
int  oh4s_transbound(int currmode, int stats);
void oh4s_exchange_border_data(void *buf, void *sbuf, void *rbuf,
                               MPI_Datatype type);
int  oh4s_map_particle_to_neighbor(struct S_particle *part, const int ps,
                                   const int s);
int  oh4s_map_particle_to_subdomain(struct S_particle *part, const int ps,
                                    const int s);
int  oh4s_inject_particle(const struct S_particle *part, const int ps);
void oh4s_remove_mapped_particle(struct S_particle *part, const int ps,
                                 const int s);
int  oh4s_remap_particle_to_neighbor(struct S_particle *part, const int ps,
                                     const int s);
int  oh4s_remap_particle_to_subdomain(struct S_particle *part, const int ps,
                                      const int s);

void oh4s_init_(int *sdid, const int *nspec, const int *maxfrac,
                const dint *npmax, const int *minmargin, const int *maxdensity,
                int *totalp, int *pbase, int *maxlocalp, int *cbufsize,
                struct S_mycommf *mycomm, int *nbor, int *pcoord, int *sdoms,
                int *scoord, const int *nbound, int *bcond, int *bounds,
                int *ftypes, int *cfields, int *ctypes, int *fsizes,
                int *zbound,
                const int *stats, const int *repiter, const int *verbose);
void oh4s_particle_buffer_(const int *maxlocalp, struct S_particle *pbuf);
void oh4s_per_grid_histogram_(int *pghgram, int *pgindex);
int  oh4s_transbound_(int *currmode, int *stats);
void oh4s_exchange_border_data_(void *buf, void *sbuf, void *rbuf, int *type);
int  oh4s_map_particle_to_neighbor_(struct S_particle *part, const int *ps,
                                    const int *s);
int  oh4s_map_particle_to_subdomain_(struct S_particle *part, const int *ps,
                                     const int *s);
int  oh4s_inject_particle_(const struct S_particle *part, const int *ps);
void oh4s_remove_mapped_particle_(struct S_particle *part, const int *ps,
                                  const int *s);
int  oh4s_remap_particle_to_neighbor_(struct S_particle *part, const int *ps,
                                      const int *s);
int  oh4s_remap_particle_to_subdomain_(struct S_particle *part, const int *ps,
                                       const int *s);
