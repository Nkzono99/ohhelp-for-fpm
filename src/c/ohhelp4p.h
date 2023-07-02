/* File: ohhelp4p.h
   Version 1.1.1 (2015/10/23)
   Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
                            (ACCMS, Kyoto University)
   This program can be freely used, redistributed and modified for non-
   commercial purpose providing that the copyright notice above remains
   unchanged.
*/
#define OH_PGRID_EXT 1
#define OH_NBR_SELF (OH_NEIGHBORS>>1)

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
EXTERN int **NOfPGridOut[2];                            /* [2][ns][z][y][x] */

struct S_griddesc {
  int x, y, z, w, d, h, dw;
};
EXTERN struct S_griddesc GridDesc[3];

EXTERN int gridOverflowLimit;
EXTERN struct S_commlist *AltSecRList;
EXTERN int SecRLIndex[OH_NEIGHBORS+1];

struct S_recvsched_context {
  int x, y, z, g, hs;
  dint nptotal, nplimit, carryover;
  struct S_commlist *cptr;
};
struct S_hotspot {
  int g, n, lev, self;
  struct S_commlist *comm;
  struct S_hotspot *next;
};
EXTERN struct S_hotspot *HotSpotList, *HotSpotTop;      /* [2*nn+2*3^D+1] */
struct S_hotspotbase {
  struct S_hotspot *head, *tail;
};
EXTERN struct S_hotspotbase HotSpot[3][OH_NEIGHBORS];

EXTERN int *HSRecv[OH_NEIGHBORS];                       /* [3^D][nn][ns] */
EXTERN int *HSSend, *HSRecvFromParent, *HSReceiver;     /* [ns] */
EXTERN MPI_Datatype T_Hgramhalf;

EXTERN int FirstNeighbor[OH_NEIGHBORS], GridOffset[2][OH_NEIGHBORS];
struct S_realneighbor {
  int n, *nbor;
};
EXTERN struct S_realneighbor RealDstNeighbors[2][2], RealSrcNeighbors[2][2];

EXTERN int BoundaryCondition[OH_DIMENSION][2];

/* Prototypes for the functions called from simulator code */
void oh4p_init(int **sdid, const int nspec, const int maxfrac, int **totalp,
               struct S_particle **pbuf, int **pbase, const int maxlocalp,
               void *mycomm, int **nbor, int *pcoord, int **sdoms, int *scoord,
               const int nbound, int *bcond, int **bounds, int *ftypes,
               int *cfields, int *ctypes, int **fsizes,
               const int stats, const int repiter, const int verbose);
int  oh4p_max_local_particles(const dint npmax, const int maxfrac,
                              const int minmargin, const int hsthresh);
void oh4p_per_grid_histogram(int **pghgram);
int  oh4p_transbound(int currmode, int stats);
int  oh4p_map_particle_to_neighbor(struct S_particle *part, const int ps,
                                   const int s);
int  oh4p_map_particle_to_subdomain(struct S_particle *part, const int ps,
                                    const int s);
int  oh4p_inject_particle(const struct S_particle *part, const int ps);
void oh4p_remove_mapped_particle(struct S_particle *part, const int ps,
                                 const int s);
int  oh4p_remap_particle_to_neighbor(struct S_particle *part, const int ps,
                                     const int s);
int  oh4p_remap_particle_to_subdomain(struct S_particle *part, const int ps,
                                      const int s);

void oh4p_init_(int *sdid, const int *nspec, const int *maxfrac, int *totalp,
                struct S_particle *pbuf, int *pbase, const int *maxlocalp,
                struct S_mycommf *mycomm, int *nbor, int *pcoord, int *sdoms,
                int *scoord, const int *nbound, int *bcond, int *bounds,
                int *ftypes, int *cfields, int *ctypes, int *fsizes,
                const int *stats, const int *repiter, const int *verbose);
int  oh4p_max_local_particles_(const dint *npmax, const int *maxfrac,
                               const int *minmargin, const int *hsthresh);
void oh4p_per_grid_histogram_(int *pghgram);
int  oh4p_transbound_(int *currmode, int *stats);
int  oh4p_map_particle_to_neighbor_(struct S_particle *part, const int *ps,
                                    const int *s);
int  oh4p_map_particle_to_subdomain_(struct S_particle *part, const int *ps,
                                     const int *s);
int  oh4p_inject_particle_(const struct S_particle *part, const int *ps);
void oh4p_remove_mapped_particle_(struct S_particle *part, const int *ps,
                                  const int *s);
int  oh4p_remap_particle_to_neighbor_(struct S_particle *part, const int *ps,
                                      const int *s);
int  oh4p_remap_particle_to_subdomain_(struct S_particle *part, const int *ps,
                                       const int *s);
