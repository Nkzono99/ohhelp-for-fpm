/* File: ohhelp2.h
   Version 1.1.1 (2015/10/23)
   Copyright (C) 2009-2015  Hiroshi Nakashima <h.nakashima@media.kyoto-u.ac.jp>
                            (ACCMS, Kyoto University)
   This program can be freely used, redistributed and modified for non-
   commercial purpose providing that the copyright notice above remains
   unchanged.
*/
#include "oh_part.h"

EXTERN int nOfLocalPLimit;
EXTERN struct S_particle *Particles;    /* [nOfLocalPLimit] */
EXTERN struct S_particle *SendBuf;      /* [nOfLocalPLimit] */
EXTERN struct S_particle **RecvBufBases;/* [2][nOfSpecies] */
EXTERN int *secondaryBase, *totalLocalParticles;
EXTERN int *SendBufDisps;               /* [nOfSpecies][nOfNodes] */
EXTERN int *RecvBufDisps;               /* [nOfNodes] */
EXTERN int nOfInjections;
EXTERN int specBase;
EXTERN MPI_Datatype T_Particle;
EXTERN MPI_Request *Requests;           /* [nOfNodes*nOfSpecies*2*2] */
EXTERN MPI_Status *Statuses;            /* [nOfNodes*nOfSpecies*2*2] */

#ifdef OH_HAS_SPEC
#define Particle_Spec(S) (S)
#else
#define Particle_Spec(S) (0)
#endif

#ifdef OH_POS_AWARE
EXTERN int gridMask, logGrid;
EXTERN int AbsNeighbors[2][OH_NEIGHBORS];
#define Decl_Grid_Info() \
  OH_nid_t nidelement;  int subdomid;\
  const int gridmask=gridMask, loggrid=logGrid
#define Subdomain_Id(ID, PS) \
  ((nidelement=(ID))<0 ? -1 :\
      ((subdomid=nidelement>>loggrid)<OH_NEIGHBORS ?\
          AbsNeighbors[PS][subdomid] : subdomid-OH_NEIGHBORS))
#define Primarize_Id(P, SD) {\
  const OH_nid_t nidelem =\
    ((P)->nid -= (OH_nid_t)(nOfNodes+OH_NEIGHBORS)<<loggrid);\
  SD = Subdomain_Id(nidelem, 1);\
}
#else
#define Decl_Grid_Info() int unusedvariable
#define Subdomain_Id(ID, PS) (ID)
#endif

/* Prototypes for the functions called from simulator code */
void oh2_set_total_particles();
int  oh2_max_local_particles(dint npmax, int maxfrac, int minmargin);
void oh2_inject_particle(struct S_particle *part);
void oh2_remap_injected_particle(struct S_particle *part);
void oh2_remove_injected_particle(struct S_particle *part);
void oh2_init(int **sdid, int nspec, int maxfrac, int **nphgram,
              int **totalp, struct S_particle **pbuf, int **pbase,
              int maxlocalp, void *mycomm, int **nbor,
              int *pcoord, int stats, int repiter, int verbose);
int  oh2_transbound(int currmode, int stats);

void oh2_set_total_particles_();
int  oh2_max_local_particles_(dint *npmax, int *maxfrac, int *minmargin);
void oh2_inject_particle_(struct S_particle *part);
void oh2_remap_injected_particle_(struct S_particle *part);
void oh2_remove_injected_particle_(struct S_particle *part);
void oh2_init_(int *sdid, int *nspec, int *maxfrac, int *nphgram,
               int *totalp, struct S_particle *pbuf, int *pbase,
               int *maxlocalp, struct S_mycommf *mycomm, int *nbor,
               int *pcoord, int *stats, int *repiter, int *verbose);
int  oh2_transbound_(int *currmode, int *stats);

/* Prototypes for the functions called from higher-level library code */
void init2(int **sdid, int nspec, int maxfrac, int **nphgram,
           int **totalp, struct S_particle **pbuf, int **pbase, int maxlocalp,
           struct S_mycommc *mycommc, struct S_mycommf *mycommf,
           int **nbor, int *pcoord, int stats, int repiter, int verbose);
int  transbound2(int currmode, int stats, int level);
void exchange_primary_particles(int currmode, int stats);
void move_to_sendbuf_primary(int secondary, int stats);
void set_sendbuf_disps(int secondary, int parent);
void exchange_particles(struct S_commlist *secrlist, int secrlsize,
                        int oldparent, int neighboring, int currmode,
                        int stats);
